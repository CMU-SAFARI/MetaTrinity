
#include "filter.h"

int filter_reads(data_manager* dm,int SIMD_WIDTH){
	printf("filter reads \n");
	struct timeval tstart, tend, pre_start, pre_end;
	double timeuse = 0;
	size_t pass_num=0;
	
	while(dm->read_index<dm->read_num){
		
		int pos_nums=find_candidates(dm,SIMD_WIDTH);
		
		
		gettimeofday(&tstart, NULL);
		
		
		omp_set_num_threads(1);
		#pragma omp parallel for schedule(dynamic, 20)
		for(int i=0;i<pos_nums;i+=SIMD_WIDTH){
			avx_shd(dm,SIMD_WIDTH,i);
		}
		
		gettimeofday(&tend, NULL);
		
		timeuse += (tend.tv_sec - tstart.tv_sec)
			+ (double) (tend.tv_usec - tstart.tv_usec) / 1000000;

	
		//check pass number
		for(int i=0;i<pos_nums;i++)
		{
			if(dm->shd_scores[i]<=dm->error_len)
				pass_num++;
		}

	}
	
	printf("pigeonhole total candidates %lu , pass %lu candidates in %f seconds\n",dm->total_candidates,pass_num,timeuse);
	
}



int avx_shd(data_manager* dm,int SIMD_WIDTH,int pos_index){
	
	int error_len = dm->error_len;
	int read_bit_len = dm->read_bit_len;
	int read_index = dm->read_maps[pos_index];
	
	WORD *readh = dm->readh;
	WORD *readl = dm->readl;
	WORD *pos_refh = dm->pos_hbits;
	WORD *pos_refl = dm->pos_lbits;
	
	WORD* shd_scores = dm->shd_scores;
	
	//local varaibles
	__m256i vzero = _mm256_set1_epi32(0);
	
	__m256i dh,dl,pre_dh,pre_dl;
	__m256i fh,fl,pre_fh,pre_fl;
	__m256i hbit,lbit,sh,sl,shift_count1;
	
	__m256i hamh,haml;
	__m256i ham_mask, fix_ham_mask , finalBV, pre_ham[2 * error_len + 1];
	__m256i sum;
	
	for(int j = 0; j < 2 * error_len + 1; j ++){
			pre_ham[j] = vzero;
	}
	
	pre_dh = vzero;
	pre_dl = vzero;
	pre_fh = vzero;
	pre_fl = vzero;
	hbit=vzero;
	lbit=vzero;
	sum=vzero;	
	
	
	for(int j = 0; j < read_bit_len; j ++){
		
		dh = _mm256_set1_epi32(readh[read_index * read_bit_len + j]);
		dl = _mm256_set1_epi32(readl[read_index * read_bit_len + j]);	
		fh = _mm256_load_si256((__m256i*)(pos_refh +  pos_index*read_bit_len + j * SIMD_WIDTH)); 
		fl = _mm256_load_si256((__m256i*)(pos_refl +  pos_index*read_bit_len + j * SIMD_WIDTH)); 	
		
		//compute ham1
		 hamh = _mm256_xor_si256(dh, fh);
		 haml = _mm256_xor_si256(dl, fl);
		 ham_mask = _mm256_or_si256(hamh, haml);
		 flip_false_zero(ham_mask, pre_ham[0],fix_ham_mask);
		 pre_ham[0] = ham_mask;
		 finalBV=fix_ham_mask;
	
		 
		 	//read shift right
			for(int k = 1; k <= error_len; k ++){

				hbit=_mm256_slli_epi32(pre_dh, BIT_PARA - k);
			    lbit=_mm256_slli_epi32(pre_dl, BIT_PARA - k);
	
				sh=_mm256_srli_epi32(dh, k);
				sl=_mm256_srli_epi32(dl, k);
				
				sh=_mm256_or_si256(hbit, sh);
				sl=_mm256_or_si256(lbit, sl);
				
				hamh = _mm256_xor_si256(sh, fh);
				haml = _mm256_xor_si256(sl, fl);
				ham_mask = _mm256_or_si256(hamh, haml);
				
				flip_false_zero(ham_mask, pre_ham[k],fix_ham_mask);
				pre_ham[k] = ham_mask;
				
				finalBV = _mm256_and_si256(finalBV, fix_ham_mask);
			
			}
			
			//ref shift right
			for(int k = 1; k <= error_len; k ++){

				hbit=_mm256_slli_epi32(pre_fh, BIT_PARA - k);
			    lbit=_mm256_slli_epi32(pre_fl, BIT_PARA - k);
				
				sh=_mm256_srli_epi32(fh, k);
				sl=_mm256_srli_epi32(fl, k);
				
				sh=_mm256_or_si256(hbit, sh);
				sl=_mm256_or_si256(lbit, sl);
				
				hamh = _mm256_xor_si256(sh, dh);
				haml = _mm256_xor_si256(sl, dl);
				ham_mask = _mm256_or_si256(hamh, haml);
				
				flip_false_zero(ham_mask, pre_ham[error_len+k],fix_ham_mask);
				pre_ham[error_len+k] = ham_mask;
				
				finalBV = _mm256_and_si256(finalBV, fix_ham_mask);
			}
		 
			pop_count(finalBV);
			sum = _mm256_add_epi32(finalBV, sum);
			pre_dh = dh;
			pre_dl = dl;
			pre_fh = fh;
			pre_fl = fl;
			
	}//core loop done
	
		
	_mm256_store_si256((__m256i*)(shd_scores + pos_index), sum);

}

void flip_false_zero(__m256i ham_mask, __m256i pre_ham, __m256i& fix_ham_mask){
	
	__m256i ham_carry;
	__m256i SL1,SL2,SR1,SR2;
	__m256i TP1_AMEND1,TP2_AMEND1,TP2_AMEND2;

	//carry bit
	ham_carry=_mm256_slli_epi32(pre_ham,BIT_PARA-1);
	SL1=_mm256_srli_epi32(ham_mask,1);
	SL1 = _mm256_or_si256(SL1, ham_carry);
	
	ham_carry=_mm256_slli_epi32(pre_ham,BIT_PARA-2);
	SL2=_mm256_srli_epi32(ham_mask,2);
	SL2 = _mm256_or_si256(SL2, ham_carry);
	
	ham_carry=_mm256_slli_epi32(pre_ham,BIT_PARA-3);
	SR1=_mm256_srli_epi32(ham_mask,3);
	SR1 = _mm256_or_si256(SR1, ham_carry);
	
	ham_carry=_mm256_slli_epi32(pre_ham,BIT_PARA-4);
	SR2=_mm256_srli_epi32(ham_mask,4);
	SR2 = _mm256_or_si256(SR2, ham_carry);
	
	TP1_AMEND1 = _mm256_and_si256(SL1, SR1);
	TP2_AMEND1 = _mm256_and_si256(SL1, SR2);
	TP2_AMEND2 = _mm256_and_si256(ham_mask, SR1);
	
	SL2 = _mm256_or_si256(SL2, TP1_AMEND1);
	SL2 = _mm256_or_si256(SL2, TP2_AMEND1);
	fix_ham_mask = _mm256_or_si256(SL2, TP2_AMEND2);

}


void pop_count(__m256i& finalBV){
	
	//constant varaibles
	__m256i mask[3];
	__m256i ltable,shift_count;
	__m256i haml,hamh;

	mask[0] = _mm256_set1_epi32(0x0f0f0f0f);
	mask[1] = _mm256_set1_epi32(0x00ff00ff);
	mask[2] = _mm256_set1_epi32(0x0000ffff);
	
		
	ltable=_mm256_set_epi8(1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 0,1, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 0);	//xin
	//低4位和高4位shuffle并求和(不会溢出)
	haml = _mm256_and_si256(finalBV, mask[0]);	
	haml = _mm256_shuffle_epi8(ltable, haml);
	
	
	hamh = _mm256_srli_epi32(finalBV, 4);
	hamh = _mm256_and_si256(hamh, mask[0]);
	hamh = _mm256_shuffle_epi8(ltable, hamh);
	
	finalBV = _mm256_add_epi8(haml, hamh);

	//低8位与高8位求和
	haml = _mm256_and_si256(finalBV, mask[1]);
	hamh = _mm256_srli_epi32(finalBV, 8);
	hamh = _mm256_and_si256(hamh, mask[1]);
	finalBV = _mm256_add_epi16(haml, hamh);
	
	//低16位与高16位求和
	haml = _mm256_and_si256(finalBV, mask[2]);
	hamh = _mm256_srli_epi32(finalBV, 16);
	hamh = _mm256_and_si256(hamh, mask[2]);
	
	finalBV = _mm256_add_epi32(haml, hamh);
}


int find_candidates(data_manager* dm,int SIMD_WIDTH){
	
	printf("find candidates \n");
	int kmer_num = dm->read_len / dm->window_size;
	int* kmer_index = (int*) malloc(kmer_num * sizeof(int));
	int* kmer_value = (int*) malloc(kmer_num * sizeof(int));
	

	int read_start=dm->read_index;
	int pos_index=0;
	
	for(int i=read_start;i<dm->read_num;i++){
		
		int single_pos = find_read_pos(dm, i, kmer_index,kmer_value,SIMD_WIDTH,pos_index);
		
		dm->read_index++;
		
		if(pos_index>MAX_POS-SLICE_MAX_POS)
		break;
		
	}
	
	#pragma omp parallel for
	for (int i = 0; i < pos_index; i++) {
		fill_regions(dm, i,SIMD_WIDTH);
	}

	dm->total_candidates+=pos_index;
	printf("pigeonhole single pass candidates %lu\n",pos_index);

	return pos_index;
}


/*
 *经过打印，是正确的
 *	padding with A
 */
void fill_regions(data_manager* dm, int pos_id,int SIMD_WIDTH){

	int pos = dm->pos_buffers[pos_id];
	int ref_pos =  dm->pos_buffers[pos_id] / BIT_PARA;
	int residues =  dm->pos_buffers[pos_id] % BIT_PARA;

	int read_bit_len = dm->read_bit_len;

	int add = pos_id / SIMD_WIDTH;
	int off = pos_id % SIMD_WIDTH;
	int start_address = add * (SIMD_WIDTH * read_bit_len) + off;

	int pad_len = dm->read_len % BIT_PARA;

	for (int i = 0; i < read_bit_len; i++) {
		
		WORD numh = dm->refh[ref_pos + i];
		WORD numl = dm->refl[ref_pos + i];
		WORD numh2 = dm->refh[ref_pos + i + 1];
		WORD numl2 = dm->refl[ref_pos + i + 1];

		numh = numh << residues;
		numl = numl << residues;
		numh2 = numh2 >> (BIT_PARA - residues);
		numl2 = numl2 >> (BIT_PARA - residues);

		numh += numh2;
		numl += numl2;

		if (i == read_bit_len - 1) {
			
			if (pad_len != 0) {
				
				//padding with A,高位(residues)设置为1,低位是0，所以无用位置就设置为0 了
				WORD pad_bit = ((WORD)1 << (pad_len + 1)) - 1;
				pad_bit = pad_bit << (BIT_PARA - pad_len);
				numh = numh & pad_bit;
				numl = numl & pad_bit;
			}
			
		}

		 dm->pos_hbits[start_address + i * SIMD_WIDTH] = numh;
		 dm->pos_lbits[start_address + i * SIMD_WIDTH] = numl;

	}

}


int find_read_pos(data_manager* dm, int read_id, int* kmer_index,int* kmer_value,int SIMD_WIDTH,int& pos_index){

	int flag = 0;
	int single_pos = 0;
	int kmer_num = dm->read_len / dm->window_size;
	int count = 0;
	int last_pos = 0;

	for (int j = 0; j < kmer_num; j++) {
		
		int key = dm->hashVal(dm->read_content + read_id * dm->read_len+ j * dm->window_size, &flag);
		
		if (flag == 1) {
			int start = dm->hv_index[key];
			int end = dm->hv_index[key + 1];
			kmer_value[j] = end - start;
			kmer_index[j] = j;
		} else {
			kmer_value[j] = 0;
			kmer_index[j] = j;
		}
	}
	

//按照每个kmer包含的位置的个数排序
	for (int j = 0; j < kmer_num; j++) {
		for (int k = j; k > 0; k--) {
			if (kmer_value[k] < kmer_value[k - 1]) {
				
				int tmp = kmer_value[k];
				kmer_value[k] = kmer_value[k - 1];
				kmer_value[k - 1] = tmp;
				
				tmp = kmer_index[k];
				kmer_index[k] = kmer_index[k - 1];
				kmer_index[k - 1] = tmp;
			}
		}
	}

	
	for (int j = 0; j < dm->error_len + 1; j++) {

		int key = dm->hashVal(
				dm->read_content + read_id * dm->read_len
						+ kmer_index[j] * dm->window_size, &flag);

		if (flag == 1) {
			//找到每一个位置
			int start = dm->hv_index[key];
			int end = dm->hv_index[key + 1];

			for (int i = start; i < end; i++) {
				last_pos = dm->hv_pos[i] - kmer_index[j] * dm->window_size;
				dm->pos_buffers[pos_index] = last_pos;
				dm->read_maps[pos_index] = read_id;
				pos_index++;
				single_pos++;
			}
		}

	}

	int residues = pos_index % SIMD_WIDTH;
	if (residues) {
		int addition = SIMD_WIDTH - residues;
		for (int i = 0; i < addition; i++) {
			dm->pos_buffers[pos_index] = last_pos;
			dm->read_maps[pos_index] = read_id;
			pos_index++;
			single_pos++;
		}
	}
	


	return single_pos;

}

