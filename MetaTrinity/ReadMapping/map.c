#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "ksw2.h"
#include "ksort.h"

#include "filter.h"
//////// filters ///////////
#include "filters/SneakySnake/SneakySnake.h" //changed from including the sneakysnake in this folder... maybe revert if this causes issues

#include "filters/adjacency-filter/AdjacencyFilter.h"
#include "filters/base-counting/Base_Counting.h"
#include "filters/magnet/MAGNET_DC.h"
#include "filters/hamming-distance/HD.h"
#include "filters/shouji/Shouji.h"
#include "filters/qgram/qgram.h"
#include "filters/shd/SHD.h"
// NEW
#include "filters/grim/grim.h"
#include "filters/edlib/edlib.h"
#include "filters/swift/swift.h"
#include "filters/pigeonhole/pigeonhole.h"
///////////////////////////

//int Edits = 0;
//filter_calls = 0;
// uint64_t (*seed_map)[2];

typedef struct {
    uint64_t seeds;
    uint64_t location;
} seed_map_entry_t;

#define ks_lt_seed_map_entry(a, b) ((a).seeds < (b).seeds)
KSORT_INIT(seed_map_sort, seed_map_entry_t, ks_lt_seed_map_entry);




struct mm_tbuf_s {
	void *km;
	int rep_len, frag_gap;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

void *mm_tbuf_get_km(mm_tbuf_t *b)
{
	return b->km;
}

char ** init_tabnt() {
    char tmp[] = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        int arr[] = {0, 0};
        tmp[i] = &arr;
    }
    return &tmp;
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && (int32_t)(dreg[u]>>32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v]>>32) < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v]>>32)? s : dreg[v]>>32;
				int ee = e < (int32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mm_match_t;

static mm_match_t *collect_matches(void *km, int *_n_m, int max_occ, const mm_idx_t *mi, const mm128_v *mv, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos)
{
	int rep_st = 0, rep_en = 0, n_m;
	size_t i;
	mm_match_t *m;
	*n_mini_pos = 0;
	*mini_pos = (uint64_t*)kmalloc(km, mv->n * sizeof(uint64_t));
	m = (mm_match_t*)kmalloc(km, mv->n * sizeof(mm_match_t));
	for (i = 0, n_m = 0, *rep_len = 0, *n_a = 0; i < mv->n; ++i) {
		const uint64_t *cr;
		mm128_t *p = &mv->a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
		int t;
		cr = mm_idx_get(mi, p->x>>8, &t);
		if (t >= max_occ) {
			int en = (q_pos >> 1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mm_match_t *q = &m[n_m++];
			q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
			q->is_tandem = 0;
			if (i > 0 && p->x>>8 == mv->a[i - 1].x>>8) q->is_tandem = 1;
			if (i < mv->n - 1 && p->x>>8 == mv->a[i + 1].x>>8) q->is_tandem = 1;
			*n_a += q->n;
			(*mini_pos)[(*n_mini_pos)++] = (uint64_t)q_span<<32 | q_pos>>1;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_m = n_m;
	return m;
}

static inline int skip_seed(int flag, uint64_t r, const mm_match_t *q, const char *qname, int qlen, const mm_idx_t *mi, int *is_self)
{
	*is_self = 0;
	if (qname && (flag & (MM_F_NO_DIAG|MM_F_NO_DUAL))) {
		const mm_idx_seq_t *s = &mi->seq[r>>32];
		int cmp;
		cmp = strcmp(qname, s->name);
		if ((flag&MM_F_NO_DIAG) && cmp == 0 && (int)s->len == qlen) {
			if ((uint32_t)r>>1 == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&1) == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&1) == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

static mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m, heap_size = 0;
	int64_t j, n_for = 0, n_rev = 0;
	mm_match_t *m;
	mm128_t *a, *heap;

	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);

	heap = (mm128_t*)kmalloc(km, n_m * sizeof(mm128_t));
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));

	for (i = 0, heap_size = 0; i < n_m; ++i) {
		if (m[i].n > 0) {
			heap[heap_size].x = m[i].cr[0];
			heap[heap_size].y = (uint64_t)i<<32;
			++heap_size;
		}
	}
	ks_heapmake_heap(heap_size, heap);
	while (heap_size > 0) {
		mm_match_t *q = &m[heap->y>>32];
		mm128_t *p;
		uint64_t r = heap->x;
		int32_t is_self, rpos = (uint32_t)r >> 1;
		if (!skip_seed(opt->flag, r, q, qname, qlen, mi, &is_self)) {
			if ((r&1) == (q->q_pos&1)) { // forward strand
				p = &a[n_for++];
				p->x = (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p = &a[(*n_a) - (++n_rev)];
				p->x = 1ULL<<63 | (r&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
		// update the heap
		if ((uint32_t)heap->y < q->n - 1) {
			++heap[0].y;
			heap[0].x = m[heap[0].y>>32].cr[(uint32_t)heap[0].y];
		} else {
			heap[0] = heap[heap_size - 1];
			--heap_size;
		}
		ks_heapdown_heap(0, heap_size, heap);
	}
	kfree(km, m);
	kfree(km, heap);

	// reverse anchors on the reverse strand, as they are in the descending order
	for (j = 0; j < n_rev>>1; ++j) {
		mm128_t t = a[(*n_a) - 1 - j];
		a[(*n_a) - 1 - j] = a[(*n_a) - (n_rev - j)];
		a[(*n_a) - (n_rev - j)] = t;
	}
	if (*n_a > n_for + n_rev) {
		memmove(a + n_for, a + (*n_a) - n_rev, n_rev * sizeof(mm128_t));
		*n_a = n_for + n_rev;
	}
	return a;
}

static mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, n_m;
	mm_match_t *m;
	mm128_t *a;
	m = collect_matches(km, &n_m, max_occ, mi, mv, n_a, rep_len, n_mini_pos, mini_pos);
	a = (mm128_t*)kmalloc(km, *n_a * sizeof(mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		uint32_t k;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self, rpos = (uint32_t)r[k] >> 1;
			mm128_t *p;
			if (skip_seed(opt->flag, r[k], q, qname, qlen, mi, &is_self)) continue;
			p = &a[(*n_a)++];
			if ((r[k]&1) == (q->q_pos&1)) { // forward strand
				p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;
			} else { // reverse strand
				p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	kfree(km, m);
	radix_sort_128x(a, a + (*n_a));
	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
		if (!(opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN))) // long join not working well without primary chains
			mm_join_long(km, opt, qlen, n_regs, regs, a);
	}
}

static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

// Comparison function for location sorting
int compare(const void *a, const void *b) {
    uint64_t *pa = (uint64_t*)a;
    uint64_t *pb = (uint64_t*)b;
    return (*pb > *pa) - (*pb < *pa);  // sort in descending order
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	int i, rep_len, qlen_sum, n_regs0, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	int64_t n_a=0;


	uint64_t *u, *mini_pos;
	mm128_t *a;
	mm128_v mv = {0,0,0};
	mm_reg1_t *regs0;
	km_stat_t kmst;

	int printSAM = (mm_dbg_flag & MM_DBG_PRINT_ALN_SEQ)?1:0;
	int Seed_Num=1;
	int32_t mappingStartingPosition=0;
	uint64_t mapStartPos = 0;
	uint64_t mapEndPos = 0;
	// uint64_t mapStartPos2 = 0;
	// uint64_t mapEndPos2 = 0;
	int ReadNo=0;
	int Accepted=0,Rejected=0, MappedReadNo=0;

	int locations_per_read = (opt->min_cnt >= 1)?opt->min_cnt:2;
	int MAX_NUM_MAPPING_LOCATION_PER_READ = 5000;
	// int MAX_NUM_MAPPING_LOCATION_PER_READ = 500000; // same workload for each filter (unclear if this affects only the filter workload... could change runtime in unpredictable ways!)
	// int SEED_NUM_PER_READ = (opt->min_cnt >= 1)?opt->min_cnt:2; // used originally
	uint64_t MIN_SEED_NUM_PER_READ = 2;
	int SEEDS_NUM_FOR_NO_FILTERING =3;
    int editWindow = 0;

	int SSWindow = qlens[0];
	int SSIteration = qlens[0];
	int SSEditThreshold = opt->bw;
	int MIN_SCORE=50;
	n_regs[0] = 0;
	int sc_mch=1;
	int sc_mis=-2;
	int gapo=2;
	int gape=1;
	int j, g = sc_mch, bb = sc_mis < 0? sc_mis : -sc_mis; // g>0 and b<0
	int8_t mat[25] = { g,bb,bb,bb,0, bb,g,bb,bb,0, bb,bb,g,bb,0, bb,bb,bb,g,0, 0,0,0,0,0 };

	int aa = 0;
	int bb2 = -1;
	int8_t score_matrix[25] = { aa,bb2,bb2,bb2,0, bb2,aa,bb2,bb2,0, bb2,bb2,aa,bb2,0, bb2,bb2,bb2,aa,0, 0,0,0,0,0 }; //make alignment score = edit distance
	int tl;
	int ql = strlen(seqs[0]);
	int Edits = 0;

	EdlibAlignResult resultEdlib;

	// only for ksw2
	ksw_extz_t ez;
	int8_t *junc;



	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;
	if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
	if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
	else a = collect_seed_hits(b->km, opt, opt->mid_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);

	if (mm_dbg_flag & MM_DBG_PRINT_SEED) { // remove if-statement to print seed hits
		fprintf(stderr, "RS\t%d\n", rep_len);
		for (i = 0; i < n_a; ++i)
			printf("SD\t%s\t%d\t%c\t%d\t%d\t%d\n", mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
					i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
					// print output: 1)RG name 2)seed 3)strand 4)read | last col important: Nr of zeros = Nr of matching seeds in same reagion
	}

	// ALLOCATE MEMORY for seed_map
	// uint64_t (*seed_map)[2];
	// seed_map = (uint64_t(*)[2])calloc(1<<30, sizeof(*seed_map));
	// seed_map = (uint64_t(*)[2])calloc(n_a, sizeof(*seed_map));
	// 	if (seed_map == NULL) {
	// 		fprintf(stderr, "Memory allocation failed\n");
	// 		exit(1);
	// 	}


/*
	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;



	a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score, opt->chain_gap_scale, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);

	if (opt->max_occ > opt->mid_occ && rep_len > 0) {
		int rechain = 0;
		if (n_regs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < n_regs0; ++i) { // find the best chain
				if (max < (int)(u[i]>>32)) max = u[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)u[i];
			}
			for (i = 1; i < (int32_t)u[max_i]; ++i) // count the number of segments in the best chain
				if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(b->km, a);
			kfree(b->km, u);
			kfree(b->km, mini_pos);
			if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			else a = collect_seed_hits(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
			a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->max_chain_iter, opt->min_cnt, opt->min_chain_score, opt->chain_gap_scale, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
		}
	}
	b->frag_gap = max_chain_gap_ref;
	b->rep_len = rep_len;
	*/
	// populate r[]
	mm_reg1_t *r;
	r = (mm_reg1_t*)calloc(n_a, sizeof(mm_reg1_t));

	// uint64_t (*seed_map)[2];
	// seed_map = malloc(n_a * sizeof(*seed_map)); // Allocate memory for the 2D array

	// uint64_t (*seed_map)[2];
	// seed_map = (uint64_t(*)[2])calloc(n_a, sizeof(*seed_map));
	
	// typedef struct {
    // uint64_t seeds;
    // uint64_t location;
	// } seed_map_entry_t;
	// #define ks_lt_seed_map_entry(a, b) ((a).seeds < (b).seeds)
	// KSORT_INIT(seed_map_sort, seed_map_entry_t, ks_lt_seed_map_entry);

	seed_map_entry_t* seed_map;
	seed_map = (seed_map_entry_t*)calloc(n_a, sizeof(seed_map_entry_t));





	if (n_segs >= 1) { // uni-segment
		//////////////////////////////////////////////////////////////////
		if (n_a>2) {
			// uint64_t (*seed_map)[2];
			// seed_map = malloc(n_a * sizeof(*seed_map)); // Allocate memory for the 2D array
			
			for (i = 1; i < n_a; ++i){
				// if (Accepted<MAX_NUM_MAPPING_LOCATION_PER_READ){
					if ((i<n_a) && (((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x) ==0)) {
						if (Seed_Num == 1) {
							mappingStartingPosition = ((int32_t)a[i-1].x - (int32_t)a[i-1].y);
							// printf("mappingStartingPosition\t%d\n", mappingStartingPosition);
						}
						Seed_Num++;
					}
					else {

						// assign to 2d array: 1st col contains number of seed hits at mapping location | 2end col contains mapping location
						seed_map[i].seeds = (uint64_t)Seed_Num;
        				// seed_map[i][1] = mi->seq[a[i-1].x<<1>>33].offset + (uint64_t)mappingStartingPosition;
						mapStartPos = mi->seq[a[i-1].x<<1>>33].offset + (uint64_t)mappingStartingPosition;
						mapEndPos = (mapStartPos + (uint64_t)qlens[0]);

						seed_map[i].location = (uint64_t)mapStartPos;
						// seed_map[i][1] = (uint64_t)mappingStartingPosition;

						if (i<n_a) {
							Seed_Num=1;
							mappingStartingPosition = ((int32_t)a[i].x - (int32_t)a[i].y);
						}
					}

					//if (mm_dbg_flag & MM_DBG_PRINT_SEED)
						//fprintf(stderr,"%d\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", Seed_Num, mappingStartingPosition, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff), i == 0? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));
				// }
				// else i=n_a;
			}
			
			// mergesort(seed_map, n_a, sizeof(*seed_map), compare);
			qsort(seed_map, n_a, sizeof(*seed_map), compare);
			// ks_heapmake_heap // ksmall

			// ks_heapmake_seed_map_sort(n_a, seed_map);
			// rs_insertsort_seed_map_sort(n_a, seed_map); 


			if (mm_dbg_flag & MM_DBG_PRINT_SEED) {
				printf("\n\nSeed Num | Map Start Loc\n");
				printf("------------------------\n");
				for (i = 1; i < n_a; ++i){
					printf("\t%d \t | \t%d\n", seed_map[i].seeds, seed_map[i].location);
				}
				printf("------------------------\n");
				printf("\n============= CURRENT READ DONE ============= \n\n\n");
			}


			for (int i = 1; i <= locations_per_read; ++i){ // consider best N mapping locations per read
				
				if (((uint64_t)seed_map[i].seeds >= (uint64_t)MIN_SEED_NUM_PER_READ)) {

					// fprintf(stderr,"mapStartPos BEFORE: %d\n",mapStartPos);
					// fprintf(stderr,"mapEndPos BEFORE: %d\n\n",mapEndPos);
			
					mapStartPos = (uint64_t)seed_map[i].location; // take first N locations from table
					// mapStartPos = mi->seq[a[i-1].x<<1>>33].offset + (uint64_t)mappingStartingPosition;
					// mapStartPos = mi->seq[a[i-1].x<<1>>33].offset + (uint64_t)seed_map[i][1];
					mapEndPos = (mapStartPos + (uint64_t)qlens[0]); // Endposition based on read length

					// fprintf(stderr,"mapStartPos AFTER: %d\n",mapStartPos);
					// fprintf(stderr,"mapEndPos AFTER: %d\n\n",mapEndPos);

					MappedReadNo=MappedReadNo+1;

					char RefSeq[qlens[0]];
						if ((a[i-1].x>>63)==0) {
							for (uint64_t mapSeqI = mapStartPos; mapSeqI < mapEndPos; ++mapSeqI)
								RefSeq[mapSeqI-(mapStartPos)]="ACGTN"[mm_seq4_get(mi->S, mapSeqI)];
						} 
						else{
							int ttttt=0;
							for (uint64_t mapSeqI = mapEndPos-1 ; mapSeqI >= mapStartPos; mapSeqI=mapSeqI-1) {
								RefSeq[ttttt]="TGCAN"[mm_seq4_get(mi->S, mapSeqI)];
								ttttt=ttttt+1;
							}
						}
						
						// ================= CALL FILTERS ==========================
						//printf("FILTER_CALLED");

						if (strcmp(filter, "adjacency-filter")==0) {
							Edits = AdjacencyFilter(qlens[0], RefSeq, seqs[0], SSEditThreshold, 5, 0);
							filter_calls++;
						}
						else if (strcmp(filter, "base-counting")==0) {
							Edits = baseCounting(qlens[0], RefSeq, seqs[0], SSEditThreshold, 0);
							filter_calls++;
						}
						else if (strcmp(filter, "magnet")==0) {
							Edits = MAGNET_DC(qlens[0], RefSeq, seqs[0], SSEditThreshold, 0);
							filter_calls++;
						}
						else if (strcmp(filter, "sneakysnake")==0) {
							Edits = SneakySnake(qlens[0], RefSeq, seqs[0], SSEditThreshold, qlens[0], 0, qlens[0]);
							filter_calls++;
						}
						else if (strcmp(filter, "hd")==0) {
							Edits = HD(qlens[0], RefSeq, seqs[0], SSEditThreshold, 0);
							filter_calls = filter_calls + 1;
						}
						else if (strcmp(filter, "shouji")==0) {
							Edits = Shouji(qlens[0], RefSeq, seqs[0], SSEditThreshold, 4, 0);
							filter_calls++;
						}
						else if (strcmp(filter, "qgram")==0) {
							Edits = qgram(qlens[0], RefSeq, seqs[0], SSEditThreshold, 5);
							filter_calls++;
						}
						else if (strcmp(filter, "shd")==0) {
							Edits = SHD(qlens[0], RefSeq, seqs[0], SSEditThreshold, 0);
							filter_calls++;
						}
						else if (strcmp(filter, "qgram_hash")==0) {
							Edits = qgram_hash(qlens[0], RefSeq, seqs[0], SSEditThreshold, 5);
							filter_calls++;
						}
						else if (strcmp(filter, "grim_original")==0) {
							Edits = grim_original(qlens[0], RefSeq, seqs[0], SSEditThreshold, 5);
							filter_calls++;
						}
						else if (strcmp(filter, "grim_original_tweak")==0) {
							Edits = grim_original_tweak(qlens[0]  , RefSeq  , seqs[0]  , SSEditThreshold, 5);
							filter_calls++;
						}
						else if (strcmp(filter, "edlib")==0) {
							resultEdlib = edlibAlign(RefSeq,  qlens[0],       seqs[0],    qlens[0],       edlibNewAlignConfig(SSEditThreshold, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0) );
							Edits = resultEdlib.editDistance;
							filter_calls++;
						}
						else if (strcmp(filter, "ksw2")==0) {
							uint8_t *ts, *qs, c[256];
							tl =strlen(RefSeq);
							memset(&ez, 0, sizeof(ksw_extz_t));
							memset(c, 4, 256);
							c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
							c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
							ts = (uint8_t*)malloc(tl);
							qs = (uint8_t*)malloc(ql);
							for (j = 0; j < tl; ++j) ts[j] = c[(uint8_t)RefSeq[j]]; // encode to 0/1/2/3
							for (j = 0; j < ql; ++j) qs[j] = c[(uint8_t)seqs[0][j]];

							//ksw_extz2_sse(void *km,   int qlen,   const uint8_t *query,   int tlen,   const uint8_t *target,  int8_t m,   const int8_t *mat,  int8_t q,   int8_t e,    int w,   int zdrop,  int end_bonus,  int flag             , ksw_extz_t *ez)
								ksw_extz2_sse(0       ,   ql      ,   qs                  ,   tl      ,   ts                   ,  5       ,   score_matrix     ,  0       ,   1       ,    -1   ,   -1        ,  0            , KSW_EZ_EXTZ_ONLY     , &ez);
							Edits = abs(ez.score);

							//free(ez.cigar);
							free(ts); free(qs);
							filter_calls++;
						}
						// ================= END CALL FILTERS ==========================


						if(Edits > -1){
							if( (Edits >= 0) && (SSEditThreshold >= 0) && (Edits <= SSEditThreshold) ) {
								
								// if (Accepted<MAX_NUM_MAPPING_LOCATION_PER_READ) {
									mm_reg1_t *ri = &r[Accepted];
									uint32_t *cigar;
									uint32_t n_cigar;
									int32_t dp_max2; //Edits in cigar in minimap.h

									// new dummy CIGAR
									cigar = (uint32_t[]){((Edits)<<4)|0x1, ((qlens[0]-(Edits))<<4)&0xFFFFFFF0};
									n_cigar = 2;
									//printf("Edits CIGAR (map.c): %d\n", Edits); //verbose
									//printf("2E%dM\n", qlens[0]-(Edits));

									mm_extra_t *p;
									if (n_cigar >= 0) {
										if (ri->p == 0) {
											uint32_t capacity = n_cigar + sizeof(mm_extra_t)/4;
											kroundup32(capacity);
											ri->p = (mm_extra_t*)calloc(capacity, 4);
											ri->p->capacity = capacity;
											ri->p->n_cigar = n_cigar;
											ri->p->dp_max2 = Edits;
										} else if (ri->p->n_cigar + ri->p->n_cigar + sizeof(mm_extra_t)/4 > ri->p->capacity) {
											ri->p->capacity = ri->p->n_cigar + ri->p->n_cigar + sizeof(mm_extra_t)/4;
											kroundup32(ri->p->capacity);
											ri->p = (mm_extra_t*)realloc(ri->p, ri->p->capacity * 4);
											ri->p->dp_max2 = Edits;
										}
										p = ri->p;
											memcpy(p->cigar, cigar, (ri->p->n_cigar )* 4);  //PROBLEM --> fixed
											//printf("===== memcpy works ===== \n");
											//////// mm_update_extra
											uint32_t k, l;
											int32_t s = 0, max = 0, qshift, tshift, toff = 0, qoff = 0;
											ri->blen = ri->mlen = 0;

											int n_ambi = 0, n_diff = 0;
											for (l = 0; l < qlens[0]; ++l) {
												if (RefSeq[l] =='N' || seqs[0][l] =='N') ++n_ambi;
												else if (RefSeq[l] != seqs[0][l]) ++n_diff;
											}
											ri->blen += qlens[0]-Edits - n_ambi, ri->mlen += qlens[0]-Edits - (n_ambi + n_diff), ri->p->n_ambi += n_ambi;
											ri->p->dp_max = opt->q + opt->e * (qlens[0]-Edits);
											ri->score = ri->score0 = ri->p->dp_max;
									}


									int count_cigar = 0;


									ri->id = Accepted;
									ri->parent = Accepted;
									ri->score = ri->score0 = ql; //TODO
									ri->hash = (uint32_t)a[i].x;
									ri->cnt = (int32_t)Seed_Num;
									ri->as = a[i].y >> 32;
									ri->div = -1.0f;
									ri->qs = 0;
									ri->qe = ql;
									ri->rs = mappingStartingPosition;
									//ri->re = mapEndPos;
									ri->rev = a[i-1].x>>63;
									ri->rid = a[i-1].x<<1>>33;
									Accepted++;
								// }
							}

						}
						else {
							Rejected++;
						}

						// if (i<n_a) {
							// Seed_Num=1;
							// mappingStartingPosition = ((int32_t)a[i].x - (int32_t)a[i].y);
						// } 	
				}			
			} 	

			
		}


		/////////
		int32_t max_score=0;
		int TMP_ID_PARENT_PRI=0;
		int n_pri = 0;
		if (Accepted >0) {
			for (i = 0; i < Accepted; ++i) {
				mm_reg1_t *ri = &r[i];
				if(ri->score > max_score) {
					max_score =ri->score;
					TMP_ID_PARENT_PRI =i;
				}
				//printf("%s\t%d\t%d\n",qname,ri->score, TMP_ID_PARENT_PRI);
			}

			for (i = 0; i < Accepted; ++i) {
				mm_reg1_t *ri = &r[i];
				ri->parent = TMP_ID_PARENT_PRI;

				if (ri->id == ri->parent) {
					++n_pri;
					ri->sam_pri = (n_pri == 1);
					//printf("%s\t%d\t%d\t%d\t%d\n",qname, ri->id, ri->parent,ri->score, ri->sam_pri);
				} else
					ri->sam_pri = 0;
				//printf("%s\t%d\t%d\t%d\t%d\n",qname, ri->id, ri->parent,ri->score, ri->sam_pri);
				/*if (r->parent == MM_PARENT_TMP_PRI)
					r->parent = i;
				else if (r->parent >= 0 && tmp[r->parent] >= 0)
					r->parent = tmp[r->parent];
				else r->parent = MM_PARENT_UNSET;
				*/
			}
		}

		// if (Accepted<MAX_NUM_MAPPING_LOCATION_PER_READ)
		    r = (mm_reg1_t*)realloc(r, (Accepted) * sizeof(mm_reg1_t));

		n_regs[0] = Accepted, regs[0] = r;




		//////////////////////////////////////////////////////////////////
		// populate r[]
		/*if (n_segs == 1) { // uni-segment

			mm_reg1_t *r;
			r = (mm_reg1_t*)calloc(n_a, sizeof(mm_reg1_t));
			for (i = 0; i < n_a; ++i) {
				mm_reg1_t *ri = &r[i];
				ri->id = i;
				ri->parent = MM_PARENT_UNSET;
				ri->score = ri->score0 = a[i].x >> 32;
				ri->hash = (uint32_t)a[i].x;
				ri->cnt = (int32_t)a[i].y;
				ri->as = a[i].y >> 32;
				ri->div = -1.0f;
			}
			n_regs[0] = n_a, regs[0] = r;
		*/
	}
	else if(n_segs>0){ // multi-segment
		n_regs0=n_a;
		/*mm_reg1_t *r;
		r = (mm_reg1_t*)calloc(n_regs0, sizeof(mm_reg1_t));
		for (i = 0; i < n_segs; ++i) {
			mm_reg1_t *regs0 = &r[i];
			regs0->id = i;
			regs0->parent = MM_PARENT_UNSET;
			regs0->score = regs0->score0 = z[i].x >> 32;
			regs0->hash = (uint32_t)z[i].x;
			regs0->cnt = (int32_t)z[i].y;
			regs0->as = z[i].y >> 32;
			regs0->div = -1.0f;
		}
		n_regs[0] = n_regs0, regs[i] = regs0;*/
	}

	/*
	// minimap2 alignment
	typedef struct {
		uint32_t capacity;                  // the capacity of cigar[]
		int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
		uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
		uint32_t n_cigar;                   // number of cigar operations in cigar[]
		uint32_t cigar[];
	} mm_extra_t;

	typedef struct {
		int32_t id;             // ID for internal uses (see also parent below)
		int32_t cnt;            // number of minimizers; if on the reverse strand
		int32_t rid;            // reference index; if this is an alignment from inversion rescue
		int32_t score;          // DP alignment score
		int32_t qs, qe, rs, re; // query start and end; reference start and end
		int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
		int32_t as;             // offset in the a[] array (for internal uses only)
		int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
		int32_t n_sub;          // number of suboptimal mappings
		int32_t score0;         // initial chaining score (before chain merging/spliting)
		uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, is_alt:1, dummy:6;
		uint32_t hash;
		float div;
		mm_extra_t *p;
	} mm_reg1_t;
	*/

	/*//regs0 = mm_gen_regs(b->km, hash, qlen_sum, n_regs0, u, a);
	if (mi->n_alt) {
		mm_mark_alt(mi, n_regs0, regs0);
		mm_hit_sort(b->km, &n_regs0, regs0, opt->alt_drop); // this step can be merged into mm_gen_regs(); will do if this shows up in profile
	}

	if (mm_dbg_flag & MM_DBG_PRINT_SEED)
		for (j = 0; j < n_regs0; ++j)
			for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
				fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
						i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

	chain_post(opt, max_chain_gap_ref, mi, b->km, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
	if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_mini_pos, mini_pos);

	if (n_segs == 1) { // uni-segment
		regs0 = align_regs(opt, mi, b->km, qlens[0], seqs[0], &n_regs0, regs0, a);
		mm_set_mapq(b->km, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
		n_regs[0] = n_regs0, regs[0] = regs0;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(b->km, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
		free(regs0);
		for (i = 0; i < n_segs; ++i) {
			mm_set_parent(b->km, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, b->km, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq(b->km, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
		}
		mm_seg_free(b->km, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(b->km, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}
	*/
	kfree(b->km, mv.a);
	kfree(b->km, a);
	kfree(b->km, u);
	kfree(b->km, mini_pos);
	free(seed_map);

	if (b->km) {
		km_stat(b->km, &kmst);
		if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
			fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
		assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
		if (kmst.largest > 1U<<28) {
			km_destroy(b->km);
			b->km = km_init();
		}
	}
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int n_processed, n_threads, n_fp;
	int64_t mini_batch_size;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;

	int n_parts;
	uint32_t *rid_shift;
	FILE *fp_split, **fp_parts;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid) // kt_for() callback
{
    step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	mm_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) { //no pairing
		for (j = 0; j < s->n_seg[i]; ++j) {
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name);
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name);
		for (j = 0; j < s->n_seg[i]; ++j) {
			s->rep_len[off + j] = b->rep_len;
			s->frag_gap[off + j] = b->frag_gap;
		}
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
			}
		}

}

static void merge_hits(step_t *s)
{
	int f, i, k0, k, max_seg = 0, *n_reg_part, *rep_len_part, *frag_gap_part, *qlens;
	void *km;
	FILE **fp = s->p->fp_parts;
	const mm_mapopt_t *opt = s->p->opt;

	km = km_init();
	for (f = 0; f < s->n_frag; ++f)
		max_seg = max_seg > s->n_seg[f]? max_seg : s->n_seg[f];
	qlens = CALLOC(int, max_seg + s->p->n_parts * 3);
	n_reg_part = qlens + max_seg;
	rep_len_part = n_reg_part + s->p->n_parts;
	frag_gap_part = rep_len_part + s->p->n_parts;
	for (f = 0, k = k0 = 0; f < s->n_frag; ++f) {
		k0 = k;
		for (i = 0; i < s->n_seg[f]; ++i, ++k) {
			int j, l, t, rep_len = 0;
			qlens[i] = s->seq[k].l_seq;
			for (j = 0, s->n_reg[k] = 0; j < s->p->n_parts; ++j) {
				mm_err_fread(&n_reg_part[j],    sizeof(int), 1, fp[j]);
				mm_err_fread(&rep_len_part[j],  sizeof(int), 1, fp[j]);
				mm_err_fread(&frag_gap_part[j], sizeof(int), 1, fp[j]);
				s->n_reg[k] += n_reg_part[j];
				if (rep_len < rep_len_part[j])
					rep_len = rep_len_part[j];
			}
			s->reg[k] = CALLOC(mm_reg1_t, s->n_reg[k]);
			for (j = 0, l = 0; j < s->p->n_parts; ++j) {
				for (t = 0; t < n_reg_part[j]; ++t, ++l) {
					mm_reg1_t *r = &s->reg[k][l];
					uint32_t capacity;
					mm_err_fread(r, sizeof(mm_reg1_t), 1, fp[j]);
					r->rid += s->p->rid_shift[j];
					if (opt->flag & MM_F_CIGAR) {
						mm_err_fread(&capacity, 4, 1, fp[j]);
						r->p = (mm_extra_t*)calloc(capacity, 4);
						r->p->capacity = capacity;
						mm_err_fread(r->p, r->p->capacity, 4, fp[j]);
					}
				}
			}
			mm_hit_sort(km, &s->n_reg[k], s->reg[k], opt->alt_drop);
			mm_set_parent(km, opt->mask_level, s->n_reg[k], s->reg[k], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL, opt->alt_drop);
			if (!(opt->flag & MM_F_ALL_CHAINS)) {
				mm_select_sub(km, opt->pri_ratio, s->p->mi->k*2, opt->best_n, &s->n_reg[k], s->reg[k]);
				mm_set_sam_pri(s->n_reg[k], s->reg[k]);
			}
			mm_set_mapq(km, s->n_reg[k], s->reg[k], opt->min_chain_score, opt->a, rep_len, !!(opt->flag & MM_F_SR));
		}
		if (s->n_seg[f] == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(km, frag_gap_part[0], opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, &s->n_reg[k0], &s->reg[k0]);
	}
	free(qlens);
	km_destroy(km);
}

static inline const mm_reg1_t *get_sam_pri(int n_regs, const mm_reg1_t *regs){
	int i;
	for (i = 0; i < n_regs; ++i)
		if (regs[i].sam_pri)
			return &regs[i];
	assert(n_regs == 0);
	return NULL;
}

static void *worker_pipeline(void *shared, int step, void *in){
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(5 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off, n_seg, rep_len and frag_gap are allocated together with n_reg
			s->n_seg = s->seg_off + s->n_seq;
			s->rep_len = s->n_seg + s->n_seq;
			s->frag_gap = s->rep_len + s->n_seq;
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
		if (p->n_parts > 0) merge_hits((step_t*)in);
		else kt_for(p->n_threads, worker_for, in, ((step_t*)in)->n_frag);
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {  //going over the reads, read by read
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				if (p->opt->split_prefix && p->n_parts == 0) { // then write to temporary files
					mm_err_fwrite(&s->n_reg[i],    sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->rep_len[i],  sizeof(int), 1, p->fp_split);
					mm_err_fwrite(&s->frag_gap[i], sizeof(int), 1, p->fp_split);
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];
						mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
						if (p->opt->flag & MM_F_CIGAR) {
							mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
							mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
						}
					}
				}

				else if (s->n_reg[i] > 0) { // the query has at least one hit
					for (j = 0; j < s->n_reg[i]; ++j) {
						mm_reg1_t *r = &s->reg[i][j];

						assert(!r->sam_pri || r->id == r->parent);
						if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
							continue;
						if (p->opt->flag & MM_F_OUT_SAM)
							mm_write_sam3(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i], 0);
						else
							mm_write_paf3(&p->str, mi, t, r, km, p->opt->flag, s->rep_len[i]);
						mm_err_puts(p->str.s);


					}

				}

				else if ((p->opt->flag & MM_F_PAF_NO_HIT) || ((p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_SAM_HIT_ONLY))) { // output an empty hit, if requested
					if (p->opt->flag & MM_F_OUT_SAM){
						mm_write_sam3(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag, s->rep_len[i], 0);
					}
					else{
						mm_write_paf3(&p->str, mi, t, 0, 0, p->opt->flag, s->rep_len[i]);
					}
					mm_err_puts(p->str.s);


				}


			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
				if (s->seq[i].comment) free(s->seq[i].comment);
			}
		}
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

static mm_bseq_file_t **open_bseqs(int n, const char **fn)
{
	mm_bseq_file_t **fp;
	int i, j;
	fp = (mm_bseq_file_t**)calloc(n, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = mm_bseq_open(fn[i])) == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j)
				mm_bseq_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	int i, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	if (opt->split_prefix)
		pl.fp_split = mm_split_init(opt->split_prefix, idx);
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);

	free(pl.str.s);
	if (pl.fp_split) fclose(pl.fp_split);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int mm_split_merge(int n_segs, const char **fn, const mm_mapopt_t *opt, int n_split_idx)
{
	int i;
	pipeline_t pl;
	mm_idx_t *mi;
	if (n_segs < 1 || n_split_idx < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt;
	pl.mini_batch_size = opt->mini_batch_size;

	pl.n_parts = n_split_idx;
	pl.fp_parts  = CALLOC(FILE*, pl.n_parts);
	pl.rid_shift = CALLOC(uint32_t, pl.n_parts);
	pl.mi = mi = mm_split_merge_prep(opt->split_prefix, n_split_idx, pl.fp_parts, pl.rid_shift);
	if (pl.mi == 0) {
		free(pl.fp_parts);
		free(pl.rid_shift);
		return -1;
	}
	for (i = n_split_idx - 1; i > 0; --i)
		pl.rid_shift[i] = pl.rid_shift[i - 1];
	for (pl.rid_shift[0] = 0, i = 1; i < n_split_idx; ++i)
		pl.rid_shift[i] += pl.rid_shift[i - 1];
	if (opt->flag & MM_F_OUT_SAM)
		for (i = 0; i < (int32_t)pl.mi->n_seq; ++i)
			printf("@SQ\tSN:%s\tLN:%d\n", pl.mi->seq[i].name, pl.mi->seq[i].len);

	kt_pipeline(2, worker_pipeline, &pl, 3);

	free(pl.str.s);
	mm_idx_destroy(mi);
	free(pl.rid_shift);
	for (i = 0; i < n_split_idx; ++i)
		fclose(pl.fp_parts[i]);
	free(pl.fp_parts);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	mm_split_rm_tmp(opt->split_prefix, n_split_idx);
	return 0;
}
void merge(void* base, size_t size, int (*cmp)(const void*, const void*), int p, int q, int r) {
    int n1 = q - p + 1;
    int n2 = r - q;
    
    void* L = malloc(n1 * size);
    void* R = malloc(n2 * size);
    
    memcpy(L, (char*)base + p * size, n1 * size);
    memcpy(R, (char*)base + (q + 1) * size, n2 * size);
    
    for (int k = p, i = 0, j = 0; k <= r; k++) {
        if (i < n1 && (j >= n2 || cmp((char*)L + i * size, (char*)R + j * size) <= 0)) {
            memcpy((char*)base + k * size, (char*)L + i * size, size);
            i++;
        } else {
            memcpy((char*)base + k * size, (char*)R + j * size, size);
            j++;
        }
    }
    
    free(L);
    free(R);
}

/* void mergesort_helper(void* base, size_t size, int (*cmp)(const void*, const void*), int p, int r) {
    if (p < r) {
        int q = (p + r) / 2;
        mergesort_helper(base, size, cmp, p, q);
        mergesort_helper(base, size, cmp, q + 1, r);
        merge(base, size, cmp, p, q, r);
    }
}

void mergesort(void* base, size_t nitems, size_t size, int (*cmp)(const void*, const void*)) {
    if (nitems > 0) {
        mergesort_helper(base, size, cmp, 0, nitems - 1);
    }
} */
