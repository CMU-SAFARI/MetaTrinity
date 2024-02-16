/*
 * Copyright (c) <2016 - 2020>, Bilkent University and ETH Zurich 
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the names of the Bilkent University, ETH Zurich,
 *   nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Authors: 
  Mohammed Alser
	  mohammedalser AT bilkent DOT edu DOT tr
  Date:
  December 3rd, 2016
*/


/* Include Files */
#include "MAGNET_DC.h"
#include "stdio.h" 

/* Function Definitions */
int MAGNET_DC(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode)
{
	int Accepted=0;
	//int *HammingMask;
	int i0;
	int ii;
	int n;
	int e;
	int d0;
	int i;
	int EarlyTermination=1;
	
	// variables for longest run detection
	int count;
	int * MAGNETMask;
	int dd;

	/* 20-Sept-2016 */
	/*  Building the Hamming masks */
	ii = ((2 * ErrorThreshold) + 1) * ReadLength;
	int * HammingMask;

	HammingMask = (int *)calloc(ii, sizeof(int));
    MAGNETMask = (int *)calloc(ReadLength, sizeof(int));
    /*
	if(!HammingMask)
	{
		printf("Out of Memory !!\n");
		exit(1);
	}
	*/	
	/*for (n = 0; n < ii; n++) {
		HammingMask[n] = 0;
	}
	
	
	//  Original Hamming Mask
	for (n = 0; n < ReadLength; n++) {
		if (ReadSeq[n]== RefSeq[n])
			HammingMask[n] = 0;
		else if (ReadSeq[n]!= RefSeq[n])
			HammingMask[n] = 1;
	}
	
	//  Incremental-step right shifted Hamming Masks (Deletion)
	for (e = 0; e < ErrorThreshold; e++) {
		// fill the shifted chars with Zeros
		for (i0 = 0; i0 <= e; i0++) {
		  HammingMask[((e +1) * ReadLength) + i0] = 0;
		}
		 ii=e+1;
		//  shift the Read chars and compare
		for (n = 0; n < (ReadLength-ii); n++) {
			if (ReadSeq[n]== RefSeq[n+ii]){
				HammingMask[((e+1) * ReadLength) + n + ii] = 0;}
			else if (ReadSeq[n]!= RefSeq[n+ii])
				HammingMask[((e+1) * ReadLength) + n + ii] = 1;
		}
	}
	
	
	//  Incremental-step left shifted Hamming Masks (Insertion)
	for (e = 0; e < ErrorThreshold; e++) {
		// fill the shifted chars with Zeros
		for (i0 = 0; i0 <= e; i0++) {
			HammingMask[((e + ErrorThreshold+1) * ReadLength) + (ReadLength-i0)-1] = 0;
		}
		ii=e+1;
		//  shift the Read chars and compare
		for (n = 0; n < ReadLength-e-1; n++) {
			//printf("%c",ReadSeq[n+ii]);
			//printf(" %c\n",RefSeq[n]);
			if (ReadSeq[n+ii]== RefSeq[n])
				HammingMask[((e + ErrorThreshold + 1) * ReadLength) + n] = 0;
			else if (ReadSeq[n+ii]!= RefSeq[n])
				HammingMask[((e + ErrorThreshold + 1) * ReadLength) + n] = 1;
		}
	}
	
	// END of Building the Hamming masks */


	// Original Hamming Mask
	for (n = (ReadLength-1); n >=0 ; n--) {
		if (ReadSeq[n]!= RefSeq[n])
			HammingMask[n] = 1;
		else
			HammingMask[n] = 0;
	}
	
	if (ErrorThreshold==0) {
		count=0;
		for (n = 0; n < ReadLength; n++) {
			if (ReadSeq[n]!= RefSeq[n])
				count++;
		}
	}
	else {

		// Shifted Hamming Masks
		for (e = 1; e <= ErrorThreshold; e++) {
			ii = e * ReadLength;
			dd = ii + (ErrorThreshold * ReadLength);
			// fill the shifted chars with Zeros
			for (i0 = e-1; i0 >=0; i0--) {
				// Deletion
				HammingMask[ii + i0] = 1;
				// Insertion
				HammingMask[dd + (ReadLength-i0)-1] = 1;
			}
			//  shift the Read chars and compare
			for (n = (ReadLength-e-1); n >=0 ; n--) {
				//  Incremental-step right shifted Hamming Masks (Deletion)
				if (ReadSeq[n]!= RefSeq[n+e])
					HammingMask[ii + n + e] = 1;
				else
					HammingMask[ii + n + e] = 0;
				// Incremental-step left shifted Hamming Masks (Insertion)
				if (ReadSeq[n+e]!= RefSeq[n])
					HammingMask[dd + n] = 1;
				else
					HammingMask[dd + n] = 0;
			}
		}
		// END of Building the Hamming masks



		/*  Looking for the longest run in each Hamming mask*/
		d0 = (2 * ErrorThreshold) + 1;
		for (i = 0; i < d0; i++) {
            ii = (i * ReadLength);
            i0 = ii + (ReadLength - 1);
            count = 0;
            for (n = ii; n <= i0; n++) {
                if (HammingMask[n] == 1)
                    count++;
            }
            if (count <= ErrorThreshold) {
                free(HammingMask);
                free(MAGNETMask);
                return count;
            }
		}
		for (n = 0; n < ReadLength; n++) {
			MAGNETMask[n] = 1;
		}
		////////////////////////////////////MAGNET here

		Extraction( ReadLength, HammingMask,  0,  ReadLength-1,  MAGNETMask,  ErrorThreshold, EarlyTermination);

		if (DebugMode==1) {
			printf("Extraction function caller= \n");
			for (n = 0; n < ReadLength; n++) {
				printf("%d",MAGNETMask[n]);
			}
			printf("\n");
		}

		count=0;
		for (n = 0; n < ReadLength; n++) {
			if (MAGNETMask[n] == 1)
				count++;
		}
	}


	Accepted = ((count-1) <= ErrorThreshold);
    free(HammingMask);
    free(MAGNETMask);
    return count;
}



int Extraction(int ReadLength, int HammingMask[], int GlobalMaxZerosStartIndex, int GlobalMaxZerosEndIndex, int MAGNETMask[], int ErrorThreshold, int EarlyTermination)
{
	int d0 = (2 * ErrorThreshold) + 1;
	int i;
	int n;
	int ii;
	int i0;
	
	// variables for longest run detection
	int LongZeros;
	int startIndex;
	int endIndex;
	int count;
	int sIndex;
	int LongestOfLongest=0;
	int endIndexMax=0;
	int startIndexMax=0;
	
	
	
	if ((GlobalMaxZerosStartIndex <= GlobalMaxZerosEndIndex) && (EarlyTermination<=ErrorThreshold+1)){
		/*if (DebugMode==1) {
			printf("\nstart %d end %d\n",GlobalMaxZerosStartIndex, GlobalMaxZerosEndIndex);
			for (n = 0; n < (((2*ErrorThreshold)+1)*ReadLength); n++) {
				if ( n % 100 == 0)
					printf("\n");
				printf("%d", HammingMask[n]  );
			}
			printf("\n");
		}*/
		LongestOfLongest=0;
		for (i = 0; i < d0; i++) {
			ii = (i*ReadLength);
			i0 = ii + GlobalMaxZerosEndIndex;

			// variables for longest run detection
			LongZeros = 0;
			startIndex = 0;
			endIndex = 0;
			count = 0;
			sIndex = GlobalMaxZerosStartIndex;

			for (n = ii+GlobalMaxZerosStartIndex; n <= i0; n++) {
				if (HammingMask[n]==0) {
					count = count+1;
					if (count > LongZeros) {
						LongZeros = count;
						startIndex = sIndex;
						endIndex = startIndex+count-1;
					}
				}
				else if (HammingMask[n]==1) {
					count = 0;
					sIndex = n+1-ii;
				}			
			}
			if (LongZeros>LongestOfLongest) {
				LongestOfLongest = LongZeros;
				endIndexMax = endIndex;
				startIndexMax = startIndex;
			}
		}
		//printf("longest: %d, from %d to %d\n",LongestOfLongest, startIndexMax, endIndexMax);
		
		//printf("MAGNETMask= \n");
		for (n = 0; n < ReadLength; n++) {
			if ( (startIndexMax<=n) && (n<=endIndexMax) ) 
				MAGNETMask[n]=0;
			//printf("%d",MAGNETMask[n]);
		}
//		printf("\n");

		EarlyTermination++;
		// Search for Longest Consecutive Zeros in the Left
		//printf("call: from %d to %d\n",GlobalMaxZerosStartIndex, startIndexMax-2);
		if ((startIndexMax-2) < GlobalMaxZerosEndIndex)
			Extraction(ReadLength, HammingMask,GlobalMaxZerosStartIndex, startIndexMax-2, MAGNETMask, ErrorThreshold, EarlyTermination);

		// Search for Longest Consecutive Zeros in the Right
		//printf("call: from %d to %d\n",endIndexMax+2, GlobalMaxZerosEndIndex);
		if ((endIndexMax+2) > GlobalMaxZerosStartIndex)
			Extraction(ReadLength, HammingMask,endIndexMax+2, GlobalMaxZerosEndIndex, MAGNETMask, ErrorThreshold, EarlyTermination);
	
		
		return count;
	}
	else
		return count;
}
