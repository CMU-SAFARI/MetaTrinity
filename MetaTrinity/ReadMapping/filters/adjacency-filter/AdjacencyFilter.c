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
#include "AdjacencyFilter.h"
#include "stdio.h"

/* Function Definitions */

/*
 * Arguments    : const char RefSeq[ReadLength]
 *                const char ReadSeq[ReadLength]
 *                int ErrorThreshold
 *                int GridSize
 * Return Type  : int
 */
int AdjacencyFilter(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int KmerSize, int DebugMode)
{
	int final_count=0;
	int n;
	int e;
	int count;
	int i;
	int dd;
	int KmerMatch;
	int LongMatch;
	/* 20-Sept-2016 */
	if (ErrorThreshold==0) {
		count=0;
		for (n = 0; n < ReadLength; n++) {
			if (ReadSeq[n]!= RefSeq[n])
				count++;
		}
		final_count = count;
	}
	else {
		// split a read into Kmer, read length and Kmer size must be divisible!!
		KmerMatch=0;
		// Go through each Kmer
		for (n = 0; n < (ReadLength/KmerSize); n++) {
			// Go through each bit of each Kmer
			dd=(n*KmerSize);
			count = 0;
			LongMatch=0;
			e=0;
			if (DebugMode==1) printf("Kmer #: %d",n+1);
			//No shift and shift the Kmer to the right-side until it matches
			while ((e<=ErrorThreshold)) {
				count = 0;
				if (DebugMode==1) printf("\n");
				for (i = 0+dd; i <(KmerSize+dd) ; i++) {
					if (i+e < ReadLength && ReadSeq[i]== RefSeq[i+e])
						count ++;
					if (DebugMode==1)  printf("%c %c %d\n",ReadSeq[i],RefSeq[i+e],count);
				}
				if (DebugMode==1) printf("\n");
				if (n==((ReadLength/KmerSize)-1)) //exclude the last Kmer from shifting to right
					e=ErrorThreshold+1;
				else
					e++;
				if (count>LongMatch)
					LongMatch=count;

			}
			if(n>0){ //exclude the first Kmer from shifting to left
				e=1;
				//shift the Kmer to the left-side until it matches
				while ((e<=ErrorThreshold)) {
					count = 0;
					if (DebugMode==1) printf("\n");
					for (i = 0+dd; i <(KmerSize+dd) ; i++) {
						if (i-e >= 0 && ReadSeq[i]==RefSeq[i-e])
							count ++;
						if (DebugMode==1) printf("%c %c %d\n",ReadSeq[i],RefSeq[i-e],count);
					}
					e++;
					if (DebugMode==1) printf("\n");
					if (count>LongMatch)
						LongMatch=count;

				}
			}
			if (LongMatch==KmerSize)
				KmerMatch++;

		}


		//printf("%d\n",MinErrors);
		final_count = (ReadLength/KmerSize)-KmerMatch;
	}

	return final_count;
}
