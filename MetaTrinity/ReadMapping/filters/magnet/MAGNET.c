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
#include "MAGNET.h"
#include "stdio.h" 

/* Function Definitions */
int MAGNET(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode)
{
	int Accepted=0;
	int i0;
	int ii;
	int n;
	int e;
	int d0;
	int i;
	
	// variables for longest run detection
	int LongZeros;
	int startIndex;
	int endIndex;
	int count;
	int sIndex;
	int LongestOfLongest;
	int endIndexMax=0;
	int startIndexMax=0;
//	int MAGNETMask[ReadLength];
	int k;
	int dd;

	/* 20-Sept-2016 */
	/*  Building the Hamming masks */
	ii = ((2 * ErrorThreshold) + 1) * ReadLength;

    int * MAGNETMask = (int *) calloc (ReadLength, sizeof(int));
	int * HammingMask = (int *) calloc (ii, sizeof(int));

    if(!HammingMask || !MAGNETMask)
	{
		printf("Out of Memory !!\n");
		exit(1);
	}

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
    count = 0;
	for (n = (ReadLength-1); n >=0 ; n--) {
		if (ReadSeq[n]!= RefSeq[n])
			HammingMask[n] = 1;
		else
			HammingMask[n] = 0;
	}

    for (n = 0; n < ReadLength; n++) {
        if (ReadSeq[n]!= RefSeq[n])
            count++;
    }

    Accepted = (count <= ErrorThreshold);

    if (Accepted == 0 && ErrorThreshold > 0) {
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
			ii = (i*ReadLength);
			i0 = ii + (ReadLength - 1);
			count = 0;
			for (n = ii; n <= i0; n++) {
				if (HammingMask[n]==1)
					count++;
			}
			if (count <= ErrorThreshold)
				return count;
		}
		for (n = 0; n < ReadLength; n++) {
			MAGNETMask[n] = 1;
		}
		for (k=0; k<= ErrorThreshold;k++) {
			if (DebugMode==1) {
				printf("\nBefore Amend: \n");
				for (n = 0; n < (((2*ErrorThreshold)+1)*ReadLength); n++) {
					if ( n % 100 == 0)
						printf("\n");
					printf("%d", HammingMask[n]  );
				}
				printf("\n\n\n\n");
			}
			LongestOfLongest=0;
			for (i = 0; i < d0; i++) {
				ii = (i*ReadLength);
				i0 = ii + (ReadLength - 1);

				// variables for longest run detection
				LongZeros = 0;
				startIndex = 0;
				endIndex = 0;
				count = 0;
				sIndex = ii;

				for (n = ii; n <= i0; n++) {
					//printf("%d", HammingMask[n]  );
					if (HammingMask[n]==0) {
						count = count+1;
						if (count > LongZeros) {
							LongZeros = count;
							startIndex = sIndex-ii;
							endIndex = startIndex+count-1;
						}
					}
					else if (HammingMask[n]==1) {
						count = 0;
						sIndex = n+1;
					}			
				}
				if (LongZeros>LongestOfLongest) {
					LongestOfLongest = LongZeros;
					endIndexMax = endIndex;
					startIndexMax = startIndex;
				}
				//printf("\n");
				//printf("\n");
				//printf("longest: %d, from %d to %d\n",LongestOfLongest, startIndexMax, endIndexMax);
			}

			//MAGNETMask= MAGNETMask+LongestOfLongest;
			//printf("Matches#: %d\n",MAGNETMask);
			// Extraction
			//printf("MAGNET\n");
			for (i = 0; i < d0; i++) {
				ii = (i*ReadLength);
				i0 = ii + (ReadLength - 1);

				for (n = ii; n <= i0; n++) {
					if ( (startIndexMax+ii-1<=n) && (n<=endIndexMax+ii+1) )
						HammingMask[n]=1;
					if ( (startIndexMax+ii<=n) && (n<=endIndexMax+ii) ) {
						if (i==0)
							MAGNETMask[n]=0;
					}
					//if (i==0)
						//printf("%d",MAGNETMask[n]);
				}
			}
			//printf("\n");

			if (DebugMode==1) {
				printf("\nAfter Amend: \n");
				for (n = 0; n < (((2*ErrorThreshold)+1)*ReadLength); n++) {
					if ( n % 100 == 0)
						printf("\n");
					printf("%d", HammingMask[n]  );
				}
				printf("\n\n\n\n");
			}
			// END of Extraction
		}
		//printf("longest Of Longest: %d, ",LongestOfLongest);
		/*  END of Amending 1001 or 101 into 1111 or 111, respectively */



	//	free(HammingMask);

		count=0;
		for (n = 0; n < ReadLength; n++) {
			if (MAGNETMask[n] == 1)
				count++;
		}
	}
	Accepted = (count <= ErrorThreshold);
    free(HammingMask);
    free(MAGNETMask);

  return count;
}
