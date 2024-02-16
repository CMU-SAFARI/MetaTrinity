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
 *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Authors: 
  Mohammed Alser
	  mohammedalser AT bilkent DOT edu DOT tr
  Date:
  December 3rd, 2016
*/


// GATEKEEPER
// Does HD first and only tries SHD if error != 0
// Get SHD from CMU-SAFARI


/* Include Files */
#include "SHD.h"
#include "stdio.h" 

//
/* Function Definitions */
int SHD(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode)
{
	int Accepted=0;
	int i0;
	int ii;
	int n;
	int e;
	int d0;
	int i;
	int MinErrors;
	int dd;
	int count;

	/* 20-Sept-2016 */
	/*  Building the Hamming masks */
	ii = ((2 * ErrorThreshold) + 1) * ReadLength;
//	int HammingMask[ii];

	int * HammingMask = (int *) calloc(ii, sizeof(int));
    int * ANDMask = (int *) calloc(ReadLength, sizeof(int));
//	if(!HammingMask)
//	{
//		printf("Out of Memory !!\n");
//		exit(1);
//	}

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
        if (ReadSeq[n] != RefSeq[n]) {
            HammingMask[n] = 1;
            count++;
        } else {
            HammingMask[n] = 0;
        }
	}

    Accepted = (count <= ErrorThreshold);

	if (Accepted == 0 && ErrorThreshold > 0) {
        // Shifted Hamming Masks
        for (e = 1; e <= ErrorThreshold; e++) {
            ii = e * ReadLength;
            dd = ii + (ErrorThreshold * ReadLength);
            // fill the shifted chars with Zeros
            for (i0 = e - 1; i0 >= 0; i0--) {
                // Deletion
                HammingMask[ii + i0] = 0;
                // Insertion
                HammingMask[dd + (ReadLength - i0) - 1] = 0;
            }
            //  shift the Read chars and compare
            for (n = (ReadLength - e - 1); n >= 0; n--) {
                //  Incremental-step right shifted Hamming Masks (Deletion)
                if (ReadSeq[n] != RefSeq[n + e])
                    HammingMask[ii + n + e] = 1;
                else
                    HammingMask[ii + n + e] = 0;
                // Incremental-step left shifted Hamming Masks (Insertion)
                if (ReadSeq[n + e] != RefSeq[n])
                    HammingMask[dd + n] = 1;
                else
                    HammingMask[dd + n] = 0;
            }
        }
        /* END of Building the Hamming masks */

        if (DebugMode == 1) {
            printf("\nBefore Amend: \n");
            for (n = 0; n < (((2 * ErrorThreshold) + 1) * ReadLength); n++) {
                if (n % 100 == 0)
                    printf("\n");
                printf("%d", HammingMask[n]);
            }
            printf("\n\n\n\n");
        }


        /*  Amending 1001 or 101 into 1111 or 111, respectively */
        d0 = (2 * ErrorThreshold) + 1;
        for (i = 0; i < d0; i++) {
            ii = (i * ReadLength);
            i0 = ii + (ReadLength - 3);
            for (n = ii; n <= i0; n++) {
                /*  pattern = '101' */
                if ((HammingMask[n] == 1) && (HammingMask[n + 1] == 0) && (HammingMask[n + 2] == 1)) {
                    HammingMask[n + 1] = 1;
                }
                /*  pattern = '1001' */
                if (n <= i0 - 1) {
                    if ((HammingMask[n] == 1) && (HammingMask[n + 1] == 0) && (HammingMask[n + 2] == 0) &&
                        (HammingMask[n + 3] == 1)) {
                        HammingMask[n + 1] = 1;
                        HammingMask[n + 2] = 1;
                    }
                }
            }
        }
        /*  END of Amending 1001 or 101 into 1111 or 111, respectively */

        if (DebugMode == 1) {
            printf("\nAfter Amend \n");
            for (n = 0; n < (((2 * ErrorThreshold) + 1) * ReadLength); n++) {
                if (n % 100 == 0)
                    printf("\n");
                printf("%d", HammingMask[n]);
            }
            printf("\n\n\n\n");
        }

        //printf("ANDing:\n");
        for (i0 = 0; i0 < ReadLength; i0++) {
            ANDMask[i0] = 1;
            for (i = 0; i < d0; i++) {
                ANDMask[i0] &= HammingMask[(i * ReadLength) + i0];
            }
            //printf("%d",ANDMask[i0]);
        }
        //printf("\n");

        /*  Speculative removal of short-matches (SRS) Count */

        MinErrors = 0;
        i = 0;

        while (i <= (ReadLength - 4)) {
            if ((ANDMask[i] == 0) && (ANDMask[i + 1] == 0) && (ANDMask[i + 2] == 0) &&
                (ANDMask[i + 3] == 0))         // '0000'
                MinErrors = MinErrors + 0;
            else if ((ANDMask[i] == 0) && (ANDMask[i + 1] == 1) && (ANDMask[i + 2] == 0) &&
                     (ANDMask[i + 3] == 1)) // '0101'
                MinErrors = MinErrors + 2;
            else if ((ANDMask[i] == 0) && (ANDMask[i + 1] == 1) && (ANDMask[i + 2] == 1) &&
                     (ANDMask[i + 3] == 0)) // '0110'
                MinErrors = MinErrors + 2;
            else if ((ANDMask[i] == 1) && (ANDMask[i + 1] == 0) && (ANDMask[i + 2] == 0) &&
                     (ANDMask[i + 3] == 1)) // '1001'
                MinErrors = MinErrors + 2;
            else if ((ANDMask[i] == 1) && (ANDMask[i + 1] == 0) && (ANDMask[i + 2] == 1) &&
                     (ANDMask[i + 3] == 0)) // '1010'
                MinErrors = MinErrors + 2;
            else if ((ANDMask[i] == 1) && (ANDMask[i + 1] == 0) && (ANDMask[i + 2] == 1) &&
                     (ANDMask[i + 3] == 1)) // '1011'
                MinErrors = MinErrors + 2;
            else if ((ANDMask[i] == 1) && (ANDMask[i + 1] == 1) && (ANDMask[i + 2] == 0) &&
                     (ANDMask[i + 3] == 1)) // '1101'
                MinErrors = MinErrors + 2;
            else
                MinErrors = MinErrors + 1;
            i = i + 4;
        }

//		free(HammingMask);

        //printf("%d\n",MinErrors);
        Accepted = (MinErrors <= ErrorThreshold);
    }

    free(HammingMask);
    free(ANDMask);
  	//return Accepted;
  	return MinErrors;
}
