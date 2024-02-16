//
// Created by mdr on 22/06/2021.
//

#include "swift.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "../../khash.h"
// SWIFT Filter as introduced by Rasmussen et al. (2006)

// The paper precomputes the a q-gram for the entire reference genome.
// We will calculate this on the fly.

static int indexFromQGram(char * gram, int qGramLength) { // calculate natural order index from QGram
    unsigned int index = 0;
    unsigned int offset = 1 << 2 * (qGramLength - 1); // is equal to 4^qGramLength / 4
    for (int i = 0; i < qGramLength; i++) {
        if (gram[i] == 'T') {
            index += offset;
        } else if (gram[i] == 'G') {
            index += 2 * offset;
        } else if (gram[i] == 'C') {
            index += 3 * offset;
        } // we do not offset the index for A... this saves a case
        offset /= 4;
    }
    return index;
}

static int getLookUpTableLength(int qGramLength) {
    return 1 << 2*qGramLength;
}

struct Parallelogram updateBin(struct Bin_swift r, int j, int diagonal, int tau, int w, int qGramLength, int refLength) { // function unused
    struct Parallelogram p;
    if (j - w + qGramLength > r.max) {
        if (r.count >= tau) {
            p.left = refLength - diagonal;
            p.top = r.max + qGramLength;
            p.bottom = r.max;
        }
        r.count = 0;
    }
    if (r.count == 0) {
        r.min = j;
    }
    if (r.max > j) {
        r.max = j;
        r.count++;
    }
    return p;
}

bool updateBinAndCheckParallelogram(struct Bin_swift * r, int j, int diagonal, int tau, int w, int qGramLength, int refLength) {
    if (j - w + qGramLength > r->max) {
        if (r->count >= tau) {
            return true;
        }
        r->count = 0;
    }
    if (r->count == 0) {
        r->min = j;
    }
    if (r->max > j) {
        r->max = j;
        r->count++;
    }
    return false;
}

bool checkAndResetBinAndCheckParallelogram(struct Bin_swift * r, int j, int diagonal, int tau, int w, int qGramLength, int refLength) {
    if (j - w + qGramLength > r->max) {
        if (r->count >= tau) {
            return true;
        }
        r->count = 0;
    }
    return false;
}

void buildRefQGramIndex(char * ref, int refLength, int qGramLength, unsigned int * lookupTable, unsigned short * occ) {
    unsigned int * occurrenceCount = (unsigned int *) calloc(getLookUpTableLength(qGramLength), sizeof(unsigned int));
    //iterate over all qGrams and count occurrences to populate lookupTable
    char gram[qGramLength];
    int qGramCount = refLength - qGramLength;
    for (int i = 0; i < qGramCount; i++) {
        memcpy(gram, &ref[i], qGramLength);
//        printf("%s\n", gram);
        occurrenceCount[indexFromQGram(gram, qGramLength)] += 1;
    }

    //build lookupTable
    int itemsInLookupTable = getLookUpTableLength(qGramLength);

    for (int i = 1; i < itemsInLookupTable; i++) { //idx 0 and final idx require special handling
        lookupTable[i] = lookupTable[i-1]+occurrenceCount[i-1];
    }
    // lookUpTable is built.

    // construct occurrences table
    unsigned int idx = 0;
    for (int i = 0; i < qGramCount; i++) {
        memcpy(gram, &ref[i], qGramLength);
        idx = indexFromQGram(gram, qGramLength);
//        printf("%i, %i, %i\n", sizeof(occurrenceCount)/sizeof(unsigned int), sizeof(lookupTable)/sizeof(unsigned int), idx);
        occ[lookupTable[idx] + occurrenceCount[idx] - 1] = (unsigned short) i;
        occurrenceCount[idx]--;
    }


    free(occurrenceCount);
    // done!
}

int exponentPowerTwo(int thresh) { // returns z s.t. 2^z > threshold
    if (thresh < 0) {
        return 0;
    }
    int z = 0;
    while (1 << z <= thresh) {
        z++;
    }
    return z;
}

int swift(char * read, int readLength, char * ref, int refLength, int errThreshold, int qGramLength) { // filter function
    // The paper assumes that qGramLength < ceil(1/errThreshold)
    // We calculate the parameters w & e that define the size of the parallelogram.
//    double epsilon = ((double) errThreshold)/((double) readLength);
//    int n0 = (int) ((double) readLength)*(1.0-epsilon);
////    n0 = 50; // TODO REMOVE
//    int n1 = (int) (ceil((floor(epsilon*((double) n0))+1.0)/epsilon));
//    if (errThreshold == 0) {
//        n1 = n0;
//        // TODO: handle case errThreshold == 0
//    }
////    n1 = 100;
//    int tau = (n0+1) - (int) qGramLength*floor((epsilon*((double) n0))+1.0);
//    if (tau > (n1+1) - qGramLength*floor((epsilon*((double) n1))+1.0)) {
//        tau = (n1+1) - (int) qGramLength*floor((epsilon*((double) n1))+1.0);
//    }
//    int e = (int) (((double) (2*(tau - 1) + (qGramLength - 1)))/(1.0/(epsilon)-qGramLength)); // review these WRONG
//    int w = (tau - 1) + qGramLength*(e+1);
////    w = 64;
////    e = 4;
//
////    int n0 =  qGramLength*ceil(((double) errThreshold + (double) qGramLength - 1.0)/(1.0/epsilon - (double) qGramLength))+errThreshold-1; //TODO: add formula
//
//    if (tau > (n0+1) - qGramLength*floor(epsilon*((double) n0)+1.0) ||
//        tau > (n1+1) - qGramLength*floor(epsilon*((double) n1)+1.0)) { //could be removed
//        printf("SWIFT called with invalid errThreshold %i\n",  tau);
//    }
//
//
//    if (qGramLength >= 1.0/epsilon) {
//        printf("SWIFT called with invalid qGramLength\n");
//    }

    double epsilon = 0.05;
    int n0 = 100;
    int w = 128;
    int e = 9;
    int tau = 59;
    int n1 = (int) ceil((floor(epsilon*n0) + 1)/epsilon);

    int z = exponentPowerTwo(e);


//    printf("epsilon = %f, e = %i, w = %i, n0 = %i, z = %i, n1 = %i, tau = %i\n", epsilon, e, w, n0, z, n1, tau);

        // We build the look up table for the reference genome here
    // This is different from the original paper as mentioned above.
    unsigned int * lookupTable = (unsigned int *) calloc(1 << 2*qGramLength, sizeof(unsigned int)); // TODO: could be replaced by a hash table
    unsigned short * occ = (unsigned short *) calloc(refLength, sizeof(unsigned short));
    buildRefQGramIndex(ref, refLength, qGramLength, lookupTable, occ);
//    printf("%i", 1 << 2*qGramLength);

    // allocate and initialize array of bin records "bins"
    int binsLength =  (refLength - e - (1 << z))/((1 << z) +1); // TODO: should be ceil() function!
    struct Bin_swift * bins = (Bin *)calloc(binsLength, sizeof(struct Bin_swift));
    for (int j = 0; j < refLength-qGramLength; j++) { // TODO: test different loops
        char gram[qGramLength]; // this may be illegal
        memcpy(gram, &read[j], qGramLength);

        int idx = indexFromQGram(gram, qGramLength);

        int noOfMatches = 0;
        if (idx < getLookUpTableLength(qGramLength)-1) {
            noOfMatches = lookupTable[idx+1]-lookupTable[idx];
        } else { //case that we are accessing last element in table
            noOfMatches = refLength - lookupTable[idx] - 1; // WARNING: Can introduce false positives
        }

        int b0 = 0;
        int bm = 0;
//        printf("alive %i \n", j);

        bool parallelogramExists = false; // if a parallelogram according to the specs in the paper exist, we should accept

        for (int count = 0; count < noOfMatches; count++) {
            int i = occ[idx + count];
            int d = refLength + j - i;
            int b0 = d >> z;
            int bm = b0 % binsLength;
//            printf("bins[%i]: max = %i, min = %i, count = %i\n", bm, bins[bm].max, bins[bm].min, bins[bm].count);
            parallelogramExists = parallelogramExists || updateBinAndCheckParallelogram(&bins[bm], j, b0 << z, tau, w, qGramLength, refLength);
            if ((d & (1 << (z - 1))) < e) {
                bm = (bm + binsLength - 1) % binsLength;
                parallelogramExists = parallelogramExists || updateBinAndCheckParallelogram(&bins[bm], j, (b0-1) << z, errThreshold, w, qGramLength, refLength);
            }
        }
        if ((j - e) % (1<<z) == 0) {
            b0 = (j - e) >> z;
            bm = b0 % binsLength;
            parallelogramExists = parallelogramExists || checkAndResetBinAndCheckParallelogram(&bins[bm], j, b0 << z, tau, w, qGramLength, refLength);
        }
        if (parallelogramExists) {
            return 1; // filter should accept
        }
    }
    free(bins);

    return 0; // no parallelogram found...

}