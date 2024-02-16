//
// Created by mdr on 06/07/2021.
//

#include "qgram.h"
#include <stdlib.h>
#include <string.h>
#include "../../khash.h"
#include <stdio.h>

KHASH_MAP_INIT_INT(occ_hash, char)

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

static int getOccurrenceTableLength(int qGramLength) {
    return 1 << 2*qGramLength;
}


// FLASH
int qgram(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength

    int occurrenceTableLength = getOccurrenceTableLength(qGramLength);
    int * occurrenceCount = (int *) calloc(occurrenceTableLength, sizeof(int));

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        occurrenceCount[indexFromQGram(refGram, qGramLength)]++;
    }

    int threshold = qGramLength*ErrorThreshold;

    int count = 0;
    for (int i = 0; i < qGramCount; i++) {
        memcpy(readGram, &ReadSeq[i], qGramLength);
        if (occurrenceCount[indexFromQGram(readGram, qGramLength)] > 0) {
            occurrenceCount[indexFromQGram(readGram, qGramLength)]--;
        } else {
            count++;
            if (count > threshold) {
                break;
            }
        }
    }

    free(occurrenceCount);

    return count/qGramLength;
}

int qgram_hash(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    int qGramCount = ReadLength - qGramLength;
    khash_t(occ_hash) *occ_table = kh_init(occ_hash);

    char refGram[qGramLength];
    char readGram[qGramLength];

    int absent, idx;

    for (int i = 0; i < qGramCount; i++) {

        memcpy(readGram, &ReadSeq[i], qGramLength);

        idx = indexFromQGram(readGram, qGramLength);

        khint_t k;

        k = kh_get(occ_hash, occ_table, idx);

        if (k == kh_end(occ_table)) { // key is not present in hash table
            khint_t entry = kh_put(occ_hash, occ_table, idx, &absent);
//            printf("absent: %i", absent);
            kh_val(occ_table, entry) = 1;
        } else {
            kh_val(occ_table, k) += 1;
        }
    }

    for (int i = 0; i < qGramCount; i++) {

        memcpy(refGram, &RefSeq[i], qGramLength);

        idx = indexFromQGram(refGram, qGramLength);

        khint_t k;

        k = kh_get(occ_hash, occ_table, idx);

        if (k == kh_end(occ_table)) { // key is not present in hash table
            khint_t entry = kh_put(occ_hash, occ_table, idx, &absent);
            kh_val(occ_table, entry) = -1;
        } else {
            kh_val(occ_table, k) -= 1;
        }
    }

    int errorTotal = 0;
    for (khint_t k = kh_begin(occ_table); k != kh_end(occ_table); ++k) {
       errorTotal += abs(kh_val(occ_table, k));
    }

    kh_destroy(occ_hash, occ_table);

    return errorTotal/(2*qGramLength);
}

int qgram_avx2(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    return 0;
}