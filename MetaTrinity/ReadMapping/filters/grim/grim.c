//
// Created by mdr on 10.01.22.
//

#include "grim.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

int grim_long(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength

    int occurrenceTableLength = getOccurrenceTableLength(qGramLength);
    char * ref_occurrence = (char *) calloc(occurrenceTableLength, sizeof(char)); // stores only binary values
    char * read_occurrence = (char *) calloc(occurrenceTableLength, sizeof(char)); // stores only binary values

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;
    int ref_one_bits = 0;
    int read_one_bits = 0;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        ref_occurrence[indexFromQGram(refGram, qGramLength)] = 1;

//        if (ref_occurrence[indexFromQGram(refGram, qGramLength)] == 0) { //read occurs in reference
//            ref_occurrence[indexFromQGram(refGram, qGramLength)] = 1;
//            ref_one_bits++;
//        }
//        if (read_occurrence[indexFromQGram(readGram, qGramLength)] == 0) {
//            read_occurrence[indexFromQGram(readGram, qGramLength)] = 1;
//            read_one_bits++;
//        }
    }

    int tmp = 0;
    int error_count = 0;
    for (int i = 0; i < qGramCount; i++) {
        memcpy(readGram, &ReadSeq[i], qGramLength);
        tmp = ref_occurrence[indexFromQGram(readGram, qGramLength)];
        if (tmp == 1) {
            read_one_bits++;
            ref_occurrence[indexFromQGram(readGram, qGramLength)] = 2;
        } else if (tmp == 0) {
            error_count++;
        }
    }

    //iterate over all qGrams and count occurrences to populate lookupTable
//    int matchTotal = 0;
//    int error_count = 0;
    int threshold = (read_one_bits - qGramLength*ErrorThreshold);

//    for (int i = 0; i < occurrenceTableLength; i++) {
////        matchTotal += (read_occurrence[i] == ref_occurrence[i]);
//        matchTotal += (read_occurrence[i] & ref_occurrence[i]);
////        error_count += (read_occurrence[i] & !ref_occurrence[i]);
//        if (matchTotal > threshold) {
//            break;
//        }
//    }

//    printf("%i == %i, ", read_one_bits - matchTotal, error_count);

    free(ref_occurrence);
    free(read_occurrence);
//    int threshold = ReadLength - qGramLength + 1 - qGramLength*ErrorThreshold;
//    int threshold = read_one_bits - qGramCount*ErrorThreshold;

//    return matchTotal >= ((1<<2*qGramLength) - 2*qGramLength*ErrorThreshold);
//    return matchTotal >= (ReadLength - (qGramLength-1) - qGramLength*ErrorThreshold); //PAPER
    return error_count <= qGramLength*ErrorThreshold;
//    return count >= threshold;
}


int grim(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength

    int occurrenceTableLength = getOccurrenceTableLength(qGramLength);
    char * ref_occurrence = (char *) calloc(occurrenceTableLength, sizeof(char)); // stores only binary values
    char * read_occurrence = (char *) calloc(occurrenceTableLength, sizeof(char)); // stores only binary values

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;
    int ref_one_bits = 0;
    int read_one_bits = 0;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        memcpy(readGram, &ReadSeq[i], qGramLength);

        if (ref_occurrence[indexFromQGram(refGram, qGramLength)] == 0) { //read occurs in reference
            ref_occurrence[indexFromQGram(refGram, qGramLength)] = 1;
            ref_one_bits++;
        }
        if (read_occurrence[indexFromQGram(readGram, qGramLength)] == 0) {
            read_occurrence[indexFromQGram(readGram, qGramLength)] = 1;
            read_one_bits++;
        }

    }



    //iterate over all qGrams and count occurrences to populate lookupTable
    int matchTotal = 0;
//    int error_count = 0;
    int threshold = (read_one_bits - qGramLength*ErrorThreshold);

    for (int i = 0; i < occurrenceTableLength; i++) {
//        matchTotal += (read_occurrence[i] == ref_occurrence[i]);
        matchTotal += (read_occurrence[i] & ref_occurrence[i]);
//        error_count += (read_occurrence[i] & !ref_occurrence[i]);
        if (matchTotal > threshold) {
            break;
        }
    }

//    printf("%i == %i, ", read_one_bits - matchTotal, error_count);

    free(ref_occurrence);
    free(read_occurrence);
//    int threshold = ReadLength - qGramLength + 1 - qGramLength*ErrorThreshold;
//    int threshold = read_one_bits - qGramCount*ErrorThreshold;

//    return matchTotal >= ((1<<2*qGramLength) - 2*qGramLength*ErrorThreshold);
//    return matchTotal >= (ReadLength - (qGramLength-1) - qGramLength*ErrorThreshold); //PAPER
    return matchTotal >= threshold;
//    return count >= threshold;
}

int grim_original(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength

    int occurrenceTableLength = getOccurrenceTableLength(qGramLength);
    char * ref_occurrence = (char *) calloc(occurrenceTableLength, sizeof(char)); // stores only binary values

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;
    int ref_one_bits = 0;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        if (ref_occurrence[indexFromQGram(refGram, qGramLength)] == 0) { //read occurs in reference
            ref_occurrence[indexFromQGram(refGram, qGramLength)] = 1;
            ref_one_bits++;
        }
    }

    // FROM GRIM-MASTER
    unsigned threshold = ReadLength - (ErrorThreshold*qGramLength); //ReadLength - 4 - (ErrorThreshold*qGramLength)
    // run the filter on that bin against the entire read sequence.
    unsigned count = 0;
    if (ErrorThreshold == 0) {
        count = 1;
    }
    int i;
    for (i = 0; i < ReadLength - qGramLength + 1; i++) {
        memcpy(readGram, &ReadSeq[i], qGramLength);

        int tmp_BV_Val = ref_occurrence[indexFromQGram(readGram, qGramLength)];

        if (ErrorThreshold == 0) {
            count &= (1 & tmp_BV_Val);
            if (count == 0) {
                break;
            }
        }
        else {
            count += (1 & tmp_BV_Val);
            if (count >= threshold) {
                break;
            }
        }
    }

    // we have to inverse output to use with our suite
    if (ErrorThreshold == 0) {
        free(ref_occurrence);
        return count/(qGramLength);
    }
    else {
        free(ref_occurrence);
        return count/(qGramLength);
    }
}

int grim_original_tweak(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength) {
    // assumes ReadLength == RefLength

    int occurrenceTableLength = getOccurrenceTableLength(qGramLength);
    char * ref_occurrence = (char *) calloc(occurrenceTableLength, sizeof(char)); // stores only binary values

    char refGram[qGramLength];
    char readGram[qGramLength];
    int qGramCount = ReadLength - qGramLength + 1;
    int ref_one_bits = 0;

    for (int i = 0; i < qGramCount; i++) {
        memcpy(refGram, &RefSeq[i], qGramLength);
        if (ref_occurrence[indexFromQGram(refGram, qGramLength)] == 0) { //read occurs in reference
            ref_occurrence[indexFromQGram(refGram, qGramLength)] = 1;
            ref_one_bits++;
        }
    }

    // FROM GRIM-MASTER
    unsigned threshold = ReadLength - (qGramLength-1) - (ErrorThreshold*qGramLength); //ReadLength - 4 - (ErrorThreshold*qGramLength)
    // run the filter on that bin against the entire read sequence.
    unsigned count = 0;
    if (ErrorThreshold == 0) {
        count = 1;
    }
    int i;
    for (i = 0; i < ReadLength - qGramLength + 1; i++) {
        memcpy(readGram, &ReadSeq[i], qGramLength);

        int tmp_BV_Val = ref_occurrence[indexFromQGram(readGram, qGramLength)];

        if (ErrorThreshold == 0) {
            count &= (1 & tmp_BV_Val);
            if (count == 0) {
                break;
            }
        }
        else {
            count += (1 & tmp_BV_Val);
            if (count >= threshold) {
                break;
            }
        }
    }

    // we have to inverse output to use with our suite
    if (ErrorThreshold == 0) {
        free(ref_occurrence);
        return count/(qGramLength);
    }
    else {
        free(ref_occurrence);
        return count/(qGramLength);
    }
}