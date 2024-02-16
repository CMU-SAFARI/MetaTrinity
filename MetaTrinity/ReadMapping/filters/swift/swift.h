//
// Created by mdr on 22/06/2021.
//

#ifndef FILTER_SWIFT_H
#define FILTER_SWIFT_H

#endif //FILTER_SWIFT_H


typedef struct Bin_swift {
    int max;
    int min;
    int count;
} Bin;

typedef struct Parallelogram {
    int left;
    int top;
    int bottom;
} Parallelogram;

//void buildRefQGramIndex(char * ref, int refLength, int qGramLength, unsigned int * lookupTable, unsigned short * occ);
int swift(char * read, int readLength, char * ref, int refLength, int errThreshold, int qGramLength);
//int indexFromQGram(char * gram, int qGramLength);
//int exponentPowerTwo(int thresh);