//
// Created by Max Rumpf on 14.05.20.
//

#ifndef FILTER_BASE_COUNTING_TEST_H
#define FILTER_BASE_COUNTING_TEST_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


/* Function Declarations */
extern int baseCounting(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode);
extern int baseCountingTest(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode);
//extern int baseCountingTest2(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode);
//extern int baseCountingTest3(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode);

#endif //FILTER_BASE_COUNTING_TEST_H
