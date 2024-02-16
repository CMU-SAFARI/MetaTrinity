//
// Created by Max Rumpf, 2020
//

#include <stdlib.h>
#include "Base_Counting.h"

int baseCountingTestArgs(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode, int argc, const char * const argv[]) {
    return baseCountingTest(ReadLength, RefSeq, ReadSeq, ErrorThreshold, DebugMode);
}

int baseCountingTest(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode) {
    return baseCounting(ReadLength, RefSeq, ReadSeq, ErrorThreshold, DebugMode);
}

int baseCounting(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode) {

    int aCount = 0;
    int tCount = 0;
    int gCount = 0;
    int cCount = 0;

    // Count bases in Reference Sequence

    for (int i = 0; i < ReadLength; i++) {
        switch (RefSeq[i]) {
            case 'A': aCount++; break;
            case 'T': tCount++; break;
            case 'G': gCount++; break;
            case 'C': cCount++; break;
        }
        switch (ReadSeq[i]) {
            case 'A': aCount--; break;
            case 'T': tCount--; break;
            case 'G': gCount--; break;
            case 'C': cCount--; break;
        }
    }

    //return (abs(aCount)+abs(tCount)+abs(gCount)+abs(cCount) <= 2*ErrorThreshold);
     return (abs(aCount)+abs(tCount)+abs(gCount)+abs(cCount)/2);// still ok??
}

//int baseCountingTest2(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode) {
//
//    int aCount;
//    int tCount;
//    int gCount;
//    int cCount;
//
//    // Count bases in Reference Sequence
//
//    for (int i = 0; i < ReadLength; i++) {
//        switch (RefSeq[i]) {
//            case 'A': aCount++; break;
//            case 'T': tCount++; break;
//            case 'G': gCount++; break;
//            case 'C': cCount++; break;
//        }
//    }
//
//    for (int i = 0; i < ReadLength; i++) {
//        switch (RefSeq[i]) {
//            case 'A':
//                aCount--;
//                break;
//            case 'T':
//                tCount--;
//                break;
//            case 'G':
//                gCount--;
//                break;
//            case 'C':
//                cCount--;
//                break;
//        }
//    }
//
//    return (abs(aCount)+abs(tCount)+abs(gCount)+abs(cCount) <= 2*ErrorThreshold);
//}
//
//int baseCountingTest3(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode) {
//
//    unsigned int aCount;
//    unsigned int tCount;
//    unsigned int gCount;
//    unsigned int cCount;
//
//    unsigned int aReadCount;
//    unsigned int tReadCount;
//    unsigned int gReadCount;
//    unsigned int cReadCount;
//
//    // Count bases in Reference Sequence
//
//    for (int i = 0; i < ReadLength; i++) {
//        switch (RefSeq[i]) {
//            case 'A': aCount++; break;
//            case 'T': tCount++; break;
//            case 'G': gCount++; break;
//            case 'C': cCount++; break;
//        }
//        switch (ReadSeq[i]) {
//            case 'A': aReadCount--; break;
//            case 'T': tReadCount--; break;
//            case 'G': gReadCount--; break;
//            case 'C': cReadCount--; break;
//        }
//    }
//    return (abs((int) (aCount-aReadCount))+abs((int) (tCount-tReadCount))+abs((int) (gCount-gReadCount))+abs((int) (cCount-cReadCount)) <= 2*ErrorThreshold);
//}
