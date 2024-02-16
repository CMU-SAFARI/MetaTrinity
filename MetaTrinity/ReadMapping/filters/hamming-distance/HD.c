//
// Created by mdr on 03.01.22.
//

#include "HD.h"

int HD(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int DebugMode)
{
    int count = 0;
    for (int i = 0; i < ReadLength; i++) {
        if (RefSeq[i] != ReadSeq[i]) {
            if (++count > ErrorThreshold) {
                break;
            }
        }
    }
    return count;
}
