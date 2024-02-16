//
// Created by mdr on 06/07/2021.
//

#ifndef FILTER_QGRAM_H
#define FILTER_QGRAM_H

#endif //FILTER_QGRAM_H


int qgram(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength);
int qgram_hash(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength);
