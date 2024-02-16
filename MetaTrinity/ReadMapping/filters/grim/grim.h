//
// Created by mdr on 10.01.22.
//

#ifndef FILTER_GRIM_H
#define FILTER_GRIM_H

int grim(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength);
int grim_original(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength);
int grim_original_tweak(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength);
int grim_long(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold, int qGramLength);

#endif //FILTER_GRIM_H
