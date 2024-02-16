

#ifndef METALIGN_MAIN_FILTER_H
#define METALIGN_MAIN_FILTER_H


#include "filters/SneakySnake/SneakySnake.h"
#include "filters/adjacency-filter/AdjacencyFilter.h"
#include "filters/base-counting/Base_Counting.h"
#include "filters/magnet/MAGNET_DC.h"
#include "filters/hamming-distance/HD.h"
#include "filters/shouji/Shouji.h"
#include "filters/qgram/qgram.h"

#include "filters/shd/SHD.h"
//#include "filters/swift/swift.h" //==> ERROR: filters/swift/swift.h:11:16: error: redefinition of ‘struct Bin_swift’

extern char* filter;

extern int filter_calls;


extern uint64_t (*seed_map)[2];

// extern int (*seed_map)[2];

#endif //METALIGN_MAIN_FILTER_H


