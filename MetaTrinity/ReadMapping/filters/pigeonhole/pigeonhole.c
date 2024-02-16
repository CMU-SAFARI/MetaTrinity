//
// Created by mdr on 18.01.22.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pigeonhole.h"
#include "../../khash.h"


#define kh_get_val(kname, hash, key, defVal) ({k=kh_get(kname, hash, key);(k!=kh_end(hash)?kh_val(hash,k):defVal);})
#define kh_set(kname, hash, key, val) ({int ret; k = kh_put(kname, hash,key,&ret); kh_value(hash,k) = val; ret;})

const int khStrInt = 33;
KHASH_MAP_INIT_STR(khStrInt, int)

int pigeonhole(int ReadLength, const char RefSeq[], const char ReadSeq[], int ErrorThreshold) {
    khiter_t k;
    khash_t(khStrInt) *h =  kh_init(khStrInt);

    //this might leave an odd tail at the end

    int n_buckets = ErrorThreshold + 1;
    int bucket_length = ReadLength / n_buckets;

//    printf("%i\n", bucket_length);

    char * subsequence = (char *) malloc(bucket_length*sizeof(char));

    //store reference to hash table
    for (int i = 0; i < ReadLength; i += bucket_length) { //i += bucket_length
        memcpy(subsequence, &RefSeq[i], bucket_length);
        kh_set(khStrInt, h, subsequence, 1);
    }

    for (int i = 0; i < ReadLength; i += bucket_length) { //i += bucket_length
        memcpy(subsequence, &ReadSeq[i], bucket_length);
        k = kh_get(khStrInt, h, subsequence);
        if (k == kh_end(h)) {
            return 1;
        }
    }
    return 0;
}