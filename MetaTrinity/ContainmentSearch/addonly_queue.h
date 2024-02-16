#ifndef ADDONLY_QUEUE_H
#define ADDONLY_QUEUE_H

#include <vector>
#include <mutex>
#include <unordered_map>

extern int cntmapsizethreshold;
typedef struct {
    std::unordered_map<std::string, unsigned long long>* map;
    std::mutex mutex;
} specific_map;
typedef struct {
    std::vector<specific_map*> maps; //TODO make array?
    int mapslen;
} ao_queue;

#endif