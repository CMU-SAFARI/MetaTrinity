#ifndef MAIN_H
#define MAIN_H

#include "addonly_queue.h"

typedef struct arg_struct {
    int argc;
    char **argv;
    ao_queue **out_queue;
} startargs;

int start(int argc, char *argv[], ao_queue **out_queue);
void* start_wrapper(void* args);

#endif
