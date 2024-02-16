#include <pthread.h>
#include <stdio.h>
#include <dirent.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <getopt.h>
#include <string.h>

#include <map>
#include <thread>
#include <chrono>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "main.h"

int str_ends_with(const char *str, const char *suffix) {
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}
bool dirfilter(const struct dirent *const entry) {
    return (strcmp(entry->d_name,".") == 0 || strcmp(entry->d_name,"..") == 0 ||
        *entry->d_name == '.' || !str_ends_with(entry->d_name, ".mmi") ||
        entry->d_type == 4); //filter directories
}
static void usage(char* myname) {
    fprintf(stderr, "Usage: %s [-n <int: minimizer-cutoff>] [-b] [-t <int: nuimber of subthreads in minimap>] <mmidir> <readsfile> <translationfile> <outfile>\n", myname);
    exit(1);
}

int main(int argc, char *argv[]) {
    int opt;
    char *n = "3";
    char *t = "1";
    char *one = "1"; //fallback for sequential flag
    bool sequential = false;

    while ((opt = getopt(argc, argv, "bn:t:")) != -1) {
        switch (opt) {
        case 'b':
            sequential = true;
            break;
        case 'n': {
            char *endptr;
            const int base = 10;
            strtol(optarg, &endptr, base);
            if (*endptr != '\0' || endptr == optarg) {
                fprintf(stderr, "-n parameter has to be a number, got: %s\n", optarg);
                return 1;
            }
            n = optarg;
            break;
        }
        case 't': {
            char *endptr;
            const int base = 10;
            strtol(optarg, &endptr, base);
            if (*endptr != '\0' || endptr == optarg) {
                fprintf(stderr, "-t parameter has to be a number, got: %s\n", optarg);
                return 1;
            }
            t = optarg;
            break;
        }
        default: /* '?' */
            usage(argv[0]);
        }
    }
    if(argc != optind + 4)
        usage(argv[0]);
    const char *mmidir = realpath(argv[optind++], NULL);
    char *reads = realpath(argv[optind++], NULL);
    const char *tr_sorted_path = realpath(argv[optind++], NULL);
    FILE *test_exist = fopen(argv[optind], "w");
    fclose(test_exist);
    const char *out_path = realpath(argv[optind], NULL);
    if(!mmidir || !reads || !tr_sorted_path || !out_path) {
        fprintf(stderr, "ERROR: Resolved inputs like this:\n");
        fprintf(stderr, "Found -n %s?, -b set %d?, -t %s?, mmi %s?, reads %s?, tr %s?, out %s?\n",
            n, sequential, t, mmidir, reads, tr_sorted_path, out_path);
        return 1;
    }
    fprintf(stderr, "Found -n %s?, -b set %d?, -t %s?, mmi %s?, reads %s?, tr %s?, out %s?\n",
            n, sequential, t, mmidir, reads, tr_sorted_path, out_path);
    const char *out_format = "taxid_%s_genomic.fna.gz,%f\n";
    const int mmidirlen = strlen(mmidir);

    if(sequential) {
        fprintf(stderr, "Sequential flag set (-b), overriding -t to 1\n");
        t = one;
    }

    DIR *inDir = opendir(mmidir);
    struct dirent *entry;
    int fCnt = 0;
    
    while((entry = readdir(inDir)) != NULL) {
        fCnt += !dirfilter(entry);
    }
    fprintf(stderr, "fCnt is %d\n", fCnt);
    inDir = opendir(mmidir);

    pthread_t threads[fCnt];
    ao_queue *table_queues[fCnt];
    
    int i = 0;
    while((entry = readdir(inDir)) != NULL) {
        if(dirfilter(entry))
            continue;
        const char* currFile = entry->d_name;
        //use this instead of array (inp[size]) since array on stack may (will) be overwritten by other loop iterations
        int inpsize = mmidirlen + strlen(currFile) + (mmidir[mmidirlen-1] != '/') + 1;
        char *inp = (char*) malloc(sizeof(*inp) * inpsize);
        snprintf(inp, inpsize, mmidir[mmidirlen-1] == '/' ? "%s%.0s%s" : "%s%s%s", mmidir, "/", currFile);
        printf("got file %s\n", inp);

        char *minimap_argv[] = { "./minimap2", "-n", n, "-t", t, &inp[0], reads, NULL };
        int minimap_argc = (int)(sizeof(minimap_argv) / sizeof(minimap_argv[0])) - 1;
        int argv_size = sizeof(char*)*(&(minimap_argv[minimap_argc])-&(minimap_argv[0]));
        char **minimap_argv_heap = (char**) malloc(argv_size);
        memcpy(minimap_argv_heap, minimap_argv, argv_size);

        startargs *arguments = (startargs*) malloc(sizeof(*arguments));
        arguments->argc = minimap_argc;
        arguments->argv = minimap_argv_heap;

        const int mapcnt = atoi(t);
        table_queues[i] = new ao_queue();
        for(int j = 0; j < mapcnt; ++j){
            table_queues[i]->maps.push_back(new specific_map());
            table_queues[i]->maps.back()->map = new std::unordered_map<std::string, unsigned long long>();
        }
        table_queues[i]->mapslen = mapcnt;
        arguments->out_queue = &(table_queues[i]);

        fprintf(stderr, "Creating thread %d\n",i);

        if(!sequential)
            pthread_create(&threads[i++], NULL, start_wrapper, (void *)arguments);
        else
            start_wrapper(arguments);
    }
    if(!sequential) {
        for(int i = 0; i < fCnt; i++) {
            if(pthread_join(threads[i], NULL) != 0)
                fprintf(stderr, "Error on join thread %d\n",i);
            else
                fprintf(stderr, "Successfully joined thread %d\n",i);
        }
    }
    fprintf(stderr, "merging %d queues (REACHED)\n", fCnt);

    std::map<std::string, unsigned long long> merged_map = {};
    for(int tbl = 0; tbl < fCnt; ++tbl) {
        for(auto mapStruct : table_queues[tbl]->maps) {
            for (auto& it: *mapStruct->map) {
                merged_map[it.first] += it.second;
            }
        }
    }
    fprintf(stderr, "len is %ld\n", merged_map.size());
    unsigned long long max_hits = 0;
    std::ifstream trsorted_file(tr_sorted_path);
    std::string line;
    std::unordered_map<std::string, unsigned long long> outmap;
    std::map<std::string, unsigned long long>::iterator it = merged_map.begin();
    while(getline(trsorted_file, line) && it != merged_map.end()) {
        std::string key;
        std::stringstream ls(line);
        getline(ls, key, ' ');
        const int cmp = it->first.compare(key);
        if(cmp < 0 && it != merged_map.begin())
            fprintf(stderr, "ERROR: SORTING IS WRONG, seen %s < %s\n", it->first, key);
        if(cmp != 0)
            continue;
        std::string newKey;
        getline(ls, newKey, ' ');
        outmap[newKey] += it->second;
        if(outmap[newKey] > max_hits)
            max_hits = outmap[newKey];
        ++it;
    }
    std::ofstream outfile(out_path);
    for(const auto& item : outmap) {
        std::string formattedKey(item.first);
        std::replace(formattedKey.begin(), formattedKey.end(), '.', '_');
        printf(out_format, formattedKey.c_str(), (double)item.second/max_hits);
        outfile << "taxid_" << formattedKey << "_genomic.fna.gz," << (double)item.second/max_hits << "\n";
    }
    printf("%s\n", out_path);
    return 0;
}