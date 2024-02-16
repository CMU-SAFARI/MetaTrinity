#include <stdio.h>
#include <stdlib.h>

#define IDX(I, J) ((I)*m + (J))
#define min(A, B) (((A)<(B))?(A):(B))

int repetitive_edit_distance(char *repetitive_pattern, int repetitive_pattern_length, char *target, int target_length){
    int n = target_length;
    int m = repetitive_pattern_length;

    int *dp = (int *) malloc(sizeof(int) * m * (n+1));

    for(int j = 0; j < m; j++){
        dp[IDX(0, j)] = j;
    }

    for(int i = 1; i < n+1; i++){
        for(int j = 1; j < m; j++){
            int from_topleft = dp[IDX(i-1, j-1)] + (repetitive_pattern[j-1] != target[i-1]);
            int from_left = dp[IDX(i-1, j)] + 1;
            dp[IDX(i, j)] = min(from_topleft, from_left);
        }

        int from_topleft = dp[IDX(i-1, m-1)] + (repetitive_pattern[m-1] != target[i-1]);
        int from_left = dp[IDX(i-1, 0)] + 1;
        dp[IDX(i, 0)] = min(from_topleft, from_left);

        for(int j = 1; j < 2*m-1; j++){
            int from_top = dp[IDX(i, (j-1)%m)] + 1;
            dp[IDX(i, j%m)] = min(dp[IDX(i, j%m)], from_top);
        }
    }

    int best = dp[IDX(n, 0)];
    for(int j = 1; j < m; j++){
        best = min(best, dp[IDX(n, j)]);
    }

    free(dp);

    return best;
}

int triangle_inequality(char *query, int query_length, char *target, int target_length, char *profile_pattern, int profile_pattern_length){
    int query_distance = repetitive_edit_distance(profile_pattern, profile_pattern_length, query, query_length);
    int target_distance = repetitive_edit_distance(profile_pattern, profile_pattern_length, target, target_length);
    int edit_distance_lower_bound = abs(query_distance - target_distance);
    return edit_distance_lower_bound;
}

int triangle_inequality_filter(char *query, int query_length, char *target, int target_length, char *profile_pattern, int profile_pattern_length, int edit_distance_threshold){
    int edit_distance_lower_bound = triangle_inequality(query, query_length, target, target_length, profile_pattern, profile_pattern_length);
    int reject = edit_distance_lower_bound > edit_distance_threshold;
    return !reject;
}

void rand_str(char *dest, size_t length) {
    char charset[] = "ACGT";

    while (length-- > 0) {
        size_t index = (double) rand() / RAND_MAX * (sizeof charset - 1);
        *dest++ = charset[index];
    }
    *dest = '\0';
}

void ith_str(char *dest, int i, size_t length) {
    //lexicographically Ith string of length LENGTH

    dest[length] = '\0';

    //no chars get set for the empty string
    if(length == 0) return;

    const char charset[] = "ACGT";
    //recursively build the first length-1 characters
    ith_str(dest, i/4, length-1);
    //add the last character
    dest[length-1] = charset[i%4];
}

int triangle_inequality_filter_random(char *query, int query_length, char *target, int target_length, int profile_pattern_length, int num_profile_patterns, int edit_distance_threshold){
    //randomly generate NUM_PROFILE_PATTERNS patterns of the given length and accept only if all of them accept

    srand(42);
    
    for(int profile_pattern_i = 0; profile_pattern_i < num_profile_patterns; profile_pattern_i++){
        char profile_pattern[profile_pattern_length + 1];
        rand_str(profile_pattern, profile_pattern_length);
        int accepted = triangle_inequality_filter(query, query_length, target, target_length, profile_pattern, profile_pattern_length, edit_distance_threshold);
        if(!accepted){
            return 0;
        }
    }
    return 1;
}

int triangle_inequality_filter_all(char *query, int query_length, char *target, int target_length, int profile_pattern_length, int edit_distance_threshold) {
    //generate all patterns of the given length and accept only if all of them accept
    int num_profile_patterns = 1<<(2*profile_pattern_length);
    for(int profile_pattern_i = 0; profile_pattern_i < num_profile_patterns; profile_pattern_i++){
        char profile_pattern[profile_pattern_length + 1];
        ith_str(profile_pattern, profile_pattern_i, profile_pattern_length);
        int accepted = triangle_inequality_filter(query, query_length, target, target_length, profile_pattern, profile_pattern_length, edit_distance_threshold);
        if(!accepted){
            return 0;
        }
    }
    return 1;
}

//int main(int argc, char ** argv){
//    char *query = "AACCGGTTAACCGGTT"; int m = 16;
//    char *target = "AAAACCCCGGGGTTTT"; int n = 16;
//    char *profile = "ACGT"; int proflen = 4;
//    //int distance = repetitive_alignment("ACG", "ACGTACGTACGTACGT", 3, 16);
//    //printf("%d\n", distance);
//    //int edit_distance_lower_bound = triangle_inequality(query, m, target, n, profile, proflen);
//    int accepted = triangle_inequality_filter_all(query, m, target, n, 4, 3);
//    //int accepted = triangle_inequality_filter_random(query, m, target, n, 4, 4, 3);
//    //printf("%d\n", edit_distance_lower_bound);
//    printf("accepted=%s\n", accepted ? "true" : "false");
//}
