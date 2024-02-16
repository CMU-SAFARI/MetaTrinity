int repetitive_edit_distance(char *repetitive_pattern, int repetitive_pattern_length, char *target, int target_length);
int triangle_inequality_filter_random(char *query, int query_length, char *target, int target_length, int profile_pattern_length, int num_profile_patterns, int edit_distance_threshold);
int triangle_inequality_filter_all(char *query, int query_length, char *target, int target_length, int profile_pattern_length, int edit_distance_threshold);

