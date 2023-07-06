#ifndef TAXI_H
#define TAXI_H

#include <quad.h>
#include <stdint.h>

/*
typedef struct {
//	quad csum;
	uint32_t a, b;
	char count;		// if another pair has the same sum, upon insertion this count will be incremented instead of adding a new item.
} ht_entry;
*/

typedef uint64_t ht_entry;

typedef struct resentry {
	quad csum;
	struct resentry *next;
	uint32_t a;
	uint32_t b;
} result_entry;

typedef struct {
	result_entry *first;
} result_list;

typedef struct {
	uint32_t thread_number;

	// intervals
	quad start, end;	// total work to do
	quad istart, iend;	// bounds of current subinterval

	// line bounds
	uint32_t *lb;	// low and high bounds for lines.
	uint32_t *hb;
	uint32_t linemax;	// highest line # used this subinterval
	uint32_t max_ll;
	char initialized;	// whether we've done the first linebound computation or not

	long istart_time;

	// hashtable
	ht_entry *hashtable;
	char *active;
	uint32_t ht_load;
	uint32_t probes;
	double target_fudge_factor;

	// results list
	result_list rlist;
	uint32_t rlist_size;
	uint32_t nresults;

} taxi_state;

double pi(double x);
void sieve_block (int *frogs, char *block, uint32_t bound, int nsmallprimes) ;
void ht_insert (taxi_state*, uint32_t, uint32_t);
result_entry *new_resentry (quad, uint32_t, uint32_t);
void add_result (taxi_state *, quad, uint32_t, uint32_t);
uint32_t get_line_idx (quad, uint32_t, uint32_t, uint32_t, uint32_t);
void clear_hashtable (taxi_state *);

taxi_state *init_taxi_state (quad start, quad end, int thread_num) ;

#endif
