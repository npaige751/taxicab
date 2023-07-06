#include "quad.h"
#include "taxi2.h"
#include <gmp.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

/* Todo:
 *
 + Redo linebounds to be vertical rather than horizontal (1/5th as many linebounds!)
 + Potentially use linear search for some linebound computations? Though this makes much more sense with horizontal lines.
 + Reduce sieve memory usage (packed sieve array, multiples of 2 not included; block sieve)
 *    - potentially store primes as gap lengths (can do 1 byte gaps since gaps will be even) + an index; perhaps do some caching and more intelligent access patterns
 * Is it possible/useful to do some cubsum computations incrementally as we add things to the table? 
 * Fix "impending" overflow problem in quad_cubsum
 * code cleanup, commenting
 * Client/server stuff? Or let ruby handle that?
 *
*/

/* What will break when?
 *
 * 2.8e28	intermediate results in quad_cubsum overflow
 * 8e28		primes no longer fit in 32 bits
 * 1e32		prime array reaches about 8GB (if it could even be in 32-bit ints anymore)
 * 1e33		prime indices no longer fit in 32 bits
 * 1e38		128-bit overflow
*/

/* Networking stuff:
 *
 * Server(s) keep track of which work units are done, allocated, and free.
 * When a client starts, it contacts the server and requests work units.
 * When a work unit finishes, client sends results back; implicitly requests a new work unit. Server responds with an acknowledgement & a new work unit.
 * Server can ask a client for a status report; client will respond with which work units it's working on, % complete, ETAs
 *
 * Work units have a time limit; if not completed within the given time interval, the server can reassign them. 
 *
 * Should there be a tcp connection between each client thread & the server, or just one per host?
 * Should the conections stay open the whole time, or should they be once-per-message? probably not persistent
 * If not persistent, can have one per thread without too much difficulty.
 *
 *
*/

uint32_t ht_bins = 80000000;
uint32_t target_load = 27000000;

uint32_t *primes;
uint32_t pbound;
uint32_t nprimes;

int nthreads = 1;
int verbose = 1;

long start_time = 0;

void sieve (int pbound) {
	printf ("Starting sieve; will allocate %d MB sieve array plus about %d MB for primes\n", pbound/1000000, (int) (pi(pbound)/250000));
	long time = clock();
	char *b = (char *) (malloc (pbound));
	b[0] = 0;
	b[1] = 0;
	memset (b+2, 1, pbound-2);
	int p = 2;
	int ct = 0;
	while (p < pbound) {
		if (b[p]) {	// if the current entry is prime, sieve with it
			ct ++;
			int pos = 2 * p;
			while (pos < pbound) {
				b[pos] = 0;
				pos += p;
			}
		}
		p++;
	}
	b[2] = 0;	// we don't want to use 2 in our searches.
	ct --;
	// extract
	primes = (uint32_t *)(malloc (4 * ct));
	int w = 0;
	for (int i=0; i < pbound; i++) {
		if (b[i]) {
			primes[w] = i;
			w++;
		}
	}
	free (b);
	nprimes = ct;
	if (verbose) printf ("Sieve (standard) took %d ms\n", (clock()-time)/1000);
}

#define BLOCKLEN (1 << 18)	// 256kb... should fit comfortably in L2 cache

int sieve2 () {
	long time = clock();
	if (pbound < 6*(BLOCKLEN*16)) {
		printf ("pbound = %d < 6*blockspan = %d; using regular sieve\n", pbound, 6*16*BLOCKLEN);
		sieve (pbound);
		return nprimes;
	}
	uint32_t bound = (int) (sqrt(pbound) + 1);
	printf ("Starting block sieve with small primes up to %d; BLOCKLEN = %d\n", bound, BLOCKLEN);
	sieve (bound);
	int nsmallprimes = nprimes;
	printf ("pbound = %d; smallprime bound = %d; found %d small primes; will do about %d blocks\n", pbound, bound, nprimes, (pbound-bound)/(16*BLOCKLEN));

	char *block = malloc (BLOCKLEN);
	int *frogs = malloc (nprimes * sizeof(int));
	int pr_alloc = (int) (pi(pbound) * 1.15);
	printf ("Allocating %d entries (%d MB) for primes\n", pr_alloc, pr_alloc/250000);
	primes = realloc (primes, pr_alloc*4);

	if (bound % 2 == 0) bound++;
	// initialize frogs
	for (int i=0; i < nprimes; i++) {
		int p = primes[i];
		frogs[i] = ((p - (bound % p)) * (p/2 + 1)) % p;	// solution to   bound + 2*frog = 0 (mod p)
	}

	while (bound < pbound) {
		sieve_block (frogs, block, bound, nsmallprimes);
		bound += BLOCKLEN * 16;
		//printf ("Done upto %d; have %d primes\n", bound, nprimes);
	}
	printf ("Finished sieve in %d ms\n", (clock()-time) / 1000);
	free(frogs);	// be free, little froggies!
	free(block);
}

void sieve_block (int *frogs, char *block, uint32_t bound, int nsmallprimes) {
	memset (block, 0, BLOCKLEN);
	/*
	// this is very slightly faster, perhaps due to smaller loop overhead or better locality
	for (int i=0; i < nsmallprimes-1; i+=2) {
		int B1 = frogs[i] / 8;
		int b1 = frogs[i] % 8;
		int B2 = frogs[i+1] / 8;
		int b2 = frogs[i+1] % 8;
		while (frogs[i] < BLOCKLEN*8 || frogs[i+1] < BLOCKLEN*8) {
			if (frogs[i] < frogs[i+1]) {
				block[B1] |= (1 << b1);
				frogs[i] += primes[i];
				B1 = frogs[i] / 8;
				b1 = frogs[i] % 8;
			} else {
				block[B2] |= (1 << b2);
				frogs[i+1] += primes[i+1];
				B2 = frogs[i+1] / 8;
				b2 = frogs[i+1] % 8;
			}
		}
		frogs[i] -= BLOCKLEN*8;	// reset for next iteration
		frogs[i+1] -= BLOCKLEN*8;	// reset for next iteration
	}
	*/
	for (int i=0; i < nsmallprimes; i++) {
		int p = primes[i];
		while (frogs[i] < BLOCKLEN*8) {
			int byte = frogs[i] / 8;
			int bit = frogs[i] % 8; 
			block[byte] |= (1 << bit);
			frogs[i] += p;	// hop forward
		}
		frogs[i] -= BLOCKLEN*8;	// reset for next iteration
	}
	// now extract
	int nonblank = 0;
	for (int i=0; i < BLOCKLEN; i++) {
		if (~block[i]) {
			nonblank ++;
			for (int j=0; j < 8; j++) {
				if ((~block[i]) & (1 << j)) {
					primes[nprimes++] = bound + (i*8+j)*2;
				}
			}
		}
	}
}

double pi (double x) {
	return x / log(x);
}

double prime_density (double x) {
	return (pi(x*1.01) - pi (x))/(0.01*x);
}
	
void get_iend (taxi_state *ts){	// computes an upper bound for a subinterval that will contain approximately target_load items and stores it in iend
	double prdensity = prime_density (quad_to_double (ts->istart));
	uint32_t ht = get_line_idx (quad_shr (ts->istart, 1), 0, 0, nprimes, nprimes/2);
	double target_linelen = ts->target_fudge_factor * (target_load / prdensity) / ht;
	double iendf = pow (cbrt (quad_to_double (ts->istart)) + target_linelen, 3);
	ts->iend = double_to_quad (iendf);

}

/* The last three parameters may be used to restrict the search and provide tighter bounds when information is known.
* They should default to 0 and nprimes */

/*Finds the largest index x such that pr[x]**3 + pr[line]**3 <= target */
uint32_t get_line_idx (quad target, uint32_t line, uint32_t lostart, uint32_t histart, uint32_t guessstart) {
	uint32_t lo, hi, guess;
	lo = lostart;
	hi = histart;
	guess = guessstart;
	while (hi-lo > 1) {
		quad v = quad_cubsum (primes[guess], primes[line]);
		int cmp = quad_cmp (v, target);
		if (cmp > 0) {	// guess is too high
			hi = guess;
			guess = (lo + hi) / 2;
		} else if (cmp < 0) {
			lo = guess;
			guess = (lo + hi) / 2;
		} else {
			return guess;
		}
	}
	return guess;
}

// searches downward from guess
uint32_t get_line_idx_linear (quad target, uint32_t line, uint32_t guess) {
	while (quad_cmp(quad_cubsum (primes[guess], primes[line]), target) > 0) {
		guess--;
	}
	return guess;
}

#define USE_LINEAR_SEARCH 0

void initialize (quad end) {
	// largest prime needed will be about cbrt(end)
	mpz_t c;
	mpz_init (c);
	mpz_set_quad (c, end);
	mpz_root (c, c, 3);
	pbound = (uint32_t) mpz_get_ui (c) + 1000;	// safety
//	sieve (pbound);
	sieve2 ();
}

taxi_state *init_taxi_state (quad start, quad end, int thread_num) {
	taxi_state *ts = calloc (1, sizeof(taxi_state));
	ts->target_fudge_factor = 0.2;

	ts->thread_number = thread_num;
	ts->start = start;
	ts->end = end;
	ts->istart = start;

	ts->hashtable = (ht_entry *) (calloc (ht_bins, sizeof(ht_entry)));
	ts->active = (char *) (calloc (ht_bins, sizeof(char)));

	return ts;
}

void clear_hashtable (taxi_state *ts) {
	memset (ts->active, 0, ht_bins * sizeof (char));
}

int compare_htentries (taxi_state *ts, quad q, ht_entry e) {
	return quad_eq (q, quad_cubsum (primes[e & 0xfffffffful], primes[e >> 32]));
}

void ht_insert (taxi_state *ts, uint32_t a, uint32_t b) {
	quad csum = quad_cubsum (primes[a], primes[b]);
	uint32_t hash = ((csum.lo/2) ^ (csum.hi*47918 - csum.lo*98721)) % ht_bins;
	ts->probes++;
	while (ts->active[hash] && !compare_htentries (ts, csum, ts->hashtable[hash])) {
		ts->probes++;
		hash++; //= (hash + 1) % ht_bins;
		if (hash == ht_bins) hash = 0;
	}
	if (ts->active[hash]) {	// we've found a 2-way (at least)
		ts->nresults ++;
		add_result (ts, csum, 0,0);
	} else {
		ts->active[hash] = 1;
		ts->hashtable[hash] = (((uint64_t) a) << 32) | (uint64_t) b;
//		ts->hashtable[hash] = (a * ts->max_ll + (b - ts->lb[a]));
	}
}

result_entry *new_resentry (quad csum, uint32_t a, uint32_t b) {
	result_entry *res = (result_entry *) (malloc (sizeof (result_entry)));
	res->csum = csum;
	res->next = NULL;
	return res;
}

void add_result (taxi_state *ts, quad csum, uint32_t a, uint32_t b) {
	ts->rlist_size++;
	if (ts->rlist.first == NULL) {
		ts->rlist.first = new_resentry (csum, a, b);
		return;
	}
	result_entry *r = new_resentry (csum, a, b);
	r->next = ts->rlist.first;
	ts->rlist.first = r;
}

int compar_entries (const void *a, const void *b) {
	result_entry *ra = (result_entry *) a;
	result_entry *rb = (result_entry *) b;
	return quad_cmp (ra->csum, rb->csum);
}

#define STEEP_THRESH 6

void run_interval (taxi_state *ts) {
	ts->istart_time = clock();
	if (ts->initialized) ts->istart = ts->iend;
	get_iend(ts);
//	printf ("Running interval from %s to %s\n", quad_get_str (istart), quad_get_str (iend));
	if (quad_cmp (ts->iend, ts->end) > 0) {
		ts->iend = ts->end;
	}
	if (ts->initialized) {
	//	recompute_linebounds (ts);
	} else {
	//	initialize_linebounds (ts);
		ts->initialized = 1;
		start_time = clock();
	}
	clear_hashtable(ts);
	ts->probes = 0;
	ts->max_ll ++;
	ts->ht_load = 0;

	uint32_t lbc_start = get_line_idx (quad_shr (ts->istart, 1), 0, 0, nprimes, nprimes/2);	// lower bound curve start
	uint32_t lbc_end   = get_line_idx (ts->istart, 0, 0, nprimes, nprimes/2);
	uint32_t ubc_start = get_line_idx (quad_shr (ts->iend, 1), 0, 0, nprimes, nprimes/2);	// lower bound curve start
	uint32_t ubc_end   = get_line_idx (ts->iend, 0, 0, nprimes, nprimes/2);
	
	uint32_t prev_bstart = nprimes;
	uint32_t prev_bend;
	bool bstart_steep = false;
	bool bend_steep= false;
	for (uint32_t a = lbc_start; a < ubc_end; a++) {
		uint32_t bstart, bend;
		if (a <= lbc_end) {
			//bstart = get_line_idx (ts->istart, a, 0, nprimes, nprimes/2) + 1;
			if (a == lbc_start || bstart_steep) {
				bstart = get_line_idx (ts->istart, a, 0, prev_bstart, prev_bstart - 20) + 1;
			} else {
				bstart = get_line_idx_linear (ts->istart, a, prev_bstart);
			}
		} else {
			bstart = 0;
		}
		if (a <= ubc_start) {
			bend = a;
		} else {
//			bend = get_line_idx (ts->iend, a, 0, nprimes, nprimes/2);
			if (a == lbc_start || bend_steep) {
				bend = get_line_idx (ts->iend, a, 0, prev_bend, prev_bend - 20);
			} else {
				bend = get_line_idx_linear (ts->iend, a, prev_bend);
			}
		}
		if (prev_bstart - bstart > STEEP_THRESH) bstart_steep = true;
		if (prev_bend - bend > STEEP_THRESH) bend_steep = true;

		prev_bstart = bstart;
		prev_bend = bend;

		for (uint32_t b = bstart; b <= bend; b++) {
			ht_insert (ts, a,b);
		}
		if (bend >= bstart) {
			ts->ht_load += bend - bstart + 1;
		}
	}
	ts->target_fudge_factor *= (double) target_load / ts->ht_load;

	// print stats
	double ilen = quad_to_double(quad_sub (ts->iend, ts->istart));
	double prog_per_entry = ilen / ts->ht_load;
	double worklen = quad_to_double (quad_sub (ts->end, ts->iend));
	double donelen = quad_to_double (quad_sub (ts->iend, ts->start));
	double completion_percent = donelen / quad_to_double(quad_sub(ts->end, ts->start)) * 100;
	double eta = 0;
	if (donelen > 0) eta = (worklen/donelen) * (clock() - start_time) / 1000000.0;

	if (verbose) printf ("\r[%.2f%%] Bound: %s. ht_load: %d; ppe: %e; ETA: %ds; probes: %d; nlines: %d; e/µs %.3f; results: %d   ", completion_percent, quad_get_str (ts->iend), ts->ht_load, prog_per_entry, (int) eta, ts->probes, ubc_end - lbc_start, ts->ht_load / (double) (clock() - ts->istart_time), ts->nresults);
	putchar('\n');
	fflush (stdout);
}

void *run_worker (void *argp) {
	taxi_state *ts = (taxi_state *) argp;
	while (!quad_eq (ts->iend, ts->end)) {
		run_interval(ts);
	}

	result_entry res[ts->rlist_size];
	result_entry *rover = ts->rlist.first;
	for (int i=0; i < ts->rlist_size; i++) {
		res[i] = *rover;
		rover = rover->next;
	}

	qsort ((void *) res, ts->rlist_size, sizeof (result_entry), compar_entries);
	/* Print results */
	char filename[8];
	sprintf(filename, "res%d", ts->thread_number);
	FILE *out = fopen(filename, "w");
	for (int i=0; i < ts->rlist_size; i++){ 
		printf("%s\n", quad_get_str(res[i].csum));
		fprintf(out, "%s\n", quad_get_str(res[i].csum));
	}
	fclose(out);
}

void run (quad *starts, quad *ends, int nthreads) {
	// find the maximum end point
	quad maxend = ends[0];
	for (int i=1; i < nthreads; i++) {
		if (quad_cmp (maxend, ends[i]) < 0) {
			maxend = ends[i];
		}
	}

	initialize (maxend);	// sieve, etc.
	pthread_t threads[nthreads];
	for (int i=0; i < nthreads; i++) {
		pthread_create (&threads[i], NULL, run_worker, init_taxi_state (starts[i], ends[i], i+1));
	}
	for (int i=0; i < nthreads; i++) {
		pthread_join (threads[i], NULL);
	}
}
			
int main (int argc, const char *argv[]) {
	quad start;
	quad interval_len;
	int a = 1;
	if (!strcmp (argv[1], "-t")) {
		nthreads = atoi(argv[2]);
		a = 3;
	}
	quad starts[nthreads], ends[nthreads];
	start = sci_to_quad(argv[a]);
	interval_len = sci_to_quad(argv[a+1]);
	for (int i=0; i < nthreads; i++) {
		starts[i] = start;
		start = quad_add(start, interval_len);
		ends[i] = start;
	}
	run (starts, ends, nthreads);
}
