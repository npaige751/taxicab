#include <gmp.h>
#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

#define USE_SEEDLISTS

typedef __int128_t quad;

mpz_t primecheck_tmp, tmp2, tmp3;

void mpz_set_quad(mpz_t x, __int128_t y) {
	bool neg = y < 0;
	if (neg) {
		y = -y;
	}
	mpz_set_ui(x, (uint64_t) (y >> 64));
	mpz_mul_2exp(x, x, 64);
	mpz_add_ui(x, x, (uint64_t) y);
	if (neg) {
		mpz_neg(x,x);
	}
}

quad gcd(quad a, quad b) {
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	while(a > 0) {
		quad t = b % a;
		b = a;
		a = t;
	}
	return b;
}

int gcdi(int a, int b) {
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	while(a > 0) {
		int t = b % a;
		b = a;
		a = t;
	}
	return b;
}

uint32_t *fb;
int fbb;
size_t fblen;
void sieve (int pbound) {
	char *b = (char *) (malloc (pbound));
	b[0] = 0;
	b[1] = 0;
	memset (b+2, 1, pbound-2);
	int p = 2;
	int powct = 0;
	while (p < pbound) {
		if (b[p]) {	// if the current entry is prime, sieve with it
			powct ++;
			int pos = 2 * p;
			while (pos < pbound) {
				b[pos] = 0;
				pos += p;
			}
			int pp = p*p;
			while (pp < pbound) {
				powct++;
				pp *= p;
			}
		}
		p++;
	}
	// extract
	fblen = powct;
	fb = (uint32_t *)(malloc (4 * powct));
	int w = 0;
	for (int i=0; i < pbound; i++) {
		if (b[i]) {
			fb[w++] = i;
			int pp = i*i;
			while (pp < pbound) {
				fb[w++] = pp;
				pp *= i;
			}
		}
	}
	free (b);
}

void print_quad_lo(quad x) {
	printf("%lld\n", (uint64_t) x);
}

void print_quad(quad x, FILE *out) {
	if (x == 0) {
		fprintf(out, "0");
		return;
	}
	char s[60];
	bool neg = x < 0;
	if (neg) x = -x;
	int w = 0;
	while (x != 0) {
		int d = (int) (x % 10);
		char c = (char) ('0' + d);
		s[w++] = (char) ('0' + d);
		x /= 10;
	}
	if (neg) fputc('-', out);
	for (int i = w-1; i >= 0; i--) {
		fputc(s[i], out);
	}
}

void print_check_result(quad *x) {
	for (int i=0; i < 4; i++) {
		print_quad(x[i], stdout);
		if (i != 3) printf(", ");
	}
	printf("\n");
}

int nresults = 0;
int res_cap = 10;
quad *results;
void add_result(quad p, quad q) {
	if (nresults == res_cap) {
		res_cap *= 2;
		results = (quad *) realloc(results, res_cap * 2 * sizeof(quad));
	}
	results[2*nresults] = p;
	results[2*nresults + 1] = q;
	nresults++;
}

int res_cmp(const void *a, const void *b) {
	quad *aq = (quad *) a;
	quad *bq = (quad *) b;
	mpz_set_quad(tmp2, aq[0]);
	mpz_pow_ui(tmp2, tmp2, 3);
	mpz_set_quad(tmp3, aq[1]);
	mpz_pow_ui(tmp3, tmp3, 3);
	mpz_add(tmp2, tmp2, tmp3);

	mpz_set_quad(primecheck_tmp, bq[0]);
	mpz_pow_ui(primecheck_tmp, primecheck_tmp, 3);
	mpz_set_quad(tmp3, bq[1]);
	mpz_pow_ui(tmp3, tmp3, 3);
	mpz_add(tmp3, tmp3, primecheck_tmp);

	int c = mpz_cmp(tmp2, tmp3);
	if (c != 0) return c;
	if (aq[0] < bq[0]) return -1;
	if (aq[0] > bq[0]) return 1;
	if (aq[1] < bq[1]) return -1;
	if (aq[1] > bq[1]) return 1;
	return 0;
}

void write_results(char *fname) {
	FILE *out = fopen(fname, "w");
	qsort(results, nresults, 2*sizeof(quad), res_cmp);
	char s[100];
	for (int i=0; i < nresults; i++) {
		if (i != 0 && res_cmp(&results[2*i], &results[2*(i-1)]) == 0) continue;
		quad p = results[2*i];
		quad q = results[2*i+1];
		mpz_set_quad(tmp2, p);
		mpz_set_quad(tmp3, q);
		mpz_pow_ui(tmp2, tmp2, 3);
		mpz_pow_ui(tmp3, tmp3, 3);
		mpz_add(tmp2, tmp2, tmp3);
		mpz_get_str(s, 10, tmp2);
		fputs(s, out);
		// print_quad(p*p*p + q*q*q, out);
		fputc(' ', out);
		print_quad(p, out);
		fputc(' ', out);
		print_quad(q, out);
		fputc('\n', out);
	}
	fclose(out);
}

FILE *diagnostic_file;
void add_diagnostic_result(int a, int b, int c, quad p, quad q, quad r, quad s) {
	fprintf(diagnostic_file, "%d %d %d ", a, b, c);
	print_quad(p, diagnostic_file);
	fputc(' ', diagnostic_file);
	print_quad(q, diagnostic_file);
	fputc(' ', diagnostic_file);
	print_quad(r, diagnostic_file);
	fputc(' ', diagnostic_file);
	print_quad(s, diagnostic_file);
	fputc('\n', diagnostic_file);
	fflush(diagnostic_file);
}

bool is_prime(quad p) {
	mpz_set_quad(primecheck_tmp, p);
	return mpz_probab_prime_p(primecheck_tmp, 30);
}

bool fermat_check(quad p, int b) {
	if (p < 0) p = -p;
	mpz_set_quad(primecheck_tmp, p);
	mpz_set_ui(tmp2, b);
	mpz_sub_ui(tmp3, primecheck_tmp, 1);
	mpz_powm(tmp2, tmp2, tmp3, primecheck_tmp);
	return mpz_cmp_ui(tmp2, 1) == 0;
}

int64_t realmod(int64_t x, int64_t p) {
	if (x >= 0) return x % p;
	int64_t z = p - ((-x) % p);
	if (z == p) return 0;
	return z;
}

int64_t calc_modular(int64_t p, int x, int64_t a, int64_t b, int64_t c) {
	int64_t v = (a*a*a*a + 2*a*a*a*b + 3*a*a*b*b + 2*a*b*b*b + b*b*b*b) % p;
	switch(x) {
		case 0:
			return realmod((v + (c*c*c)*(2*a+b)), p);
		case 1:
			return realmod(-(v - (a-b)*(c*c*c)), p);
		case 2:
			return realmod(c*(-a*a*a + b*b*b + c*c*c), p);
		case 3:
			return realmod(-c*(2*a*a*a + 3*a*a*b + 3*a*b*b + b*b*b) - c*c*c*c, p);
	}
	return -1;
}

int quad_cmp(const void *x, const void *y) {
	quad a = *((quad *) x);
	quad b = *((quad *) y);
	if (a < b) return -1;
	if (a > b) return 1;
	return 0;
}

int headroom = 128;
bool check(quad *x, int64_t a, int64_t b, int64_t c) {
	quad a2 = a*a;
	quad b2 = b*b;
	quad c2 = c*c;
	quad a3 = a2 * a;
	quad b3 = b2 * b;
	quad c3 = c2 * c;

	quad v = a2*a2 + a3*(2*b) + a2*(3*b*b) + b3*(2*a) + b2*b2;
	x[0] = v + c3*(2*a+b); 
	x[1] = -(v - (a-b)*c3);
	x[2] = c * (-a3 + b3 + c3);
	x[3] = -c * (2*a3 + 3*a2*b + 3*a*b2 + b3) - c2*c2;

	for (int i=0; i < 4; i++) {
		if (x[i] == 0) {
			return false;
		}
		quad z = x[i] < 0 ? -x[i] : x[i];
		uint64_t hi = (uint64_t) (z >> 64);
		if (hi != 0) {
			int h = __builtin_clzll(hi);
			if (h < headroom) headroom = h;
		} else {
			uint64_t lo = (uint64_t) z;
			if (lo != 0) {
				int h = __builtin_clzll(lo) + 64;
				if (h < headroom) headroom = h;
			}
		}
	}
	qsort(x, 4, sizeof(quad), quad_cmp);
	
	if (x[0] + x[3] == 0 || x[1] + x[2] == 0) {
		return false;
	}

	quad d = gcd(x[0], gcd(x[1], gcd(x[2], x[3])));
	for (int i=0; i < 4; i++) {
		x[i] /= d;
		if (x[i] == 1 || x[i] == -1) return false;
		if (!fermat_check(x[i], 2)) {
			return false;
		}
	}

	for (int i=0; i < 4; i++) {
		if (!is_prime(x[i])) return false;
	}
	// normalize: largest in absolute value is positive.
	if (x[0] < 0 && -x[0] > x[3]) {
		quad t = x[0];
		x[0] = x[3];
		x[3] = t;
		t = x[1];
		x[1] = x[2];
		x[2] = t;
		for (int i=0; i < 4; i++) x[i] = -x[i];
	}
	add_result(-x[0], -x[1]);
	add_result(x[2], x[3]);
	add_diagnostic_result(a,b,c,-x[0],-x[1],x[2],x[3]);
	return true;
}

uint64_t continue_mask(uint64_t mask, int p) {
	int ti = p;
	uint64_t basemask = mask;
	while (ti < 64) {
		mask = (mask << p) | basemask;
		ti += p;
	}
	return mask;
}

uint64_t rotate(uint64_t x, int c, int p) {
	if (c == 0) return x;
	uint64_t next = x >> (64%p);
	x >>= c;
	return x | (next << (64-c));
}

uint64_t **seed_masks;
uint16_t ***seed_lists;

void init_seed(int a) {
	int16_t templist[fbb];
	for (int i=0; i < fblen; i++) {
		int p = fb[i];
		for (int b=0; b < p; b++) {
			for (int i=0; i < fbb; i++) templist[i] = 0;
			uint64_t mask = 0;
			int w = 0;
			for (int c = 0; c < p; c++) {
				int ct = 0;
				for (int x=0; x < 4; x++) {
					if (calc_modular(p, x, a, b, c) == 0) {
						ct++;
					}
				}
				if (ct >= 1 && ct <= 3) {
					if (p < 64) {
						mask |= 1L << c;
					} else {
#ifdef USE_SEEDLISTS
						templist[w++] = c;
#endif
					}
				}
			}
			if (p < 64) {
				seed_masks[i][b] = continue_mask(mask, p);
			} else {
#ifdef USE_SEEDLISTS
				int idx = b;
				seed_lists[i][idx] = (uint16_t *) malloc(2 * (w+1));
				memcpy(seed_lists[i][idx] + 1, templist, 2*w);
				seed_lists[i][idx][0] = w;
#endif
			}
		}
	}
}

void clear_seed_lists() {
#ifdef USE_SEEDLISTS
	for (int i=0; i < fblen; i++) {
		if (fb[i] >= 64) {
			for (int j=0; j < fb[i]; j++) {
				free(seed_lists[i][j]);
				seed_lists[i][j] = NULL;
			}
		}
	}
#endif
}

void lattice_sieve(int a, int bstart, int bend, int cstart, int cend) {
	long time = clock();
	int gcd_cull = 0;
	int sl_cull = 0;
	int ncheck = 0;
	long npos = 0;
	int nc = cend - cstart;
	int nchunks = (nc + 63)/64;
	for (int b = bstart; b < bend; b++) {
		uint64_t s[nchunks];
		memset(s, 0, nchunks * sizeof(uint64_t));
		int cend_altbound = cstart + a + b - bstart;
		int local_cend = cend_altbound < cend ? cend_altbound : cend;
		npos += local_cend - cstart;
		if (local_cend < cstart) continue;
		int final_chunk = (local_cend - cstart + 63)/64; 
		for (int i=0; i < fblen; i++) {
			int p = fb[i];
			if (p >= 64) continue;
			int sb = (int) realmod(b, p);
			int sc = (int) realmod(cstart, p);
			uint64_t masks[p];
			int advance = 64%p;
			masks[0] = seed_masks[i][sb];
			for (int sh = 1; sh < p; sh++) {
				masks[sh] = rotate(masks[sh-1], 1, p);
			}
			int mask_idx = sc;
			for (int chunk = 0; chunk < final_chunk; chunk++) {
				s[chunk] |= masks[mask_idx];
				mask_idx += advance;
				if (mask_idx >= p) mask_idx -= p;
			}
		}
		int abgcd = gcd(a,b);
		for (int chunk = 0; chunk < final_chunk; chunk++) {
			uint64_t acc = ~s[chunk];
			// acc now has 1s where we should check.
			if (acc != 0) {
				for (int bit = 0; bit < 64; bit++) {
					if ((acc & 1) == 1) {
						int c = cstart + chunk*64 + bit;
						if (c >= local_cend) break;
						if (abgcd > 1 && gcd(c, abgcd) > 1) {
							gcd_cull++;
						} else {
							bool ok = true;
#ifdef USE_SEEDLISTS
							for (int i=0; i < fblen; i++) {
								int p = fb[i];
								if (p >= 64) {
									int sb = (int) realmod(b, p);
									int sc = (int) realmod(c, p);
									uint16_t *list = seed_lists[i][sb];
									int len = list[0];
									for (int j=0; j < len; j++) {
										if (list[1+j] == sc) {
											ok = false;
											break;
										}
									}
									if (ok == false) break;
								}
							}
#endif
							if (ok) {
								ncheck++;
								quad x[4];
								if (check(&x[0], a, b, c)) {
									printf("(a,b,c) = %d, %d, %d: ", a, b, c);
									print_check_result(&x[0]);
								}
							} else {
								sl_cull++;
							}
						}
					}
					acc = acc >> 1;
				}
			}
		}
	}
	printf("Sieve took %ldms\n", (clock()-time)/1000);
	printf("Lattice sieve searched %ld entries, GCD culled %d, SL culled %d, checked %d\n", npos, gcd_cull, sl_cull, ncheck);
}


void init(int astart, int amax, int bound) {
	mpz_inits(primecheck_tmp, tmp2, tmp3, NULL);
	seed_masks = (uint64_t **) malloc(sizeof(uint64_t*) * fblen);
	seed_lists = (uint16_t ***) malloc(sizeof(uint16_t**) * fblen);
	for (int i=0; i < fblen; i++) {
		seed_masks[i] = (uint64_t *) malloc(8 * fb[i]);
		seed_lists[i] = (uint16_t **) malloc(sizeof(uint16_t*) * fb[i]);
	}
    results	= (quad *) malloc(res_cap*2*sizeof(quad));
	char res_file_name[50];
	sprintf(&res_file_name[0], "presults/diagnostic_%d_%d_%d", astart, amax, bound);
 	diagnostic_file = fopen(res_file_name, "w");
}

int main(int argc, const char *argv[]) {
	fbb = 256;
	sieve(fbb);
	printf("Sieve complete, fblen = %ld\n", fblen);
	int astart = atoi(argv[1]);
	int amax = atoi(argv[2]);
	int bound = atoi(argv[3]);
	int a = astart;
	init(astart, amax, bound);
	printf("Init complete\n");

	while (a < amax) {
		printf("==== Starting a = %d ====\n", a);
		long time = clock();
		init_seed(a);
		printf("Init seed complete in %ldms\n", (clock()-time)/1000);
		lattice_sieve(a, -bound, bound, -bound, bound);
		printf("Headroom: %d\n", headroom);
		a += 4;
		clear_seed_lists();
	}
	fclose(diagnostic_file);
	char res_file_name[40];
	sprintf(&res_file_name[0], "presults/results_%d_%d_%d", astart, amax, bound);
	write_results(res_file_name);
}
