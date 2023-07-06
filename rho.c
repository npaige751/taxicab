#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>

void tdiv (mpz_t n, int bound){
	while (mpz_divisible_ui_p (n, 2)){
		mpz_divexact_ui (n, n, 2);
		printf("2 ");
	}
	long t = 3;
	while (t < bound){
		while (mpz_divisible_ui_p (n, t)){
			mpz_divexact_ui (n, n, t);
			printf("%ld ", t);
		}
		t += 2;
	}
}

#define GCD_BLOCK_SIZE 128

void rho (mpz_t n, unsigned int steps, int seed){	// printrem controls whether a composite cofactor is printed.
	mpz_t x, y, g, temp;
	mpz_inits (x, y, g, temp, NULL);
	mpz_set_ui (x, 2);
	mpz_set_ui (y, 2);

	int ct = 0;
	while (ct < steps){
		mpz_set_ui (g, 1);
		for (int i=0; i < GCD_BLOCK_SIZE; i++){
			mpz_mul    (temp, x, x);
			mpz_add_ui (temp, temp, seed);
			mpz_mod    (x, temp, n);

			mpz_mul	   (temp, y, y);
			mpz_add_ui (temp, temp, seed);
			mpz_mul    (y, temp, temp);
			mpz_add_ui (temp, y, seed);
			mpz_mod    (y, temp, n);

			mpz_sub    (temp, x, y);
			mpz_mul    (g, g, temp);
			mpz_mod    (g, g, n);
		}
		mpz_gcd (g, g, n);
		if (mpz_cmp_ui (g, 1) != 0){
			// found a factor
			if (mpz_probab_prime_p (g, 32)){	// a prime one!
				mpz_out_str (stdout, 10, g);
				printf(" ");
				mpz_divexact (n, n, g);
			} else {
				mpz_set(temp, g);
				rho(temp, 0xffffffff, seed+1);
				/*
				mpz_out_str (stdout, 10, g);
				printf(" ");
//				printf (" (composite)\n");
				*/
				mpz_divexact (n, n, g);
			}
			if (mpz_cmp_ui (n, 1) == 0){
				mpz_clears (x, y, g, temp, NULL);
				return;
			}
			if (mpz_probab_prime_p (n, 32)){
				mpz_out_str (stdout, 10, n);
				printf(" ");
				mpz_clears (x, y, g, temp, NULL);
				return;
			}
		}
		ct += GCD_BLOCK_SIZE;
	}
}

int main(int argc, const char *argv[]) {
	mpz_t n;
	mpz_init(n);
	while (!feof(stdin)) {
		int e = mpz_inp_str(n, stdin, 10);
		if (e == 0) break;
		tdiv(n, 50000);
		if (mpz_cmp_ui(n, 1) != 0) {
			if (mpz_probab_prime_p(n, 32)) {
				mpz_out_str(stdout, 10, n);
			} else {
				rho(n, 0xffffffff, 0);
			}
		}
		printf("\n");
	}
}

