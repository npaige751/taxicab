#include "quad.h"
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

quad quad_init (uint64_t lo) {
	quad res;
	res.lo = lo;
	res.hi = 0;	
	return res;
}

quad quad_init2 (uint64_t hi, uint64_t lo) {
	quad res;
	res.lo = lo;
	res.hi = hi;
	return res;
}

void mpz_set_quad (mpz_t res, quad x) {
	mpz_set_ui (res, x.hi);
	mpz_mul_2exp (res, res, 64);
	mpz_add_ui (res, res, x.lo);
}

char *quad_get_str (quad q) {
	mpz_t x;
	mpz_init (x);
	mpz_set_quad (x, q);
	char *res = mpz_get_str (NULL, 10, x);
	mpz_clear (x);
	return res;
}

quad quad_init_str (char *x) {
	mpz_t n;
	mpz_init (n);
	mpz_set_str (n, x, 10);
	quad q;
	q.lo = mpz_getlimbn (n, 0);
	q.hi = mpz_getlimbn (n, 1);
	mpz_clear (n);
	return q;
}

double quad_to_double (quad q) {
	double res = (double) q.hi;
	res *= (double) ULONG_MAX;
	res += (double) q.lo;
	return res;
}

quad double_to_quad (double x) {
	return quad_init2 ((uint64_t) (x / pow (2, 64)), (uint64_t) fmod(x, pow(2,64)));
}

quad sci_to_quad (char *sci) {	// string of form ####e##
	uint64_t mult = 0;
	while (*sci != '\0' && *sci != 'e') {
		mult = mult * 10 + (*sci) - '0';
		sci++;
	}
	if (*sci == '\0') {	// no exponent
		return quad_init (mult);
	}
	sci++;
	int exp = 0;
	while (*sci != '\0') {
		exp = exp * 10 + (*sci) - '0';
		sci++;
	}
	quad ten = quad_init (10);
	quad res = quad_init (1);
	for (int i=0; i < exp; i++) {
		res = quad_mul (res, ten);
	}
	return quad_mul (res, quad_init (mult));
}

/*
quad quad_cube (unsigned int a) {
	unsigned long al = (unsigned long) a;
	al *= al;
	quad x, y;
	x.hi = 0;
	y.hi = 0;
	x.lo = a;
	y.lo = al;
	return quad_mul (x, y);
}
*/

#ifndef ASM_CUBSUM
quad quad_cubsum (unsigned int a, unsigned int b) {
	unsigned long al = (unsigned long) a;
	unsigned long bl = (unsigned long) b;
	quad x, y;
	x.hi = 0;
	y.hi = 0;
	x.lo = al + bl;
	y.lo = al*al - al*bl + bl*bl;
	return quad_mul (x, y);
}
#endif

int quad_cmp (quad a, quad b) {
	if (a.hi < b.hi) return -1;
	if (a.hi > b.hi) return 1;
	if (a.lo < b.lo) return -1;
	if (a.lo > b.lo) return 1;
	return 0;
}

quad quad_shr (quad a, int amt) {
	quad res;
	res.lo = (a.lo >> amt) | (a.hi << (64 - amt));
	res.hi = a.hi >> amt;
	return res;
}

quad quad_shl (quad a, int amt);

/*
int main () {
	quad a, b;
	a.lo = 18;
	a.hi = 10;
	b.lo = 19;
	b.hi = 10;
	int eq = c_quad_cmp (a,b);
	int aeq = quad_cmp (a,b);
	printf ("C:   %d\n", eq);
	printf ("asm: %d\n", aeq);
}
*/
