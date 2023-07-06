#ifndef QUAD_H
#define QUAD_H

#include <stdint.h>
#include <gmp.h>

typedef struct {
	uint64_t hi;
	uint64_t lo;
} quad;

quad quad_init  (uint64_t lo);
quad quad_init2 (uint64_t hi, uint64_t lo);
quad quad_init_str (char *);

quad sci_to_quad (char *);
quad double_to_quad (double);
double quad_to_double (quad);
char *quad_get_str (quad);
void mpz_set_quad (mpz_t, quad);

extern quad quad_add (quad a, quad b);
extern quad quad_sub (quad a, quad b);
extern quad quad_mul (quad a, quad b);
extern int quad_eq (quad a, quad b);
int quad_cmp (quad a, quad b);

quad quad_cube (unsigned int a);
#ifdef ASM_CUBSUM
extern quad quad_cubsum (unsigned int a, unsigned int b);
#else
quad quad_cubsum (unsigned int a, unsigned int b);
#endif

quad quad_shr (quad a, int amt);
quad quad_shl (quad a, int amt);

#endif
