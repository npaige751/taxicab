.globl _quad_add
.globl quad_add
_quad_add:
quad_add:
	addq	%rsi, %rcx
	adcq	%rdi, %rdx
	movq	%rdx, %rax
	movq	%rcx, %rdx
	ret

.globl _quad_sub
.globl quad_sub
_quad_sub:
quad_sub:
	subq	%rcx, %rsi
	sbbq	%rdx, %rdi
	movq	%rdi, %rax
	movq	%rsi, %rdx
	ret

.globl _quad_mul
.globl quad_mul
_quad_mul:
quad_mul:
/*	rdi, rsi  is a.hi, a.lo
	rdx, rcx  is b.hi, b.lo
	want to return res.hi in rax, res.lo in rdx

	we have:
	   ah  al
	 * bh  bl
	 ---------  
	   |al*bl|
	 |ah*bl|
	 |al*bh|
   |ah*bh|
   
   We only want the low-order 128 bits of this. So we need to compute
   al*bl in 128-bit precision, but only need the low-order 64 bits of the
   al*bh and ah*bl multiplies.
   */

   /* We're going to do al*bh first since rdx (b.hi) will get clobbered by mul */
	movq	%rsi, %rax
	mulq	%rdx	/* RDX:RAX now has the 128 bit result of al*bh. We only care about RAX. */
	movq	%rax, %r8

	/* Now multiply ah*bl */
	movq	%rcx, %rax
	mulq	%rdi
	addq	%rax, %r8

	/* Now multiply al*bl */
	movq	%rcx, %rax
	mulq	%rsi
	addq	%rdx, %r8	/* High order bits get added to r8 this time */

	movq	%rax, %rdx
	movq	%r8, %rax
	ret

.globl _quad_cube
_quad_cube:
	/* Get value to cube in %rdi 
	   result.hi in rax
	   result.lo in rdx

		The multiplication to do is:
		- get [h l] = x^2
		Then:
		
		h   l
	   *    x
	   ------
	   [ xl ]
    [ xh ] 
		 
	*/
	
	movq	%rdi, %rax
	mulq	%rax	/* RDX:RAX now has x^2  = [h l] */

	movq	%rdx, %rsi	/* Copy H into %rsi */

	mulq	%rdi	/* L is in rax already */
					/* RDX:RAX now contains XL */
	xchg	%rax, %rsi	/* Save low 64 bits of XL in rsi, and get H into rax */
	movq	%rdx, %rcx	/* Save high 64 bits of XL in rcx */
	
	mulq	%rdi	/* Now RDX:RAX has XH */
	/* We will ignore the top 64 bits of XH. Conveniently, the low 64 bits just
	   need to have the high 64 of XL added to them and they're in the right place. */
	addq	%rcx, %rax
	movq	%rsi, %rdx
	ret
	
#ifdef ASM_CUBSUM

.globl _quad_cubsum
.globl quad_cubsum
_quad_cubsum:
quad_cubsum:
	/*	Inputs are in %rdi (a)  and %rsi (b). 
		result.hi in %rax
		result.lo in %rdx
		
		Compute (a+b)(a^2-ab+b^2), assuming each factor fits in 64 bits.
		Do this by first computing a+b, then making (a+b)^2 = a^2 + 2ab + b^2. Then compute ab, shift/add to multiply it by 3, and subtract from the (a+b)^2.
		Then multiply the terms together.
	*/

	movq	%rdi, %rax
	mulq	%rsi		/* RDX:RAX now contains a*b. Ignore rdx. */
	movq	%rax, %rcx	/* duplicate */
	shlq	$1, %rcx	/* Can we use leaq here instead? */
	addq	%rax, %rcx	/* rcx now has 3ab */
		
	addq	%rdi, %rsi	/* rsi now contains a+b. We don't need the individual values anymore, so it's OK that b got overwritten */
	movq	%rsi, %rax
	mulq	%rax		/* RDX:RAX now contains (a+b)^2 */
	subq	%rcx, %rax	/* rax now contains a^2-ab+b^2 */
	mulq	%rsi
	xchg	%rax, %rdx	/* get the order right */
	ret

#endif

.globl _quad_incr_cubsum
.globl quad_incr_cubsum
_quad_incr_cubsum:
quad_incr_cubsum:
	/* 	Inputs are  (%rdi, %rsi) --- prev value of x^3 + y^3
					(%rdx, %rcx) --- value of x^2
					%r8			 --- value of x
					%r9			 --- value of k

		Want to compute (x+k)^3 + y^3 --> (%rax, %rdx).
		Add 3x^2k + 3xk^2 + k^3 to prev; assume 3x^2k does not necessarily fit in 64 bits, but the other terms do.

	*/
.globl _quad_eq
.globl quad_eq
_quad_eq:
quad_eq:
	xorq	%rax, %rax
	cmpq	%rdi, %rdx
	sete	%al
	cmpq	%rsi, %rcx
	sete	%cl
	andb	%cl, %al
	ret

/*

.globl _quad_cmp
_quad_cmp:
	xorq	%r8, %r8
	movq	$2, %rax
	cmpq	%rdi, %rdx
	
	cmovb	%rax, %r8
	movq	$-2, %rax
	cmova	%rax, %r8

	xorq	%rax, %rax
	xorq	%rdx, %rdx

	cmpq	%rsi, %rcx
	seta	%al
	setb	%dl
	cmovb	%rdx, %rax
	
	addq	%r8, %rax
	ret
*/
