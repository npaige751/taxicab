taxi2:	taxi2.c taxi2.h quad.c quad.s quad.h
	gcc -O3 -std=c99 -o taxi2 taxi2.c quad.c quad.s -D ASM_CUBSUM -I ./ -lgmp -lm -lpthread

taxi:	taxi.c taxi.h quad.c quad.s quad.h
	gcc -O3 -std=c99 -o taxi taxi.c quad.c quad.s -D ASM_CUBSUM -I ./ -lgmp -lm

parametric: parametric.c
	gcc -O3 -std=c99 -o parametric parametric.c -I ./ -lgmp -lm

rho: rho.c
	gcc -O3 -std=c99 -o rho rho.c -I ./ -lgmp -lm

clean:
	rm -f taxi taxi2 parametric rho
