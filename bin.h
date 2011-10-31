#ifndef BIN_H
#define BIN_H

#include "particle.h"

typedef struct bin_t {
	float len;						/* len of an edge	*/
	int n;							/* #bins in a row	*/
	particle_t** restrict bins;		/* points to bins	*/
	float* restrict center;			/* center of bins	*/
	int* restrict size;				/* #particles in bin*/
} bin_t;

bin_t* alloc_bin(float h);
void free_bin(bin_t* b);
void clear_bins(bin_t* b);
void assign_bin(bin_t* b, particle_t* p);
particle_t* find_vertical_bin(bin_t* b, particle_t* p);
particle_t* find_horizontal_bin(bin_t* b, particle_t* p);
particle_t* find_diagonal_bin(bin_t* b, particle_t* p);

#endif /* BIN_H */
