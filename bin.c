#include <stdlib.h>

#include "dump.h"
#include "particle.h"
#include "bin.h"

bin_t* alloc_bin(float h)
{
	int n = 1/h/2;
	bin_t* b = (bin_t*) calloc(1, sizeof(bin_t));
	b->bins = (particle_t**) calloc(n*n, sizeof(particle_t*));
	b->len = 2*h;
	b->n = n; 
	b->center = (float*) calloc(2*n*n, sizeof(float));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			b->center[2*(i*n+j)+0] = h*(2*j+1); 
			b->center[2*(i*n+j)+1] = h*(2*i+1); 
		}
	}
	b->size = (int*) calloc(n*n, sizeof(int));
	return b;
}

void free_bin(bin_t* b)
{
	free(b->size);
	free(b->center);
	free(b);
}

void clear_bins(bin_t* b)
{
	for (int i = 0; i < b->n * b->n; ++i) {
		b->bins[i] = NULL;
		b->size[i] = 0;
	}
}

void assign_bin(bin_t* b, particle_t* p)
{
	int bx = p->x[0] / b->len;
	int by = p->x[1] / b->len;
	if (bx < 0) bx = 0;
	if (bx > b->n - 1) bx = b->n - 1;
	if (by < 0) by = 0;
	if (by > b->n - 1) by = b->n - 1;
	int bn = by * b->n + bx;
	p->bn = bn;
	p->bnext = b->bins[bn];
	b->bins[bn] = p;
	//TODO: use ++
	b->size[bn] = b->size[bn] + 1;
}

particle_t* find_vertical_bin(bin_t* b, particle_t* p)
{
	const float e = -1e-9;
	int bn = p->bn;
	int by = bn / b->n;
	if (p->x[1] - b->center[bn*2+1] >= e) {	
		if (by != b->n - 1)
			return b->bins[bn + b->n];
		else
			return NULL;
	}
	if (b->center[bn*2+1] - p->x[1] >= e) {
		if (by != 0)
			return b->bins[bn - b->n];
		else
			return NULL;
	}
	return NULL;
}
		
particle_t* find_horizontal_bin(bin_t* b, particle_t* p)
{
	const float e = -1e-9;
	int bn = p->bn;
	int bx = bn % b->n;
	if (p->x[0] - b->center[bn*2+0] >= e) {	
		if (bx != b->n - 1)
			return b->bins[bn + 1];
		else
			return NULL;
	}
	if (b->center[bn*2+0] - p->x[0] >= e) {
		if (bx != 0) 
			return b->bins[bn - 1];
		else
			return NULL;
	}
	return NULL;
}
		
particle_t* find_diagonal_bin(bin_t* b, particle_t* p)
{
	const float e = -1e-9;
	int bn = p->bn;
	int bx = bn % b->n;
	int by = bn / b->n;
	if (p->x[0] - b->center[bn*2+0] >= e && p->x[1] - b->center[bn*2+1] >= e) {
		if (bx == b->n - 1)
			return NULL;
		if (by == b->n - 1)
			return NULL;
		return b->bins[bn + b->n + 1];
	}
	if (p->x[0] - b->center[bn*2+0] >= e && b->center[bn*2+1] - p->x[1] >= e) {
		if (bx == b->n - 1) 
			return NULL;
		if (by == 0)
			return NULL;
		return b->bins[bn - b->n + 1];
	}
	if (b->center[bn*2+0] - p->x[0] >= e && p->x[1] - b->center[bn*2+1] >= e) {
		if (bx == 0)
			return NULL;
		if (by == b->n - 1)
			return NULL;
		return b->bins[bn + b->n - 1];
	}
	if (b->center[bn*2+0] - p->x[0] >= e && b->center[bn*2+1] - p->x[1] >= e) {
		if (bx == 0)
			return NULL;
		if (by == 0)
			return NULL;
		return b->bins[bn - b->n - 1];
	}
	return NULL;
}
