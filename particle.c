#include <stdlib.h>

#include "particle.h"

particle_t* alloc_particle(int n)
{
    particle_t* p = (particle_t*) calloc(1, sizeof(particle_t));
    p->rho =  (float*) calloc(1, sizeof(float));
    p->x =    (float*) calloc(2, sizeof(float));
    p->vh =   (float*) calloc(2, sizeof(float));
    p->v =    (float*) calloc(2, sizeof(float));
    p->a =    (float*) calloc(2, sizeof(float));
    if (n > 1) {
        p->next = alloc_particle(n - 1);
		p->bnext = NULL;
	}
    return p;
}

void free_particle(particle_t* p)
{
    if (p->next != NULL)
    	free_particle(p->next);
    free(p->a);
    free(p->v);
    free(p->vh);
    free(p->x);
    free(p->rho);
    free(p);
}
