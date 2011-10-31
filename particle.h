#ifndef PARTICLE_T_H
#define PARTICLE_T_H 

typedef struct particle_t {
	int bn;		/* bin index	*/
    float* rho;		/* Density		*/
    float* x;		/* Position		*/
    float* vh;		/* Velocity (half step)	*/
    float* v;		/* Velocity (full step)	*/
    float* a;		/* Acceleration		*/

    struct particle_t* bnext;	/* Next particle in bin			*/
	struct particle_t* next;	/* Next particle in iterator	*/
} particle_t;

particle_t* alloc_particle(int n);
void free_particle(particle_t* p);

#endif /* PARTICLE_T_H */ 
