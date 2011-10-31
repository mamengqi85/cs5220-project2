#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "bin.h"
#include "state.h"
#include "particle.h"

void leapfrog_start(sim_state_t* s, bin_t* b, particle_t* p, double dt);
void leapfrog_step(sim_state_t* s, bin_t* b, particle_t* p, double dt);

#endif /* LEAPFROG_H */
