#ifndef INTERACT_H
#define INTERACT_H

#include "params.h"
#include "state.h"
#include "particle.h"

void compute_density(sim_state_t* s, bin_t* b, particle_t* p,
		sim_param_t* params);
void compute_accel(sim_state_t* state, bin_t* b, particle_t* p,
		sim_param_t* params);

#endif /* INTERACT_H */
