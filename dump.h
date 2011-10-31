#ifndef DUMP_H
#define DUMP_H

#include "bin.h"
#include "particle.h"
#include "state.h"

void dump(char* fname, sim_state_t* s, particle_t* p, int t);
void dump_bins(bin_t* b, int t);

#endif /* DUMP_H */
