#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "dump.h"
#include "bin.h"
#include "params.h"
#include "state.h"
#include "particle.h"
#include "interact.h"

/*@q
 * ====================================================================
 */

/*@T
 * \subsection{Density computations}
 * 
 * The formula for density is
 * \[
 *   \rho_i = \frac{4m}{\pi h^8} \sum_{j \in N_i} (h^2 - r^2)^3.
 * \]
 * We search for neighbors of node $i$ by checking every particle,
 * which is not very efficient.  We do at least take advange of
 * the symmetry of the update ($i$ contributes to $j$ in the same
 * way that $j$ contributes to $i$).
 *@c*/

void compute_density(sim_state_t* s, bin_t* b, particle_t* p,
		sim_param_t* params)
{
	//TODO
    int n = s->n;
    float h  = params->h;
    float h2 = h*h;
    float h8 = ( h2*h2 )*( h2*h2 );
    float C  = 4 * s->mass / M_PI / h8;

    particle_t* p_tmp = p;
	for (int i = 0; i < n; ++i) {
		p_tmp->rho[0] = 0;
		p_tmp = p_tmp->next;
	}

	p_tmp = p;
    particle_t* p_tmp2;
    for (int i = 0; i < n; ++i) {
        p_tmp->rho[0] += 4 * s->mass / M_PI / h2;
		int bn = p_tmp->bn;
        p_tmp2 = b->bins[bn];
		int size = b->size[bn];
        for (int j = 0; j < size; ++j) {
			if (p_tmp != p_tmp2) {
            	float dx = p_tmp->x[0]-p_tmp2->x[0];
            	float dy = p_tmp->x[1]-p_tmp2->x[1];
            	float r2 = dx*dx + dy*dy;
            	float z  = h2-r2;
            	if (z > 0) {
            	    float rho_ij = C*z*z*z;
            	    p_tmp->rho[0] += rho_ij;
            	}
			}
            p_tmp2 = p_tmp2->bnext;
        }

        p_tmp2 = find_vertical_bin(b, p_tmp);
		if (p_tmp2 != NULL) {
			int size = b->size[p_tmp2->bn];
        	for (int j = 0; j < size; ++j) {
        	    float dx = p_tmp->x[0]-p_tmp2->x[0];
        	    float dy = p_tmp->x[1]-p_tmp2->x[1];
        	    float r2 = dx*dx + dy*dy;
        	    float z  = h2-r2;
        	    if (z > 0) {
        	        float rho_ij = C*z*z*z;
        	        p_tmp->rho[0] += rho_ij;
        	    }
        	    p_tmp2 = p_tmp2->bnext;
        	}
		}

		p_tmp2 = find_horizontal_bin(b, p_tmp);
		if (p_tmp2 != NULL) {
			int size = b->size[p_tmp2->bn];
        	for (int j = 0; j < size; ++j) {
        	    float dx = p_tmp->x[0]-p_tmp2->x[0];
        	    float dy = p_tmp->x[1]-p_tmp2->x[1];
        	    float r2 = dx*dx + dy*dy;
        	    float z  = h2-r2;
        	    if (z > 0) {
        	        float rho_ij = C*z*z*z;
        	        p_tmp->rho[0] += rho_ij;
        	    }
        	    p_tmp2 = p_tmp2->bnext;
        	}
		}

		p_tmp2 = find_diagonal_bin(b, p_tmp);
		if (p_tmp2 != NULL) {
			size = b->size[p_tmp2->bn];
        	for (int j = 0; j < size; ++j) {
        	    float dx = p_tmp->x[0]-p_tmp2->x[0];
        	    float dy = p_tmp->x[1]-p_tmp2->x[1];
        	    float r2 = dx*dx + dy*dy;
        	    float z  = h2-r2;
        	    if (z > 0) {
        	        float rho_ij = C*z*z*z;
        	        p_tmp->rho[0] += rho_ij;
        	    }
        	    p_tmp2 = p_tmp2->bnext;
        	}
		}
        p_tmp = p_tmp->next;
    }
}

/*@T
 * \subsection{Computing forces}
 * 
 * The acceleration is computed by the rule
 * \[
 *   \bfa_i = \frac{1}{\rho_i} \sum_{j \in N_i} 
 *     \bff_{ij}^{\mathrm{interact}} + \bfg,
 * \]
 * where the pair interaction formula is as previously described.
 * Like [[compute_density]], the [[compute_accel]] routine takes
 * advantage of the symmetry of the interaction forces
 * ($\bff_{ij}^{\mathrm{interact}} = -\bff_{ji}^{\mathrm{interact}}$)
 * but it does a very expensive brute force search for neighbors.
 *@c*/
void compute_accel(sim_state_t* state, bin_t* b, particle_t* part,
		sim_param_t* params)
{
    // Unpack basic parameters
    const float h    = params->h;
    const float rho0 = params->rho0;
    const float k    = params->k;
    const float mu   = params->mu;
    const float g    = params->g;
    const float mass = state->mass;
    const float h2   = h*h;

    // Unpack system state
    int n = state->n;
    particle_t* p_tmp;

    // Compute density and color
    compute_density(state, b, part, params);

    // Start with gravity and surface forces
    p_tmp = part;
    for (int i = 0; i < n; ++i) {
        p_tmp->a[0] = 0;
        p_tmp->a[1] = -g;
        p_tmp = p_tmp->next;
    }

    // Constants for interaction term
    float C0 = mass / M_PI / ( (h2)*(h2) );
    float Cp =  15*k;
    float Cv = -40*mu;

    // Now compute interaction forces
    p_tmp = part;
    particle_t* p_tmp2;

    for (int i = 0; i < n; ++i) {
	//dump("dump_hash.out", state, part, 5);
        const float rhoi = p_tmp->rho[0];
		int bn = p_tmp->bn;
        p_tmp2 = b->bins[bn];
		int size = b->size[bn];
        for (int j = 0; j < size; ++j) {
			if (p_tmp != p_tmp2) {
            	float dx = p_tmp->x[0]-p_tmp2->x[0];
            	float dy = p_tmp->x[1]-p_tmp2->x[1];
            	float r2 = dx*dx + dy*dy;
            	if (r2 < h2) {
            	    const float rhoj = p_tmp2->rho[0];
            	    float q = sqrt(r2)/h;
            	    float u = 1-q;
            	    float w0 = C0 * u/rhoi/rhoj;
            	    float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
            	    float wv = w0 * Cv;
            	    float dvx = p_tmp->v[0]-p_tmp2->v[0];
            	    float dvy = p_tmp->v[1]-p_tmp2->v[1];
            	    p_tmp->a[0] += (wp*dx + wv*dvx);
            	    p_tmp->a[1] += (wp*dy + wv*dvy);
            	}
			}
            p_tmp2 = p_tmp2->bnext;
        }

		p_tmp2 = find_vertical_bin(b, p_tmp);
		if (p_tmp2 != NULL) {
			size = b->size[p_tmp2->bn];
			for (int j = 0; j < size; ++j) {
            	float dx = p_tmp->x[0]-p_tmp2->x[0];
            	float dy = p_tmp->x[1]-p_tmp2->x[1];
            	float r2 = dx*dx + dy*dy;
            	if (r2 < h2) {
            	    const float rhoj = p_tmp2->rho[0];
            	    float q = sqrt(r2)/h;
            	    float u = 1-q;
            	    float w0 = C0 * u/rhoi/rhoj;
            	    float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
            	    float wv = w0 * Cv;
            	    float dvx = p_tmp->v[0]-p_tmp2->v[0];
            	    float dvy = p_tmp->v[1]-p_tmp2->v[1];
            	    p_tmp->a[0] += (wp*dx + wv*dvx);
            	    p_tmp->a[1] += (wp*dy + wv*dvy);
            	    //p_tmp2->a[0] -= (wp*dx + wv*dvx);
            	    //p_tmp2->a[1] -= (wp*dy + wv*dvy);
            	}
				p_tmp2 = p_tmp2->bnext;
			}
		}

		p_tmp2 = find_horizontal_bin(b, p_tmp);
		if (p_tmp2 != NULL) {
			size = b->size[p_tmp2->bn];
			for (int j = 0; j < size; ++j) {
            	float dx = p_tmp->x[0]-p_tmp2->x[0];
            	float dy = p_tmp->x[1]-p_tmp2->x[1];
            	float r2 = dx*dx + dy*dy;
            	if (r2 < h2) {
            	    const float rhoj = p_tmp2->rho[0];
            	    float q = sqrt(r2)/h;
            	    float u = 1-q;
            	    float w0 = C0 * u/rhoi/rhoj;
            	    float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
            	    float wv = w0 * Cv;
            	    float dvx = p_tmp->v[0]-p_tmp2->v[0];
            	    float dvy = p_tmp->v[1]-p_tmp2->v[1];
            	    p_tmp->a[0] += (wp*dx + wv*dvx);
            	    p_tmp->a[1] += (wp*dy + wv*dvy);
 
            	    //p_tmp2->a[0] -= (wp*dx + wv*dvx);
            	    //p_tmp2->a[1] -= (wp*dy + wv*dvy);
            	}
				p_tmp2 = p_tmp2->bnext;
			}
		}

		p_tmp2 = find_diagonal_bin(b, p_tmp);
		if (p_tmp2 != NULL) {
			size = b->size[p_tmp2->bn];
			for (int j = 0; j < size; ++j) {
            	float dx = p_tmp->x[0]-p_tmp2->x[0];
            	float dy = p_tmp->x[1]-p_tmp2->x[1];
            	float r2 = dx*dx + dy*dy;
            	if (r2 < h2) {
            	    const float rhoj = p_tmp2->rho[0];
            	    float q = sqrt(r2)/h;
            	    float u = 1-q;
            	    float w0 = C0 * u/rhoi/rhoj;
            	    float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
            	    float wv = w0 * Cv;
            	    float dvx = p_tmp->v[0]-p_tmp2->v[0];
            	    float dvy = p_tmp->v[1]-p_tmp2->v[1];
            	    p_tmp->a[0] += (wp*dx + wv*dvx);
            	    p_tmp->a[1] += (wp*dy + wv*dvy);
            	    //p_tmp2->a[0] -= (wp*dx + wv*dvx);
            	    //p_tmp2->a[1] -= (wp*dy + wv*dvy);
            	}
				p_tmp2 = p_tmp2->bnext;
			}
		}

        p_tmp = p_tmp->next;
    }
	//dump("dump_hash_fin.out", state, part, 601);
}

