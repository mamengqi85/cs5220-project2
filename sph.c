#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "io.h"
#include "dump.h"
#include "bin.h"
#include "params.h"
#include "state.h"
#include "particle.h"
#include "interact.h"
#include "leapfrog.h"
#include "timing.h"

#define N_THREADS 8

/*@q
 * ====================================================================
 */

/*@T
 * \section{Initialization}
 *
 * We've hard coded the computational domain to a unit box, but we'd prefer
 * to do something more flexible for the initial distribution of fluid.
 * In particular, we define the initial geometry of the fluid in terms of an
 * {\em indicator function} that is one for points in the domain occupied
 * by fluid and zero elsewhere.  A [[domain_fun_t]] is a pointer to an
 * indicator for a domain, which is a function that takes two floats and
 * returns 0 or 1.  Two examples of indicator functions are a little box
 * of fluid in the corner of the domain and a circular drop.
 *@c*/
typedef int (*domain_fun_t)(float, float);

int box_indicator(float x, float y)
{
    return (x < 0.5) && (y < 0.5);
}

int circ_indicator(float x, float y)
{
    float dx = (x-0.5);
    float dy = (y-0.3);
    float r2 = dx*dx + dy*dy;
    return (r2 < 0.25*0.25);
}

/*@T
 *
 * The [[place_particles]] routine fills a region (indicated by the
 * [[indicatef]] argument) with fluid particles.  The fluid particles
 * are placed at points inside the domain that lie on a regular mesh
 * with cell sizes of $h/1.3$.  This is close enough to allow the
 * particles to overlap somewhat, but not too much.
 *@c*/
particle_t* place_particles(sim_state_t* s, bin_t* b, sim_param_t* param, 
                             domain_fun_t indicatef)
{
    float h  = param->h;
    float hh = h/1.3;

    // Count mesh points that fall in indicated region.
    int count = 0;
    //TODO:omp does not accept float
//    #pragma omp parallel	\
	shared(hh, indicatef)	\
	private(x, y)		\
	reduction(+:count)
    {
    	for (float x = 0; x < 1; x += hh)
    	    for (float y = 0; y < 1; y += hh)
    	        count += indicatef(x,y);
    }

    // Populate the particle data structure
    s->n = count;
    particle_t* p = alloc_particle(count);
    particle_t* p_tmp = p;
    //TODO:omp:float?
    for (float x = 0; x < 1; x += hh) {
        for (float y = 0; y < 1; y += hh) {
            if (indicatef(x,y)) {
                p_tmp->x[0] = x;
                p_tmp->x[1] = y;
                p_tmp->v[0] = 0;
                p_tmp->v[1] = 0;
				assign_bin(b, p_tmp);
                p_tmp = p_tmp->next;
            }
        }
    }
    return p;
}

/*@T
 *
 * The [[place_particle]] routine determines the initial particle
 * placement, but not the desired mass.  We want the fluid in the
 * initial configuration to exist roughly at the reference density.
 * One way to do this is to take the volume in the indicated body of
 * fluid, multiply by the mass density, and divide by the number of
 * particles; but that requires that we be able to compute the volume
 * of the fluid region.  Alternately, we can simply compute the
 * average mass density assuming each particle has mass one, then use
 * that to compute the particle mass necessary in order to achieve the
 * desired reference density.  We do this with [[normalize_mass]].
 * @c*/
void normalize_mass(sim_state_t* s, bin_t* b, particle_t* p, sim_param_t* param)
{
    s->mass = 1;
    compute_density(s, b, p, param);
    float rho0 = param->rho0;
    float rho2s = 0;
    float rhos  = 0;
    int i;
    particle_t* p_tmp = p;
//    #pragma omp parallel for num_threads(N_THREADS)	\
	shared(s, p)	\
	private(i)	\
        reduction(+:rho2s, rhos)
    for (i = 0; i < s->n; ++i) {
        rho2s += (p_tmp->rho[0])*(p_tmp->rho[0]);
        rhos  += p_tmp->rho[0];
        p_tmp = p_tmp->next;
    }
    s->mass *= ( rho0*rhos / rho2s );
}

particle_t* init_particles(sim_state_t* s, bin_t* b, sim_param_t* param)
{
    particle_t* p = place_particles(s, b, param, box_indicator);
    normalize_mass(s, b, p, param);
    return p;
}

/*@T
 * \section{The [[main]] event}
 *
 * The [[main]] routine actually runs the time step loop, writing
 * out files for visualization every few steps.  For debugging
 * convenience, we use [[check_state]] before writing out frames,
 * just so that we don't spend a lot of time on a simulation that
 * has gone berserk.
 *@c*/

void check_state(sim_state_t* s, particle_t* p)
{
    int i;
    float xi, yi;
    particle_t* p_tmp = p;
    //TODO:become slower if uncommented.
    //#pragma omp parallel for num_threads(N_THREADS)	\
	shared(s)	\
	private(i, xi, yi)
    for (i = 0; i < s->n; ++i) {
        xi = p_tmp->x[0];
        yi = p_tmp->x[1];
        assert( xi >= 0 || xi <= 1 );
        assert( yi >= 0 || yi <= 1 );
        p_tmp = p_tmp->next;
    }
}

int main(int argc, char** argv)
{
    sim_param_t params;
    sim_state_t s_instance;
    if (get_params(argc, argv, &params) != 0)
        exit(-1);
	bin_t* bin = alloc_bin(params.h);
    sim_state_t* state = &s_instance;
    particle_t* part = init_particles(state, bin, &params);
    FILE* fp    = fopen(params.fname, "w");
    int nframes = params.nframes;
    int npframe = params.npframe;
    float dt    = params.dt;
    int n       = state->n;

    tic(0);
    write_header(fp, n);
    write_frame_data(fp, n, part, NULL);
    compute_accel(state, bin, part, &params);
    leapfrog_start(state, bin, part, dt);
	//dump_bins(bin,1);
    check_state(state, part);
    for (int frame = 1; frame < nframes; ++frame) {
        for (int i = 0; i < npframe; ++i) {
			if (frame == 71 && i == 0) {
            	compute_accel(state, bin, part, &params);
			}
			else
            compute_accel(state, bin, part, &params);
			//dump("dump_hash.out", state, part, 7000);
            leapfrog_step(state, bin, part, dt);
            check_state(state, part);
        }
        write_frame_data(fp, n, part, NULL);
    }
    printf("Ran in %g seconds\n", toc(0));

    fclose(fp);
    free_particle(part);
	free_bin(bin);
}
