#ifndef STATE_H
#define STATE_H

/*@T
 * \section{System state}
 * 
 * The [[sim_state_t]] structure holds the information for the current
 * state of the system and of the integration algorithm.  The array
 * [[x]] has length $2n$, with [[ x[2*i+0] ]] and [[ x[2*i+1] ]] representing
 * the $x$ and $y$ coordinates of the particle positions.  The layout
 * for [[v]], [[vh]], and [[a]] is similar, while [[rho]] only has one
 * entry per particle.
 * 
 * The [[alloc_state]] and [[free_state]] functions take care of storage
 * for the local simulation state.
 *@c*/
typedef struct sim_state_t {
    int n;                /* Number of particles    */
    float mass;           /* Particle mass          */
} sim_state_t;

/*@q*/
#endif /* STATE_H */
