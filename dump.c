#include <stdio.h>

#include "bin.h"
#include "particle.h"
#include "state.h"
#include "dump.h"

static int toDump = 1;
static int toDump2 = 1;

/* when the function is called at the tst time, outout the particle state to the file */
void dump(char* fname, sim_state_t* s, particle_t* p, int t)
{
	if (toDump == t) {
		FILE* fp = fopen(fname, "w");
		int n = s->n;
		fprintf(fp, "n = %d, mass = %f\n", s->n, s->mass);
		particle_t* p_tmp = p;
		for (int i = 0; i < n; ++i) {
			float rhoi = p_tmp->rho[0];
			float xi = p_tmp->x[0];
			float yi = p_tmp->x[1];
			float vxi = p_tmp->v[0];
			float vyi = p_tmp->v[1];
			float vhxi = p_tmp->vh[0];
			float vhyi = p_tmp->vh[1];
			float axi = p_tmp->a[0];
			float ayi = p_tmp->a[1];
			int bni = p_tmp->bn;

			fprintf(fp, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n", i, rhoi, xi, yi, vxi, vyi, vhxi, vhyi, axi, ayi, bni);
			p_tmp = p_tmp->next;
		}
		fclose(fp);
	}
	toDump++;
}

/* when the function is called at the tst time, outout the bin state to the file */
void dump_bins(bin_t* b, int t)
{
	if (t == toDump2) {
		FILE* fp = fopen("dumpbin.out", "w");
		for (int i = 0; i < b->n; ++i) {
			for (int j = 0; j < b->n; ++j) {
				fprintf(fp, "%d,%d:%d\t", i, j, b->size[i*b->n+j]);
				particle_t* p = b->bins[i*b->n+j];
				for (int k = 0; k < b->size[i*b->n+j]; ++k) {
					fprintf(fp, "%f,%f\t", p->x[0], p->x[1]);
					p = p->bnext;
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);
	}
	toDump2++;
}
