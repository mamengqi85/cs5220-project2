#include <stdio.h>

#include "particle.h"
#include "io.h"

#define VERSION_TAG "SPHView00 "

void write_header(FILE* fp, int n)
{
    fprintf(fp, "%s%d 1\n", VERSION_TAG, n);
}


void write_frame_data(FILE* fp, int n, particle_t* p, int* c)
{
    particle_t* p_tmp = p;
    for (int i = 0; i < n; ++i) {
        float xi = p_tmp->x[0];
        float yi = p_tmp->x[1];
        int   ci = c ? *c++ : 0;
        fprintf(fp, "%e %e %d\n", xi, yi, ci);
        p_tmp = p_tmp->next;
    }
}
