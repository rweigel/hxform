#include <stdio.h>
#include "cxform.h"

int cxformv(const void *indatav, const void *timev, const char *from, const char *to, void *outdatav, int Nv, int Nt) {

    int i, es, ret;
    int debug = 0;
    double v_in[3];
    double v_out[3];

    const double *indata = (double *) indatav;
    const int *time = (int *) timev;
    double *outdata = (double *) outdatav;

    if (debug) printf("%s to %s\n", from, to);
    if (debug) printf("es = %d\n", es);
    for (i = 0; i < Nv; i++) {
        if (debug) printf("%d %d %d %d %d %d\n", time[6*i], time[6*i+1], time[6*i+2], time[6*i+3], time[6*i+4], time[6*i+5]);
        es = date2es(time[6*i], time[6*i+1], time[6*i+2], time[6*i+3], time[6*i+4], time[6*i+5]);
        v_in[0] = indata[3*i];
        v_in[1] = indata[3*i+1];
        v_in[2] = indata[3*i+2];
        if (debug) printf("%f %f %f\n",v_in[0],v_in[1],v_in[2]);
        ret = cxform(from, to, es, v_in, v_out);
        if (ret != 0) {
            return ret; /* return of 1 or 2 indicates error in from or to strings */
        }
        if (debug) printf("%f %f %f\n",v_out[0],v_out[1],v_out[2]);
        outdata[3*i] = v_out[0];
        outdata[3*i+1] = v_out[1];
        outdata[3*i+2] = v_out[2];
    }
    return ret;
}