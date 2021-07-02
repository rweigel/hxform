#include <stdio.h>
#include "cxform.h"

int cxform_wrapper(const void *indatav, const void *timev, const char *from, const char *to, void *outdatav, int Nv, int Nt) {

    int i, es, ret;
    int debug = 0;
    double v_in[3];
    double v_out[3];

    const double *indata = (double *) indatav;
    const int *time = (int *) timev;
    double *outdata = (double *) outdatav;

    if (debug) printf("cxform_wrapper(): %s to %s\n", from, to);
    if (debug) printf("cxform_wrapper(): es = %d\n", es);
    if (Nv == Nt) {
        for (i = 0; i < Nv; i++) {
            if (debug) printf("cxform_wrapper(): time = [%d %d %d %d %d %d]\n", time[6*i], time[6*i+1], time[6*i+2], time[6*i+3], time[6*i+4], time[6*i+5]);
            es = date2es(time[6*i], time[6*i+1], time[6*i+2], time[6*i+3], time[6*i+4], time[6*i+5]);
            v_in[0] = indata[3*i];
            v_in[1] = indata[3*i+1];
            v_in[2] = indata[3*i+2];
            if (debug) printf("cxform_wrapper(): v_in  = [%f %f %f]\n",v_in[0],v_in[1],v_in[2]);
            ret = cxform(from, to, es, v_in, v_out);
            if (ret != 0) {
                return ret;
                // return of 1 or 2 indicates error in 'from' or 'to' strings.
                // TODO: Check for valid from/to strings at start of function.
                // https://github.com/edsantiago/cxform/blob/added-interfaces/cxform-auto.c#L33
                // J2000, GEI, GEO, MAG, GSE, GSM, SM, RTN, GSEQ, HEE, HAE, HEEQ
            }
            if (debug) printf("cxform_wrapper(): v_out = [%f %f %f]\n",v_out[0],v_out[1],v_out[2]);
            outdata[3*i] = v_out[0];
            outdata[3*i+1] = v_out[1];
            outdata[3*i+2] = v_out[2];
        }
    } else if (Nv == 1) {
        v_in[0] = indata[0];
        v_in[1] = indata[1];
        v_in[2] = indata[2];
        for (i = 0; i < Nt; i++) {
            es = date2es(time[6*i], time[6*i+1], time[6*i+2], time[6*i+3], time[6*i+4], time[6*i+5]);
            ret = cxform(from, to, es, indata, v_out);
            if (ret != 0) {return ret;}
            outdata[3*i] = v_out[0];
            outdata[3*i+1] = v_out[1];
            outdata[3*i+2] = v_out[2];
        }
    } else {
        es = date2es(time[0], time[1], time[2], time[3], time[4], time[5]);
        for (i = 0; i < Nv; i++) {
            v_in[0] = indata[3*i];
            v_in[1] = indata[3*i+1];
            v_in[2] = indata[3*i+2];
            ret = cxform(from, to, es, v_in, v_out);
            if (ret != 0) {return ret;}
            outdata[3*i] = v_out[0];
            outdata[3*i+1] = v_out[1];
            outdata[3*i+2] = v_out[2];
        }
    }
    return ret;
}