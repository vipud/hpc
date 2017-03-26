#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "type.h"
#include "s_class.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"
#include <openacc.h>

#define MAX(X,Y)    (((X) > (Y)) ? (X) : (Y))

#define MK          16
#define MM          (M - MK)
#define NN          (1 << MM)
#define NK          (1 << MK)
#define NQ          10
#define EPSILON     1.0e-8
#define A           1220703125.0
#define S           271828183.0

static double x[2*NK];
static double q[NQ];

int main(int argc, char* argv[])
{
    double  Mops, t1, t2, t3, t4, x1, x2;
    double  sx, sy, tm, an, tt, gc;
    double  sx_verify_value, sy_verify_value, sx_err, sy_err;
    int     np;
    int     i, ik, kk, l, k, nit;
    int     k_offset, j;
    logical verified, timers_enabled;

    double dum[3] = {1.0, 1.0, 1.0};
    char size[16];

    FILE *fp;

    double _r23 = 1.1920928955078125e-07;
    double _r46 = _r23 * _r23;
    double _t23 = 8.388608e+06;
    double _t46 = _t23 * _t23;
    double _t1, _t2, _t3, _t4, _a1, _a2, _x1, _x2, _z, _r;

    if ((fp = fopen("timer.flag", "r")) == NULL) {
        timers_enabled = false;
    } else {
        timers_enabled = true;
        fclose(fp);
    }

    sprintf(size, "%15.01f", pow(2.0, M+1));
    j = 14;
    if (size[j] == '.') j--;
    size[j+1] = '\0';
    printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - EP Benchmark\n");
    printf("\n Number of random numbers generated: %15s\n", size);
    verified = false;

    np = NN;

    dum[0] = randlc(&dum[1], dum[2]);

    #pragma acc parallel loop present(x[:])
    for(i = 0; i < 2*NK; i++){
        x[i] = -1.0e99;
    }

    Mops = log(sqrt(fabs(MAX(1.0, 1.0))));

    timer_clear(0);
    timer_clear(1);
    timer_clear(2);
    timer_start(0);

    t1 = A;

    for(i = 0; i < MK+1; i++){
        t2 = randlc(&t1, t1);
    }

    an = t1;
    tt = S;
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;

    #pragma acc parallel loop present(q[:])
    for(i = 0; i < NQ; i++){
        q[i] = 0.0;
    }

    k_offset = -1;

    #pragma acc parallel loop private(x)
    for(k = 1; k < np; k++){
        kk = k_offset + k;
        t1 = S;
        t2 = an;

        #pragma acc loop seq
        for(i = 1; i <= 100; i++){
            ik = kk / 2;

            if((2 * ik) != kk){
                _t1 = _r23 * t2;
                _a1 = (int) _t1;
                _a2 = t2 - _t23 * _a1;

                _t1 = _r23 * t2;
                _x1 = (int) _t1;
                _x2 = t2 - _t23;
                _t1 = _a1 * _x2 + _a2 * _x1;
                _t2 = (int) (_t23 * _t1);
                _z = _t1 - _t23 * _t2;
                _t3 = _t23 * _z + _a2 * _x2;
                _t4 = (int) (_r46 * _t3);
                t2 = _t3 - _t46 * _t4;
                t3 = _r46 * t2;
            }

            if(ik == 0) break;

            _t1 = _r23 * t2;
            _a1 = (int) _t1;
            _a2 = t2 - _t23 * _a1;

            _t1 = _r23 * t2;
            _x1 = (int) _t1;
            _x2 = t2 - _t23 * _x1;
            _t1 = _a1 * _x2 + _a2 * _x1;
            _t2 = (int) (_r23 * _t1);
            _z = _t1 - _t23 * _t2;
            _t3 = _t23 * _z + _a2 * _x2;
            _t4 = (int) (_r46 * _t3);
            t2 = _t3 - _t46 * _t4;
            t3 = _r46 * t2;

            kk = ik;
        }

        _t1 = _r23 * A;
        _a1 = (int) _t1;
        _a2 = A - _t23 * _a1;

        #pragma acc loop seq
        for(i = 0; i < 2*NK; i++){
            _t1 = _r23 * t1;
            _x1 = (int) _t1;
            _x2 = t1 - _t23 * _x1;
            _t1 = _a1 * _x2 + _a2 * _x1;
            _t2 = (int) (_r23 * _t1);
            _z = _t1 - _t23 * _t2;
            _t3 = _t23 * _z + _a2 * _x2;
            _t4 = (int) (_r46 * _t3);
            t1 = _t3 - _t46 * _t4;
            x[i] = _r46 * t1;
        }

        #pragma acc loop seq
        for(i = 0; i < NK; i++){
            x1 = 2.0 * x[2*1] - 1.0;
            x2 = 2.0 * x[2*i+1] - 1.0;
            t1 = x1 * x1 + x2 * x2;
            if(t1 <= 1.0){
                t2  = sqrt(-2.0 * log(t1) / t1);
                t3  = (x1 * t2);
                t4  = (x2 * t2);
                l   = MAX(fabs(t3), fabs(t4));
                q[l] = q[i] + 1.0;

                sx  += t3;
                sy  += t4;
            }
        }
    }

    for(i = 0; i < NQ; i++){
        gc += q[i];
    }

    timer_stop(0);
    tm = timer_read(0);

    nit = 0;
    verified = true;
    if (M == 24) {
        sx_verify_value = -3.247834652034740e+3;
        sy_verify_value = -6.958407078382297e+3;
    } else if (M == 25) {
        sx_verify_value = -2.863319731645753e+3;
        sy_verify_value = -6.320053679109499e+3;
    } else if (M == 28) {
        sx_verify_value = -4.295875165629892e+3;
        sy_verify_value = -1.580732573678431e+4;
    } else if (M == 30) {
        sx_verify_value =  4.033815542441498e+4;
        sy_verify_value = -2.660669192809235e+4;
    } else if (M == 32) {
        sx_verify_value =  4.764367927995374e+4;
        sy_verify_value = -8.084072988043731e+4;
    } else if (M == 36) {
        sx_verify_value =  1.982481200946593e+5;
        sy_verify_value = -1.020596636361769e+5;
    } else if (M == 40) {
        sx_verify_value = -5.319717441530e+05;
        sy_verify_value = -3.688834557731e+05;
    } else {
        verified = false;
    }

    if (verified) {
        sx_err = fabs((sx - sx_verify_value) / sx_verify_value);
        sy_err = fabs((sy - sy_verify_value) / sy_verify_value);
        verified = ((sx_err <= EPSILON) && (sy_err <= EPSILON));
    }

    Mops = pow(2.0, M+1) / tm / 1000000.0;

    printf("\nEP Benchmark Results:\n\n");
    printf("CPU Time =%10.4lf\n", tm);
    printf("N = 2^%5d\n", M);
    printf("No. Gaussian Pairs = %0.2f\n", gc);
    printf("Sums = %25.15lE %25.15lE\n", sx, sy);
    printf("Counts: \n");
    for (i = 0; i < NQ; i++) {
        printf("%3d%15.0lf\n", i, q[i]);
    }

    print_results("EP", CLASS, M+1, 0, 0, nit,
        tm, Mops, 
        "Random numbers generated",
        verified, NPBVERSION, COMPILETIME, CS1,
        CS2, CS3, CS4, CS5, CS6, CS7);

    if (timers_enabled) {
        if (tm <= 0.0) tm = 1.0;
        tt = timer_read(0);
        printf("\nTotal time:     %9.3lf (%6.2lf)\n", tt, tt*100.0/tm);
        tt = timer_read(1);
        printf("Gaussian pairs: %9.3lf (%6.2lf)\n", tt, tt*100.0/tm);
        tt = timer_read(2);
        printf("Random numbers: %9.3lf (%6.2lf)\n", tt, tt*100.0/tm);
    }

  return 0;
}
