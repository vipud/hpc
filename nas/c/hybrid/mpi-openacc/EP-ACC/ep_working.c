//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB EP code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

//--------------------------------------------------------------------
//      program EMBAR
//--------------------------------------------------------------------
//  This is the serial version of the APP Benchmark 1,
//  the "embarassingly parallel" benchmark.
//
//
//  M is the Log_2 of the number of complex pairs of uniform (0, 1) random
//  numbers.  MK is the Log_2 of the size of each batch of uniform random
//  numbers.  MK can be set for convenience on a given system, since it does
//  not affect the results.
//--------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "type.h"
#include "a_class.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"

#define MAX(X,Y)  (((X) > (Y)) ? (X) : (Y))

#define MK        16
#define MM        (M - MK)
#define NN        (1 << MM)
#define NK        (1 << MK)
#define NQ        10
#define EPSILON   1.0e-8
#define A         1220703125.0
#define S         271828183.0
#define r23       1.1920928955078125e-07
#define r46       r23 * r23
#define t23       8.388608e+06
#define t46       t23 * t23

static double x[2*NK];
static double q[NQ];


int main(int argc, char* argv[])
{
    // Print Values
    printf("\nM   = %d\n", M);
    printf("MK  = %d\n", MK);
    printf("MM  = %d\n", MM);
    printf("NN  = %d\n", NN);
    printf("NK  = %d\n", NK);
    printf("NQ  = %d\n", NQ);

    double Mops, t1, t2, t3, t4, x1, x2;
    double sx, sy, tm, an, tt, gc;
    double sx_verify_value, sy_verify_value, sx_err, sy_err;
    int    np;
    int    i, ik, kk, l, k, nit;
    int    k_offset, j;
    logical verified, timers_enabled;

    double dum[3] = {1.0, 1.0, 1.0};
    char   size[16];

    double _t1, _t2, _t3, _t4, _a1, _a2, _x1, _x2, _z;
    double q0, q1, q2, q3, q4, q5, q6, q7, q8, q9;
    double q_global[NQ];

    FILE *fp;

    // Read from file for additional timer info
    if ((fp = fopen("timer.flag", "r")) == NULL) {
        timers_enabled = false;
    } else {
        timers_enabled = true;
        fclose(fp);
    }

    //--------------------------------------------------------------------
    //  Because the size of the problem is too large to store in a 32-bit
    //  integer for some classes, we put it into a string (for printing).
    //  Have to strip off the decimal point put in there by the floating
    //  point print statement (internal file)
    //--------------------------------------------------------------------
    // Print out the header and the size of the problem
    sprintf(size, "%15.0lf", pow(2.0, M+1));
    j = 14;
    if (size[j] == '.') j--;
    size[j+1] = '\0';
    printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - EP Benchmark\n");
    printf("\n Number of random numbers generated: %15s\n", size);
    verified = false;

    //--------------------------------------------------------------------
    //  Compute the number of "batches" of random number pairs generated
    //  per processor. Adjust if the number of processors does not evenly
    //  divide the total number
    //--------------------------------------------------------------------

    np = NN; // NN is the second largest number
    // In class S it is 255

    //--------------------------------------------------------------------
    //  Call the random number generator functions and initialize
    //  the x-array to reduce the effects of paging on the timings.
    //  Also, call all mathematical functions that are used. Make
    //  sure these initializations cannot be eliminated as dead code.
    //--------------------------------------------------------------------

    // vranlc(0, &dum[0], dum[1], &dum[2]); // This actually does nothing.

    // dum[0] is now a random number between [0, 1]
    // dum[1] is now a seed for random number generation
    dum[0] = randlc(&dum[1], dum[2]);


    // Initialize every value in x to be some really small number
    // Can be parallized with shared memory
    #pragma acc kernels
    {
        // NK is the largest number, 65000 in class S
        // So x is the size of our problem
        for (i = 0; i < 2 * NK; i++) {
            x[i] = -1.0e99; // each value in x starts at -1.0
        }
    } // end acc kernel

    // this line of code also seems redundant...
    // Mops = log(1) basically
    Mops = log(sqrt(fabs(MAX(1.0, 1.0))));

    timer_clear(0);
    timer_clear(1);
    timer_clear(2);
    timer_start(0);

    t1 = A; // A is some generically long double value
    // vranlc(0, &t1, A, x); // Also does nothing

    //--------------------------------------------------------------------
    //  Compute AN = A ^ (2 * NK) (mod 2^46).
    //--------------------------------------------------------------------

    t1 = A; // Again, A is some generically long double value

    // Repeatedly set t2 to be a random value between [0, 1]
    // Also, t1 is constantly set as a new seed
    // This is not parallizable, data dependency
    for (i = 0; i < MK + 1; i++) {
        t2 = randlc(&t1, t1);
    }

    an = t1; // seed value
    tt = S; // S is a generically large number
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;

    // Set every value in q to 0.0
    // Parallizable
    #pragma acc kernels
    {
        // NQ is 10 values, and q is printed at the end
        for (i = 0; i < NQ; i++) {
            q[i] = 0.0;
        }
    }

    //--------------------------------------------------------------------
    //  Each instance of this loop may be performed independently. We compute
    //  the k offsets separately to take into account the fact that some nodes
    //  have more numbers to generate than others
    //--------------------------------------------------------------------

    k_offset = -1; // some offset

    // Loop should be parallizable
    // The only reason it is paralizable is because ti and t2 get set to the same
    // value each iteration
    // And it seems that like every value in this loop will have to be private
//    #pragma acc data copy(q)
    #pragma acc parallel loop gang private(q, x, x1, x2, l, k, kk, ik, t1, t2, t3, t4, _t1, _t2, _t3, _t4, _a1, _a2, _x1, _x2, _z) reduction(+:sx) reduction(+:sy)
    for (k = 1; k <= np; k++) { // np is a large value, 65000 in class s
        // kk will be 1 less than k, k being the loop variable
        kk = k_offset + k;
        t1 = S; // S is some generically large number
        t2 = an; // an at this point is some sort of seed value

        // Find starting seed t1 for this kk.

        #pragma acc loop vector
        for (i = 1; i <= 100; i++) {
            ik = kk / 2; // some interger kk >= 0
            // basically if k is even
            if ((2 * ik) != kk){ // t3 = randlc(&t1, t2); // t3 is a random number
                // t3 = randlc(&t1, t2)
                _t1 = r23 * t2;
                _a1 = (int) _t1;
                _a2 = t2 - t23 * _a1;

                _t1 = r23 * t1;
                _x1 = (int) _t1;
                _x2 = t1 - t23 * _x1;
                _t1 = _a1 * _x2 + _a2 * _x1;
                _t2 = (int) (r23 * _t1);
                _z = _t1 - t23 * _t2;
                _t3 = t23 * _z + _a2 * _x2;
                _t4 = (int) (r46 * _t3);
                t1 = _t3 - t46 * _t4;
                t3 = r46 * t1;
                // t3 = randlc(&t1, t2)
            }
            // t1 is used as the seed, t1 started at S
            if (ik == 0) break; // this should only be true when k = 1 or 2
            // t3 = randlc(&t2, t2)
            _t1 = r23 * t2;
            _a1 = (int) _t1;
            _a2 = t2 - t23 * _a1;

            _t1 = r23 * t2;
            _x1 = (int) _t1;
            _x2 = t2 - t23 * _x1;
            _t1 = _a1 * _x2 + _a2 * _x1;
            _t2 = (int) (r23 * _t1);
            _z = _t1 - t23 * _t2;
            _t3 = t23 * _z + _a2 * _x2;
            _t4 = (int) (r46 * _t3);
            t2 = _t3 - t46 * _t4;
            t3 = r46 * t2;
            // t3 = randlc(&t2, t2)
            // where t2 starts as the old seed
            kk = ik; // kk is now half of what it used to be
        }

        //--------------------------------------------------------------------
        //  Compute uniform pseudorandom numbers.
        //--------------------------------------------------------------------
        //if (timers_enabled) timer_start(2);
        // vranlc(2 * NK, &t1, A, x)
        _t1 = r23 * A;
        _a1 = (int) _t1;
        _a2 = A - t23 * _a1;

        #pragma acc loop vector
        for(i = 0; i < 2*NK; i++){
            _t1 = r23 * t1;
            _x1 = (int) _t1;
            _x2 = t1 - t23 * _x1;
            _t1 = _a1 * _x2 + _a2 * _x1;
            _t2 = (int) (r23 * _t1);
            _z = _t1 - t23 * _t2;
            _t3 = t23 * _z + _a2 * _x2;
            _t4 = (int) (r46 * _t3);
            t1 = _t3 - t46 * _t4;
            x[i] = r46 * t1;
        }
        // vranlc(2 * NK, &t1, A, x)
        //if (timers_enabled) timer_stop(2);

        //--------------------------------------------------------------------
        //  Compute Gaussian deviates by acceptance-rejection method and 
        //  tally counts in concentri//square annuli.  This loop is not 
        //  vectorizable. 
        //--------------------------------------------------------------------
        //if (timers_enabled) timer_start(1);

        for (i = 0; i < NK; i++) {
            x1 = 2.0 * x[2*i] - 1.0;
            x2 = 2.0 * x[2*i+1] - 1.0;
            t1 = x1 * x1 + x2 * x2;
            if (t1 <= 1.0) {
                t2   = sqrt(-2.0 * log(t1) / t1);
                t3   = (x1 * t2);
                t4   = (x2 * t2);
                l    = MAX(fabs(t3), fabs(t4));
                q[l] = q[l] + 1.0;
                // add one to q at wherever the random value falls into
                sx   = sx + t3;
                sy   = sy + t4;
                // sx and sy keep a sum
            }
        }

        #pragma acc loop vector
        for(i = 0; i < NQ; i++){
            q_global[i] += q[i];
            q[i] = 0.0;
        }
/*
        q0 += q[0];
        q1 += q[1];
        q2 += q[2];
        q3 += q[3];
        q4 += q[4];
        q5 += q[5];
        q6 += q[6];
        q7 += q[7];
        q8 += q[8];
        q9 += q[9];
        */
        //if (timers_enabled) timer_stop(1);
    }
    // Parallel using a sum reduction
	#pragma acc kernels
	{
    for (i = 0; i < NQ; i++) {
        gc = gc + q_global[i];
    }
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
    printf("No. Gaussian Pairs = %15.0lf\n", gc);
    printf("Sums = %25.15lE %25.15lE\n", sx, sy);
    printf("Counts: \n");
    for (i = 0; i < NQ; i++) {
        printf("%3d%15.0lf\n", i, q_global[i]);
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
