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
#include "npbparams.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"
#include <mpi.h>
#include <omp.h>

#define MAX(X,Y)  (((X) > (Y)) ? (X) : (Y))

#define MK        16
#define MM        (M - MK)
#define NN        (1 << MM)
#define NK        (1 << MK)
#define NQ        10
#define EPSILON   1.0e-8
#define A         1220703125.0
#define S         271828183.0

static double x[2*NK];
static double qq[NQ];
static double global_q[NQ];
#pragma omp threadprivate(x,qq)
static double q[NQ];


int main(int argc, char* argv[])
{
  double Mops, t1, t2, t3, t4, x1, x2;
  double sx, sy, tm, an, tt, gc;
  double sx_verify_value, sy_verify_value, sx_err, sy_err;
  int    np;
  int    i, ik, kk, l, k, nit;
  int    k_offset, j;
  logical verified, timers_enabled;

  double dum[3] = {1.0, 1.0, 1.0};
  char   size[16], timers_enabled_char;

  FILE *fp;
  int required = MPI_THREAD_SERIALIZED;
  int provided;
  int rank, npes;
  MPI_Init_thread(&argc, &argv, required, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &npes);

  if(rank == 0){
    if ((fp = fopen("timer.flag", "r")) == NULL) {
        timers_enabled = false;
        timers_enabled_char = 0;
    } else {
        timers_enabled = true;
        fclose(fp);
        timers_enabled_char = 1;
    }
  }

  MPI_Bcast(&timers_enabled_char, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
  timers_enabled = timers_enabled ? true : false;

  //--------------------------------------------------------------------
  //  Because the size of the problem is too large to store in a 32-bit
  //  integer for some classes, we put it into a string (for printing).
  //  Have to strip off the decimal point put in there by the floating
  //  point print statement (internal file)
  //--------------------------------------------------------------------
  if(rank == 0){
    sprintf(size, "%15.0lf", pow(2.0, M+1));
    j = 14;
    if (size[j] == '.') j--;
    size[j+1] = '\0';
    printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - EP Benchmark\n");
    printf("\n Number of random numbers generated: %15s\n", size);
  }
  verified = false;

  //--------------------------------------------------------------------
  //  Compute the number of "batches" of random number pairs generated
  //  per processor. Adjust if the number of processors does not evenly
  //  divide the total number
  //--------------------------------------------------------------------

  np = NN; // np is just some integer value

  //--------------------------------------------------------------------
  //  Call the random number generator functions and initialize
  //  the x-array to reduce the effects of paging on the timings.
  //  Also, call all mathematical functions that are used. Make
  //  sure these initializations cannot be eliminated as dead code.
  //--------------------------------------------------------------------

  vranlc(0, &dum[0], dum[1], &dum[2]); // This actually does nothing.

  // dum[0] is now a random number between [0, 1]
  // dum[1] is now a seed for random number generation

  if(rank == 0)
    dum[0] = randlc(&dum[1], dum[2]);

  MPI_Bcast(dum, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Initialize every value in x to be some really small number
  // Can be parallized with shared memory
  #pragma omp parallel default(shared) private(i)
  {
      for (i = 0; i < 2 * NK; i++) {
          x[i] = -1.0e99;
      }
  }

  // this line of code also seems redundant...
  // Mops = log(1) basically
  Mops = log(sqrt(fabs(MAX(1.0, 1.0))));

  timer_clear(0);
  timer_clear(1);
  timer_clear(2);
  timer_start(0);

  t1 = A; // A is some generically long double value
  vranlc(0, &t1, A, x); // Also does nothing

  //--------------------------------------------------------------------
  //  Compute AN = A ^ (2 * NK) (mod 2^46).
  //--------------------------------------------------------------------

  t1 = A; // Again, A is some generically long double value

  // Repeatedly set t2 to be a random value between [0, 1]
  // Also, t1 is constantly set as a new seed
  // This is not parallizable, data dependency
  if(rank == 0){
    for (i = 0; i < MK + 1; i++) {
        t2 = randlc(&t1, t1);
    }
  }
  MPI_Bcast(&t1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  an = t1; // seed value
  tt = S; // S is a generically large number
  gc = 0.0;
  sx = 0.0;
  sy = 0.0;

  // Set every value in q to 0.0
  // Parallizable

  for (i = 0; i < NQ; i++) {
    q[i] = 0.0;
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

#pragma omp parallel default(shared) private(k,kk,t1,t2,t3,t4,i,ik,x1,x2,l)
{
    #pragma omp parallel default(shared) private(i)
    {
        for (i=0; i<NQ;i++){
            qq[i]=0.0;
        }
    }
    #pragma omp for reduction(+:sx,sy) nowait
    for (k = rank+1; k <= np; k+=npes) {
        // kk will be 1 less than k, k being the loop variable
        kk = k_offset + k;
        t1 = S; // S is some generically large number
        t2 = an; // an at this point is some sort of seed value

        // Find starting seed t1 for this kk.

        for (i = 1; i <= 100; i++) {
            ik = kk / 2; // some interger kk >= 0
            // basically if k is even
            if ((2 * ik) != kk) t3 = randlc(&t1, t2); // t3 is a random number
            // t1 is used as the seed, t1 started at S
            if (ik == 0) break; // this should only be true when k = 1 or 2
            t3 = randlc(&t2, t2); // t3 is a random number using t2 as seed
            // where t2 starts as the old seed
            kk = ik; // kk is now half of what it used to be
        }

    //--------------------------------------------------------------------
    //  Compute uniform pseudorandom numbers.
    //--------------------------------------------------------------------
        if (timers_enabled) timer_start(2);
        vranlc(2 * NK, &t1, A, x); // fill x with a bunch of random numbers
        if (timers_enabled) timer_stop(2);

    //--------------------------------------------------------------------
    //  Compute Gaussian deviates by acceptance-rejection method and
    //  tally counts in concentri//square annuli.  This loop is not
    //  vectorizable.
    //--------------------------------------------------------------------
        if (timers_enabled) timer_start(1);

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

        if (timers_enabled) timer_stop(1);
    }
    for(i=0;i<NQ; i++){
        #pragma omp atomic
        q[i]+=qq[i];
    }
}

  // Parallel using a sum reduction

  for (i = 0; i < NQ; i++) {
    gc = gc + q[i];
  }
  double gc_global;
  double sx_global;
  double sy_global;
  MPI_Reduce(&gc, &gc_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sx, &sx_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&sy, &sy_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(q, global_q, 10, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  gc = gc_global;
  sx = sx_global;
  sy = sy_global;

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

  if(rank == 0){
    printf("\nEP Benchmark Results:\n\n");
    printf("CPU Time =%10.4lf\n", tm);
    printf("N = 2^%5d\n", M);
    printf("No. Gaussian Pairs = %15.0lf\n", gc);
    printf("Sums = %25.15lE %25.15lE\n", sx, sy);
    printf("Counts: \n");
    for (i = 0; i < NQ; i++) {
        printf("%3d%15.0lf\n", i, global_q[i]);
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
  }

  MPI_Finalize();

  return 0;
}
