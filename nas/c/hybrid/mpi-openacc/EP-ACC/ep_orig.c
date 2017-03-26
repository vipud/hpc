#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "type.h"
#include "a_class.h"
#include "timers.h"
#include "print_results.h"
#include <openacc.h>

#define MAX(X,Y)  (((X) > (Y)) ? (X) : (Y))

int MK;
int MM;
int NN;
double EPSILON;
double A;
double S;
int NK;
int NQ;

int BLKSIZE;

double r23;
double r46;
double t23;
double t46;

inline double randlc_ep( double *x, double a )
{
  double t1, t2, t3, t4, a1, a2, x1, x2, z;
  double r;

  t1 = r23 * a;
  a1 = (int) t1;
  a2 = a - t23 * a1;

  t1 = r23 * (*x);
  x1 = (int) t1;
  x2 = *x - t23 * x1;
  t1 = a1 * x2 + a2 * x1;
  t2 = (int) (r23 * t1);
  z = t1 - t23 * t2;
  t3 = t23 * z + a2 * x2;
  t4 = (int) (r46 * t3);
  *x = t3 - t46 * t4;
  r = r46 * (*x);

  return r;
}

int main()
{
  double Mops, t1, t2, t3, t4, x1, x2;
  double sx, sy, tm, an, tt, gc;
  double sx_verify_value, sy_verify_value, sx_err, sy_err;
  int    np;
  int    i, ik, kk, l, k, nit;
  int    k_offset, j;
  int verified, timers_enabled;
  double q0, q1, q2, q3, q4, q5, q6, q7, q8, q9;

    MK =       16;
    MM =       (M - MK);
    NN =       (1 << MM);
    EPSILON =  1.0e-8;
    A =        1220703125.0;
    S =        271828183.0;
    NK =       1 << MK;
    NQ =       10;

    BLKSIZE = 2048;

    r23 = 1.1920928955078125e-07;
    r46 = r23 * r23;
    t23 = 8.388608e+06;
    t46 = t23 * t23;

  double x[2*(1<<16)];
  double q[10];
  double *xx, *qq;

  double in_t1, in_t2, in_t3, in_t4;
  double in_a1, in_a2, in_x1, in_x2, in_z;

  double tmp_sx, tmp_sy;
  double dum[3] = {1.0, 1.0, 1.0};
  char   size[16];

  int blksize = BLKSIZE;
  int blk, koff, numblks;

  FILE *fp;

  acc_init(acc_device_default);

  if ((fp = fopen("timer.flag", "r")) == NULL) {
    timers_enabled = 0;
  } else {
    timers_enabled = 1;
    fclose(fp);
  }

  if (NN < blksize) {
     blksize = NN;
  }
  numblks = ceil( (double)NN / (double) blksize);

  xx = (double*)malloc(blksize*2*NK*sizeof(double));
  qq = (double*)malloc(blksize*NQ*sizeof(double));

  sprintf(size, "%15.0lf", pow(2.0, M+1));
  j = 14;
  if (size[j] == '.') j--;
  size[j+1] = '\0';
  printf("\n\n NAS Parallel Benchmarks (NPB3.3-ACC-C) - EP Benchmark\n");
  printf("\n Number of random numbers generated: %15s\n", size);

  verified = 0;

  np = NN;


  #pragma acc data create(xx[0:blksize*2*NK],qq[0:blksize*NQ]) copyout(q[0:NQ])
  {
    vranlc(0, &dum[0], dum[1], &dum[2]); // I don't think this function call does anything
    dum[0] = randlc_ep(&dum[1], dum[2]);

    #pragma acc parallel num_gangs((NQ+127)/128) vector_length(128) present(q[0:NQ])
    {
      #pragma acc loop gang vector
      for (i = 0; i < NQ; i++) {
        q[i] = 0.0;
      } // end acc loop
    } // end acc parallel
    Mops = log(sqrt(fabs(MAX(1.0, 1.0))));

    timer_clear(0);
    timer_clear(1);
    timer_clear(2);
    timer_start(0);

    t1 = A;

    for (i = 0; i < MK + 1; i++) {
      t2 = randlc_ep(&t1, t1);
    }

    an = t1;
    tt = S;
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;
    k_offset = -1;

    // Break up the work into blocks so that multiple blocks can be run at once
    for (blk=0; blk < numblks; ++blk) {

      // The actual k_offset with respect to block
      koff = blk*blksize;

      // If there is left over after a block, adjust
      if (koff + blksize > np) {
        blksize = np - (blk*blksize);
      }

      // Parallel block of code, gangs = block size
      #pragma acc parallel num_gangs(blksize) vector_length(128) present(qq[0:blksize*NQ])
      {
        // Each gang will execute this loop
        #pragma acc loop gang
        for(k=0; k<blksize; k++) // Loop through and initialize qq
        {
          #pragma acc loop vector
  	      for(i=0; i<NQ; i++) // Not entirely sure what vector does yet
		        qq[k*NQ + i] = 0.0;
          }
        }

        // Parallel block of code, the gangs are split up more
        #pragma acc parallel num_gangs((blksize+255)/256) num_workers(1) vector_length(256) \
          present(xx[0:blksize*2*NK],qq[0:blksize*NQ])
        {
          // Generate random numbers based on blocksize
          #pragma acc loop gang worker vector reduction(+:sx,sy)
          for (k = 1; k <= blksize; k++) {
            kk = k_offset + k + koff;
            t1 = S;
            t2 = an;

            for (i = 1; i <= 100; i++) {
              ik = kk / 2;
              if ((2 * ik) != kk)
              {
                in_t1 = r23 * t2;
                in_a1 = (int)in_t1;
                in_a2 = t2 - t23 * in_a1;

                in_t1 = r23 * t1;
                in_x1 = (int)in_t1;
                in_x2 = t1 - t23 * in_x1;
                in_t1 = in_a1 * in_x2 + in_a2 * in_x1;
                in_t2 = (int)(r23 * in_t1);
                in_z = in_t1 - t23 * in_t2;
                in_t3 = t23 * in_z + in_a2 * in_x2;
                in_t4 = (int)(r46 * in_t3);
                t1 = in_t3 - t46 * in_t4;
                t3 = r46 * t1;
              }
              if (ik == 0) break;
              in_t1 = r23 * t2;
              in_a1 = (int)in_t1;
              in_a2 = t2 - t23 * in_a1;

              in_t1 = r23 * t2;
              in_x1 = (int)in_t1;
              in_x2 = t2 - t23 * in_x1;
              in_t1 = in_a1 * in_x2 + in_a2 * in_x1;
              in_t2 = (int)(r23 * in_t1);
              in_z = in_t1 - t23 * in_t2;
              in_t3 = t23 * in_z + in_a2 * in_x2;
              in_t4 = (int)(r46 * in_t3);
              t2 = in_t3 - t46 * in_t4;
              t3 = r46 * t2;
              kk = ik;
            }

            in_t1 = r23 * A;
            in_a1 = (int)in_t1;
            in_a2 = A - t23 * in_a1;

            for(i=0; i<2*NK; i++)
            {
		          in_t1 = r23 * t1;
		          in_x1 = (int)in_t1;
		          in_x2 = t1 - t23 * in_x1;
		          in_t1 = in_a1 * in_x2 + in_a2 * in_x1;
		          in_t2 = (int)(r23 * in_t1);
		          in_z = in_t1 - t23 * in_t2;
		          in_t3 = t23*in_z + in_a2 *in_x2;
		          in_t4 = (int)(r46 * in_t3);
		          t1 = in_t3 - t46 * in_t4;
              xx[i*blksize + (k-1)] = r46 * t1;
            }

	          tmp_sx = 0.0;
	          tmp_sy = 0.0;

            for (i = 0; i < NK; i++) {
              x1 = 2.0 * xx[2*i*blksize + (k-1)] - 1.0;
              x2 = 2.0 * xx[(2*i+1)*blksize + (k-1)] - 1.0;
              t1 = x1 * x1 + x2 * x2;
              if (t1 <= 1.0) {
              t2   = sqrt(-2.0 * log(t1) / t1);
              t3   = (x1 * t2);
              t4   = (x2 * t2);
              l    = MAX(fabs(t3), fabs(t4));
              qq[l*blksize + (k-1)] += 1.0;
              tmp_sx   = tmp_sx + t3;
              tmp_sy   = tmp_sy + t4;
            }
          }

          sx += tmp_sx;
          sy += tmp_sy;

        } // end pragma
      } // end block loop

  #pragma acc parallel num_gangs(NQ) num_workers(4) vector_length(32) \
                       present(q[0:NQ],qq[0:blksize*NQ])
  {
	#pragma acc loop gang reduction(+:gc)
	for(i=0; i<NQ; i++)
	{
		double sum_qi = 0.0;
		#pragma acc loop worker vector reduction(+:sum_qi)
		for(k=0; k<blksize; k++)
			sum_qi = sum_qi + qq[i*blksize + k];
		/*sum of each column of qq/q[i] */
		q[i] += sum_qi;
		/*final sum of q*/
		gc += sum_qi;
	}
   }
 }

}/*end acc data*/

  timer_stop(0);
  tm = timer_read(0);

  nit = 0;
  verified = 1;
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
    verified = 0;
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

    // Free up our allocated memory
	free(xx);
	free(qq);

  return 0;
}
