#include <stdio.h>
#include <stdlib.h>
#include <openacc.h>

int main(int arg, char* argv[])
{
    double a, b, c;
    int i, j, k;
    double x[100000], q[100000];
    double sum;

    #pragma acc parallel loop present(x[:])
    for(i = 0; i < 100000; i++)
    {
        x[i] = 0.0;
    }

    #pragma acc parallel loop present(q[:])
    for(j = 0; j < 100000; j++)
    {
        q[j] = 0.0;
    }

    sum = 0.0;

    #pragma acc parallel loop independent
    for(k = 1; k < 99999; k++)
    {
        a = 1.0 + k;
        b = a + a;
        c = b + b;
        sum += 1.0;
        x[k-1] += 1.0;
        x[k] += 1.0;
        x[k+1] += 1.0;
        q[k-1] += 1.0;
        q[k] += 1.0;
        q[k+1] += 1.0;
    }





}
