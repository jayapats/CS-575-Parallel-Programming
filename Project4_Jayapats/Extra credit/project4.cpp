#define _USE_MATH_DEFINES
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <xmmintrin.h>
#include <time.h>
#include <stdlib.h>

#define SSE_WIDTH		4

#ifndef NUMT
#define NUMT	4
#endif

//int NUMT=4;

//#define SIZE       	(8*1024*1024)
#ifndef ArraySize
#define ArraySize	2048
#endif

//#define SIZE       	2097152
#define NUMTRIES     30
#define NUM_ELEMENTS_PER_CORE ArraySize/NUMT

float SimdMulSum(float*, float*, int);
float SIMDMul(float* a, float* b, int len);

float A[ArraySize];
float B[ArraySize];
float C[ArraySize];

int
main(int argc, char* argv[])
{
#ifndef _OPENMP
    fprintf(stderr, "OpenMP is not supported here -- sorry.\n");
    return 1;
#endif
    /*if (argc >= 2)
    {
        NUMT = atoi(argv[1]);
    }

    if (argc >= 3)

    {
        DSIZE = atoi(argv[2]);
    }*/

    // inialize the arrays:
    for (int i = 0; i < ArraySize; i++)
    {
        A[i] = 1.;
        B[i] = 2.;
    }

    omp_set_num_threads(NUMT);
    //fprintf(stderr, "Using %d threads\n", NUMT);

    double maxMegaMults1 = 0.;
    double maxMegaMults2 = 0.;

    for (int t = 0; t < NUMTRIES; t++)
    {
        double time0 = omp_get_wtime();

#pragma omp parallel 
        {
            int first = omp_get_thread_num() * NUM_ELEMENTS_PER_CORE;
            SimdMulSum(&A[first], &B[first], NUM_ELEMENTS_PER_CORE);
        }
        double time1 = omp_get_wtime();
        double megaMults1 = (double)ArraySize / (time1 - time0) / 1000000.;

        if (megaMults1 > maxMegaMults1)
            maxMegaMults1 = megaMults1;

        //printf("Peak Performance = %8.2lf MegaMults/Sec\n", maxMegaMults1);
        //printf("Average Performance = %8.2lf MFLOPS\n", summflops / (double)NUMTRIES);
    }

    for (int t = 0; t < NUMTRIES; t++)
    {
        double time01 = omp_get_wtime();

#pragma omp parallel
        /*for (int i = 0; i < ArraySize; i++)
        {
            C[i] = A[i] * B[i];
        }*/
        

        SIMDMul(A, B, ArraySize);

        double time11 = omp_get_wtime();
        double megaMults2 = (double)ArraySize / (time11 - time01) / 1000000.;

        if (megaMults2 > maxMegaMults2)
            maxMegaMults2 = megaMults2;
    }


    //printf("Peak Performance = %8.2lf MegaMults/Sec\n", maxMegaMults2);
    //printf("Average Performance = %8.2lf MFLOPS\n", summflops / (double)NUMTRIES);

    //double sp = maxMegaMults1 / maxMegaMults2;
    double sp =  maxMegaMults1/maxMegaMults2;
    printf("SSE Performance=%8.2lf \t Non-SSE Performance=%8.2lf \t NUMT= \t Array Size =%d \t SpeedUp =%8.2lf ", maxMegaMults1, maxMegaMults2,NUMT, ArraySize, sp);

}
float
SimdMulSum(float* a, float* b, int len)
{
    float sum[4] = { 0., 0., 0., 0. };
    int limit = (len / SSE_WIDTH) * SSE_WIDTH;
    register float* pa = a;
    register float* pb = b;

    __m128 ss = _mm_loadu_ps(&sum[0]);
    for (int i = 0; i < limit; i += SSE_WIDTH)
    {
        ss = _mm_add_ps(ss, _mm_mul_ps(_mm_loadu_ps(pa), _mm_loadu_ps(pb)));
        pa += SSE_WIDTH;
        pb += SSE_WIDTH;
    }
    _mm_storeu_ps(&sum[0], ss);

    for (int i = limit; i < len; i++)
    {
        sum[0] += a[i] * b[i];
    }

    return sum[0] + sum[1] + sum[2] + sum[3];
}


float SIMDMul(float* a, float* b, int len) {

    float sum = 0.;

    for (int i = 0; i < len; i++) {

        sum += a[i] * b[i];
    }

    return sum;
}
