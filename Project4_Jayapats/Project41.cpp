#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <xmmintrin.h>
#define SSE_WIDTH		4

#ifndef NUMT
#define NUMT	1
#endif

#ifndef ArraySize
#define ArraySize	2048
#endif

#define NUMTRIES     30

float SimdMulSum(float* , float* , int);
float SIMDMul(float* a, float* b, int len);

float A[ArraySize];
float B[ArraySize];
float C[ArraySize];

int
main()
{
#ifndef _OPENMP
    fprintf(stderr, "OpenMP is not supported here -- sorry.\n");
    return 1;
#endif

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
    //double summflops = 0.;

    for (int t = 0; t < NUMTRIES; t++)
    {
        double time0 = omp_get_wtime();
        SimdMulSum(A, B, ArraySize);
        double time1 = omp_get_wtime();
        double megaMults1 = (double)ArraySize / (time1 - time0) / 1000000.;

        if (megaMults1 > maxMegaMults1)
            maxMegaMults1 = megaMults1;

    }

    for (int t = 0; t < NUMTRIES; t++)
    {
        double time01 = omp_get_wtime();
        
        #pragma omp parallel
              
                SIMDMul(A, B, ArraySize);
                double time11 = omp_get_wtime();
                double megaMults2 = (double)ArraySize / (time11 - time01) / 1000000.;
                //summflops += megaMults;

                if (megaMults2 > maxMegaMults2)
                    maxMegaMults2 = megaMults2;
    }     

    double speedUp = maxMegaMults1/maxMegaMults2;
    //printf("SSE Performance = %8.2lf  Non-SSE Performance = %8.2lf Array Size =%d SpeedUp = %8.2lf \n", maxMegaMults1, maxMegaMults2, ArraySize,sp);
	printf("SSE=%8.2lf Non SSE=%8.2lf ArraySize=%d SpeedUp=%8.2lf\n",maxMegaMults1,maxMegaMults2,ArraySize,speedUp);
    
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

   