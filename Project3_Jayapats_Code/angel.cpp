#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <cmath>

// system state variables
int NowYear = 2020;          // 2020 - 2025
int NowMonth = 0;         // 0 - 11

float NowPrecip;      // inches of rain per month
float NowTemp;        // temperature this month
float NowHeight = 1.0;      // grain height in inches
int NowNumDeer = 1;       // number of deer in the current population
float GrainQty = 1.0;	// How much grain was grown in the past year
int NumDeerSaved = 0;			// Number of Deers saved by the angel
int NumDeerBef = 1;

// constants
const float GRAIN_GROWS_PER_MONTH = 9.0;
const float ONE_DEER_EATS_PER_MONTH = 1.0;
const float AVG_PRECIP_PER_MONTH = 7.0;    // average
const float AMP_PRECIP_PER_MONTH = 6.0;    // plus or minus
const float RANDOM_PRECIP = 2.0;           // plus or minus noise
const float SQUATCH_CHANCE = 10.0;         
const float AVG_TEMP = 60.0;               // average
const float AMP_TEMP = 20.0;               // plus or minus
const float RANDOM_TEMP = 10.0;            // plus or minus noise
const float MIDTEMP = 40.0;
const float MIDPRECIP = 10.0;

float Ranf(unsigned int* seedp, float low, float high) {
	float r = (float)rand_r(seedp);
	return(low + r * (high - low) / (float)RAND_MAX);
}

int Ranf(unsigned int* seedp, int ilow, int ihigh) {
	float low = (float)ilow;
	float high = (float)ihigh + 0.9999f;
	return (int)(Ranf(seedp, low, high));
}

float sqr(float n) {
	return n * n;
}

void GrainDeer() {
	while (NowYear < 2026) {
		NumDeerBef = NowNumDeer;
		if (NowNumDeer < NowHeight) {
			NowNumDeer++;
		}
		else if (NowNumDeer > NowHeight) {
			NowNumDeer--;
			if (NowNumDeer < 0) { NowNumDeer = 0; }
		}

#pragma omp barrier
		//NowNumDeer = NowNumDeer - NumDeerSaved;
		NowNumDeer = NowNumDeer + NumDeerSaved;
		if (NowNumDeer < 0) { NowNumDeer = 0; }


#pragma omp barrier

#pragma omp barrier

	}
}

void Grain() {
	while (NowYear < 2026) {
		// Using a temporary local var so the grain or the deer don't get confused
		int tempHeight = NowHeight; // "Temp" as in "temporary" not "temperature"
		float tempFactor = exp(-sqr((NowTemp - MIDTEMP) / 10.0));
		float precipFactor = exp(-sqr((NowPrecip - MIDPRECIP) / 10.0));

		tempHeight += tempFactor * precipFactor * GRAIN_GROWS_PER_MONTH;
		tempHeight -= (float)NowNumDeer * ONE_DEER_EATS_PER_MONTH;
		if (tempHeight < 0.0) { tempHeight = 0.0; }
		if (NowMonth == 0) {
			GrainQty = tempHeight;
		}
		else {
			GrainQty += tempHeight;
		}

#pragma omp barrier

		NowHeight = tempHeight;

#pragma omp barrier

#pragma omp barrier

	}
}

void Watcher() {
	unsigned int seed = time(NULL);
	printf("year,month,preciptation,temp,grain height,deer,Deers saved by Angel,Grain Quantity\n");
	while (NowYear < 2026) {

#pragma omp barrier

#pragma omp barrier
		// Everything's done running by now
		printf("%d,%d,%f,%f,%f,%d,%d,%f\n", NowYear, NowMonth, NowPrecip, (5.0 / 9.0) * (NowTemp - 32.0), (NowHeight * 2.54), NowNumDeer, NumDeerSaved, (GrainQty * 2.54));

		NowMonth++;
		if (NowMonth == 12) {
			NowYear += 1;
			NowMonth = 0;
		}

		// Start prepping for the next cycle
		if (NowMonth != 11) {
			NumDeerSaved = 0;
		}
		float ang = (30.0 * (float)NowMonth + 15.0) * (M_PI / 180.0);

		float temp = AVG_TEMP - AMP_TEMP * cos(ang);
		NowTemp = temp + Ranf(&seed, -RANDOM_TEMP, RANDOM_TEMP);

		float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin(ang);
		NowPrecip = precip + Ranf(&seed, -RANDOM_PRECIP, RANDOM_PRECIP);
		if (NowPrecip < 0.0) {
			NowPrecip = 0.0;
		}
#pragma omp barrier

	}
}

void Angel() {
	unsigned int seed = time(NULL);
	while (NowYear < 2026) {
		
		//Angels are sent to planet earth to save the Deer from getting extinct.
		//So at the end of each year if the grain quantity goes low, angels appear on earth and save the deers

		if (NowMonth > 10) {
			if (GrainQty > 80.0) {
				NumDeerSaved = 0;
			}
			if ((GrainQty > 40.0) && (GrainQty < 50.0)) {
				NumDeerSaved = Ranf(&seed, 0, 2);
			}
			if ((GrainQty > 30.0) && (GrainQty < 40.0)) {
				NumDeerSaved = Ranf(&seed, 1, 3);
			}
			if (GrainQty < 40.0) {
				NumDeerSaved = Ranf(&seed, 0, 4);
			}

		}
#pragma omp barrier

#pragma omp barrier

#pragma omp barrier

	}
}

int main(int argc, char* argv[]) {
	// This stuff is also repeated every month
	float ang = (30.0 * (float)NowMonth + 15.0) * (M_PI / 180.0);

	float temp = AVG_TEMP - AMP_TEMP * cos(ang);

	unsigned int seed = time(NULL);
	NowTemp = temp + Ranf(&seed, -RANDOM_TEMP, RANDOM_TEMP);

	float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin(ang);
	NowPrecip = precip + Ranf(&seed, -RANDOM_PRECIP, RANDOM_PRECIP);
	if (NowPrecip < 0.0) {
		NowPrecip = 0.0;
	}

	omp_set_num_threads(4);	// same as # of sections
#pragma omp parallel sections
	{
#pragma omp section
		{
			GrainDeer();
		}

#pragma omp section
		{
			Grain();
		}

#pragma omp section
		{
			Watcher();
		}

#pragma omp section
		{
			Angel();
		}
	}       // implied barrier -- all functions must return in order
			// to allow any of them to get past here

	return 0;
}
