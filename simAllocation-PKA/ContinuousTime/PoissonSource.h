#pragma once
#include "../common.h"

class PoissonSource
{
public:
	PoissonSource();
	PoissonSource(double lamda,double timeStep,double noiseStdev,gsl_rng* prng);
	bool drawSpikeEvent(void);
	double randGauss( double min, double max, double sigma, double centre);
	double getTimeUntilNextEvent ();
	double getDelayUntilNextEvent();
	unsigned int getTimestepsUntilNextEvent();
	double getRate();
	~PoissonSource(void);

private:
	double mlamda; //Rate of event Arrival given in number of timesteps - So realtime is mlamda*h
	unsigned int mlamdaInTs; //The rate of Events in number of timesteps
	double h; //Simulation Timestep
	double sigma; //Gaussian Noise StdDev in Sec
	gsl_rng * rng_r; //Used by GSL Rand Num Generator

	time_t t; //Used for random num generation

};
