/*
############## POISSON SOURCE #########################################
# Implements a spike event source timing trigger, where on each simulation step 
# the source returns if a spike occured stochastically according to a poisson distribution
# parametrised by lamda (�) the rate of spike occurance in Hz (no of spikes per second).
# The probability of each event is exponentially distributed in time.
# Kostantinos Lagogiannis 07/2007
########################################################################
*/

#include "../common.h"
#include "../util.h"
#include "../../synapseModels/common.h" //For the Global GSL Instance
#include "PoissonSource.h"
/*
 * lamda is the rate of the occurance of events per second.
 * timestep is the simulation timestep used - This will be used to obtain a poisson distibuted time of next event by counting the number of timesteps in lamda
 */
PoissonSource::PoissonSource(double lamda,double timeStep,double noiseStdev)
{
	
	mlamda	= (lamda < 0.0)?0.0:lamda;
	h		= timeStep;
	sigma	= noiseStdev;
	mlamdaInTs = mlamda/h; ////The rate of Events in number of timesteps - The number of timesteps in the Event rate

	rng_r = g_getRandGeneratorInstance(false);
	if (!rng_r) throw "GSL RNG Init Failed. Out Of Memory?";

	 //unsigned int seed = unsigned(time(&t)) + rand()*100;
	 // gsl_rng_set(rng_r,seed);
	 //Seed random number generator
	 //srand(seed);
}

bool PoissonSource::drawSpikeEvent()
{
	float noise =0; // PoissonSource::randGauss(0,4.0f*sigma,sigma,2.0f*sigma);

	//poisson spike draw
	//double r = (rand()/(double)RAND_MAX);
	double r = gsl_rng_uniform(rng_r);
	if ((mlamda*h + noise)< r) return false;
	else		return true;
}

double PoissonSource::getRate()
{
	return mlamda;
}

//gaussian distribution
double PoissonSource::randGauss( double min, double max, double sigma, double centre)
{
/*double random = (min + (max-min) * (double)rand()/RAND_MAX); //create random domain between [min,max]

double tmp = (random-centre)/sigma; 
double gauss = exp(-tmp*tmp/2); //gaussian formula
*/
//Use Chi Square distribution with 2 degrees of freedom
double r1 = (double)gsl_rng_uniform_pos (rng_r);
double r2 = (double)gsl_rng_uniform_pos (rng_r);

//Note Log =ln
double gauss = sqrt(-2.0f*log(r1))*cos(2*M_PI*r2);

return gauss*sigma;
}

//Returns the time of the next spike under an exponential distribution
double PoissonSource::getTimeUntilNextEvent ()
{
	//Scnd param should be 1/Encoding rate
	return gsl_ran_gamma(rng_r,1.0,1.0);  //The return value is the time until next event occurs

	//return gsl_ran_gamma(rng_r,1.0,1.0/dEncodingRate);  //The return value is the time until next event occurs

}

unsigned int PoissonSource::getTimestepsUntilNextEvent()
{
	//return gsl_ran_exponential(rng_r,mlamdaInTs); // This Does not Work Properly The Variance is not Big enough!!
	return gsl_ran_gamma(rng_r,1.0,1.0)/h; //Last Param 1/dEncodingRate assumed to be 1


	//->Draw a random from Poisson With a mean given as the number of timesteps until the next Event
	//  double mu = (double)(1.0/(mlamda*h)); //Probability of spike event
	//  double u = gsl_rng_uniform_pos (rng_r); //This returns a non-zero probability less that 1
	//  return -mu * log (u);
}
///Previous method:
double PoissonSource::getDelayUntilNextEvent()
{
  double mu = (double)(1.0/(mlamda)); //Probability of spike event
  double u = gsl_rng_uniform_pos (rng_r); //This returns a non-zero probability less that 1
  return -mu * log (u);
 }

PoissonSource::~PoissonSource(void)
{
}
