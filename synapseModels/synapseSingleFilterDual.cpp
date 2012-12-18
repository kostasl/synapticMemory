/*
 * synapseFilterUnified.cpp
 *
 * Threshold LookUp tables define a set for {q,p} thresholds for transitions
 *	The threshold used for p or q depends on the strength state of the synapse
 *
 *  Created on: 27 Mar 2011
 *      Author: kostasl
 *
 *These filters are a stochastic updaters equivalent. They do not grow or follow a cascade
 *They always remain in the same cascade state and do not use a reflecting boundary formulation
 *But  once a threshold is reached -> re-injection to zero state.
 */
#include "common.h"
#include "synapseSingleFilterDual.h"

//Default ConsTructor - Never Actually Used
synapseSingleFilterDual::synapseSingleFilterDual() {

	liberrexit(500,"Constructor Not Implemented");
}

//Init With StartStrength Used by SynapticConnection
synapseSingleFilterDual::synapseSingleFilterDual(int piCascadeSize,float startStrength,gsl_rng * rng_r)
{
	liberrexit(500,"Constructor Not Implemented");
}

synapseSingleFilterDual::synapseSingleFilterDual(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r)
{
	liberrexit(500,"Constructor Not Implemented");

}

//Initialise with cascadeSize
synapseSingleFilterDual::synapseSingleFilterDual(int piCascadeSize,gsl_rng * rng_r)
{
	liberrexit(500,"Constructor Not Implemented");

}

//Both Thresholds Are Positive
synapseSingleFilterDual::synapseSingleFilterDual(int piLThres,int piHThres,double pdDecay)
{

	if (!mprng)
		mprng = g_getRandGeneratorInstance(false);

	miCascadeIndex = miStartIndex = 0;
	miTerminalIndex = 0;// miCascadeIndex + 1;

	dDDecayRate = dPDecayRate = pdDecay;


	miHThres = abs(piHThres);
	miLThres = abs(piLThres);
	//miRPFilterValue = gsl_ran_binomial (mprng,0.5,miHThres-1);
	//miRDFilterValue = gsl_ran_binomial (mprng,0.5,miLThres-1);


	double r = gsl_rng_uniform(mprng);
	if (r<0.5)
		penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
	else
		penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;

	mbIsMonitored = false;
	mbNoPlasticity = false;

	initialiseFilterState();
	assert(miRPFilterValue < miHThres && miRDFilterValue < miLThres);

}

//Used to Start From A fixed point and test how distribution evolves
synapseSingleFilterDual::synapseSingleFilterDual(int piCascadeSize,int piStartIndex,gsl_rng * rng_r)
{
	mprng = rng_r;
	piCascadeSize = 1; //Fixed
	init( piCascadeSize,piStartIndex,SYN_STRENGTH_NOTSET,0);
}

synapseSingleFilterDual::synapseSingleFilterDual(int piCascadeSize,int piStartIndex,
		ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
		gsl_rng *  rng_r,
		int iRateDependentParameterSet) //Default Value is rate=1.0

{
	piCascadeSize = 1; //Fixed
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,penumStartStrength,iRateDependentParameterSet);
}


synapseSingleFilterDual::synapseSingleFilterDual(int piCascadeSize,int piStartIndex,
					ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
					int iStartFilterState,
					gsl_rng *  rng_r,
					int iRateDependentParameterSet) //Default Value is rate=1.0
{
	mprng = rng_r;
	piCascadeSize = 1;//Fixed
	init( piCascadeSize,piStartIndex,penumStartStrength,iRateDependentParameterSet);
	miRPFilterValue= miRDFilterValue = iStartFilterState; //Some simple init
	assert (miRPFilterValue >= miLThres && miRPFilterValue <= miHThres);
}


void synapseSingleFilterDual::init(int piCascadeSize,int iThresholdsIndex,
								SYN_STRENGTH_STATE enumCurrStrengthState,int ParameterSetIndex)
{
	miParameterSetUsed = ParameterSetIndex;
	miCascadeIndex = miStartIndex = iThresholdsIndex;
	miTerminalIndex = miCascadeIndex + 1; //Make it look like a non terminal State
	//miTerminalIndex = piCascadeSize-1;
	mbNoPlasticity = false;
	mbIsMonitored = false;

	//g_rng_r = getRandGeneratorInstance();
	if (!mprng)
		mprng = g_getRandGeneratorInstance(false);

	if (enumCurrStrengthState == SYN_STRENGTH_NOTSET)
	{
		double r = gsl_rng_uniform(mprng);
		if (r < 0.5)
			penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
		else
			penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;
	}else
		penumStartStrength = penumStrength = (uint8_t)enumCurrStrengthState;

	setFilterThresholds();
	initialiseFilterState();
}

int synapseSingleFilterDual::getCascadeIndex() const
{

	return miCascadeIndex;
}

int synapseSingleFilterDual::getCascadeSize() const
{
	return (miTerminalIndex+1);
}

bool synapseSingleFilterDual::hasIndexModified() const
{
	//return mbCascadeIndexChanged;
	return (miCascadeIndex	!= miStartIndex);
}


bool synapseSingleFilterDual::hasStrengthModified() const
{
	return (penumStrength != penumStartStrength);
}

//Induces an LTP event depending on current state either switches cascades or moves down
int8_t synapseSingleFilterDual::handlePOT()
{
	return addSample(+1);
}

//Induces an LTD event depending on current state either switches cascades or moves down
int8_t synapseSingleFilterDual::handleDEP()
{
	return addSample(-1);
}

//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
void synapseSingleFilterDual::handleNOP()
{
	addSample(0);
}


void synapseSingleFilterDual::doStochasticDecay()
{
	//Do Decay on both running sums
	unsigned int rDecaySteps;
	//double dPDecayRate,dDDecayRate;

	double pH;
	if (miRDFilterValue > 0) //Not required to CHeck but saves time
	{
		pH = 1.0-exp(-dDDecayRate*miTimeSinceLastInduction);
		rDecaySteps = gsl_ran_binomial(mprng,pH,miRDFilterValue);
		miRDFilterValue = miRDFilterValue -  rDecaySteps;
	}
	//Now Do the decay on each filter
	if (miRPFilterValue > 0)
	{
		pH = 1.0-exp(-dPDecayRate*miTimeSinceLastInduction); //Probability of Particle Decay
		rDecaySteps = gsl_ran_binomial(mprng,pH,miRPFilterValue);
		miRPFilterValue = miRPFilterValue -  rDecaySteps;
	}

}

///@Returns -1 For a Plastic Transition, +1 For Metaplastic and 0 For No Transition
int synapseSingleFilterDual::addSample(int iValue)
{
	int iRet = 0;
	miTimeSinceLastInduction++; //Increment time since last induction

	if (iValue != 0) // WHY iValue != 0 ??
	{//First draw the number of decay steps that have occurred up to current time t
		//Do Decay on both running sums
		doStochasticDecay();

		assert (miRDFilterValue >= 0 && miRPFilterValue >= 0); //Check no -Ve values emerge
		miTimeSinceLastInduction = 0; //Reset Counter of No-induction timesteps -
	}

	if (iValue < 0) //-ve value added so add to low Running sum miRValue_1
	{	miRDFilterValue -= iValue; //iVal is +Ve even for low threshold!
		if (miRDFilterValue >= miLThres ) 	//Check the modified sum if it has exceeded threshold
			iRet = lThresReached();
	}

	if (iValue > 0)
	{
		miRPFilterValue += iValue; //iVal is +Ve even for low threshold!
		if (miRPFilterValue >= miHThres ) 	//Check the modified sum if it has exceeded threshold
			iRet = hThresReached(); //Return -1 if the cascade is to increase
	}
	assert(miRPFilterValue < miHThres && miRDFilterValue < miLThres);
	return iRet;
}

int synapseSingleFilterDual::lThresReached()
{
	int Ret = 0;
	//Add new sample to the correct internal sum
	if (penumStrength == SYN_STRENGTH_STRONG)
	{ //Then lThres is a q transition
		switchReset();
		Ret = -1;
	}
	else
	{
		Ret = 0;//No Transition
		reInjectFilterStateToCascadeState(); //Reset running Sum
		miStartIndex = miCascadeIndex+1;//Makes Index Changed Flag true
		uiSameThresholdTransitionCounter++;// Counter - Used for Statistics Or Allocation
	}

	return Ret;
}

int synapseSingleFilterDual::hThresReached()
{
	int Ret = 0;
	//Add new sample to the correct internal sum
	if (penumStrength == SYN_STRENGTH_WEAK)
	{ //Then lThres is a q transition
		switchReset();
		Ret = -1;
	}else
	{
		Ret = 0;//No Transition
		reInjectFilterStateToCascadeState(); //Reset running Sum - No reflecting Boundary
		miStartIndex = miCascadeIndex+1; //Makes Index Changed Flag true
		uiSameThresholdTransitionCounter++; //Metaplasticity Counter - Used for Statistics Or Allocation
	}


	return Ret;
}

int synapseSingleFilterDual::getRunningValueH() const
{
	return miRPFilterValue;
}


int synapseSingleFilterDual::getRunningValueL() const
{
	return miRDFilterValue;
}

int synapseSingleFilterDual::getStrength() const
{
	return (int)penumStrength;
}

int synapseSingleFilterDual:: getStartStrength() const
{
return (int)penumStartStrength;
}

int synapseSingleFilterDual::getLThres() const
{
	return miLThres;
}
int synapseSingleFilterDual::getHThres() const
{
	return miHThres;
}

double synapseSingleFilterDual::getHDecay() const
{
	return 0.0;
}
double synapseSingleFilterDual::getLDecay() const{
	return 0.0;
}

double synapseSingleFilterDual::getDecay() const{

	return dDDecayRate; //(*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
}

//Returns the number of times the Filter Is allowed to Grow Threshold (The Max Threshold Index)
int synapseSingleFilterDual::getMaxSteps() const
{
	return miTerminalIndex;
}

int synapseSingleFilterDual::getThresholdIndex() const
{
  return miCascadeIndex;
}

//Save Current State As Default
void synapseSingleFilterDual::startMonitoring()
{
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	miStartIndex = miCascadeIndex;
	penumStartStrength = penumStrength;
	mbIsMonitored = true;
}

void synapseSingleFilterDual::stopMonitoring()
{
	mbIsMonitored = false;
}

bool synapseSingleFilterDual::isMonitored() const
{
	return mbIsMonitored;
}

//Method called by testEscapeTime
void synapseSingleFilterDual::reset()
{
	ICascadeSynapse::reset(); //Randomizes Start Strength, Unfreezeplasticity

	//penumStrength = (SYN_STRENGTH_STATE)penumStartStrength;
	//miStartIndex = miCascadeIndex; //Index is frozen in this CLass - So this is used to detect changes in P Transitions
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	mbNoPlasticity = false;
	mbIsMonitored = false;

	//setFilterThresholds();
	//reInjectFilterStateToCascadeState();
	//initialiseFilterState(); //Like Starting Over - BUT If we Re-init from Non zero position Then The escape time through any boundary Reduces
	 //The reset is only called by the testEscapetime Function - And thus it affects the result of this function
}

void synapseSingleFilterDual::switchReset()
{
	if (penumStrength == SYN_STRENGTH_STRONG){
		penumStrength = SYN_STRENGTH_WEAK;
	}else{
		penumStrength = SYN_STRENGTH_STRONG;
		}


		//setFilterThresholds();
		reInjectFilterStateToCascadeState(); //Reset To Zero On New State
}


//Use Index / strength and the Lookup arrays to set new L/H thresholds
//Called after changes to cascade index
void synapseSingleFilterDual::setFilterThresholds()
{
	//Since we are using only terminals and n_nonterm = (n+1)-term then use the index at +1
	if (penumStrength == SYN_STRENGTH_STRONG)
	{
		miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0]; //P
		miLThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
		//WEAK SYNAPSE
	}
else
	{
		miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
		miLThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
	} //ENDOF IF STRONG

	//SET DECAY
	//Add 1 to index as Terminal States match againste n-1
	if (penumStrength == SYN_STRENGTH_STRONG){
			dPDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0]; //P is 1st values
			dDDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][1]; //Q is 2nd Values
		}else {
			dPDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][1];
			dDDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
		}

	//throw "Set Threshold Called";
}

//Sets where the running sum should be initialized after every change to Cascade state -
void synapseSingleFilterDual::reInjectFilterStateToCascadeState()
{
	//Inject Back to 0
	miRPFilterValue = 0;
	miRDFilterValue = 0;
}

////Sets where the running sum should be initialise When the Object is first created
void synapseSingleFilterDual::initialiseFilterState()
{
	//Use injection PDF - Here bothe Running Sums Are +Ve and so are the thresholds
	//They act as counters of +ve and -ve stimuli separetely
	double p = 0.0;

	mbIsMonitored = false;

	//When Not Using a particular Cascade state we set Terminal to 0 - The injection PDF for no decay rate can then be simply calculated
	//This Injection Assumes HThres is always attached to the PRunning Sum
	if (dPDecayRate == 0.0 && dDDecayRate == 0.0) //We know The distribution of the No decay case
		{
			miRPFilterValue = miRDFilterValue = 0;//Set To Invalid Value Initially
			 //The PDF is a rising straight line with a fixed step increase at each step toward threshold
			const double probPStep = 2.0/(double)(miHThres*miHThres);
			//const double probDStep = 1/miLThres;

		//for (int i=1;i<miHThres;i++)
		//		p+=(i+1)*probPStep;

			//Inject On one filter
			assert (miHThres == miLThres);
			double r = gsl_rng_uniform(mprng);
			double q = gsl_rng_uniform(mprng);
			p=0.0;
			for (int i=0;i<miHThres;i++)
			{
				p+=(i+1)*probPStep;
				if (p >= r)
				{
					miRPFilterValue = i;
					break;
				}
			}

			p = 0.0;
			for (int i=0;i<miHThres;i++)
			{
				p+=(i+1)*probPStep;
				if (p >= q)
				{
					miRDFilterValue = i;
					break;
				}
			}

		} //If no Decay

	//miRDFilterValue = 0;
	//miRPFilterValue = 0;
	return;

}

void  synapseSingleFilterDual::getTypeAsString(char* buff)
{
	getTypeName(buff);
}


void synapseSingleFilterDual::getTypeName(char* buff)
{
	strcpy(buff,"_synapseSingleFilterDual");
}

synapseSingleFilterDual::~synapseSingleFilterDual() {
	//  Auto-generated destructor stub
}
