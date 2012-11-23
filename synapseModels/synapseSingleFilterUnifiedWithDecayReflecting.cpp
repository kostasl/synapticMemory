/*
 * synapseFilterUnifiedWithDecay.cpp
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 *
 *      This is A single unified filter with one absorbing and One Reflecting Threshold
 *      - Re-Injection to Zero even after crossing q-threshold only
 */

#include "synapseSingleFilterUnifiedWithDecayReflecting.h"


synapseSingleFilterUnifiedWithDecayReflecting::synapseSingleFilterUnifiedWithDecayReflecting():super() {
	//  Auto-generated constructor stub
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miTerminalIndex = miCascadeIndex+1; //Make It look as non terminal so the injection Uses the non Terminal PDF
	setFilterThresholds();
	initialiseFilterState();
	mbIsMonitored = false;
}
synapseSingleFilterUnifiedWithDecayReflecting::synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,gsl_rng * rng_r):super(piCascadeSize,rng_r)
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miStartIndex = miCascadeIndex;
	miTerminalIndex = miCascadeIndex+ 1; //Make it Look as non Terminal State

	setFilterThresholds();
	initialiseFilterState();
}

synapseSingleFilterUnifiedWithDecayReflecting::synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,
															ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
															gsl_rng * rng_r):super(piCascadeSize,penumStartStrength,rng_r)
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miCascadeIndex = miStartIndex;
	miTerminalIndex = miCascadeIndex+ 1; //Make it Look as non Terminal State

	setFilterThresholds();
	initialiseFilterState();
}
//Used to Start From A fixed point and test how distribution evolves
synapseSingleFilterUnifiedWithDecayReflecting::synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,int piStartIndex,gsl_rng * rng_r):super(piStartIndex+1, piStartIndex, rng_r) //Call Base Constructor
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miCascadeIndex = miStartIndex = piStartIndex;
	miTerminalIndex = miCascadeIndex+ 1; //Make it Look as non Terminal State


	setFilterThresholds();
	initialiseFilterState();

}
//A Generic Constructor
synapseSingleFilterUnifiedWithDecayReflecting::synapseSingleFilterUnifiedWithDecayReflecting(int piLThres,int piHThres,double pdDecay)
{

	miTerminalIndex = miCascadeIndex = miStartIndex = 0;
	miLThres = piLThres;
	miHThres = piHThres;
	mdDecayRate = pdDecay;
	assert(miLThres <  miHThres);

	double r = gsl_rng_uniform(mprng);
	if (r<0.5)
		penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
	else
		penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;

	//Inject Somewhere
	initialiseFilterState();
	assert(miRFilterValue > piLThres && miRFilterValue < piHThres);

}

//Default Value is rate=1.0
synapseSingleFilterUnifiedWithDecayReflecting::synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,int piStartIndex,
															ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
															gsl_rng *  prng_r,
															int iRateDependentParameterSet)
{
	miCascadeIndex = miStartIndex = piStartIndex;
	miTerminalIndex = miCascadeIndex + 1; //Make it Look as non Terminal State

	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	setFilterThresholds();
	initialiseFilterState();

}


//Works Only if Theta_m = Theta_p - Initialiases
//This Init Code Is for a Filter With Two Absorbing Boundaries - Init Period Required
//Or Import Tables Using Math. code double-zero-Assym-Reflecting.nb
void synapseSingleFilterUnifiedWithDecayReflecting::initialiseFilterState()
{

	const int ciMaxInternalStates = miHThres-miLThres-1; //Remove the Threshold states as the are absorbing
	double dPDF[ciMaxInternalStates];
	double dProbTransition = 1.0/(miHThres*miHThres);

	//Make The triangular Distribution
	for (int i=1;i<=miHThres;i++)
	{
		dPDF[i-1] = dProbTransition*i;
		dPDF[ciMaxInternalStates-i] = dPDF[i-1]; //Symmetric
	}
	double p =  gsl_rng_uniform(mprng);
	double cp = 0.0;
	for (int i=0;i<ciMaxInternalStates;i++)
	{
		cp += dPDF[i];
		if (cp > p)
		{
			miRFilterValue = i + miLThres+1; //miLThres is -ve
			break;
		}
	}
			//if (miRFilterValue >= miHThres || miRFilterValue <= miLThres)
	//	miRFilterValue=0;
}

int synapseSingleFilterUnifiedWithDecayReflecting::addSample(int iValue)
{
	int Ret = 0;

	miTimeSinceLastInduction++; //Increment time since last induction
	if (iValue != 0) //Some Induction step - So Calculate Decay Steps in between
	{
		doStochasticDecay();
		miTimeSinceLastInduction = 0; //Reset time since last event
	}

	//NOP - No decay on this synapse so we don't count time
	if (iValue == 0) return 0;

	miRFilterValue +=iValue;

	 //Now Check Threshold Condition;
	if (miRFilterValue <= miLThres)
	{
		//p Threshold Reached
		if (penumStrength == SYN_STRENGTH_WEAK)
		{
			Ret = 0;//No Transition
			//miRFilterValue = miLThres; //Holding Barrier
			miRFilterValue -=iValue; //Reflecting Barrier
			//reInjectFilterStateToCascadeState(); //Reset running Sum - No reflecting Boundary
			miStartIndex = miCascadeIndex+1; //Makes Index Changed Flag true
		}else{//q thres Reached
			switchReset(); //Change to STRONG And Reset Cascade Index - Reset Runnning Sum
			Ret = -1; //Plastic Transition Occured
		}
	}else //Make So As Both Cannot Occur
		if (miRFilterValue >= miHThres)
		{ //lOW tRHES
			if (penumStrength == SYN_STRENGTH_STRONG)
			{
				Ret = 0;//No Transition
				//miRFilterValue = miHThres; //Holding Barrier
				miRFilterValue -=iValue; //Reflecting Barrier
				//reInjectFilterStateToCascadeState(); //Reset running Sum
				miStartIndex = miCascadeIndex+1;//Makes Index Changed Flag true
			}else{ //WEAK SYNAPSE -> SWITCH TO STRONG
				switchReset(); //Change to weak And Reset Cascade Index Reset Runnning Sum
				Ret = -1; //Plastic Transition Occurred
			}
		}

	return Ret;
}


void synapseSingleFilterUnifiedWithDecayReflecting::switchReset()
{
	if (penumStrength == SYN_STRENGTH_STRONG){
		penumStrength = SYN_STRENGTH_WEAK;
	}else{
		penumStrength = SYN_STRENGTH_STRONG;
		}

		//setFilterThresholds(); //Mo need to re-check thresholds
		reInjectFilterStateToCascadeState(); //Reset To Zero On New State
}
////////DECAY FUNCTIONS///////////
/// There is a single decay rate depending on the sign of the running sum
//Decrement running sum based on decay probability eta * filter state (running sum)
//As we approach the threshold the decay increases proportionally
int synapseSingleFilterUnifiedWithDecayReflecting::doStochasticDecay()
{

	double p;
	unsigned int rDecaySteps;

	//g_rng_r = getRandGeneratorInstance();
	if (miRFilterValue == 0) return 0;

	if (miRFilterValue > 0) //p side decay
	{
		p = 1.0-exp(-mdDecayRate*miTimeSinceLastInduction);
		rDecaySteps = gsl_ran_binomial(mprng,p,miRFilterValue);
		miRFilterValue -= rDecaySteps; //Subtract to move to zero
	}
	else //Increment Sum - q side is decaying back to zero
	{
		p = 1.0-exp(-mdDecayRate*miTimeSinceLastInduction);
		rDecaySteps = gsl_ran_binomial(mprng,p,-miRFilterValue);
		miRFilterValue += rDecaySteps; //Add To move to zero
	}


	return rDecaySteps;
}


//Sets where the running sum should be initialized after every change to Cascade state -
void synapseSingleFilterUnifiedWithDecayReflecting::reInjectFilterStateToCascadeState()
{
	miRFilterValue = 0;
}


//Use Index / strength and the Lookup arrays to set new L/H thresholds
//Called after changes to cascade index
//TODO : Add the decay
void  synapseSingleFilterUnifiedWithDecayReflecting::setFilterThresholds()
{

	miLThres = -(*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
	miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
	mdDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];

	//g_rng_r = getRandGeneratorInstance(false);
	//Shuffle Thresholds
	double r = gsl_rng_uniform(mprng);
	if (r < 0.5)
	{
		dummy4 = miLThres;
		miLThres = -miHThres;
		miHThres = -dummy4;
	}
}

double synapseSingleFilterUnifiedWithDecayReflecting::getHDecay() const
{
	return mdDecayRate;
}

double synapseSingleFilterUnifiedWithDecayReflecting::getLDecay() const
{
	return mdDecayRate;
}

double synapseSingleFilterUnifiedWithDecayReflecting::getDecay() const
{
	return mdDecayRate;
}

void synapseSingleFilterUnifiedWithDecayReflecting::setDecay(double pdNewDecayRate)
{
	mdDecayRate = pdNewDecayRate;
}

int8_t synapseSingleFilterUnifiedWithDecayReflecting::handlePOT()
{
	return super::handlePOT();
}

int8_t synapseSingleFilterUnifiedWithDecayReflecting::handleDEP()
{
	return super::handleDEP();
}

void synapseSingleFilterUnifiedWithDecayReflecting::handleNOP()
{
	super::handleNOP();
}

void synapseSingleFilterUnifiedWithDecayReflecting::reset()
{
	//super::reset(); //If we Re-init from Non zero position Then The escape time through any boundary Reduces
	 //The reset is only called by the testEscapetime Function - And thus it affects the result of this function
	//initialiseFilterState(); //Reset running Sum - Do not - Let it be as it has been Inited by previous experience
	miStartIndex = miCascadeIndex; //Reset Index
	penumStrength = penumStartStrength;
}

void synapseSingleFilterUnifiedWithDecayReflecting::getTypeAsString(char* buff)
{
	getTypeName(buff);
}



void synapseSingleFilterUnifiedWithDecayReflecting::getTypeName(char* buff)
{
	strcpy(buff,"_singleFilterUnifiedWDecayRefl");
}

synapseSingleFilterUnifiedWithDecayReflecting::~synapseSingleFilterUnifiedWithDecayReflecting() {

}
