/*
 * synapseFilterUnifiedWithDecay.cpp
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 *
 *      This is A single unified filter with two absorbing Thresholds - Re-Injection to Zero even after crossing p-thresholds
 */

#include "synapseSingleFilterUnifiedWithDecay.h"

synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay():super() {

	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miTerminalIndex = miCascadeIndex+1; //Make It look as non terminal so the injection Uses the non Terminal PDF


	miLThres = -6;
	miHThres = 6;

	//setFilterThresholds();
	initialiseFilterState();
	mbIsMonitored = false;
	mdDecayRate = 0.0;
	miTimeSinceLastInduction = 0;

}


///Used by Threshold Cycle Measurements - Simulation To start all synapses from the same Internal State
//StartState: Asssume Desired Threshold Crossing is from Upper thres, for synapses whose target threshold is the lower one, setTracked State will mirror
//			miRFilterValue to the relevant distance from desired threshold
synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piLThres,int piHThres,double pdDecay,int StartState,int CycleSamples, gsl_rng* prng)
{
	 //Enable Metaplastic Transition Allocations
	mbPlasticAlloc = false;
	mbMetaplasticAlloc = false;
	mprng = prng;
	miTerminalIndex = miCascadeIndex = miStartIndex = 0;
	mdDecayRate 			= pdDecay;
	miLThres 				= piLThres;
	miHThres 				= piHThres;
	miTimeSinceLastInduction = 0;

	setAllocationThreshold(0); //No refraction Period
	double r=gsl_rng_uniform(mprng);

	if (r<0.5)
	{
		penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
		miHThres = piHThres;
		miLThres = piLThres;
	}
	else
	{
		penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;
		miHThres = piLThres;
		miLThres = piHThres;
	}
	//Inject Somewhere
	//initialiseFilterState();
	miRFilterValue = StartState;
	iCycleSamplesRemaining = CycleSamples; //The CountDown Starts from CycleSamples->On 0 No more Distribution Sampling occurs
	assert(miLThres <  miHThres);
	assert(miRFilterValue > piLThres && miRFilterValue < piHThres);

}

///Used by Threshold Cycle Measurements - Simulation To start all synapses from the same Internal State
//StartState: Asssume Desired Threshold Crossing is from Upper thres, for synapses whose target threshold is the lower one, setTracked State will mirror
//			miRFilterValue to the relevant distance from desired threshold
synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piLThres,int piHThres,double pdDecay,int CycleSamples, gsl_rng* prng)
{
	 //Enable Metaplastic Transition Allocations
	mbPlasticAlloc = false;
	mbMetaplasticAlloc = false;
	mprng = prng;
	miTerminalIndex = miCascadeIndex = miStartIndex = 0;
	mdDecayRate 			= pdDecay;
	miLThres 				= piLThres;
	miHThres 				= piHThres;
	miTimeSinceLastInduction = 0;


	double r=gsl_rng_uniform(mprng);

	if (r<0.5)
		penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
	else
		penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;

	uiSameThresholdTransitionCounter = gsl_ran_geometric(mprng,0.5);

	//Inject Somewhere
	initialiseFilterState();
	iCycleSamplesRemaining = CycleSamples; //The CountDown Starts from CycleSamples->On 0 No more Distribution Sampling occurs


	assert(miLThres <  miHThres);
	assert(miRFilterValue > piLThres && miRFilterValue < piHThres);

}
/*
 * This is the constructor Called by the allocation Function for Single Filters.
 */
//A Generic Constructor
synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piLThres,int piHThres,double pdDecay)
{
	 //Enable Metaplastic Transition Allocations
	mbPlasticAlloc = false;
	mbMetaplasticAlloc = false;

	miTerminalIndex = miCascadeIndex = miStartIndex = 0;
	mdDecayRate = pdDecay;
	miTimeSinceLastInduction = 0;

	setAllocationThreshold(0); //No refraction Period
	mprng = g_getRandGeneratorInstance(false);
	double r = gsl_rng_uniform(mprng);
	if (r<0.5)
	{
		penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
		miHThres = piHThres;
		miLThres = piLThres;
	}
	else
	{
		penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;
		miHThres = -piLThres;
		miLThres = -piHThres;
	}

	assert(miLThres <=  miHThres);
	//Inject Somewhere
	initialiseFilterState();
	assert(miRFilterValue > piLThres && miRFilterValue < piHThres);

}

/*
 * This is the constructor Called by the allocation Function for Single Filters.
 */
//A Generic Constructor
synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piLThres,int piHThres,double pdDecay,uint uiAllocationRefraction,double dAllocDecayRate)
{
	 //Disable all Transition Allocations
	mbPlasticAlloc = false;
	mbMetaplasticAlloc = false;
	mbStabilityAlloc = false;
	miTimeSinceLastInduction = 0;
	mdAllocationDecayRate = dAllocDecayRate; //Set The allocation Signal Stochastic Decay Rate
	miTerminalIndex = miCascadeIndex = miStartIndex = 0;
	miLThres = piLThres;
	miHThres = piHThres;
	mdDecayRate = pdDecay;
	assert(miLThres <  miHThres);
	setAllocationThreshold(uiAllocationRefraction); //The time a strength state needs to be stable for Before allocation is allowed
	assert(uiThresholdForAllocation >= 0);
	mprng = g_getRandGeneratorInstance(false);

	double r = gsl_rng_uniform(mprng); //Random Start State
	if (r < 0.5)
		penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
	else
		penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;

	//Inject Somewhere
	initialiseFilterState();
	assert(miRFilterValue > piLThres && miRFilterValue < piHThres);

}


//Works Only if Theta_m = Theta_p - Initialiases
void synapseSingleFilterUnifiedWithDecay::initialiseFilterState()
{

	const int ciMaxInternalStates = miHThres-miLThres-1; //Remove the Threshold states as the are absorbing
	double dPDF[ciMaxInternalStates];
	double dProbTransition = 1.0/(miHThres*miHThres);

	uiStateLifetime = 0;//Reset The Strength State Allocation Counter

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


int synapseSingleFilterUnifiedWithDecay::getRunningValue() const
{
	return miRFilterValue;
}

int synapseSingleFilterUnifiedWithDecay::addSample(int iValue)
{
	int Ret = 0;
	freezePlasticity(); //Lock Synapse If any of the Allocation Criteria was met at the previous memory storage event

	miTimeSinceLastInduction++; //Increment time since last induction - For Decay

	uiStateLifetime++; /////Strength Stability Allocation Mechanism Increment Induction time in the Same Strength State

	///DECAY MECHANISM
	if (iValue != 0) //Some Induction step - So Calculate Decay Steps in between
	{
		doStochasticDecay();
		miTimeSinceLastInduction = 0; //Reset time since last event
	}

	//NOP - No decay on this synapse so we don't count time
	if (iValue == 0) return 0;

	miRFilterValue +=iValue; //Integrate Stimulus Value

	//Now Check Threshold Condition;
	if (miRFilterValue <= miLThres) //Check Lower Threshold
	{
		if (mbStopRecordingOfMHistogramAtNextThresholdEvent) //Stop Measuring - THis is the quickest way to stop saving into the distribution
		{ //resetMetaplasticCounter();
		  mbSaveMetaplasticHistogram = false;
		}

		//p Threshold Reached
		if (penumStrength == SYN_STRENGTH_WEAK)
		{
			Ret = 0;//No Transition
			reInjectFilterStateToCascadeState(); //Reset running Sum - No reflecting Boundary
			miStartIndex = miCascadeIndex+1; //Makes Index Changed Flag true
			uiSameThresholdTransitionCounter++; //Metaplasticity Counter - Used for Statistics Or Allocation

		}else{//q thres Reached
			switchReset(); //Change to STRONG And Reset Cascade Index - Reset Runnning Sum
			Ret = -1; //Plastic Transition Occurred
		}

	}else{ //Make So As Both Cannot Occur

		if (miRFilterValue >= miHThres) //Check Higher Threshold
		{ //lOW tRHES
			//Threshold Event - Check if Recording should stop
			if (mbStopRecordingOfMHistogramAtNextThresholdEvent) //Stop Measuring - THis is the quickest way to stop saving into the distribution
				{ //resetMetaplasticCounter();
				  mbSaveMetaplasticHistogram = false;
				}

			if (penumStrength == SYN_STRENGTH_STRONG)
			{
				Ret = 0;//No Transition
				reInjectFilterStateToCascadeState(); //Reset running Sum
				miStartIndex = miCascadeIndex+1;//Makes Index Changed Flag true
				uiSameThresholdTransitionCounter++;// Counter - Used for Statistics Or Allocation
			}else{ //WEAK SYNAPSE -> SWITCH TO STRONG
				switchReset(); //Change to weak And Reset Cascade Index Reset Running Sum
				Ret = -1; //Plastic Transition Occurred
			}

		}
	}

	return Ret;
}

/*
 * Changes the strength of a synapse
 * If However the Plasticity is locked --The Internal state is reset without changing the Strength
 */
void synapseSingleFilterUnifiedWithDecay::switchReset()
{

	reInjectFilterStateToCascadeState(); //Reset To Zero On New State
	resetAllocationRefraction(); //uiStateLifetime = 0;//Reset The Strength State Allocation Counter
	//Do not Update Metaplastic Statistics Is Synapse Plasticity Is Locked
	if (mbNoPlasticity)//If Plasticity Is Frozen and NOT During Allocation Signal Then Do not Switch - Implements The Locking
			return;

	if (mbPlasticAlloc)
		mbNoPlasticity = true; //Freeze This Synapse from now on if Flag Of Alloc Through Plasticity is set

	//-Update Max Metaplastic Transitions and then Reset The Counter
	if (uiSameThresholdTransitionCounter > uiMaxMetaplasticTransitions)
			uiMaxMetaplasticTransitions = uiSameThresholdTransitionCounter;
	//If Strength Changes Towards correct target Strength Increment the counter - Simplifies the effect of 1st encoding

	resetMetaplasticCounter();//Sets uiMetaplasticTransitionCounter = 0;
	uiSameThresholdTransitionCounter=1;//Count 1st Threshold Crossing due to this plasticity event

	//uiSameThresholdTransitionCounter++;
	//Finally Update Synaptic Strength
	if (penumStrength == SYN_STRENGTH_STRONG){
		//Incremend Correct Threshold Crossing

		penumStrength = SYN_STRENGTH_WEAK;
	}else{
		//Incremend Correct Threshold Crossing
//		if (penumTrackedStrength == SYN_STRENGTH_STRONG) //This is used at 1st memory encoding
//			uiSameThresholdTransitionCounter++;
		penumStrength = SYN_STRENGTH_STRONG;
	}

		//setFilterThresholds(); //Mo need to re-check thresholds
}

////////DECAY FUNCTIONS///////////
/// There is a single decay rate depending on the sign of the running sum
//Decrement running sum based on decay probability eta * filter state (running sum)
//As we approach the threshold the decay increases proportionally
int synapseSingleFilterUnifiedWithDecay::doStochasticDecay()
{

	//miRFilterValue = 0;
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
void synapseSingleFilterUnifiedWithDecay::reInjectFilterStateToCascadeState()
{
	//mbNoPlasticity = true; //Lock Plasticity
	miRFilterValue = 0;
}


//Use Index / strength and the Lookup arrays to set new L/H thresholds
//Called after changes to cascade index
//TODO : Add the decay
void  synapseSingleFilterUnifiedWithDecay::setFilterThresholds()
{

//Legacy From Cascade Init Times
//	miLThres = -(*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
//	miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
//	mdDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];

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


//Overload From ICascade So I can revert the Filter State
//This Function Must be called Only Once Per trial
void  synapseSingleFilterUnifiedWithDecay::setTrackedState(uint8_t TargetStrength)
{
	//Setting To zero allows counting of Threshold Events and not just metaplastic cycles, because after then next threshold The counter will be set to 1 and will close at n
	//at which time it will be labelled as Correct or Wrong cycle. The counter will increment To 1 (at reset) whether the next threshold is same/opposite.
	// This will sync the counter to count the number of threshold crossings at the end of each cycle
	penumTrackedStrength = TargetStrength;

	//Revert Set Internal Filter State To mirror Distance From "Correct Threshold"

	if (penumTrackedStrength == SYN_STRENGTH_WEAK)
		miRFilterValue = -miRFilterValue;//Now the Filter state is set relevant to its desired Threshold State for both Strong And Weak Synapses
	//mbStopRecordingOfMHistogramAtNextThresholdEvent = true; //Testing
}



double synapseSingleFilterUnifiedWithDecay::getDecay() const
{
	return mdDecayRate;
}

double synapseSingleFilterUnifiedWithDecay::getHDecay() const
{
	return mdDecayRate;
}

double synapseSingleFilterUnifiedWithDecay::getLDecay() const
{
	return mdDecayRate;
}

int synapseSingleFilterUnifiedWithDecay::getLThres() const
{
	return miLThres;
}


int synapseSingleFilterUnifiedWithDecay::getHThres() const
{
	return miHThres;
}


void synapseSingleFilterUnifiedWithDecay::setDecay(double pdNewDecayRate)
{
	mdDecayRate = pdNewDecayRate;
}

int8_t synapseSingleFilterUnifiedWithDecay::handlePOT()
{
	return addSample(+1);
}

int8_t synapseSingleFilterUnifiedWithDecay::handleDEP()
{
	return addSample(-1);
}

void synapseSingleFilterUnifiedWithDecay::handleNOP()
{
	addSample(0);
}

//Reset is called at the beginning of each trial in simLifetime -ThresholdCycle Experiments
void synapseSingleFilterUnifiedWithDecay::reset()
{
	ICascadeSynapse::reset(); //Randomizes Start Strength, Unfreezeplasticity

	//super::reset();//If we Re-init from Non zero position Then The escape time through any boundary Reduces
		//initialiseFilterState(); //Reset running Sum - Do not - Let it be as it has been Inited by previous experience-Not the case if starting over fixed position

	//miStartIndex = miCascadeIndex; //Reset Index ??What is the poinht of this????
	//penumStrength = penumStartStrength;


	//this->unfreezePlasticity();//Base Class:: reset Does this
	//disableMetaplasticAllocation(); //Base Class : Switch off the Allocation Signal

	///Moved to ICascadeSynapse
	//uiSameThresholdTransitionCounter = gsl_ran_geometric(mprng,0.5);
}

void synapseSingleFilterUnifiedWithDecay::getTypeAsString(char* buff)
{
	getTypeName(buff);
}

void synapseSingleFilterUnifiedWithDecay::getTypeName(char* buff)
{
	strcpy(buff,"_synapseSingleFilterUnifiedWithDecay");
}


synapseSingleFilterUnifiedWithDecay::~synapseSingleFilterUnifiedWithDecay() {

}





/*
 * Removed Constructors- If you need to load a specific index Then just pass the required Values Thresh and Decay from lookup table

synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piCascadeSize,gsl_rng * rng_r):super(piCascadeSize,rng_r)
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miStartIndex = miCascadeIndex;
	miTerminalIndex = miCascadeIndex+ 1; //Make it Look as non Terminal State

	setFilterThresholds();
	initialiseFilterState();
}

synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piCascadeSize,
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
synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,gsl_rng * rng_r):super(piStartIndex+1, piStartIndex, rng_r) //Call Base Constructor
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	miCascadeIndex = miStartIndex = piStartIndex;
	miTerminalIndex = miCascadeIndex+ 1; //Make it Look as non Terminal State

	setFilterThresholds();
	initialiseFilterState();
}


//Default Value is rate=1.0
synapseSingleFilterUnifiedWithDecay::synapseSingleFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,
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
*/
