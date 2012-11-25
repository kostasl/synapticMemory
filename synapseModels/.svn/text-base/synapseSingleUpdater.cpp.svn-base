/*
 * synapseCascade.cpp
 *
 *  Created on: 22 Mar 2011
 *      Author: kostasl
 *
 *   No Need to initialize synapse with allocation thresh.
 *     setAllocationThreshold is called during simulation to set uiThresholdForAllocation
 */

#include "synapseSingleUpdater.h"

//Static Variables - Accessible By All Instances
float synapseSingleUpdater::mfQprob[MAX_CASCADE_SIZE];
float synapseSingleUpdater::mfPprob[MAX_CASCADE_SIZE];

synapseSingleUpdater::synapseSingleUpdater() {
	//  Auto-generated constructor stub


	miCascadeSize = DEFAULT_CASCADE_SIZE;
	initDefaultTransitionProb();

}


//Initialize As Stochastic Updater - A cascade of Size 1 with a given Q probability
synapseSingleUpdater::synapseSingleUpdater(double pfQ, gsl_rng * rng_r, int MetaCycleSamples)
{
	miCascadeSize = 1; //Because S.U Synapse
	mfQprob[0] = pfQ;
	mfPprob[0] = pfQ; //Symmetric Metaplastic Probability
	iCycleSamplesRemaining = MetaCycleSamples;
	init(miCascadeSize 	,mfQprob,mfPprob,0,SYN_STRENGTH_NOTSET,rng_r);
}


//Initialize with cascadeSize - So transition probability is set to cascade State n
synapseSingleUpdater::synapseSingleUpdater(int piCascadeSize,gsl_rng * rng_r){

	miCascadeSize 	= piCascadeSize;
	initDefaultTransitionProb();

	//Set to NOT-SET so INIT function will set the random start values - Randomly Set Cascade Index
	init(piCascadeSize 	,mfQprob,mfPprob,-1,SYN_STRENGTH_NOTSET,rng_r);

}


//Used to Start From A fixed point and test how distribution evolves
synapseSingleUpdater::synapseSingleUpdater(int piCascadeSize,int piStartIndex,gsl_rng * rng_r)
{
	mprng = rng_r;
	//mpFilter = 0; //Null pointer for Filter
	miCascadeSize 	= piCascadeSize;
	miCascadeStartIndex = piStartIndex;
	if ((piCascadeSize-1) < piStartIndex) //Then Single Stochastic Updater - Or Filter
	{
		miCascadeSize = 1; //Just being Explicit
		//Take the n=6 q and p values as the transition probs of a single Stoch. Updater
		mfQprob[0] = pow(DEFAULT_CASCADE_X,piStartIndex); //q_i=x^(i-1)
		mfPprob[0] = pow(DEFAULT_CASCADE_X,piStartIndex+1)/(1-DEFAULT_CASCADE_X); //As given in Fusi '05 p=x^i/(1-x) x=1/2
		miCascadeIndex = miCascadeStartIndex = 0; //Now Change Running Cascade index to 0  - So as to have a stochUpdater
	}
	else
		initDefaultTransitionProb();;

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,miCascadeStartIndex,SYN_STRENGTH_NOTSET,rng_r);

}

synapseSingleUpdater::synapseSingleUpdater(int piCascadeSize,int piStartIndex,SYN_STRENGTH_STATE setStrength,gsl_rng * rng_r)
{
		miCascadeSize 	= piCascadeSize; //Have to be set before calling InitDefault

		if ((piCascadeSize-1) < piStartIndex) //Then Single Stochastic Updater - Or Filter
		{
			miCascadeSize = 1; //Just being Explicit
			mfQprob[0] = pow(DEFAULT_CASCADE_X,piStartIndex); //q_i=x^(i-1)
			mfPprob[0] = pow(DEFAULT_CASCADE_X,piStartIndex+1)/(1-DEFAULT_CASCADE_X); //As given in Fusi '05 p=x^i/(1-x) x=1/2
		}
		else
			initDefaultTransitionProb();

		//Set to NOT-SET so INIT function will set the random start values
		init(piCascadeSize 	,mfQprob,mfPprob,piStartIndex,setStrength,rng_r);

		penumStrength = penumStartStrength = setStrength;
		miCascadeIndex = miCascadeStartIndex = piStartIndex; //Fix Start position
}
//Init With StartStrength Used by SynapticConnection
synapseSingleUpdater::synapseSingleUpdater(int piCascadeSize,float startStrength,gsl_rng * rng_r)
{
///
	if (piCascadeSize < 1)
		liberrexit(500,"Invalid CascadeSize Parameter");

	initDefaultTransitionProb();
	//Set to NOT-SET so INIT function will set the random start values
	if (startStrength > 0.0)
	{
		init(miCascadeSize 	,mfQprob,mfPprob,0,SYN_STRENGTH_STRONG,rng_r);
	}
	else
	{
		init(miCascadeSize 	,mfQprob,mfPprob,0,SYN_STRENGTH_WEAK,rng_r);
	}



}

synapseSingleUpdater::synapseSingleUpdater(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r)
{


	miCascadeSize 	= piCascadeSize; //Have to be set before calling InitDefault
	initDefaultTransitionProb();

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,-1,penumStartStrength,rng_r);

	}


//Initializes the Transition probabilities P's and Qs in a geometric progression as defined in the 05 Fusi Paper
void synapseSingleUpdater::initDefaultTransitionProb()
{

	assert(miCascadeSize > 0); //Check for errors

	//double sumQ = 0;
	//BUG FIX Using CascadeIndex As Index! 8/3/10
	for (int i =0;i<miCascadeSize;i++)
	{
		mfQprob[i] = pow(DEFAULT_CASCADE_X,i); //q_i=x^(i-1)
		mfPprob[i] = pow(DEFAULT_CASCADE_X,i+1)/(1-DEFAULT_CASCADE_X); //As given in Fusi '05 p=x^i/(1-x) x=1/2
		//sumQ +=		mfQprob[miCascadeIndex] ; //For Verification
	}
	//For Boundary effects last q is divided by (1-x) --Keep as 1-0.5 for clarity of the value X
	//As the Arrays are static and global the boundary condition is taken care of in the HandlePot/DEP functions

}



void synapseSingleUpdater::init(int piStatesCount, //Number of Cascade states n
						float pfPCascadeSwitch[], //Probability array of q's
						float pfPCascadeTrans[],  //Probability array of p's
						int piCurrCascadeIndex, //Cascade Current State
						SYN_STRENGTH_STATE pfCurrStrengthState,
						gsl_rng * prng_r
						)
{

	double r; //random num
	///Check Cascade Size Parameter
	if (piStatesCount > MAX_CASCADE_SIZE)
		liberrexit(900,"Max Number of Cascade Size Violation. Decrease or Change constant.");
	//mbNoPlasticity is Inited in ISynapse
	mprng = prng_r;
	mbIsMonitored = false;
	mbStrengthChanged = false; //Reset Monitoring Flag
	//Copy to Member Variables
	miCascadeSize = piStatesCount;
	miCascadeStartIndex = miCascadeIndex = piCurrCascadeIndex; //if -1 then it will be randomly assigned later
	penumStrength = penumStartStrength = pfCurrStrengthState;
	uiSameThresholdTransitionCounter = 0; //Should Be initialized Geometrically But For sufficient Trials This will quickly converge Given at least 300 mems per trial.
	mbNoPlasticity = false;

	if (mfQprob != pfPCascadeSwitch )//Check If pointer to Same Data
	{
		memcpy(mfQprob,pfPCascadeSwitch,sizeof(int)*piStatesCount);
	}

	if (mfPprob != pfPCascadeTrans) //Check If pointer to Same Data
	{
		memcpy(mfPprob,pfPCascadeTrans,sizeof(int)*piStatesCount);
	}

	//Init Random Members if required
	if (pfCurrStrengthState == SYN_STRENGTH_NOTSET)
	{
		//Save Start State So We can count Synapse change in lifetime
		r = gsl_rng_uniform(mprng);
		penumStrength = penumStartStrength = (r < 0.5)?SYN_STRENGTH_WEAK:SYN_STRENGTH_STRONG;
	}

//	//Set the Random Index if required
//	if (piCurrCascadeIndex == -1)
//	{
//		//Set Random INdex Uniformly Distributed among Cascade States
//		r = gsl_rng_uniform(mprng)*0.999;
//		miCascadeIndex = miCascadeStartIndex = floor(miCascadeSize*r);//Take Cascade Index as concat double to give 0-Count
//
//	}

	assert(miCascadeIndex >= 0 && miCascadeIndex < miCascadeSize);
	reset();
	///Finished Initialisation
}

//Induces an LTP event depending on current state either switches cascades or moves down
int8_t synapseSingleUpdater::handlePOT()
{
	int8_t Ds = 0;
	float fQProb;
	float fPProb;
	assert(miCascadeIndex == 0);
	fQProb = mfQprob[miCascadeIndex];
	fPProb = mfPprob[miCascadeIndex];
	freezePlasticity(); //Lock Synapse If any of the Allocation Criteria was met at the previous memory storage event

	//Draw random p
	double lfp = gsl_rng_uniform(mprng);
    //Switch operational In POT Mode, so potentiate
    ///Check Appropriate Transition Probability
    if (penumStrength == SYN_STRENGTH_WEAK)
    {
		//Check Q for transition to Opposite Cascade
		if (lfp < fQProb)
		{
			//Ds = SYN_STRENGTH_STRONG-SYN_STRENGTH_WEAK;
			//miCascadeIndex = 0;

			if (!mbNoPlasticity) //If Plasticity Is Frozen Then Do not Switch
			{
				if (uiSameThresholdTransitionCounter > uiMaxMetaplasticTransitions)
					uiMaxMetaplasticTransitions = uiSameThresholdTransitionCounter;

				resetAllocationRefraction(); //uiStateLifetime = 0;//Reset The Strength State Allocation Counter
				resetMetaplasticCounter();//Sets uiMetaplasticTransitionCounter = 0;
				uiSameThresholdTransitionCounter = 1;//Count 1st Threshold Crossing due to this plasticity event

				///Change Strength After Meta.Cycles Distrib Has been saved
				penumStrength = SYN_STRENGTH_STRONG;
				mbStrengthChanged = true;

				Ds = -1; //Signal a switch Change
			}
			//mbNoPlasticity = true; //Freeze This Synapse from now on
			mbNoPlasticity = mbPlasticAlloc;//If Allocation through Plastic Transition is set Lock Plasticity
		}
    }
    else
    { //Already in STRONG  - Do Nothing
    	if (penumStrength != SYN_STRENGTH_STRONG)
    		liberrexit(666,"Error:S.Updater Synapse not in STRONG state as expected!");
		//mbCascadeIndexChanged = true; //Flag Index Change So Synapse Moves into tracked Group - But for Stoch.

    	///LOCK METAPLASTIC TRANS
		 if ((fPProb) > lfp){ //Accept Locking with p rate
			 uiSameThresholdTransitionCounter++;
			 mbNoPlasticity = mbMetaplasticAlloc; //Freeze This Metaplastic Synapse from now on if Flag Is Set
		 }
	}

    //mbNoPlasticity = false;
return Ds;
}

//Induces an LTD event depending on current state either switches cascades or moves down
int8_t synapseSingleUpdater::handleDEP()
{
	int8_t Ds = 0;
	float fQProb; //The Probability of Switching Strength
	float fPProb; //There is no P transition Unless it is used for allocation-Locking Plasticity

	fQProb = mfQprob[miCascadeIndex];
	fPProb = mfPprob[miCascadeIndex];
	freezePlasticity(); //Lock Synapse If any of the Allocation Criteria was met at the previous memory storage event

		//Draw random p
	double lfp = gsl_rng_uniform(mprng);

    //Switch operational In POT Mode, so potentate
    ///Check Appropriate Transition Probability
    if (penumStrength == SYN_STRENGTH_STRONG)
    {
		//Check Q for transition to Opposite Cascade
		if (lfp < fQProb)
		{
			//Ds = SYN_STRENGTH_WEAK-SYN_STRENGTH_STRONG;

			if (!mbNoPlasticity) //If Plasticity Is Frozen Then Do not Switch
			{

				//-Update Max Metaplastic Transitions and then Reset The Counter
				if (uiSameThresholdTransitionCounter > uiMaxMetaplasticTransitions)	uiMaxMetaplasticTransitions = uiSameThresholdTransitionCounter;

				resetAllocationRefraction(); //uiStateLifetime = 0;//Reset The Strength State Allocation Counter
				resetMetaplasticCounter();//Sets uiMetaplasticTransitionCounter = 0;
				uiSameThresholdTransitionCounter = 1;//Count 1st Threshold Crossing due to this plasticity event

				penumStrength = SYN_STRENGTH_WEAK;
				mbStrengthChanged = true;
				Ds = -1;
			}
			mbNoPlasticity = mbPlasticAlloc;//If Allocation through Plastic Transition is set Lock Plasticity

		}
    }
    else
    { //Already in WEAK cascade - Move down in Cascade
    	if (penumStrength != SYN_STRENGTH_WEAK)
    		liberrexit(100,"Error:S.Updater Synapse not in WEAK state as expected!");
			//mbCascadeIndexChanged = true; //Flag Index Change So Synapse Moves into tracked Group

    	///LOCK METAPLASTIC TRANS
		 if ((fPProb) > lfp)
		 { //Accept Locking with p rate
			uiSameThresholdTransitionCounter++;
			mbNoPlasticity = mbMetaplasticAlloc; //Freeze This Metaplastic Synapse from now on if Flag Is Set
		 }
    }

   // mbNoPlasticity = false;
    return Ds;
}

//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
void synapseSingleUpdater::handleNOP()
{
////Nothing to see here
}


//Overloaded Methods
int synapseSingleUpdater::getCascadeIndex() const
{
	return miCascadeIndex;
}
int synapseSingleUpdater::getCascadeSize() const
{
return miCascadeSize;
}
int synapseSingleUpdater::getStrength() const
{
return penumStrength;
}

//Returns whether a change of state has occured
int synapseSingleUpdater::getStartStrength() const
{
return penumStartStrength;
}
//Returns whether a change of state has occured
bool synapseSingleUpdater::hasStrengthModified() const
{

	return (penumStrength != penumStartStrength); //It appears to be so for the '05 Paper
	//return mbStrengthChanged;
}

// Returns true if Index has Been Incremented
bool synapseSingleUpdater::hasIndexModified() const
{
	return mbCascadeIndexChanged;
}


//Resets Flag which holds if switchSynapse has changed.
void synapseSingleUpdater::startMonitoring()
{
	mbStrengthChanged = false;
	mbCascadeIndexChanged = false;
	penumStartStrength = penumStrength; //Reset
	mbIsMonitored = true;
}
void synapseSingleUpdater::stopMonitoring()
{
	mbIsMonitored = false;
}

bool synapseSingleUpdater::isMonitored() const
{
	return mbIsMonitored;
}

void synapseSingleUpdater::reset()
{
	ICascadeSynapse::reset(); //Randomizes Start Strength, Unfreezes plasticity

	//menumStrength_State = menumStartStrength_State;

	mbStrengthChanged = false;
	mbCascadeIndexChanged = false;
	mbNoPlasticity = false;
	miCascadeIndex =  miCascadeStartIndex;
	mbIsMonitored = false;

}

void synapseSingleUpdater::getTypeAsString(char* buff)
{
	synapseSingleUpdater::getTypeName(buff);
}

void synapseSingleUpdater::getTypeName(char* buff)
{
	strcpy(buff,"_synapseSingleUpdater");
}
//Filter Baggage -Inherited from the ICASCADESYnapse
int synapseSingleUpdater::getLThres() const
{
	return 0;
}

int synapseSingleUpdater::getHThres() const
{
	return 0;
}
//The Decay Currently running
double synapseSingleUpdater::getDecay() const
{
	return 0;
}



synapseSingleUpdater::~synapseSingleUpdater() {
	//  Auto-generated destructor stub
}
