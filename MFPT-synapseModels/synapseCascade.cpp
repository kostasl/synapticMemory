/*
 * synapseCascade.cpp
 *
 *  Created on: 22 Mar 2011
 *      Author: kostasl
 */

#include "synapseCascade.h"

//Static Variables - Accessible By All Instances
float synapseCascade::mfQprob[MAX_CASCADE_SIZE];
float synapseCascade::mfPprob[MAX_CASCADE_SIZE];

synapseCascade::synapseCascade() {



	miCascadeSize = DEFAULT_CASCADE_SIZE;
	initDefaultTransitionProb();

}
//Initialize As Stochastic Updater - A cascade of Size 1 with a given Q probability
synapseCascade::synapseCascade(double pfQ,gsl_rng * rng_r)
{
	miCascadeSize = 1;
	mfQprob[0] = pfQ;
	init(1 	,mfQprob,mfPprob,0,SYN_STRENGTH_NOTSET,rng_r);
}

//Initialize with cascadeSize
synapseCascade::synapseCascade(int piCascadeSize,gsl_rng * rng_r){

	miCascadeSize 	= piCascadeSize;
	initDefaultTransitionProb();


	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,-1,SYN_STRENGTH_NOTSET,rng_r);

}
// Custom Transition Probability Initialization
synapseCascade::synapseCascade(float pQprob[],float pPprob[],int piCascadeSize,gsl_rng * rng_r)
{
	miCascadeSize 	= piCascadeSize;

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize,pQprob,pPprob,-1,SYN_STRENGTH_NOTSET,rng_r);

}
//Used to Start From A fixed point and test how distribution evolves
synapseCascade::synapseCascade(int piCascadeSize,int piStartIndex,gsl_rng * rng_r)
{
	mprng = rng_r;
	//mpFilter = 0; //Null pointer for Filter
	miCascadeSize 	= piCascadeSize;
	miCascadeStartIndex = piStartIndex;

	initDefaultTransitionProb();;

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,miCascadeStartIndex,SYN_STRENGTH_NOTSET,rng_r);

}

synapseCascade::synapseCascade(int piCascadeSize,int piStartIndex,SYN_STRENGTH_STATE penumStrength,gsl_rng * rng_r)
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
		init(piCascadeSize 	,mfQprob,mfPprob,piStartIndex,penumStrength,rng_r);

		menumStrength_State = menumStartStrength_State = penumStrength;
		miCascadeIndex = miCascadeStartIndex = piStartIndex; //Fix Start position
}
//Init With StartStrength Used by SynapticConnection
synapseCascade::synapseCascade(int piCascadeSize,float startStrength,gsl_rng * rng_r)
{

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

synapseCascade::synapseCascade(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r)
{

	miCascadeSize 	= piCascadeSize; //Have to be set before calling InitDefault
	initDefaultTransitionProb();

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,-1,penumStartStrength,rng_r);

}


//Initializes the Transition probabilities P's and Qs in a geometric progression as defined in the 05 Fusi Paper
void synapseCascade::initDefaultTransitionProb()
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



void synapseCascade::init(int piStatesCount, //Number of Cascade states n
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
	uiSameThresholdTransitionCounter = miCascadeStartIndex; //Set Same threshold Count the same as index

	menumStrength_State = menumStartStrength_State = pfCurrStrengthState;

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
		menumStrength_State = menumStartStrength_State = (r < 0.5)?SYN_STRENGTH_WEAK:SYN_STRENGTH_STRONG;
	}

	//Set the Random Index if required
	if (piCurrCascadeIndex == -1)
	{
		//Set Random INdex Uniformly Distributed among Cascade States
		r = gsl_rng_uniform(mprng)*0.999;
		miCascadeIndex = miCascadeStartIndex = floor(miCascadeSize*r);//Take Cascade Index as concat double to give 0-Count

	}

	assert(miCascadeIndex >= 0 && miCascadeIndex < miCascadeSize);
	reset();
	///Finished Initialisation
}

//Induces an LTP event depending on current state either switches cascades or moves down
int8_t synapseCascade::handlePOT()
{
	int8_t Ds = 0;
	float fQProb,fPProb;

	if (miCascadeIndex < (miCascadeSize-1)) //Not Terminal State
	{
		fQProb = mfQprob[miCascadeIndex];
		fPProb = mfPprob[miCascadeIndex];
	}
	else
	{
		fQProb = mfQprob[miCascadeIndex]/(1-DEFAULT_CASCADE_X);
		fPProb = 0;
	}
	//Draw random p
	double lfp = gsl_rng_uniform(mprng);
    //Switch operational In POT Mode, so potentiate
    ///Check Appropriate Transition Probability
    if (menumStrength_State == SYN_STRENGTH_WEAK)
    {
		//Check Q for transition to Opposite Cascade
		if (fQProb > lfp)
		{
			//Ds = SYN_STRENGTH_STRONG-SYN_STRENGTH_WEAK;

			mbCascadeIndexChanged = (miCascadeIndex !=0);
			miCascadeIndex = 0;
			//Update Max Metaplastic Transitions and then Reset The Counter
			if (uiSameThresholdTransitionCounter > uiMaxMetaplasticTransitions)
					uiMaxMetaplasticTransitions = uiSameThresholdTransitionCounter;
			uiSameThresholdTransitionCounter = 0; //Reset Same threshold Count
			Ds = -1; //Signal a switch Change

			menumStrength_State = SYN_STRENGTH_STRONG;
			mbStrengthChanged = true;
		}
    }
    else
    { //Already in STRONG cascade - Move down in Cascade
    	if (menumStrength_State != SYN_STRENGTH_STRONG)
    		liberrexit(666,"Error:Cascade Synapse not in STRONG state as expected!");

    	//Check p for transition Down the Cascade
		if (fPProb > lfp) //fPProb=0 At Terminal So No metaplastic Transition Occurs
		{
			uiSameThresholdTransitionCounter++; //Increment On all P transitions
			assert(miCascadeSize > miCascadeIndex+1); //Check for end of Cascade (No need cause fPProb=0 at terminal)

			miCascadeIndex++;//Move down Cascade
			mbCascadeIndexChanged = true;
			Ds = 1; //Signal A segregation change
		}


    }

return Ds;
}

//Induces an LTD event depending on current state either switches cascades or moves down
int8_t synapseCascade::handleDEP()
{
	int8_t Ds = 0;
	float fQProb,fPProb;

	if (miCascadeIndex < (miCascadeSize-1)) //Not Terminal State
	{
		fQProb = mfQprob[miCascadeIndex];
		fPProb = mfPprob[miCascadeIndex];
	}
	else
	{
		fQProb = mfQprob[miCascadeIndex]/(1-DEFAULT_CASCADE_X);
		fPProb = 0;
	}
//Draw random p
	double lfp = gsl_rng_uniform(mprng);


    //Switch operational In POT Mode, so potentate
    ///Check Appropriate Transition Probability
    if (menumStrength_State == SYN_STRENGTH_STRONG)
    {
		//Check Q for transition to Opposite Cascade
		if (fQProb > lfp)
		{
			mbCascadeIndexChanged = (miCascadeIndex !=0);
			miCascadeIndex = 0;
			//Update Max Metaplastic Transitions and then Reset The Counter
			if (uiSameThresholdTransitionCounter > uiMaxMetaplasticTransitions)
					uiMaxMetaplasticTransitions = uiSameThresholdTransitionCounter;
			uiSameThresholdTransitionCounter = 0; //Reset Same threshold Count

			//Ds = SYN_STRENGTH_WEAK-SYN_STRENGTH_STRONG;
			Ds = -1;
			menumStrength_State = SYN_STRENGTH_WEAK;
			mbStrengthChanged = true;
		}
    }
    else
    { //Already in WEAK cascade - Move down in Cascade
    	if (menumStrength_State != SYN_STRENGTH_WEAK)
    		liberrexit(100,"Error:Cascade Synapse not in WEAK state as expected!");

    	//Check p for transition Down the Cascade
		if (fPProb > lfp)
		{
			uiSameThresholdTransitionCounter++; //Increment Regardless Of Cascade Size
			//menumStrength_State = SYN_STRENGTH_STRONG; //Same State
			assert(miCascadeSize > miCascadeIndex+1); //Check for end of Cascade

			miCascadeIndex++;//Move down Cascade
			mbCascadeIndexChanged = true;
			Ds = 1; //Segregation signal

		}

    }

    return Ds;

}

//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
void synapseCascade::handleNOP()
{
////Nothing to see here
}


//Overloaded Methods
int synapseCascade::getCascadeIndex() const
{
	return miCascadeIndex;
}
int synapseCascade::getCascadeSize() const
{
return miCascadeSize;
}
int synapseCascade::getStrength() const
{
return menumStrength_State;
}

//Returns whether a change of state has occured
int synapseCascade::getStartStrength() const
{
return menumStartStrength_State;
}
//Returns whether a change of state has occured
bool synapseCascade::hasStrengthModified() const
{

	return (menumStrength_State != menumStartStrength_State); //It appears to be so for the '05 Paper
	//return mbStrengthChanged;
}

// Returns true if Index has Been Incremented
bool synapseCascade::hasIndexModified() const
{
	return mbCascadeIndexChanged;
}


//Resets Flag which holds if switchSynapse has changed.
void synapseCascade::startMonitoring()
{
	mbStrengthChanged = false;
	mbCascadeIndexChanged = false;
	menumStartStrength_State = menumStrength_State; //Reset
	mbIsMonitored = true;
}
void synapseCascade::stopMonitoring()
{
	mbIsMonitored = false;
}

bool synapseCascade::isMonitored() const
{
	return mbIsMonitored;
}

void synapseCascade::reset()
{
	menumStrength_State = menumStartStrength_State;
	mbStrengthChanged = false;
	mbCascadeIndexChanged = false;
	mbNoPlasticity = false;
	miCascadeIndex =  miCascadeStartIndex;
	mbIsMonitored = false;
}

void synapseCascade::getTypeAsString(char* buff)
{
	synapseCascade::getTypeName(buff);
}

void synapseCascade::getTypeName(char* buff)
{
	strcpy(buff,"_synapseCascade");
}
//Filter Baggage -Inherited from the ICASCADESYnapse
int synapseCascade::getLThres() const
{
	return 0;
}

int synapseCascade::getHThres() const
{
	return 0;
}
//The Decay Currently running
double synapseCascade::getDecay() const
{
	return 0;
}



synapseCascade::~synapseCascade() {
	//  Auto-generated destructor stub
}
