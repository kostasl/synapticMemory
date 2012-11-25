/*
 * synapseCascade.cpp
 *
 * This is a transfer of a the Fusi cascade Synapse To a serial state cascade model of Leibolt and Kempter
 * S<->S<->S<->S<->W<->W<->W<->W<->W  Synapses move serially between states-Cascade Indexes
 *
 * Unlike the Cascade - The serial Synapse Moves Between states deterministically
 *  There are no cross transitions so Q represents a movement towards the opposite Strength
 * Cascade Index Can be = Cascade Size - Because each strength state has Cascade Indexes from 1 to Casc.Size
 *
 * - NEEDS VERIFICATION -
 *  Created on: 5 Oct 2011
 *      Author: kostasl
 */

#include "synapseSerialCascade.h"

//Static Variables - Accessible By All Instances
float synapseSerialCascade::mfQprob[MAX_CASCADE_SIZE];
float synapseSerialCascade::mfPprob[MAX_CASCADE_SIZE];

synapseSerialCascade::synapseSerialCascade() {
	//  Auto-generated constructor stub


	miCascadeSize = DEFAULT_CASCADE_SIZE;
	initDefaultTransitionProb();

}
//Initialize As Stochastic Updater - A cascade of Size 1 with a given Q probability
synapseSerialCascade::synapseSerialCascade(double pfQ,gsl_rng * rng_r)
{
	miCascadeSize = 1;
	mfQprob[0] = pfQ;
	init(1 	,mfQprob,mfPprob,0,SYN_STRENGTH_NOTSET,rng_r);
}

//Initialize with cascadeSize
synapseSerialCascade::synapseSerialCascade(int piCascadeSize,gsl_rng * rng_r){

	miCascadeSize 	= piCascadeSize;
	initDefaultTransitionProb();


	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,-1,SYN_STRENGTH_NOTSET,rng_r);

}


//Used to Start From A fixed point and test how distribution evolves
synapseSerialCascade::synapseSerialCascade(int piCascadeSize,int piStartIndex,gsl_rng * rng_r)
{
	mprng = rng_r;
	//mpFilter = 0; //Null pointer for Filter
	miCascadeSize 	= piCascadeSize;
	miCascadeStartIndex = piStartIndex;
	initDefaultTransitionProb();

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,miCascadeStartIndex,SYN_STRENGTH_NOTSET,rng_r);

}

synapseSerialCascade::synapseSerialCascade(int piCascadeSize,int piStartIndex,SYN_STRENGTH_STATE penumStrength,gsl_rng * rng_r)
{
		miCascadeSize 	= piCascadeSize; //Have to be set before calling InitDefault

		initDefaultTransitionProb();

		//Set to NOT-SET so INIT function will set the random start values
		init(piCascadeSize 	,mfQprob,mfPprob,piStartIndex,penumStrength,rng_r);

		menumStrength_State = menumStartStrength_State = penumStrength;
		miCascadeIndex = miCascadeStartIndex = piStartIndex; //Fix Start position
}
//Init With StartStrength Used by SynapticConnection
synapseSerialCascade::synapseSerialCascade(int piCascadeSize,float startStrength,gsl_rng * rng_r)
{

	if (piCascadeSize < 1)
		liberrexit(500,"Invalid CascadeSize Parameter");

	initDefaultTransitionProb();
	//Set to NOT-SET so INIT function will set the random start values
	if (startStrength > 0.0)
	{
		init(miCascadeSize,mfQprob, mfPprob,0,SYN_STRENGTH_STRONG,rng_r);
	}
	else
	{
		init(miCascadeSize,mfQprob, mfPprob,0,SYN_STRENGTH_WEAK,rng_r);
	}



}

synapseSerialCascade::synapseSerialCascade(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r)
{


	miCascadeSize 	= piCascadeSize; //Have to be set before calling InitDefault
	initDefaultTransitionProb();

	//Set to NOT-SET so INIT function will set the random start values
	init(piCascadeSize 	,mfQprob,mfPprob,-1,penumStartStrength,rng_r);

	}


//Initializes the Transition probabilities P's and Qs in a geometric progression as defined in the 05 Fusi Paper
void synapseSerialCascade::initDefaultTransitionProb()
{
	assert(miCascadeSize > 0); //Check for errors
	//double sumQ = 0;
	//BUG FIX Using CascadeIndex As Index! 8/3/10
	for (int i =0;i<miCascadeSize;i++)
	{
		mfQprob[i] = 1.0; //Unlike the Cascade - The serial Synapse Moves Between states deterministically
		mfPprob[i] = 1.0; //There are no cross transitions so Q represents a movement towards the opposite Strength
		//sumQ +=		mfQprob[miCascadeIndex] ; //For Verification
	}
	//For Boundary effects last q is divided by (1-x) --Keep as 1-0.5 for clarity of the value X
	//As the Arrays are static and global the boundary condition is taken care of in the HandlePot/DEP functions

}


//Cascade Index Can be = Cascade Size - Because each strength state has Cascade Indexes from 1 to Casc.Size
void synapseSerialCascade::init(int piStatesCount, //Number of Cascade states n
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

		miCascadeIndex = miCascadeStartIndex = gsl_rng_uniform_int(mprng,miCascadeSize)+1;

	}

	assert(miCascadeIndex >= 0 && miCascadeIndex <= miCascadeSize);
	reset();
	///Finished Initialisation
}

//Induces an LTP event depending on current state either switches cascades or moves down
int8_t synapseSerialCascade::handlePOT()
{
	int8_t Ds = 0;
	float fQProb,fPProb;

	if (miCascadeIndex < (miCascadeSize)) //Not Terminal State
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
	double lfp = 0.5 ; //gsl_rng_uniform(mprng);
    //Switch operational In POT Mode, so potentiate
    ///Check Appropriate Transition Probability
    if (menumStrength_State == SYN_STRENGTH_WEAK)
    {
		//Check Q for transition to Opposite Cascade
		if (fQProb > lfp)
		{
			//Ds = SYN_STRENGTH_STRONG-SYN_STRENGTH_WEAK;

			mbCascadeIndexChanged = true;
			miCascadeIndex--;

			Ds = -1; //Signal a switch Change
			if (miCascadeIndex==0) //If reached the Top Then We switch from Weak to STRONG
			{
				miCascadeIndex = 1; //1st state in strong
				menumStrength_State = SYN_STRENGTH_STRONG;
				mbStrengthChanged = true;
			}
		}
    }
    else
    { //Already in STRONG cascade - Move down in Cascade

    	//Check p for transition Down the Cascade
		if (fPProb > lfp)
		{
			//menumStrength_State = SYN_STRENGTH_STRONG; //Same State
			if (miCascadeSize > miCascadeIndex)//Check for end of Cascade condition
			{
					miCascadeIndex++;//Move down Cascade
					mbCascadeIndexChanged = true;
			}
			Ds = 1; //Signal A segregation change
		}


    }

return Ds;
}

//Induces an LTD event depending on current state either switches cascades or moves down
int8_t synapseSerialCascade::handleDEP()
{
	int8_t Ds = 0;
	float fQProb,fPProb;

	if (miCascadeIndex < (miCascadeSize)) //Not Terminal State
	{
		fQProb = mfQprob[miCascadeIndex];
		fPProb = mfPprob[miCascadeIndex];
	}
	else
	{
		fQProb = mfQprob[miCascadeIndex];
		fPProb = 0;
	}
//Draw random p
	double lfp = 0.5;//gsl_rng_uniform(mprng);


    //Switch operational In POT Mode, so potentate
    ///Check Appropriate Transition Probability
    if (menumStrength_State == SYN_STRENGTH_STRONG)
    {
		//Check Q for transition to Opposite Cascade
		if (fQProb > lfp)
		{
			mbCascadeIndexChanged = true;
			miCascadeIndex--;
			//Ds = SYN_STRENGTH_WEAK-SYN_STRENGTH_STRONG;
			if (miCascadeIndex == 0)
			{
				Ds = -1;
				menumStrength_State = SYN_STRENGTH_WEAK;
				mbStrengthChanged = true;
				miCascadeIndex = 1;
			}
		}
    }
    else
    { //Already in WEAK cascade - Move down in Cascade

    	//Check p for transition Down the Cascade
		if (fPProb > lfp)
		{
			//menumStrength_State = SYN_STRENGTH_STRONG; //Same State
			if (miCascadeSize > miCascadeIndex) //Check for end of Cascade (fP will be 0 anyway)
			{
				miCascadeIndex++;//Move down Cascade
				mbCascadeIndexChanged = true;
				Ds = 1; //Segregation signal
			}
		}

    }

    return Ds;

}

//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
void synapseSerialCascade::handleNOP()
{
////Nothing to see here
}


//Overloaded Methods
int synapseSerialCascade::getCascadeIndex() const
{
	return miCascadeIndex;
}
int synapseSerialCascade::getCascadeSize() const
{
return miCascadeSize;
}
int synapseSerialCascade::getStrength() const
{
return menumStrength_State;
}

//Returns whether a change of state has occured
int synapseSerialCascade::getStartStrength() const
{
return menumStartStrength_State;
}
//Returns whether a change of state has occured
bool synapseSerialCascade::hasStrengthModified() const
{

	return (menumStrength_State != menumStartStrength_State); //It appears to be so for the '05 Paper
	//return mbStrengthChanged;
}

// Returns true if Index has Been Incremented
bool synapseSerialCascade::hasIndexModified() const
{
	return mbCascadeIndexChanged;
}


//Resets Flag which holds if switchSynapse has changed.
void synapseSerialCascade::startMonitoring()
{
	mbStrengthChanged = false;
	mbCascadeIndexChanged = false;
	menumStartStrength_State = menumStrength_State; //Reset
	mbIsMonitored = true;
}
void synapseSerialCascade::stopMonitoring()
{
	mbIsMonitored = false;
}

bool synapseSerialCascade::isMonitored() const
{
	return mbIsMonitored;
}

void synapseSerialCascade::reset()
{
	menumStrength_State = menumStartStrength_State;
	mbStrengthChanged = false;
	mbCascadeIndexChanged = false;
	mbNoPlasticity = false;
	miCascadeIndex =  miCascadeStartIndex;
	mbIsMonitored = false;
}

void synapseSerialCascade::getTypeAsString(char* buff)
{
	synapseSerialCascade::getTypeName(buff);
}

void synapseSerialCascade::getTypeName(char* buff)
{
	strcpy(buff,"_synapseSerialCascade");
}
//Filter Baggage -Inherited from the ICASCADESYnapse
int synapseSerialCascade::getLThres() const
{
	return 0;
}

int synapseSerialCascade::getHThres() const
{
	return 0;
}
//The Decay Currently running
double synapseSerialCascade::getDecay() const
{
	return 0;
}

synapseSerialCascade::~synapseSerialCascade() {
	//  Auto-generated destructor stub
}
