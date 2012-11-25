/*
 * synapseFilterUnified.cpp
 *
 * Single Unified Filter No Decay. Can be initialized to a set of thresholds for p/q that match one of the cascades states.
 * Threshold LookUp tables define a set for {q,p} thresholds for transitions
 *	The threshold used for p or q depends on the strength state of the synapse
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 *
 */
#include "common.h"
#include "synapseCascadeFilterUnified.h"

//alpha = 0.5
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
const int synapseCascadeFilterUnified::miThreshold_r100[][2] =
{
  {6,6},   {1,2},		//{1,1},   {1,2},
  {2,2},    {2,4},
  {4,4},    {4,8},
  {8,8},    {8,16},
  {16,16},  {16,32}, //n=10
  {32,32},  {32,64},
  {64,64},  {64,128},
  {128,128},{128,128},{128,256},{256,256},{256,512},{512,512}};
 /*,{128,256}, n=16  {256,256},{256,512},  {512,512}}; //n=19*/

 //TODO: UNIFIED FILTER WITH DECAY: NEED TO MAKE THE REFLECTING VALUES FOR R=0.1 and r=0.01
const int synapseCascadeFilterUnified::miThresholdTerminal_r100[][2] =
{{8,8},     {8,8},
 {2,1},     {4,1},
 {8,1},     {16,1},
 {32,1},    {64,1},
 {128,1},   {256,1},
 {512,1},   {1024,1},
 {2048,1},  {4096,1},
 {8192,1},{16384,1},{32768,1},{65536,1},{131072,1},{262144,1}}; //n=15
 /*{1,16384}, {1,32768}, {1,65536}, {1,131072},{1,262144}};*/

/* My previous Values
{ {0,1},{0,1},
{11,1},{3,1},
{2,2},{3,3}, //n=6
{2,2},{3,3},
{3,3},{3,3}, //n=10
{3,3},{3,3}, //n=12
{3,3},{3,3}, //n=14
{3,3} }; ////n=15
*/

//alpha = 0.05
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter

const int synapseCascadeFilterUnified::miThreshold_r010[][2] ={ {1,1},{2,1},//n=2
														{2,2},{2,2},//n=4
														{4,4},{4,4}, //n=6
														{4,4},{4,4},//n=8
														{4,4},{4,4}, //n=10
														{4,4},{4,4}, //n=12
														{4,4},{4,4}, //n=14
														{4,4},{4,4},{4,4},{4,4},{4,4},{4,4} }; ////n=15

 //TODO MISSING:NOT CALCULATED - MAKE THE REFLECTING VALUES FOR R=0.1
const int synapseCascadeFilterUnified::miThresholdTerminal_r010[][2] ={ {0,1},{0,1},
  												{11,1},{3,1},
  												{2,2},{3,3}, //n=6
  												{2,2},{3,3},
  												{3,3},{3,3}, //n=10
  												{3,3},{3,3}, //n=12
  												{3,3},{3,3}, //n=14
  												{3,3},{3,3},{3,3},{3,3},{3,3},{3,3} }; ////n=15

//alpha = 0.005
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
const int synapseCascadeFilterUnified::miThreshold_r001[][2] ={  {1,1},{2,1}, //n=2
														{2,2},{2,2}, //n=4
														{2,2},{2,2}, //n=6
														{2,2},{2,2}, //n=8
														{3,3},{4,4}, //n=10
														{3,3},{3,3}, //n=12
														{4,4},{4,4}, //n=14
														{4,4},{4,4},{4,4},{4,4},{4,4},{4,4} };	 //n=15

 //TODO MISSING:NOT CALCULATED - MAKE THE REFLECTING VALUES FOR R=0.1
const int synapseCascadeFilterUnified::miThresholdTerminal_r001[][2] ={ {0,1},{0,1},
   												{11,1},{3,1},
   												{2,2},{3,3}, //n=6
   												{2,2},{3,3},
   												{3,3},{3,3}, //n=10
   												{3,3},{3,3}, //n=12
   												{3,3},{3,3}, //n=14
   												{3,3},{3,3},{3,3},{3,3},{3,3},{3,3} }; ////n=15


//Set Up the Array with the sets indexed by rate r
//This is a pointer to an array - Not An array of pointers
const PtrThresholdSet synapseCascadeFilterUnified::miThresholdsSet[] = {&synapseCascadeFilterUnified::miThreshold_r100,
																 &synapseCascadeFilterUnified::miThreshold_r010,
																 &synapseCascadeFilterUnified::miThreshold_r001};

//This is a pointer to an array of thresholds- Not An array of pointers
const PtrThresholdSet synapseCascadeFilterUnified::miTerminalThresholdsSet[] = {&synapseCascadeFilterUnified::miThresholdTerminal_r100,
																		 &synapseCascadeFilterUnified::miThresholdTerminal_r010,
																		 &synapseCascadeFilterUnified::miThresholdTerminal_r001};

//Default ConsTructor - Never Actually Used
synapseCascadeFilterUnified::synapseCascadeFilterUnified() {

//	//RAND INDEX
//	mprng  = 0; //Init Will Handle It
//	miCascadeIndex = 0;
//	init( DEFAULT_CASCADE_SIZE,miCascadeIndex,SYN_STRENGTH_NOTSET,0);
//	//Now Rng Is init From Global INstance
//	double r = gsl_rng_uniform(mprng);
//	miCascadeIndex = miStartIndex = round(miTerminalIndex*r);
//	mbIsMonitored = false;
	uiStateLifetime = 0;
	miRFilterValue = 0;
}

//Init With StartStrength Used by SynapticConnection
synapseCascadeFilterUnified::synapseCascadeFilterUnified(int piCascadeSize,float startStrength,gsl_rng * rng_r)
{

	mprng = rng_r;
	miCascadeIndex = 0;
	ICascadeSynapse::SYN_STRENGTH_STATE enumStartStrength;
	if (startStrength>0) //Convert Float to Enum  ....
		enumStartStrength = SYN_STRENGTH_STRONG;
	else
		enumStartStrength = SYN_STRENGTH_WEAK;

	init( piCascadeSize,miCascadeIndex,enumStartStrength,0);

	//RAND INDEX
	double r = gsl_rng_uniform(mprng);
	miCascadeIndex = miStartIndex = round(miTerminalIndex*r);
	mbIsMonitored = false;

}

synapseCascadeFilterUnified::synapseCascadeFilterUnified(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r)
{
	miCascadeIndex = 0;
	mprng = rng_r;
	init( piCascadeSize,miCascadeIndex,penumStartStrength,0);

	//RAND INDEX
	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(piCascadeSize*r);

}

//Initialise with cascadeSize //Do Not INIT so Super Class Has Faster Init Period
synapseCascadeFilterUnified::synapseCascadeFilterUnified(int piCascadeSize,gsl_rng * rng_r)
{


	mprng = rng_r;
	miTerminalIndex = piCascadeSize - 1;
	//RAND INDEX
	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(piCascadeSize*r);//Take Cascade Index as concat double to give 0-Count

	init( piCascadeSize,miCascadeIndex,SYN_STRENGTH_NOTSET,0);

}

//Used to Start From A fixed point and test how distribution evolves
synapseCascadeFilterUnified::synapseCascadeFilterUnified(int piCascadeSize,int piStartIndex,gsl_rng * rng_r)
{
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,SYN_STRENGTH_NOTSET,0);
}

synapseCascadeFilterUnified::synapseCascadeFilterUnified(int piCascadeSize,int piStartIndex,
		ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
		gsl_rng *  rng_r,
		int iRateDependentParameterSet) //Default Value is rate=1.0

{
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,penumStartStrength,iRateDependentParameterSet);
}


synapseCascadeFilterUnified::synapseCascadeFilterUnified(int piCascadeSize,int piStartIndex,
					ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
					int iStartFilterState,
					gsl_rng *  rng_r,
					int iRateDependentParameterSet) //Default Value is rate=1.0
{
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,penumStartStrength,iRateDependentParameterSet);
	miRFilterValue = iStartFilterState;
	assert (miRFilterValue >= miLThres && miRFilterValue <= miHThres);
}


void synapseCascadeFilterUnified::init(int piCascadeSize,int iThresholdsIndex,
								SYN_STRENGTH_STATE enumCurrStrengthState,int ParameterSetIndex)
{

	uiStateLifetime = 0;//Reset The Strength State Allocation Counter
	miParameterSetUsed = ParameterSetIndex;
	miCascadeIndex = miStartIndex = iThresholdsIndex;
	miTerminalIndex = piCascadeSize-1;

	assert(miCascadeIndex <= miTerminalIndex);
	//mbStrengthChanged =false;
	//mbCascadeIndexChanged = false; //When A switchSynapse Cascade Index is increased
	mbNoPlasticity = false;
	mbInternallyManagedGSL = false;
	mbIsMonitored = false;

	//g_rng_r = getRandGeneratorInstance();
	if (!mprng)
		mprng = g_getRandGeneratorInstance(false);

	if (enumCurrStrengthState == SYN_STRENGTH_NOTSET)
	{
		double r = gsl_rng_uniform(mprng);
		if (r<0.5)
			penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
		else
			penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;
	}else
	penumStartStrength = penumStrength = (uint8_t)enumCurrStrengthState;

	setFilterThresholds(); //These Calls do not Call the Virtual Functions But rather this functions implementations
	initialiseFilterState();
}

int synapseCascadeFilterUnified::getCascadeIndex() const
{

	return miCascadeIndex;
}

int synapseCascadeFilterUnified::getCascadeSize() const
{
	return (miTerminalIndex+1);
}

bool synapseCascadeFilterUnified::hasIndexModified() const
{
	return (miStartIndex !=miCascadeIndex);
}


bool synapseCascadeFilterUnified::hasStrengthModified() const
{
	return (penumStrength != penumStartStrength);
}

//Induces an LTP event depending on current state either switches cascades or moves down
int8_t synapseCascadeFilterUnified::handlePOT()
{
	return addSample(+1);
}

//Induces an LTD event depending on current state either switches cascades or moves down
int8_t synapseCascadeFilterUnified::handleDEP()
{
	return addSample(-1);
}

//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
void synapseCascadeFilterUnified::handleNOP()
{
	addSample(0);
}


///THE IFILTER INTERFACE

///@Returns -1 For a Plastic Transition, +1 For Metaplastic and 0 For No Transition
int synapseCascadeFilterUnified::addSample(int iValue)
{
	int Ret = 0;
	//NOP - No decay on this synapse so we don't count time
	if (iValue == 0) return 0;

	miRFilterValue +=iValue;
	uiStateLifetime+=1;

	 //Now Check Threshold Condition;
	if (miRFilterValue <= miLThres)
	{
		//p Threshold Reached
		if (penumStrength == SYN_STRENGTH_WEAK)
		{
			if (miCascadeIndex < miTerminalIndex){
				miCascadeIndex++; //Do A p Transitions
				uiSameThresholdTransitionCounter++; //Metaplasticity Counter - Used for Statistics Or Allocation
				//mbCascadeIndexChanged = true;
				Ret = 1; //Signal METAplastic Transition Occurred
				setFilterThresholds(); //Set New Thresholds
				reInjectFilterStateToCascadeState(); //Reset To Zero On New State
			}else{//Terminal State //No p Transitions
				Ret = 0;//No Transition
				//miRFilterValue = miLThres; //Holding Barrier
				miRFilterValue -=iValue; //Reflecting Barrier
			}
		}else{//q thres Reached
			switchReset(); //Change to STRONG And Reset Cascade Index - Reset Runnning Sum
			Ret = -1; //Plastic Transition Occured
			}
	}else // Both Transitions Cannot Occur
		if (miRFilterValue >= miHThres)
		{ //lOW tRHES
			if (penumStrength == SYN_STRENGTH_STRONG)
			{
				if (miCascadeIndex < miTerminalIndex){
					miCascadeIndex++; //Do A p Transitions
					uiSameThresholdTransitionCounter++; //Metaplasticity Counter - Used for Statistics Or Allocation
					//mbCascadeIndexChanged = true;
					setFilterThresholds(); //Set New Thresholds On this Index
					reInjectFilterStateToCascadeState(); //Reset To Zero On New State
					Ret = 1; //METAplastic Transition Occurred
				}else{//Terminal State //No p Transitions
					Ret = 0;//No Transition
					//miRFilterValue = miHThres; //Holding Barrier
					miRFilterValue -=iValue; //Reflecting Barrier
				}
			}else{ //WEAK SYNAPSE -> SWITCH TO STRONG
				switchReset(); //Change to weak And Reset Cascade Index Reset Runnning Sum
				Ret = -1; //Plastic Transition Occurred
				}
		}

	return Ret;
}

int synapseCascadeFilterUnified::getRunningValue() const
{
	return miRFilterValue;
}

int synapseCascadeFilterUnified::getStrength() const
{
	return (int)penumStrength;
}

int synapseCascadeFilterUnified:: getStartStrength() const
{
return (int)penumStartStrength;
}

int synapseCascadeFilterUnified::getLThres() const
{
	return miLThres;
}
int synapseCascadeFilterUnified::getHThres() const
{
	return miHThres;
}

double synapseCascadeFilterUnified::getHDecay() const
{
	return 0.0;
}
double synapseCascadeFilterUnified::getLDecay() const{
	return 0.0;
}

double synapseCascadeFilterUnified::getDecay() const{
	return 0.0;
}

//Returns the number of times the Filter Is allowed to Grow Threshold (The Max Threshold Index)
int synapseCascadeFilterUnified::getMaxSteps() const
{
	return miTerminalIndex;
}

int synapseCascadeFilterUnified::getThresholdIndex() const
{
  return miCascadeIndex;
}

//Save Current State As Default
void synapseCascadeFilterUnified::startMonitoring()
{
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	miStartIndex = miCascadeIndex;
	penumStartStrength = penumStrength;
	mbIsMonitored = true;

}

void synapseCascadeFilterUnified::stopMonitoring()
{
	mbIsMonitored = false;
}

bool synapseCascadeFilterUnified::isMonitored() const
{
	return mbIsMonitored;
}

//Method called by testEscapeTime and MLT Lifetime Tests
void synapseCascadeFilterUnified::reset()
{

	penumStrength = (SYN_STRENGTH_STATE)penumStartStrength;
	miCascadeIndex	= miStartIndex;
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	mbNoPlasticity = false;
	mbIsMonitored = false;
	uiStateLifetime = 0; //At Each Trial This Counter Needs to Be reset
	setFilterThresholds();
	reInjectFilterStateToCascadeState();
	//initialiseFilterState(); //Like Starting Over - Problem When Measuring Escape time BUT If we Re-init from Non zero position Then The escape time through any boundary Reduces
	 //The reset is only called by the testEscapetime Function - And thus it affects the result of this function
}

void synapseCascadeFilterUnified::switchReset()
{
	if (penumStrength == SYN_STRENGTH_STRONG){
		penumStrength = SYN_STRENGTH_WEAK;
	}else{
		penumStrength = SYN_STRENGTH_STRONG;
		}

		miCascadeIndex = 0;
		uiStateLifetime = 0;///Reset Lifetime Counter
		//mbCascadeIndexChanged = true;
		//mbStrengthChanged 	= true;
		setFilterThresholds();
		reInjectFilterStateToCascadeState(); //Reset To Zero On New State
}


//Use Index / strength and the Lookup arrays to set new L/H thresholds
//Called after changes to cascade index -This code Does not like Changes - Binary Is not Updating after modifications!
void synapseCascadeFilterUnified::setFilterThresholds()
{
	if (penumStrength == SYN_STRENGTH_STRONG)
	{
		if (miCascadeIndex < miTerminalIndex)
		{
			miLThres = -(*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1] + 0;
			miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
		}
		else //Terminal States
		{
			miLThres = -(*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
			miHThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
		}
//WEAK SYNAPSE   sss
	}else{
		if (miCascadeIndex < miTerminalIndex)
		{
			miLThres = -(*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
			miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
		}
		else //Terminal States
		{
			miLThres = -(*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
			miHThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
		}

	} //ENDOF IF STRONG

		//implement Symmetric Assymetric Idea By Randomly Switching Between thres
		//Should We Randomly Switch Between Thresholds ?
	if (miCascadeIndex != miTerminalIndex)
	{
	 	//g_rng_r = getRandGeneratorInstance(false);
		double r = gsl_rng_uniform(mprng);
		if (r < 0.5)
		{
			int temp = miLThres;
			miLThres = -miHThres;
			miHThres = -temp;
		}
	}

	assert(miLThres*miHThres < 0); //Sign Should mismatch
}

//Sets where the running sum should be initialized after every change to Cascade state -
void synapseCascadeFilterUnified::reInjectFilterStateToCascadeState()
{
	miRFilterValue = 0;
}

////Sets where the running sum should be initialise When the Object is first created
void synapseCascadeFilterUnified::initialiseFilterState()
{
	miRFilterValue = 0;
//	return;

	uint iUpper,iLower;
	iUpper = miHThres-1;
	iLower = -(miLThres +1); //Do not inject on to Absorbing Boundary!
	const int LENGTH = 1023;
	int i, zero = (LENGTH-1)/2;

	double cdf[LENGTH],norm[LENGTH];
	double p,total = 0;

	  //CHECKTERMINAL STATES
	if (miCascadeIndex == miTerminalIndex)
	{//Terminal State Draw Uniform
		p = gsl_rng_uniform(mprng);
		miRFilterValue = round(p*(iUpper+iLower+1))-iLower;
		return;
	}

	  memset(norm,0,sizeof(double)*LENGTH);
	  // First make the distribution

	  for(i=(int)(zero-iLower);i<=(int)(zero+iUpper);i++){
	    if (iLower == iUpper)
	    {
	      if (i<=zero) norm[i] = 1 + i - (zero - iLower);
	      else         norm[i] = 1 + (zero + iUpper) - i;
	    }
	    if (iLower < iUpper)
	    {
	      if (i<=zero) norm[i] = 2 * (1 + i - (zero - iLower));
	      else         norm[i] = 1 * (1 + (zero + iUpper) - i);
	    }
	    if (iLower > iUpper)
	    {
	      if (i<=zero) norm[i] = 1 * (1 + i - (zero - iLower));
	      else         norm[i] = 2 * (1 + (zero + iUpper) - i);
	    }
	    total += norm[i];
	  }

	  for(i=0;i<LENGTH;i++) norm[i] = norm[i] / total;

	  // Now we can draw from the distribution

	  cdf[0] = norm[0];

	  for(i=1;i<LENGTH;i++) cdf[i] = cdf[i-1] + norm[i];

	  p = gsl_rng_uniform(mprng);

	  for(i=0;i<LENGTH-1;i++) if (p < cdf[i]) break;

	  // Return the picked entry
	  miRFilterValue =i - zero;

	  assert(miRFilterValue < miHThres && miRFilterValue > miLThres);
	  return;

}
void  synapseCascadeFilterUnified::getTypeAsString(char* buff)
{
	getTypeName(buff); ///Doesn like changes
}


void synapseCascadeFilterUnified::getTypeName(char* buff)
{
	strcpy(buff,"_synapseCascadeFilterUnified");
}

synapseCascadeFilterUnified::~synapseCascadeFilterUnified() {
	//  Auto-generated destructor stub
}
