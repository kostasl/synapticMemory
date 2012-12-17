/*
 * synapseFilterUnified.cpp
 *
 * Threshold LookUp tables define a set for {q,p} thresholds for transitions
 *	The threshold used for p or q depends on the strength state of the synapse
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 *
 */
#include "common.h"
#include "synapseFilterDual.h"

//alpha = 0.5
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
/*
 * VALUES VALID ONLY FOR UP TO n=15 --- THE 16-20 Are artificial And have not been calculated.
 */

const int synapseFilterDual::miThreshold_r100[20][2] =
{ {1,1},  {2,2},
  {2,2},  {3,3},
  {3,3},  {4,4},
  {4,4},  {4,4},
  {4,4},  {4,4}, //n=10
  {4,4},  {4,4},
  {4,4},  {4,4},
  {4,4}};
 /*//Matched to T.E values given by email {1,2,2,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4};*/

 //TODO: UNIFIED FILTER WITH DECAY: NEED TO MAKE THE REFLECTING VALUES FOR R=0.1 and r=0.01
const int synapseFilterDual::miThresholdTerminal_r100[20][2] =
{	  {1,1}, {1,1}, //1st and 2nd Terminal State matches the n-1
	  {2,2}, {2,2},
	  {3,3}, {3,3},
	  {4,4}, {4,4},
	  {4,4}, {4,4},
	  {4,4}, {4,4}, //n=10
	  {4,4}, {4,4},
	  {4,4} }; //n=15
 /*T.E : {1,1,2,2,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4};*/

//Single Decay Value for both p&q filters
const double synapseFilterDual::mdDecay_r100[20][2] =
		  {
				  {0.0000000000,0.0000000000 },
				  {0.0000000000,0.0000000000},
				  {0.6369496463,0.6369496463},
				  {0.3555279044,0.3555279044},
				  {0.6021765555,0.6021765555},
				  {0.4029658475,0.4029658475},
				  {0.5068588714,0.5068588714}, //n=7
				  {0.6072944823,0.6072944823},
				  {0.7072748450,0.7072748450},
				  {0.8079752083,0.8079752083}, //n=10
				  {0.9098652084,0.9098652084}, //n=11
				  {1.0130996579,1.0130996579},
				  {1.1176817270,1.1176817270},
				  {1.2235406356,1.2235406356},
				  {1.3305718238,1.3305718238} //n=15
		  };

const double synapseFilterDual::mdDecayTerminal_r100[20][2] =
		  {
				  {0.0000000000,0.0000000000}, //Same As non -terminal but shifted 1 >>
				  {0.0000000000,0.0000000000},
				  {0.0000000000,0.0000000000},
				  {0.6931471806,0.6931471806},
				  {0.38515754842,0.38515754842},//n=6
				  {0.6086804892,0.6086804892},
				  {0.4060507262,0.4060507262},
				  {0.5077448731,0.5077448731}, //n=8
				  {0.6075621254,0.6075621254},
				  {0.70735835949,0.70735835949},
				  {0.8080018018,0.8080018018}, //n=11
				  {0.9098737784,0.9098737784}, //Terminal State n=12
				  {1.0131024374,1.0131024374},
				  {1.1176826309,1.1176826309},
				  {1.2235409297,1.2235409297} //n=116
		  };

///RE-Injection At ZERO - PDFs for each of the filters (Assumes matching Parameters for p&q filters)
const double synapseFilterDual::mdPDF_r100[][4] = {
		   {1.000000,0.000000,0.000000,0.000000},
		   {0.599754,0.400246,0.000000,0.000000},
		   {0.527021,0.472979,0.000000,0.000000},
		   {0.333482,0.412148,0.254370,0.000000}, //n=4
		   {0.330825,0.461110,0.208065,0.000000}, //n=5
		   {0.211369,0.389231,0.294212,0.105188},
		   {0.236352,0.425755,0.266321,0.071572},
		   {0.267351,0.451809,0.232556,0.048284},
		   {0.296800,0.467404,0.203341,0.032455},
		   {0.325199,0.477696,0.175126,0.021979}, //n=10
		   {0.348149,0.485032,0.151700,0.015118},
		   {0.367877,0.490223,0.131665,0.010236},
		   {0.385957,0.493280,0.113807,0.006955},
		   {0.400556,0.495453,0.099179,0.004813},
		   {0.413272,0.496758,0.086718,0.003252} //n=15
};


const double synapseFilterDual::mdPDFTerminal_r100[][4]=
{
	   {1.000000,0.000000,0.000000,0.000000},
	   {1.000000,0.000000,0.000000,0.000000},
	   {0.500398,0.499602,0.000000,0.000000},
	   {0.499645,0.500355,0.000000,0.000000},
	   {0.292058,0.437918,0.270024,0.000000},
	   {0.320513,0.468204,0.211283,0.000000},
	   {0.200237,0.394079,0.299931,0.105754},
	   {0.232919,0.428285,0.267106,0.071690},
	   {0.265864,0.451856,0.233608,0.048672},
	   {0.298212,0.466363,0.202994,0.032430},//n=10
	   {0.323904,0.478648,0.175479,0.021969},
	   {0.348132,0.485024,0.151834,0.015010},//n=12
	   {0.368276,0.490410,0.131190,0.010124},
	   {0.385858,0.492938,0.114224,0.006980},
	   {0.401359,0.494360,0.099398,0.004884} //n=15
};

//alpha = 0.05
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
//NOT CALCULATED - MISSING
const int synapseFilterDual::miThreshold_r010[20][2] ={ {1,1},{2,1},//n=2
														{2,2},{2,2},//n=4
														{4,4},{4,4}, //n=6
														{4,4},{4,4},//n=8
														{4,4},{4,4}, //n=10
														{4,4},{4,4}, //n=12
														{4,4},{4,4}, //n=14
														{4,4} }; ////n=15

 //TODO MISSING:NOT CALCULATED - MAKE THE REFLECTING VALUES FOR R=0.1
const int synapseFilterDual::miThresholdTerminal_r010[20][2] ={ {0,1},{0,1},
  												{11,1},{3,1},
  												{2,2},{3,3}, //n=6
  												{2,2},{3,3},
  												{3,3},{3,3}, //n=10
  												{3,3},{3,3}, //n=12
  												{3,3},{3,3}, //n=14
  												{3,3} }; ////n=15




//alpha = 0.005
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
//NOT CALCULATED - MISSING
const int synapseFilterDual::miThreshold_r001[20][2] ={  {1,1},{2,1}, //n=2
														{2,2},{2,2}, //n=4
														{2,2},{2,2}, //n=6
														{2,2},{2,2}, //n=8
														{3,3},{4,4}, //n=10
														{3,3},{3,3}, //n=12
														{4,4},{4,4}, //n=14
														{4,4} };	 //n=15

 //TODO MISSING:NOT CALCULATED - MAKE THE REFLECTING VALUES FOR R=0.1
const int synapseFilterDual::miThresholdTerminal_r001[20][2] ={ {0,1},{0,1},
   												{0,0},{0,0},
   												{0,0},{0,0}, //n=6
   												{2,2},{3,3},
   												{3,3},{3,3}, //n=10
   												{3,3},{3,3}, //n=12
   												{3,3},{3,3}, //n=14
   												{3,3} }; ////n=15




//Set Up the Array with the sets indexed by rate r
//This is a pointer to an array - Not An array of pointers
const PtrThresholdSet synapseFilterDual::miThresholdsSet[] = {&synapseFilterDual::miThreshold_r100,
															  &synapseFilterDual::miThreshold_r010,
															  &synapseFilterDual::miThreshold_r001};

//This is a pointer to an array of thresholds- Not An array of pointers
const PtrThresholdSet synapseFilterDual::miTerminalThresholdsSet[] = {&synapseFilterDual::miThresholdTerminal_r100,
																      &synapseFilterDual::miThresholdTerminal_r010,
																	  &synapseFilterDual::miThresholdTerminal_r001};


//Set Up the Array with the sets indexed by rate r For Decay and Terminal States Decay
const PtrDecaySet synapseFilterDual::mdDecaysSet[] = {&synapseFilterDual::mdDecay_r100,
													  &synapseFilterDual::mdDecay_r100, //Should be mdDecay_r100 but they do not exist yet
													  &synapseFilterDual::mdDecay_r100};
//
////This is a pointer to an array of thresholds- Not An array of pointers
const PtrDecaySet synapseFilterDual::mdTermDecaySet[] = {&synapseFilterDual::mdDecayTerminal_r100,
														 &synapseFilterDual::mdDecayTerminal_r100,
														 &synapseFilterDual::mdDecayTerminal_r100};

//Default ConsTructor - Never Actually Used
synapseFilterDual::synapseFilterDual() {

	//RAND INDEX
	mprng  = 0; //Init Will Handle It
	miCascadeIndex = 0;
	init( DEFAULT_CASCADE_SIZE,miCascadeIndex,SYN_STRENGTH_NOTSET,0);
	//Now Rng Is init From Global INstance

	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(miTerminalIndex*r);//Take Cascade Index as concat double to give 0-Count

	mbIsMonitored = false;

}

//Init With StartStrength Used by SynapticConnection
synapseFilterDual::synapseFilterDual(int piCascadeSize,float startStrength,gsl_rng * rng_r)
{

	mprng = rng_r;
	miCascadeIndex = 0;
	ICascadeSynapse::SYN_STRENGTH_STATE enumStartStrength;
	if (startStrength>0) //Convert Float to Enum
		enumStartStrength = SYN_STRENGTH_STRONG;
	else
		enumStartStrength = SYN_STRENGTH_WEAK;

	init( piCascadeSize,miCascadeIndex,enumStartStrength,0);

	//RAND INDEX
	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(piCascadeSize*r);//Take Cascade Index as concat double to give 0-Count


	mbIsMonitored = false;

}

synapseFilterDual::synapseFilterDual(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r)
{
	miCascadeIndex = 0;
	mprng = rng_r;
	init( piCascadeSize,miCascadeIndex,penumStartStrength,0);

	//RAND INDEX
	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(piCascadeSize*r);//Take Cascade Index as concat double to give 0-Count

}

//Initialise with cascadeSize
synapseFilterDual::synapseFilterDual(int piCascadeSize,gsl_rng * rng_r)
{

	mprng = rng_r;
	miTerminalIndex = piCascadeSize - 1;
	//RAND INDEX
	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(piCascadeSize*r);//Take Cascade Index as concat double to give 0-Count

	init( piCascadeSize,miCascadeIndex,SYN_STRENGTH_NOTSET,0);

}

//Used to Start From A fixed point and test how distribution evolves
synapseFilterDual::synapseFilterDual(int piCascadeSize,int piStartIndex,gsl_rng * rng_r)
{
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,SYN_STRENGTH_NOTSET,0);
}

synapseFilterDual::synapseFilterDual(int piCascadeSize,int piStartIndex,
		ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
		gsl_rng *  rng_r,
		int iRateDependentParameterSet) //Default Value is rate=1.0

{
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,penumStartStrength,iRateDependentParameterSet);
}


synapseFilterDual::synapseFilterDual(int piCascadeSize,int piStartIndex,
					ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
					int iStartFilterState,
					gsl_rng *  rng_r,
					int iRateDependentParameterSet) //Default Value is rate=1.0
{
	mprng = rng_r;
	init( piCascadeSize,piStartIndex,penumStartStrength,iRateDependentParameterSet);
	miRPFilterValue= miRDFilterValue = iStartFilterState; //Some simple init
	assert (miRPFilterValue >= miLThres && miRPFilterValue <= miHThres);
}


void synapseFilterDual::init(int piCascadeSize,int iThresholdsIndex,
								SYN_STRENGTH_STATE enumCurrStrengthState,int ParameterSetIndex)
{
	assert(iThresholdsIndex < piCascadeSize);

	miParameterSetUsed = ParameterSetIndex;
	miCascadeIndex = miStartIndex = iThresholdsIndex;
	miTerminalIndex = piCascadeSize-1;
//	mbStrengthChanged =false;
//	mbCascadeIndexChanged = false; //When A switchSynapse Cascade Index is increased
	mbNoPlasticity = false;
	//mbInternallyManagedGSL = false;
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

	setFilterThresholds();
	initialiseFilterState();
}

int synapseFilterDual::getCascadeIndex() const
{

	return miCascadeIndex;
}

int synapseFilterDual::getCascadeSize() const
{
	return (miTerminalIndex+1);
}

bool synapseFilterDual::hasIndexModified() const
{
	//return mbCascadeIndexChanged;
	return (miCascadeIndex	!= miStartIndex);
}


bool synapseFilterDual::hasStrengthModified() const
{
	return (penumStrength != penumStartStrength);
}

//Induces an LTP event depending on current state either switches cascades or moves down
int8_t synapseFilterDual::handlePOT()
{
	return addSample(+1);
}

//Induces an LTD event depending on current state either switches cascades or moves down
int8_t synapseFilterDual::handleDEP()
{
	return addSample(-1);
}

//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
void synapseFilterDual::handleNOP()
{
	addSample(0);
}


void synapseFilterDual::doStochasticDecay()
{
	//Do Decay on both running sums
	unsigned int rDecaySteps;
	double dPDecayRate,dDDecayRate;

	if (miCascadeIndex != miTerminalIndex)
	{
		if (penumStrength == SYN_STRENGTH_STRONG){
			dPDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
			dDDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][1];
		}
		else ///When Symmetric it doesnt really matter - But now different decay for p&q can be set
		{
			dPDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][1];
			dDDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
		}
	}
	else//Terminal INdex
	{
		if (penumStrength == SYN_STRENGTH_STRONG){
			dPDecayRate = (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][0]; //P is 1st values
			dDDecayRate = (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][1]; //Q is 2nd Values
		}else {
			dPDecayRate = (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][1];
			dDDecayRate = (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][0];
		}
	}

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
int synapseFilterDual::addSample(int iValue)
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
	return iRet;
}

int synapseFilterDual::lThresReached()
{
	int Ret = 0;
	//Add new sample to the correct internal sum
	if (penumStrength == SYN_STRENGTH_STRONG)
	{ //Then lThres is a q transition
		switchReset();
		Ret = -1;
	}
	else //P - Transitions
	{
		if (miCascadeIndex < miTerminalIndex)
		{
			miCascadeIndex++;
			setFilterThresholds();
			reInjectFilterStateToCascadeState();
			Ret = 1;
			//mbCascadeIndexChanged = true;
		}
		else //REFLECT - Terminal State
			//Reached Terminal Then Just reflect
			miRDFilterValue = miLThres; //Not required as we can ignore this filters state
	}

	return Ret;
}

int synapseFilterDual::hThresReached()
{
	int Ret = 0;
	//Add new sample to the correct internal sum
	if (penumStrength == SYN_STRENGTH_WEAK)
	{ //Then lThres is a q transition
		switchReset();
		Ret = -1;
	}else
	{//This is a p transition
		if (miCascadeIndex < miTerminalIndex)
		{
			miCascadeIndex++;
			setFilterThresholds();
			reInjectFilterStateToCascadeState();
			//mbCascadeIndexChanged = true;
			Ret = 1;
		}
		else //REFLECT - Terminal State
			//Reached Terminal Then Just reflect
			miRPFilterValue = miHThres; //Not required as we can ignore this filters state
	}

	return Ret;
}

int synapseFilterDual::getRunningValueH() const
{
	return miRPFilterValue;
}


int synapseFilterDual::getRunningValueL() const
{
	return miRDFilterValue;
}

int synapseFilterDual::getStrength() const
{
	return (int)penumStrength;
}

int synapseFilterDual:: getStartStrength() const
{
return (int)penumStartStrength;
}

int synapseFilterDual::getLThres() const
{
	return miLThres;
}
int synapseFilterDual::getHThres() const
{
	return miHThres;
}

double synapseFilterDual::getHDecay() const
{
	return 0.0;
}
double synapseFilterDual::getLDecay() const{
	return 0.0;
}

double synapseFilterDual::getDecay() const{

	if (miCascadeIndex != miTerminalIndex)
	{
		return (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
	}
	else
	{
		return (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][0];
	}
}

//Returns the number of times the Filter Is allowed to Grow Threshold (The Max Threshold Index)
int synapseFilterDual::getMaxSteps() const
{
	return miTerminalIndex;
}

int synapseFilterDual::getThresholdIndex() const
{
  return miCascadeIndex;
}

//Save Current State As Default
void synapseFilterDual::startMonitoring()
{
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	miStartIndex = miCascadeIndex;
	penumStartStrength = penumStrength;
	mbIsMonitored = true;
}

void synapseFilterDual::stopMonitoring()
{
	mbIsMonitored = false;
}

bool synapseFilterDual::isMonitored() const
{
	return mbIsMonitored;
}

//Method called by testEscapeTime
void synapseFilterDual::reset()
{

	ICascadeSynapse::reset(); //Randomizes Start Strength, Unfreezeplasticity

	//super::reset();//If we Re-init from Non zero position Then The escape time through any boundary Reduces
		//initialiseFilterState(); //Reset running Sum - Do not - Let it be as it has been Inited by previous experience-Not the case if starting over fixed position
	miStartIndex = miCascadeIndex; //Reset Index

	penumStrength = (SYN_STRENGTH_STATE)penumStartStrength;

	miCascadeIndex	= miStartIndex;
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	mbNoPlasticity = false;
	mbIsMonitored = false;

	//setFilterThresholds();
	//reInjectFilterStateToCascadeState();
	//initialiseFilterState(); //Like Starting Over - BUT If we Re-init from Non zero position Then The escape time through any boundary Reduces
	 //The reset is only called by the testEscapetime Function - And thus it affects the result of this function
}

void synapseFilterDual::switchReset()
{
	if (penumStrength == SYN_STRENGTH_STRONG){
		penumStrength = SYN_STRENGTH_WEAK;
	}else{
		penumStrength = SYN_STRENGTH_STRONG;
		}

		miCascadeIndex = 0;

//		mbCascadeIndexChanged = true;
//		mbStrengthChanged 	= true;
		setFilterThresholds();
		reInjectFilterStateToCascadeState(); //Reset To Zero On New State
}


//Use Index / strength and the Lookup arrays to set new L/H thresholds
//Called after changes to cascade index
void synapseFilterDual::setFilterThresholds()
{
	if (penumStrength == SYN_STRENGTH_STRONG)
	{
		if (miCascadeIndex < miTerminalIndex)
		{
			miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0]; //P
			miLThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];//Q
		}
		else //Terminal States
		{
			miHThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][0]; //P
			miLThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
		}
//WEAK SYNAPSE
	}else{
		if (miCascadeIndex < miTerminalIndex)
		{
			miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1]; //P
			miLThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0]; //Q
		}
		else //Terminal States
		{
			miHThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
			miLThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
		}

	} //ENDOF IF STRONG

}

//Sets where the running sum should be initialized after every change to Cascade state -
void synapseFilterDual::reInjectFilterStateToCascadeState()
{
	//Inject Back to 0
	miRPFilterValue = 0;
	miRDFilterValue = 0;
}

////Sets where the running sum should be initialise When the Object is first created
void synapseFilterDual::initialiseFilterState()
{
	//Use injection PDF - Here bothe Running Sums Are +Ve and so are the thresholds
	//They act as counters of +ve and -ve stimuli separetely
	double p = 0.0;
	int i; //Iterator
	const double (*dPDFUsed)[15][ciMaxInternalStates]; //Pointer to PDF. If at terminal State then Use The reflecting Boundary PDF

	mbIsMonitored = false;

	//return;
	if (((miTerminalIndex) < 1)) return; //Not For 1,2 thresholds filter

	//TODO: Add PDF Sets For All Rates -
	if (miCascadeIndex < miTerminalIndex )
		dPDFUsed = (&mdPDF_r100);
	else
		dPDFUsed = (&mdPDFTerminal_r100);

	//DO P-Filter State
	double r = gsl_rng_uniform(mprng);
	for ( i = 0; i< ciMaxInternalStates;i++)
	{
		p += (*dPDFUsed)[miCascadeIndex][i]; //Accumulate the Pdf
		if (p > r){
			//Found the spot since r was just exceeded
			miRPFilterValue = i; //Remove Offset so i=0 becomes state -3 floor(ciMaxInternalStates/2)
			break;
		}
	}

	//Do-D Filter State Init - Same As Above
	r = gsl_rng_uniform(mprng); //For the other filter now
	for ( i = 0; i< ciMaxInternalStates;i++)
	{
		p += (*dPDFUsed)[miCascadeIndex][i]; //Accumulate the Pdf
		if (p > r){
			//Found the spot since r was just exceeded
			miRDFilterValue = i; //Remove Offset so i=0 becomes state -3 floor(ciMaxInternalStates/2)
			break;
		}
	}

	//Catch stupid Errors
	if(( (miRDFilterValue >= miLThres) || (miRPFilterValue >= miHThres) ))
	{
		assert((miRDFilterValue < miLThres) && ((miRPFilterValue < miHThres)));
	}

}
void  synapseFilterDual::getTypeAsString(char* buff)
{
	getTypeName(buff);
}


void synapseFilterDual::getTypeName(char* buff)
{
	strcpy(buff,"_synapseCascadeDualFilter");
}

synapseFilterDual::~synapseFilterDual() {
	//  Auto-generated destructor stub
}
