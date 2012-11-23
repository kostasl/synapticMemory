/*
 * ICascadeSynapse.cpp
 *
 *  Created on: 9 Mar 2011
 *      Author: kostasl
 */

#include "ICascadeSynapse.h"

ICascadeSynapse::ICascadeSynapse() {

	mbNoPlasticity = false;
	mbMetaplasticAlloc = false;
	mbPlasticAlloc		= false;
	mbStabilityAlloc =false;
	mbStopRecordingOfMHistogramAtNextThresholdEvent = false;
	mbSaveMetaplasticHistogram = false;
	iCycleSamplesRemaining = 0;
	uiStateLifetime = 0;

	setAllocationThreshold(0); //No refraction Period -- uiThresholdForAllocation = 0;
	uiSameThresholdTransitionCounter = 0;

	uiMaxMetaplasticTransitions = 0;
	mpMDistribinTime = NULL;
	penumTrackedStrength = penumStartStrength = penumStrength  = SYN_STRENGTH_NOTSET;

	mdAllocationDecayRate = 0.0;
}

///Constr. With Allocation Switches
ICascadeSynapse::ICascadeSynapse(bool PlasticAlloc, bool MetaplasticAlloc,bool bStabilityAlloc)
{
	mbStopRecordingOfMHistogramAtNextThresholdEvent = false;
	mbMetaplasticAlloc 	= MetaplasticAlloc;
	mbPlasticAlloc		= PlasticAlloc;
	mbStabilityAlloc	= bStabilityAlloc;
	mbNoPlasticity 		= false;
	uiStateLifetime 	= 0;
	uiThresholdForAllocation 		= 0;
	uiSameThresholdTransitionCounter 	= 0;
	mpMDistribinTime 				= NULL;
	penumTrackedStrength = penumStartStrength = penumStrength  =SYN_STRENGTH_NOTSET;
	iCycleSamplesRemaining = 0;
}

void ICascadeSynapse::freezePlasticity()
{
	if (mbStabilityAlloc && (uiStateLifetime > uiThresholdForAllocation)) //This criterion counts the number of memories stored without this synapse changing strength
		mbNoPlasticity = true;

	//If Metaplastic -- Then Compare Metaplastic Cycles against threshold
	if (mbMetaplasticAlloc)
	{
		if(uiSameThresholdTransitionCounter >= uiThresholdForAllocation) //This criterion counts the number of memories stored without this synapse changing strength
		{
			mbNoPlasticity = true;
			resetMetaplasticCounter(); //Once Allocation Period has expired we reset the counters--
			//Without this line allocation Freezes synapses and a metaplastic sample May never be reached
		}

/*		//Set Lifetime of Metaplastic Alloc Signal
		if (mdAllocationDecayRate > 0.0)
		{
			double r = gsl_rng_uniform(mprng);
			if (r < mdAllocationDecayRate) //Allocation Signal Stochatically Switches off
				mbMetaplasticAlloc = false;//Allocation Stops at this Synapse
			//SameThresholdCounters Were Reset upon Allocation
		}*/

	}
}

bool ICascadeSynapse::isPlastic()
{
 return !mbNoPlasticity;
}

//Pass pointer of Simulation Histogram
void ICascadeSynapse::setMetaplasticDistribution(map<uint,uint>* pMDistribinTime)
{
	mpMDistribinTime = pMDistribinTime;
}

//Called by encodeMemory To set the strength over which this synapses transitions is compared against
//When saving the distribution of threshold Crossings
//Tracked State Is set to the initial Synaptic strength in the constructor
void ICascadeSynapse::setTrackedState(uint8_t TargetStrength)
{
	//Setting To zero allows counting of Threshold Events and not just metaplastic cycles, because after then next threshold The counter will be set to 1 and will close at n
	//at which time it will be labelled as Correct or Wrong cycle. The counter will increment To 1 (at reset) whether the next threshold is same/opposite.
	// This will sync the counter to count the number of threshold crossings at the end of each cycle
	//uiSameThresholdTransitionCounter = 0; //Setting To zero allows counting of Threshold Events
	penumTrackedStrength = TargetStrength;
	//mbStopRecordingOfMHistogramAtNextThresholdEvent = true; //Testing
}


/*
 * The value uiThresholdForAllocation is a threshold that can either be used
 * to allocate based on strength stability or Same threshold crossings.
 */
void ICascadeSynapse::setAllocationThreshold(uint uiThres)
{
	uiThresholdForAllocation = uiThres;
}

uint ICascadeSynapse::getAllocationThreshold()
{
	return uiThresholdForAllocation;
}

uint ICascadeSynapse::getRefractionCounter()
{
	return uiStateLifetime;
}

void ICascadeSynapse::resetAllocationRefraction()
{
	uiStateLifetime = 0; //Reset The counteruiStateLifetime

	//**The calling Code Should Explicitly set StabilityAlloc to True -
//	if (uiThresholdForAllocation > 0)
//		mbStabilityAlloc = true; //When Set, the Synapse Will lock After An induction stimulus causes it to cross the stability Threshold
}

/*
Count Same Threshold Crossings - Save distribution at the closing of a cycle
Ignore the cycle closing just after the 1st encoding event t
*/
void ICascadeSynapse::saveMetaplasticDistribution()
{
	assert(mpMDistribinTime!=NULL);
	if (uiSameThresholdTransitionCounter == 0) //This is the First threshold Crossing So Do not save - Counters of threshold crossings get incremented
		return; // Ignore the effects of start state - Count only correct threshold crossing-No Zero cycles possible Except after Init(ignore these)-

	if (iCycleSamplesRemaining == 0)
	{
		mbSaveMetaplasticHistogram = false;
		return;
	}
	iCycleSamplesRemaining--; //Decrement Counter of Number of samples Left

	assert(penumTrackedStrength != SYN_STRENGTH_NOTSET );

	//Separate Correct Strength State Cycles from the Others
	if (penumTrackedStrength != penumStrength) //Wrong State Cycle
	{
		(*mpMDistribinTime)[uiSameThresholdTransitionCounter]+=1; //Increment the occupancy of the current Counter on Global Histogram
	}
	else
	{//Correct Strength Cycle Ended
		(*mpMDistribinTime)[1000+uiSameThresholdTransitionCounter]+=1; //Increment the occupancy of the current Counter on Global Histogram
	}

	(*mpMDistribinTime)[0] +=1; //Increment total Number of samples in the distribution-
}


void ICascadeSynapse::disableDistributionSampleLimit()
{
	iCycleSamplesRemaining = -1; //Make -1 so saveDistribution Can Carry on
}

uint ICascadeSynapse::getMetaplasticDistributionSampleSize()
{
	return (*mpMDistribinTime)[0];
}
//-- Called Everytime the strength of a synapse Changes -
void ICascadeSynapse::resetMetaplasticCounter()
{

	if (mbSaveMetaplasticHistogram)
	{
		saveMetaplasticDistribution();
	}
	//This Flag is only Used for One Cycle To Set mbSaveMetaplasticHistogram = false
	//It has been removed from reset() because reset is called at first encoding and thus we cannot tell the system to sample at t=0
	//because reset would reset this flag
	mbStopRecordingOfMHistogramAtNextThresholdEvent = false;

	//26/3 - Changed To match T.E Counting - THe 0 histState Holds the total number of samples
	//Set to 0 here, but if called after a plasticity step switchreset it is Immediatelly set 1 to count the threshold crossing
	uiSameThresholdTransitionCounter = 0;

	//uiMaxMetaplasticTransitions = 0; //Do not Reset MAX TOO - Only at beginning of new trial
}

//Save Current State As Default
void ICascadeSynapse::startMonitoring()
{
	//mbStrengthChanged = false;
	//mbCascadeIndexChanged = false;
	miStartIndex = miCascadeIndex;
	penumStartStrength = penumStrength;
	mbIsMonitored = true;

}

void ICascadeSynapse::stopMonitoring()
{
	mbIsMonitored = false;
}

bool ICascadeSynapse::isMonitored() const
{
	return mbIsMonitored;
}
//Called At the beginning of a trial -
void ICascadeSynapse::reset()
{
	//mbStopRecordingOfMHistogramAtNextThresholdEvent = false; //This is now set to false After a the threshold is resetMetaplasticCounter
	//uiSameThresholdTransitionCounter	= 0; --Let Asymptotic Value

	uiMaxMetaplasticTransitions 		= 0; //Reset MAX TOO
	uiStateLifetime 					= 0;

	double r = gsl_rng_uniform(mprng);
	if (r < 0.5)
		penumTrackedStrength = penumStartStrength = penumStrength = SYN_STRENGTH_STRONG; //Set Tracked too So we can Obtain a distribution of an unencoded memory
		//penumStartStrength = penumStrength = SYN_STRENGTH_STRONG;
	else
		//penumStartStrength = penumStrength = SYN_STRENGTH_WEAK;
		penumTrackedStrength = penumStartStrength = penumStrength = SYN_STRENGTH_WEAK; //Set Tracked too So we can Obtain a distribution of an unencoded memory

	disableMetaplasticCounting(); //Stop the histogram  Updating
	disableMetaplasticAllocation(); //Removes the allocation Signal That may have been enabled on previous trial
	disableStabilityAllocation(); //Unset Flag For Automatic Allocation
	unfreezePlasticity(); //De-Allocate Synapse
}

//NOTE HACK: Stops the Histogram of Metaplastic Counters - Until Next Reset
uint ICascadeSynapse::getMetaplasticCount()
{
	 ///UNCOMMENT TO MAKE SIMULATION RECORD ONLY COMPLETE CYCLES Do not Record Any More Cycles
	//mbSaveMetaplasticHistogram = false;

	//THE FOLLOWING FLAG Allows for truncating the cycles up to the next threshold event.
	//mbStopRecordingOfMHistogramAtNextThresholdEvent = true; ///This flag tells addSample to set the mbSaveMetaplasticHistogram=false when next theshold is reached

	return uiSameThresholdTransitionCounter;
}

void ICascadeSynapse::disableMetaplasticCounting()
{
	mbSaveMetaplasticHistogram = false;
}

void ICascadeSynapse::enableMetaplasticCounting()
{
	mbSaveMetaplasticHistogram = true;
}
uint ICascadeSynapse::getMaxMetaplasticCount()
{
	return uiMaxMetaplasticTransitions;
}


void ICascadeSynapse::unfreezePlasticity()
{
	mbNoPlasticity = false;
}

void ICascadeSynapse::enableMetaplasticAllocation()
{
	if (uiThresholdForAllocation > 0)
		mbMetaplasticAlloc = true;
}

void ICascadeSynapse::disableMetaplasticAllocation()
{
	mbMetaplasticAlloc = false;
}

//If True Then The synapse May Lock Strength After A plastic Transition
void ICascadeSynapse::enablePlasticAllocation()
{
	mbPlasticAlloc = true;
}


void ICascadeSynapse::disablePlasticAllocation()
{
	mbPlasticAlloc = false;
}

void ICascadeSynapse::disableStabilityAllocation()
{
	mbStabilityAlloc = false;
}


void ICascadeSynapse::enableStabilityAllocation()
{
	if (uiThresholdForAllocation > 0)
		mbStabilityAlloc = true;
}
//
//void ICascadeSynapse::setDecay(double n) {
////Empty
//}

bool ICascadeSynapse::hasStrengthModified() const
{
	return (penumStrength != penumStartStrength);
}

bool ICascadeSynapse::hasIndexModified() const
{
	return (miStartIndex !=miCascadeIndex);
}

void ICascadeSynapse::handleNOP()
{
return;
}
int8_t ICascadeSynapse::handleDEP()
{
	return 0;
}
int8_t ICascadeSynapse::handlePOT()
{
	return 0;
}
int ICascadeSynapse::getCascadeIndex() const
{
	return miCascadeIndex;
}

int ICascadeSynapse::getCascadeSize() const
{
	return (miTerminalIndex+1);
}

int ICascadeSynapse::getStrength() const
{
	return (int)penumStrength;
}

int ICascadeSynapse::getStartStrength() const
{
	return (int)penumStartStrength;
}


void ICascadeSynapse::getTypeAsString(char* buff)
{
	getTypeName(buff);
}

void ICascadeSynapse::getTypeName(char* buff)
{
	strcpy(buff,"_ICascadeSynapse");
}

ICascadeSynapse::~ICascadeSynapse() {
	//  Auto-generated destructor stub
}
