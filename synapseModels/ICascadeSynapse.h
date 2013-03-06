/*
 * ICascadeSynapse.h
 *
 *  Created on: 9 Mar 2011
 *      Author: kostasl
 */

#ifndef ICASCADESYNAPSE_H_
#define ICASCADESYNAPSE_H_

#define MAX_CASCADE_SIZE 20
#define DEFAULT_CASCADE_SIZE  1
#define DEFAULT_CASCADE_X  0.5f

#include "common.h" //For Gsl
#include <string>


//Define A datatype for the Threshold Arrays So as To use it to define a super Set of all Sets for different Rates
typedef const int (*PtrThresholdSet)[20][2];
typedef const int8_t (*PtrThresholdSetSmall)[20][2];

typedef const double (*PtrDecaySet)[20][2];
typedef const double (*PtrPDFSet)[20][2];


class ICascadeSynapse {
public:
	enum SYN_STRENGTH_STATE {SYN_STRENGTH_WEAK_ALLOC = -2, SYN_STRENGTH_WEAK = -1, SYN_STRENGTH_STRONG = 1,SYN_STRENGTH_STRONG_ALLOC = 2,SYN_STRENGTH_NOTSET = -10};
	ICascadeSynapse();
	ICascadeSynapse(bool PlasticAlloc, bool MetaplasticAlloc,bool bStabilityAlloc);
	ICascadeSynapse(int piCascadeSize,float startStrength,gsl_rng * rng_r); //Init With StartStrength Used by SynapticConnection
	ICascadeSynapse(bool bNoPlasticity);
	ICascadeSynapse(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize
	ICascadeSynapse(int piCascadeSize,int piStartIndex,gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	ICascadeSynapse(int piCascadeSize,	int piStartIndex, SYN_STRENGTH_STATE penumStrength, gsl_rng * rng_r);

//TODO: THESE CAN BE MADE NON-VIRTUAL AND REDUCE MEMORY FOOTPRINT
	//By using the getTypeAsString to identify the object
	int getCascadeIndex() const; // Returns the current position in Cascade - Either Weak or Strong
	int getCascadeSize() const;
	int getStrength() const;
	int getStartStrength() const;

	bool hasStrengthModified() const; //Returns whether a change of state has occured
	bool hasIndexModified() const; // Returns true if Index has Been Incremented

	//Temporary
///BEGIN OF PURE VIRTUAL
	//virtual void setDecay(double newDecayRate);
	//virtual int getLThres() const=0;
	//virtual int getHThres() const=0;
	//virtual double getDecay() const;
	int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	void handleNOP();//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
	void getTypeName(char* buff);
	virtual void getTypeAsString(char*);
///END OF PURE VIRTUAL

	void startMonitoring(); //Resets Flag which holds if switchSynapse has changed.
	void stopMonitoring();
	bool isMonitored() const; //Pure Virtual

	void reset(); // Reset Of Synapse - Beginning of new trial - Reset Strength and State

	//virtual void getTypeName(char* buff) =0;
	void freezePlasticity();
	void unfreezePlasticity();
	void setTrackedState(uint8_t StrengthState);//Set the penumTrackedStrength

	void setAllocationThreshold(uint); //Set the Threshold of stability over which plasticity is locked automatically
	void setMetaplasticDistribution(map<uint,uint>* mpMDistribinTime); //Pass pointer of Simulation Histogram

	uint getAllocationThreshold(); //Return the stability threshold set
	uint getRefractionCounter();
	uint getMetaplasticCount(); //Get the number of metaplastic steps since counter reset
	uint getMaxMetaplasticCount(); //Returns the statistic of Max Metaplastic Transitions Since Last Reset.
	uint getMetaplasticDistributionSampleSize();

	void resetAllocationRefraction();
	void resetMetaplasticCounter();
	void enableMetaplasticAllocation(); //Set a flag that will enable the locking of plasticity if number of metaplastic steps exceeds a threshold
	void disableMetaplasticAllocation();

	void enableMetaplasticCounting(); //Set Flag To Save to Histogram
	void disableMetaplasticCounting();

	void enablePlasticAllocation(); //Allocate the next plastic transition
	void disablePlasticAllocation();
	void enableStabilityAllocation(); //Counts the number of memories that this synapse has remain in the same strength - Allocation If threshold Is exceeded
	void disableStabilityAllocation(); ///Stable Synapses Will not be automatically Allocated
	void saveMetaplasticDistribution();

	void disableDistributionSampleLimit(); //Remove the limitation on the number of samples obtained for thres-cycle sampling

	bool isPlastic(); //Returns true if synapse is not Locked (Allocated)

	virtual ~ICascadeSynapse();

protected:
	gsl_rng* mprng;
	map<uint,uint>* mpMDistribinTime; //Pointer to Main simulation Histogram - Updated everytime Counters Are reset

	double mdAllocationDecayRate; //Probability of Disabling Metaplastic Allocation On the next Step - Decay Rate of AllocState
	bool mbNoPlasticity;
	bool mbMetaplasticAlloc; //Lock Plasticity After A metaplastic Transition
	bool mbPlasticAlloc; //Lock Plasticity After A plastic Transition
	bool mbStabilityAlloc; //If Set Then A synapse Locks when automatically it has exceeded the stability Criteria
	bool mbAllocationTag; //A SynapseMark Set for when local allocation conditions have been met. Here Set when a synapse is Strong and Receives a Strong Input +1*+1

	bool mbSaveMetaplasticHistogram; //If Set Then each reset of the Metaplastic counter is saved in this histogram
	bool mbStopRecordingOfMHistogramAtNextThresholdEvent; //Change mbSaveMetaplasticHistogram to true on next threshold event
	bool mbIsMonitored; //true Means this synapse is part of the Tracked group - As fusi defines it
	uint uiStateLifetime; //Used for Allocation - Timer for how long a strength state has been alive Measured In Induction stimuli

	int iCycleSamplesRemaining; //Once it reaches 0, Cycle Sampling Ends

	uint8_t miTerminalIndex; // Reflects the Cascade Size - But Index is Size-1
	uint8_t miCascadeIndex; //As from The CascadeSynapse it counts how deep we go down each cascade.
	uint8_t miStartIndex; //The Current Running Index, And the Index Obtained at Initialization
	int8_t penumStartStrength;

	int8_t penumStrength; //Not Typed to Enum  SYN_STRENGTH_STATE but Using it
	int8_t penumTrackedStrength; //The strength against which metaplastic cycles are evaluated(Must be set to the tracked memory state)
	int8_t	dummy1;
	int8_t dummy2;

	int8_t dummy4;
	uint8_t uiThresholdForAllocation; //The time spent in this state used as a threshold for allocation
	uint8_t uiSameThresholdTransitionCounter; //Count The number of correct threshold crossing transitions. This may match CascadeIndex
	uint8_t uiMaxMetaplasticTransitions; //Holds the statistic of Maximum value obtained in the metaplastic counter

};

#endif /* ICASCADESYNAPSE_H_ */
