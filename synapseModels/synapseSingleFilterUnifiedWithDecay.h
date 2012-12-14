/*
 * synapseSingleFilterUnifiedWithDecay.h
 *
 *  Created on: 30 Mar 2011
 *      Author: kostasl
 *      TODO: Remove The dependency on CascadeFilter Class - Simplify by using the ICascade Interface directly
 */

#ifndef SYNAPSESINGLEFILTERUNIFIEDWITHDECAY_H_
#define SYNAPSESINGLEFILTERUNIFIEDWITHDECAY_H_

//#include "synapseCascadeFilterUnifiedWithDecay.h"
#include "ICascadeSynapse.h"

class synapseSingleFilterUnifiedWithDecay: public ICascadeSynapse {

public:
	typedef ICascadeSynapse super;//Creates the keyword super to refer to base class
	synapseSingleFilterUnifiedWithDecay();
	synapseSingleFilterUnifiedWithDecay(int piFilterSize,gsl_rng * rng_r); //Initialize with cascadeSize
	synapseSingleFilterUnifiedWithDecay(int piFilterSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //Initialize with Filter And Strength
	//synapseSingleFilterUnifiedWithDecay(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //
	//synapseSingleFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	//synapseSingleFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng *  prng_r,int iRateDependentParameterSet=0); //Default Value is rate=1.0

	synapseSingleFilterUnifiedWithDecay(int iLThres,int iHThres,double pdDecay); //A Generic Constructor
	synapseSingleFilterUnifiedWithDecay(int iLThres, int iHThres, double pdDecay, int iStartState, int iSampleCycles, gsl_rng* prng); //Constructor sets Running Value - Filter State
	synapseSingleFilterUnifiedWithDecay(int iLThres, int iHThres, double pdDecay, int iSampleCycles, gsl_rng* prng); //Constructor sets Running Value - Filter State
	synapseSingleFilterUnifiedWithDecay(int iLThres,int iHThres,double pdDecay,uint uiRefractionThreshold,double dAllocDecayRate); //A Generic Constructor

	///Overloaded Methods
	int doStochasticDecay();
	int addSample(int iValue);
	double getHDecay() const;
	double getLDecay() const;
	int getLThres() const;
	int getHThres() const;
	double getDecay() const; //The Decay Currently running
	void setDecay(double newDecay);
	int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	void handleNOP();
	void reset();
	void switchReset();
	void getTypeAsString(char*);
	void setTrackedState(uint8_t TargetStrength);///Overload From ICascade So I can revert the set internal Filter State--Used for Sampling the distribution
	static void getTypeName(char*);

	//int getCascadeIndex() const;
	//int getCascadeSize() const;
	//int getStrength() const;
	//int getStartStrength() const;
	//bool hasStrengthModified() const; //Returns whether a change of state has occured
	//bool hasIndexModified() const; // Returns true if Index has Been Incremented
	//void startMonitoring(); //Resets Flag which holds if switchSynapse has changed.
	//void stopMonitoring();
	//bool isMonitored() const;


	int getRunningValue() const ;
	virtual ~synapseSingleFilterUnifiedWithDecay();


	//Overloaded Methods
	//int getCascadeIndex() const;
	//int getCascadeSize() const;
	//int getStrength() const;
	//int getStartStrength() const;
	//bool hasStrengthModified() const; //Returns whether a change of state has occured
	//bool hasIndexModified() const; // Returns true if Index has Been Incremented
	//void startMonitoring(); //Resets Flag which holds if switchSynapse has changed.
	//void stopMonitoring();
	//bool isMonitored() const;





protected:
	void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	void initialiseFilterState(); //Where the Running Sum should be set when the object is created

	void setFilterThresholds();
	bool mbInternallyManagedGSL; //Padding For Word Boundary
	double	mdDecayRate; //The Decay rate at the current Step
	//int8_t miParameterSetUsed; //uint8_t;
	int8_t miLThres;
	int8_t miHThres;
	int8_t miTimeSinceLastInduction;
	int8_t  miRFilterValue; //int8_t; //The Internal Running Sum At the current State



private:
	//<LOOKUP TABLES>
	//These are found in the base class
	//</LOOKUP TABLES>
};

#endif /* synapseSingleFilterUnifiedWithDecay_H_ */
