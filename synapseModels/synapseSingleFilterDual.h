/*
 * synapseFilterUnified.h
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 */

#ifndef SYNAPSESINGLEFILTERDUAL_H_
#define SYNAPSESINGLEFILTERDUAL_H_
#include "common.h"
#include <string>
//#include "IFilter.h"
#include "ICascadeSynapse.h"
#include "synapseFilterDual.h"

class synapseSingleFilterDual:public synapseFilterDual{

public:
	synapseSingleFilterDual();
	synapseSingleFilterDual(bool bNoPlasticity);
	synapseSingleFilterDual(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize
	synapseSingleFilterDual(int piLThres,int piHThres,double pdDecay); //Initialize with Fixed Threshold And Decay
	synapseSingleFilterDual(int piCascadeSize,float startStrength,gsl_rng * rng_r); //Init With StartStrength Used by SynapticConnection
	synapseSingleFilterDual(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //

	//These Constructors Handles the construction of Stochastic Updaters by setting CascadeSize = 1 and StartIndex > 1
	synapseSingleFilterDual(int piCascadeSize,int piStartIndex, gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	synapseSingleFilterDual(int piCascadeSize,int piStartIndex,
						ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
						gsl_rng *  prng_r,
						int iRateDependentParameterSet=0); //Default Value is rate=1.0

	synapseSingleFilterDual(int piCascadeSize,int piStartIndex,
						ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
						int iStartFilterState,
						gsl_rng *  prng_r,
						int iRateDependentParameterSet=0); //Default Value is rate=1.0

	virtual ~synapseSingleFilterDual();

//Overloaded Methods
	int getCascadeIndex() const;
	int getCascadeSize() const;
	int getStrength() const;
	int getStartStrength() const;
	bool hasStrengthModified() const; //Returns whether a change of state has occured
	bool hasIndexModified() const; // Returns true if Index has Been Incremented
	void startMonitoring(); //Resets Flag which holds if switchSynapse has changed.
	void stopMonitoring();
	bool isMonitored() const;

	//FILTER INTERFACE
	virtual int addSample(int iValue); //Interacts with Filter
	virtual void setFilterThresholds(); //Uses STRENGTH and Index to set thresholds


	int getRunningValueH() const ;
	int getRunningValueL() const ;
	int getLThres() const; //OVERLOADED
	int getHThres() const;
	int getMaxSteps() const; //Returns the number of times the Filter Is allowed to Grow Threshold (The Max Threshold Index)
	int getThresholdIndex() const;
	int lThresReached();
	int hThresReached();
	void doStochasticDecay();

	virtual double getHDecay() const; //No Decay on this class make it Virtual
	virtual double getLDecay() const;
	virtual double getDecay() const; //The Decay Currently running

	virtual int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	virtual int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	virtual void handleNOP();//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
	virtual void reset();
	virtual void switchReset();
	virtual void getTypeAsString(char*);
	static void getTypeName(char*);

protected:
	virtual void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	virtual void initialiseFilterState(); //Where the Running Sum should be set when the object is created

//Footprint of Class Up to Here is 32 Bytes -- OverHead Due to Virtual Functions of Base class and GSL Pointer
	//With GSL Ptr :Unused Bytes 32 - (12+8)+8(vptr) = 4 Dummy Variables To be used
	//Without GSL Ptr :Unused Bytes 24 - 12(data)+8(vptr) = 4 Dummy Variables To be used
private:
	 //Called By Constructors to Initialise common object variables
	void init(int piCascadeSize,int iThresholdsIndex,ICascadeSynapse::SYN_STRENGTH_STATE enumCurrStrengthState,int ParameterSetIndex);

	double dPDecayRate;
	double dDDecayRate;
};

#endif /* synapseSingleFilterDual_H_ */
