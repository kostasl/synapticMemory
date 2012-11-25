/*
 * synapseSingleFilterUnifiedWithDecay.h
 *
 *  Created on: 30 Mar 2011
 *      Author: kostasl
 */

#ifndef SYNAPSESINGLEFILTERUNIFIEDWITHDECAYREFLECTING_H_
#define SYNAPSESINGLEFILTERUNIFIEDWITHDECAYREFLECTING_H_

#include "synapseCascadeFilterUnifiedWithDecay.h"

class synapseSingleFilterUnifiedWithDecayReflecting: public synapseCascadeFilterUnifiedWithDecay {

public:
	typedef synapseCascadeFilterUnifiedWithDecay super;//Creates the keyword super to refer to base class
	synapseSingleFilterUnifiedWithDecayReflecting();
	synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize
	synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //
	synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,int piStartIndex,gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	synapseSingleFilterUnifiedWithDecayReflecting(int piCascadeSize,int piStartIndex,
								ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
								gsl_rng *  prng_r,
								int iRateDependentParameterSet=0); //Default Value is rate=1.0

	synapseSingleFilterUnifiedWithDecayReflecting(int piLThres,int piHThres,double pdDecay); //A Generic Constructor

///Overloaded Methods
	int doStochasticDecay();
	int addSample(int iValue);
	double getHDecay() const;
	double getLDecay() const;
	double getDecay() const; //The Decay Currently running
	void setDecay(double newDecay);
	int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	void handleNOP();
	void reset();
	void switchReset();
	//void switchReset();
	void getTypeAsString(char*);
	static void getTypeName(char*);

	virtual ~synapseSingleFilterUnifiedWithDecayReflecting();

protected:
	void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	void initialiseFilterState(); //Where the Running Sum should be set when the object is created

	void setFilterThresholds();

private:
	//<LOOKUP TABLES>
	//These are found in the base class
	//</LOOKUP TABLES>
};

#endif /* synapseSingleFilterUnifiedWithDecayReflecting_H_ */
