/*
 * synapseFilterUnifiedWithDecay.h
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 */

#ifndef SYNAPSECASCADEFILTERUNIFIEDWITHDECAY_H_
#define SYNAPSECASCADEFILTERUNIFIEDWITHDECAY_H_

#include "synapseCascadeFilterUnified.h"

class synapseCascadeFilterUnifiedWithDecay: public synapseCascadeFilterUnified {

public:
	typedef synapseCascadeFilterUnified super;//Creates the keyword super to refer to base class
	synapseCascadeFilterUnifiedWithDecay();
	synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize

	synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //
	synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,
								ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
								gsl_rng *  prng_r,
								int iRateDependentParameterSet=0); //Default Value is rate=1.0

///Overloaded Methods
	int doStochasticDecay();
	int addSample(int iValue);
	double getHDecay() const;
	double getLDecay() const;
	double getDecay() const; //The Decay Currently running
	int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	void handleNOP();
	void reset();
	//void switchReset();
	void getTypeAsString(char*);
	static void getTypeName(char*);

	virtual ~synapseCascadeFilterUnifiedWithDecay();
	virtual void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	virtual void initialiseFilterState(); //Where the Running Sum should be set when the object is created
	virtual void setFilterThresholds();


	//PDF Lookup tables for the PDF or re-injection to zero - For different rates, cascade states and for the cascade terminal states
	const static double mdPDF_r100[][ciMaxInternalStates],mdPDF_r010[][ciMaxInternalStates],mdPDF_r001[][ciMaxInternalStates];
	const static double mdPDFTerminal_r100[][ciMaxInternalStates],mdPDFTerminal_r010[][ciMaxInternalStates],mdPDFTerminal_r001[][ciMaxInternalStates];


	const static uint8_t ciMaxInternalStates = 7; //Maximum Number of Filter Internal States LThres+HThres+1

	//Publish Decay value set
	const static PtrDecaySet mdDecaysSet[]; //SuperSet Variable Holding All Thresholds For the 3 rate values
	const static PtrDecaySet mdTermDecaySet[];
	//Publish Threshold values set
	const static PtrThresholdSetSmall miThresholdsSet[]; //SuperSet Variable Holding All Thresholds For the 3 rate values
	const static PtrThresholdSetSmall miTerminalThresholdsSet[];

protected:
//	void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	//void initialiseFilterState(); //Where the Running Sum should be set when the object is created
//	virtual void initialiseFilterState(); //Where the Running Sum should be set when the object is created
	//void setFilterThresholds();
	uint32_t miTimeSinceLastInduction;
	double	mdDecayRate; //The Decay rate at the current Step



	//<LOOKUP TABLES> //Thresholds & PDFs for Re injections
	//THresholds

	const static  int8_t  miThreshold_r100[][2],miThreshold_r010[][2],miThreshold_r001[][2];
	const static  int8_t  miThresholdTerminal_r100[][2],miThresholdTerminal_r010[][2],miThresholdTerminal_r001[][2];

	//DECAY
	const static double mdDecayTerminal_r100[][2],mdDecayTerminal_r010[][2],mdDecayTerminal_r001[][2];
	const static  double mdDecay_r100[][2],mdDecay_r010[][2],mdDecay_r001[][2];



	//</LOOKUP TABLES>
private:
};

#endif /* SYNAPSEFILTERUNIFIEDWITHDECAY_H_ */
