/*
 * synapseFilterUnified.h
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 */

#ifndef SYNAPSECASCADEFILTERUNIFIED_H_
#define SYNAPSECASCADEFILTERUNIFIED_H_
#include "common.h"
#include <string>
//#include "IFilter.h"
#include "ICascadeSynapse.h"

//Removed public IFilter as Unnecessary
class synapseCascadeFilterUnified:public ICascadeSynapse{

public:
	synapseCascadeFilterUnified();
	synapseCascadeFilterUnified(bool bNoPlasticity);
	synapseCascadeFilterUnified(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize
	synapseCascadeFilterUnified(int piCascadeSize,float startStrength,gsl_rng * rng_r); //Init With StartStrength Used by SynapticConnection
	synapseCascadeFilterUnified(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //
	synapseCascadeFilterUnified(int piCascadeSize,int piStartIndex,gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	synapseCascadeFilterUnified(int piCascadeSize,int piStartIndex,
						ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
						gsl_rng *  prng_r,
						int iRateDependentParameterSet=0); //Default Value is rate=1.0

	synapseCascadeFilterUnified(int piCascadeSize,int piStartIndex,
						ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
						int iStartFilterState,
						gsl_rng *  prng_r,
						int iRateDependentParameterSet=0); //Default Value is rate=1.0

	virtual ~synapseCascadeFilterUnified();

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
	virtual int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	virtual int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	virtual void handleNOP();//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
	virtual void reset();
	virtual void getTypeAsString(char*);

	//FILTER INTERFACE
	virtual int addSample(int iValue); //Interacts with Filter
	virtual void setFilterThresholds(); //Uses STRENGTH and Index to set thresholds
	virtual void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	virtual void initialiseFilterState(); //Where the Running Sum should be set when the object is created

	int getRunningValue() const ;
	int getLThres() const; //OVERLOADED
	int getHThres() const;
	int getMaxSteps() const; //Returns the number of times the Filter Is allowed to Grow Threshold (The Max Threshold Index)
	int getThresholdIndex() const;

	virtual double getHDecay() const; //No Decay on this class make it Virtual
	virtual double getLDecay() const;
	virtual double getDecay() const; //The Decay Currently running

	virtual void switchReset();
	static void getTypeName(char*);

	//<LOOKUP TABLES> //Thresholds & PDFs for Re injections
	const static uint8_t ciMaxInternalStates = 7; //Maximum Number of Filter Internal States LThres+HThres+1
	const static  int  miThreshold_r100[][2],miThreshold_r010[][2],miThreshold_r001[][2];
	const static  int  miThresholdTerminal_r100[][2],miThresholdTerminal_r010[][2],miThresholdTerminal_r001[][2];
	const static PtrThresholdSet miThresholdsSet[]; //SuperSet Variable Holding All Thresholds For the 3 rate values
	const static PtrThresholdSet miTerminalThresholdsSet[];
	//Lookup tables for the PDF or re-injection to zero - For different rates, cascade states and for the cascade terminal states
	const static double mdPDF_r100[][ciMaxInternalStates],mdPDF_r010[][ciMaxInternalStates],mdPDF_r001[][ciMaxInternalStates];
	const static double mdPDFTerminal_r100[][ciMaxInternalStates],mdPDFTerminal_r010[][ciMaxInternalStates],mdPDFTerminal_r001[][ciMaxInternalStates];;
	//</LOOKUP TABLES>



protected:
	 //Called By Constructors to Initialize common object variables
	void init(int piCascadeSize,int iThresholdsIndex,ICascadeSynapse::SYN_STRENGTH_STATE enumCurrStrengthState,int ParameterSetIndex);


	//Filter State Variables // THESE WERE ALIGNED TO WORDSIZE=64 BUT NOT ANYMORE - Moved 4x8bit Variables To ICascadeIndex

	int8_t  miRFilterValue; //The Internal Running Sum At the current State
	int16_t miLThres;
	int16_t miHThres;
	bool mbInternallyManagedGSL; //Padding For Word Boundary
	uint8_t miParameterSetUsed;

	int8_t dummy3;
	int8_t dummy4;

//--CHANGED Most vars moved to ICascadeSynapse needs recalc:Footprint of Class Up to Here is 32 Bytes -- OverHead Due to Virtual Functions of Base class
	//With GSL Ptr :Unused Bytes 32 - (12+8)+8(vptr) = 4 Dummy Variables To be used
	//Without GSL Ptr :Unused Bytes 24 - 12(data)+8(vptr) = 4 Dummy Variables To be used
private:
};

#endif /* SYNAPSEFILTERUNIFIED_H_ */
