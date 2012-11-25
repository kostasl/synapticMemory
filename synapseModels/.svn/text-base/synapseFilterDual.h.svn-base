/*
 * synapseFilterUnified.h
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 */

#ifndef SYNAPSEFILTERDUAL_H_
#define SYNAPSEFILTERDUAL_H_
#include "common.h"
#include <string>
//#include "IFilter.h"
#include "ICascadeSynapse.h"


//Removed public IFilter as Unnecessary
class synapseFilterDual:public ICascadeSynapse{

public:
	synapseFilterDual();
	synapseFilterDual(bool bNoPlasticity);
	synapseFilterDual(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize
	synapseFilterDual(int piCascadeSize,float startStrength,gsl_rng * rng_r); //Init With StartStrength Used by SynapticConnection
	synapseFilterDual(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //

	//These Constructors Handles the construction of Stochastic Updaters by setting CascadeSize = 1 and StartIndex > 1
	synapseFilterDual(int piCascadeSize,int piStartIndex, gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	synapseFilterDual(int piCascadeSize,int piStartIndex,
						ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
						gsl_rng *  prng_r,
						int iRateDependentParameterSet=0); //Default Value is rate=1.0

	synapseFilterDual(int piCascadeSize,int piStartIndex,
						ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
						int iStartFilterState,
						gsl_rng *  prng_r,
						int iRateDependentParameterSet=0); //Default Value is rate=1.0

	virtual ~synapseFilterDual();

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

	//<LOOKUP TABLES> //Thresholds & PDFs for Re injections
	const static uint8_t ciMaxInternalStates = 4; //Maximum Number of Filter Internal States LThres+HThres+1
	const static  int  miThreshold_r100[][2],miThreshold_r010[][2],miThreshold_r001[][2];
	const static  int  miThresholdTerminal_r100[][2],miThresholdTerminal_r010[][2],miThresholdTerminal_r001[][2];
	const static double mdDecayTerminal_r100[][2],mdDecayTerminal_r010[][2],mdDecayTerminal_r001[][2];
	const static double mdDecay_r100[][2],mdDecay_r010[][2],mdDecay_r001[][2];

	const static PtrThresholdSet miThresholdsSet[]; //SuperSet Variable Holding All Thresholds For the 3 rate values
	const static PtrThresholdSet miTerminalThresholdsSet[];
	const static PtrDecaySet mdDecaysSet[]; //SuperSet Variable Holding All Thresholds For the 3 rate values
	const static PtrDecaySet mdTermDecaySet[];
	//Lookup tables for the PDF or re-injection to zero - For different rates, cascade states and for the cascade terminal states
	const static double mdPDF_r100[][ciMaxInternalStates],mdPDF_r010[][ciMaxInternalStates],mdPDF_r001[][ciMaxInternalStates];
	const static double mdPDFTerminal_r100[][ciMaxInternalStates],mdPDFTerminal_r010[][ciMaxInternalStates],mdPDFTerminal_r001[][ciMaxInternalStates];;
	//</LOOKUP TABLES>

protected:
	virtual void reInjectFilterStateToCascadeState(); //Resets The Filter Internal Sum According to PDF
	virtual void initialiseFilterState(); //Where the Running Sum should be set when the object is created

	//Filter State Variables //ALIGN TO WORDSIZE=64 On this Maching
	uint8_t miTerminalIndex; // Reflects the Cascade Size - But Index is Size-1
	uint8_t miCascadeIndex;
	uint8_t miStartIndex; //The Current Running Index, And the Index Obtained at Initialization
	uint8_t miLThres;
	uint8_t miHThres;
	uint16_t  miRPFilterValue; //The Internal Running Sum POT Signals
	bool mbNoPlasticity;
	///WORD FULL 8 bytes


	uint16_t  miRDFilterValue; //The Internal Running Sum DEP Signals
	uint16_t miTimeSinceLastInduction;
	int8_t penumStartStrength;
	int8_t penumStrength; //Not Using Enum
	uint8_t miParameterSetUsed;
	bool mbIsMonitored; //true Means this synapse is part of the Tracked group - As fusi defines it
	///2nd WORD FULL 8 bytes

	//bool mbStrengthChanged;
	//bool mbCascadeIndexChanged; //When A switchSynapse Cascade Index is increased

	//SYN_STRENGTH_STATE penumStrength; //Current Synapse State and Initial Strength On Construction

	//bool mbInternallyManagedGSL; //Padding For Word Boundary
	//int8_t dummy3,dummy4,dummy5,dummy6,dummy7,dummy8;
	///3rd WORD FULL 8 bytes

	//Pointer held In ICascadeSynapse +1 WORD
	//gsl_rng * mprng_r; //Pointer To RandGen


//Footprint of Class Up to Here is 32 Bytes -- OverHead Due to Virtual Functions of Base class and GSL Pointer
	//With GSL Ptr :Unused Bytes 32 - (12+8)+8(vptr) = 4 Dummy Variables To be used
	//Without GSL Ptr :Unused Bytes 24 - 12(data)+8(vptr) = 4 Dummy Variables To be used
private:
	 //Called By Constructors to Initialise common object variables
	void init(int piCascadeSize,int iThresholdsIndex,ICascadeSynapse::SYN_STRENGTH_STATE enumCurrStrengthState,int ParameterSetIndex);
};

#endif /* synapseFilterDual_H_ */
