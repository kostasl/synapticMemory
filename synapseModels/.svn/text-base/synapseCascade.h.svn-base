/*
 * synapseCascade.h
 *
 *  Created on: 22 Mar 2011
 *      Author: kostasl
 */

#ifndef SYNAPSECASCADE_H_
#define SYNAPSECASCADE_H_
//#include "stdafx.h"
#include "ICascadeSynapse.h"

class synapseCascade: public ICascadeSynapse {
public:
	synapseCascade();
	synapseCascade(int piCascadeSize,gsl_rng * rng_r); //Initialize with cascadeSize
	synapseCascade(int piCascadeSize,float startStrength,gsl_rng * rng_r); //Init With StartStrength Used by SynapticConnection
	synapseCascade(int piCascadeSize,ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,gsl_rng * rng_r); //
	synapseCascade(int piCascadeSize,int piStartIndex,gsl_rng * rng_r); //Used to Start From A fixed point and test how distribution evolves
	synapseCascade(int piCascadeSize,int piStartIndex,SYN_STRENGTH_STATE penumStrength,gsl_rng * rng_r);
	synapseCascade(double pQprob,gsl_rng * rng_r); // A single Updater INstantiation
	synapseCascade(float pQprob[],float pPprob[],int piCascadeSize,gsl_rng * rng_r); // Custom Transition Probability Initialization

//Overloaded Methods
	int getCascadeIndex() const;
	int getCascadeSize() const;
	int getStrength() const;
	int getStartStrength() const;
	bool hasStrengthModified() const; //Returns whether a change of state has occured
	bool hasIndexModified() const; // Returns true if Index has Been Incremented
	void startMonitoring(); //Resets Flag which holds if switchSynapse has changed.
	void stopMonitoring(); //Resets Flag which holds if switchSynapse has changed.
	bool isMonitored() const;
	int8_t handlePOT(); //Induces an LTP event depending on current state either switches cascades or moves down
	int8_t handleDEP();//Induces an LTD event depending on current state either switches cascades or moves down
	void handleNOP();//Time Increment where No Induction of Stimulus on Synapse - Can be empty - But used with the Filter Case
	void reset();

	void getTypeAsString(char*);
	static void getTypeName(char*);
	//static void getTypeName(string& buff);

//Filter Baggage -Inherited from the ICASCADESYnapse
	int getLThres() const;
	int getHThres() const;
	double getDecay() const; //The Decay Currently running

	virtual ~synapseCascade();

protected:
	float doStochasticDEP();
	float doStochasticPOT();
	void initDefaultTransitionProb();//Initialises the Transition probabilities P's and Qs in a geometric progression as defined in the 05 Fusi Paper
	void init(int piStatesCount, //Number of Cascade states n
							float pfPCascadeSwitch[], //Probability array of q's
							float pfPCascadeTrans[],  //Probability array of p's
							int piCurrCascadeIndex, //Cascade Current State
							SYN_STRENGTH_STATE pfCurrStrengthState,
							gsl_rng * prng_r
							);

	uint8_t miCascadeIndex;
	uint8_t miCascadeStartIndex; //The Current position in the Cascade
	uint8_t miCascadeSize; // The Number of States in each cascade
	bool mbStrengthChanged;
	bool mbCascadeIndexChanged; //When A switchSynapse Cascade Index is increased
	bool mbNoPlasticity;
	int8_t menumStrength_State;
	int8_t menumStartStrength_State; //Not Using Enum - But Will be set using SYN_STRENGTH_STATE
	//One Word//

	bool mbIsMonitored; //true Means this synapse is part of the Tracked group - As fusi defines it
	static float mfQprob[MAX_CASCADE_SIZE];
	static float mfPprob[MAX_CASCADE_SIZE];

};



///THESE PARAMETERS ARE USED FOR CUSTOM CASCADE PROB INITIALIZATION WHEN TRYING TO MATCH A FILTER
//N=10 Global Transition Probabilities for Cascade Synapses Initialization up to n=15												//12	   13			14			15
/*
float gfQCascTransition[] = {1,		 0.25,		 0.25, 	 0.11111,
							 0.0625, 0.0277778, 0.015625, 0.0082644,
							 0.0039062, 0.00384615, 0.000976563, 0.000493827,
							 0.000244141,0.000123457,0.000121197};

float gfPCascTransition[] = {1,		 0.25,		 0.25, 	 0.11111,
							0.0625, 0.0277778, 0.015625, 0.0082644,
							0.0039062, 0.001890359, 0.000976563, 0.000493827,
							0.000244141,0.000123457,0.000061035};
*/

//N=6 Global Transition Probabilities for Cascade Synapses Initialization up to n=15												//12	   13			14			15
/*
float gfQCascTransition[] = {1,		 0.25,		 0.25, 	 0.11111,
							 0.0625, 0.0666667, 0.015625, 0.0082644,
							 0.0039062, 0.00384615, 0.000976563, 0.000493827,
							 0.000244141,0.000123457,0.000121197};

float gfPCascTransition[] = {1,		 0.25,		 0.25, 	 0.11111,
							0.0625, 0.0277778, 0.015625, 0.0082644,
							0.0039062, 0.001890359, 0.000976563, 0.000493827,
							0.000244141,0.000123457,0.000061035};

*/


#endif /* SYNAPSECASCADE_H_ */
