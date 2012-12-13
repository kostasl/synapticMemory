/*
 * HopfieldMemoryTests.h
 *
 *  Created on: 17 Dec 2010
 *      Author: kostasl
 */

#ifndef HOPFIELDMEMORYTESTS_H_
#define HOPFIELDMEMORYTESTS_H_

//Function Prototypes
//Hopfield - One Function for All Models - Generic
//void doHopfieldCapacityTest(pAllocationFunct pF,uint initPatterns,int maxCascSize,int startIndex);


int recallHopfieldNetPattern(uint _uiNeuronCount,uint StartNeuron,int* tX, int* X,float** W,
							  float& AvgSignal, uint trials ); //Used for Demo

void doHopfieldCapacityTest(int modelType,string modelName,uint iSynCount,uint trials, uint initPatterns,int maxCascSize,int startIndex);
//float** makeWeightMatrixCascadeSyns(int NetSize,float** Xin,uint patCount,uint iTrackedIndex,int _CascadeSize,char**& buffer,float**& w,vector<ICascadeSynapse*>& vSyn,pAllocationFunct2 pF,vector<string>& slogfiles);
float** makeWeightMatrix(int NetSize,int* Xin,float** W);


template<class T>
void deleteMemoryBuffer(uint NetSize,T**& buffer);

//Constants And Global Vars
//Simulation Time step
//const float h=0.0002f;
const float h=0.0001f;
extern uint g_time; //Define a global Discrete Time

//bool bRecallInProgress; //Used by Camera Demo

#endif /* HOPFIELDMEMORYTESTS_H_ */
