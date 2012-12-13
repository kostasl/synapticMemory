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
void doHopfieldCapacityTest(pAllocationFunct2 pF,uint initPatterns,int maxCascSize,int startIndex);

int testHopfieldBinarySyns(int _iPatCount,uint _uiNeuronCount,float** X,float* tX,
						 float _fProbeNoiseLevel,int _iStoredPatIndex,int _iCascadeSize,
						 pAllocationFunct2 pf,
						 float& AvgSignal, uint trials ,vector<string>& slogFiles,
						 vector<ICascadeSynapse*>& vSyn,char**& mem_buffer); //Runs 100 Trials of storing and recalling a tracked memory

int recallHopfieldNetPattern(uint _uiNeuronCount,uint StartNeuron,int* tX, int* X,float** W,
							  float& AvgSignal, uint trials ); //Used for Demo

float** makeWeightMatrixCascadeSyns(int NetSize,float** Xin,uint patCount,uint iTrackedIndex,int _CascadeSize,char**& buffer,float**& w,vector<ICascadeSynapse*>& vSyn,pAllocationFunct2 pF,vector<string>& slogfiles);
float** makeWeightMatrix(int NetSize,int* Xin,float** W);
int SearchForNetCapacity(uint _uiNeuronCount,uint iTrackedMemIndex,float _fProbeNoiseLevel,int _iCascadeSize,pAllocationFunct2 pF,uint itrials);


template<class T>
void deleteMemoryBuffer(uint NetSize,T**& buffer);

//Constants And Global Vars
//Simulation Time step
//const float h=0.0002f;
const float h=0.0001f;
extern uint g_time; //Define a global Discrete Time

//bool bRecallInProgress; //Used by Camera Demo

#endif /* HOPFIELDMEMORYTESTS_H_ */
