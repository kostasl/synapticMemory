/*
 * InputVectorHandling.h
 *
 *  Created on: 10 Mar 2011
 *      Author: kostasl
 */

#ifndef INPUTVECTORHANDLING_H_
#define INPUTVECTORHANDLING_H_


///VECTOR FILE FUNCTIONS
void makeTestVectorsInFile(uint PatCount,uint VecSize,float fbitBalance);
void makeDriftlessTestVectorsInFile(uint PatCount,uint VecSize,float fbitBalance);
uint makeGSTestVectorsInFile(int PatCount,uint VecSize,float fbitBalance);
void makeGNTestVectorsInFile(int PatCount,uint VecSize,float fbitBalance);
void makeGN2TestVectorsInFile(int PatCount,uint VecSize,float fbitBalance); //Population Mutations
void makeHadamardShuffledSet(uint PatCount,uint VecSize);//Opens a Hadamard Input - Shuffles n Vectors and Writes to Output
void readTestVectorsFromFile(const char* fname,int** X,uint VecSize, uint& iPatsToLoad,uint indexToInsertVectors, bool bShuffleVectors = true);
void writeTestVectorsToFile(char* cfOutName,t_inVal** X,uint VecSize,uint PatCount);
void writeTestVectorsToFileTxt(char* cfOutName,t_inVal**X,uint VecSize,uint PatCount);

void initPatternMemory(int **X,uint PatCount,uint _uiNeuronCount,uint iTrackedMemIndex,float Fp,gsl_rng* mprng,bool bUseRandomPatterns);

void convertVectorFileToTxt(uint uiPatCount,uint iSynCount,char* pInputFile); //Made for T.E
//MEASUREMENT AND UTIL FUNCTIONS
double dotprod(t_inVal* X1,t_inVal* X2,uint VecSize);

#endif /* INPUTVECTORHANDLING_H_ */
