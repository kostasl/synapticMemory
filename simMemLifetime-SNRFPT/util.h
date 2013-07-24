/*
 * util.h
 *
 *  Created on: Dec 19, 2011
 *      Author: kostasl
 */

#ifndef UTIL_H_
#define UTIL_H_

extern int g_FilterTh;

#include "common.h"
#include <sstream> //Used in MakeLogFileNames
#include "../synapseModels/ICascadeSynapse.h"
#include "../synapseModels/synapseCascade.h"


gsl_rng* g_getRandGeneratorInstance(bool newInstance=false); ///GLOBAL Instance and FunctioN - Called by some Default Constructors even Synapse Classes

void MakeListOfFiles(vector<string>& vLogFiles,unsigned int ciInitPeriod,int iCascadeSize,double mdRate,double dFp, int iSynCount);

pair<string,string> SplitFilename (const std::string& str);

//For Rev 1 Objects of ICascadeSynapse Interface
template <class T>
void makeLogFileNames(vector<string>& vLogFiles,unsigned int ciInitPeriod,int iCascadeSize,double mdRate,double dSignalThres,uint trials, int iSynCount, synapseAllocator<ICascadeSynapse>::pFunct pallocF)
{
	if (vLogFiles.empty())
		ERREXIT(500,"Base name missing from Logfile Vector");

	string fOutName(vLogFiles[0].c_str()); //Append Output Directory Which should be stored at location 0 of the vLogFiles Vector

	vLogFiles.clear(); //Clear And Start again
	char buff[200];
	strcpy(buff,"\n");
	T::getTypeName(buff);
	//getModelNameFromFunctionPointer(pallocF,buff);

	fOutName.append(buff);
	fOutName.append("//Capacity");

	std::ostringstream oss;
	oss << "-" << dSignalThres << "-T" << trials;
	//oss << "_NS" << iSynCount << "_I" << ciInitPeriod << "_n" <<  iCascadeSize << "_r" << mdRate;
	oss << "_NS" << iSynCount << "_I%d" << "_n" <<  iCascadeSize << "_r" << mdRate;

	//"-" << dSignalThres
	if (iCascadeSize <= 1)
		oss << "_Th" << g_FilterTh;

	fOutName.append(oss.str());

	vLogFiles.push_back(fOutName); //Corpus
	//HardCode Induction Rate is Fp=0.5

	MakeListOfFiles(vLogFiles,ciInitPeriod,iCascadeSize,mdRate,0.5,iSynCount);
}



//This Is called By Memory Test - Hopfield
template <class T>
void reportStateDistribution(vector<T*>& vCS,int piCount,const char* pcFilename)
{
	if (vCS.empty())
		ERREXIT(500,"Synapse Vector is empty!");

	ofstream *ofile=0;
	int iStateCount = vCS[0]->getCascadeSize();
	int is,js = 0;
	//2 Strength States
	int iOccupancy[2][iStateCount]; //Store counters for each Strength-State Pair
	int iThreshold[2][iStateCount];
	float fDecay[iStateCount];

	memset(iOccupancy,0,sizeof(int)*2*iStateCount); //Empty the memory Buffer
	memset(iThreshold,0,sizeof(int)*2*iStateCount); //Empty the memory Buffer
	memset(fDecay,0,sizeof(float)*iStateCount); //Empty the memory Buffer

	vector<ICascadeSynapse*>::iterator	it;

	 for (typename vector<T*>::iterator it =vCS.begin();it != vCS.end();it++)
	 {
		 is = ((int)(*it)->getStrength() == ICascadeSynapse::SYN_STRENGTH_STRONG)?1:0;
		 js = (int)(*it)->getCascadeIndex();
		 //js = 2;
		 iOccupancy[is][js] +=1;

			if ((*it)->getStrength() == ICascadeSynapse::SYN_STRENGTH_STRONG)
			{
				if (iThreshold[0][js] == 0)
				{
					iThreshold[0][js] 	= (*it)->getLThres();
					iThreshold[1][js] 	= (*it)->getHThres();
					fDecay[js]		 	= (*it)->getDecay();
				}
				else
				{
					//Cant Test Anymore Cause We are switching Randomly - Symmetric Assymetric
				//	assert(iThreshold[0][js] == (*it)->getFilter()->getLThres());
				}
			}

	 }


	cout << "State \t Weak \t Strong \t LTh \t HTh \t Decay" << endl;
	 if (pcFilename)
		 ofile =  new ofstream(pcFilename, ios::out ); //OPEN OUTPUT FILE

	 int isumWeak=0;
	 int isumStrong = 0;

	 //Report Weak Distribution
	 for(int k=(iStateCount-1);k>-1;k--)
	 {
		  isumWeak += iOccupancy[0][k];
		 if (pcFilename)
		  *ofile << (-(k+1)) << "\t" <<  iOccupancy[0][k] << endl;
	 }

	 for(int k=0;k<(iStateCount);k++)
	 {
			 ///Show On Screen
			 cout << (k+1) << " \t " << iOccupancy[0][k] << " \t " << iOccupancy[1][k] << "\t \t" << iThreshold[0][k] << "\t" << iThreshold[1][k] << "\t" << fDecay[k] << endl;
			 isumStrong += iOccupancy[1][k];

			 //Save to File Also
		 if (pcFilename)
		  *ofile << ((k+1)) << "\t" <<  iOccupancy[1][k] << endl;
	 } //For Each State

	 cout << "----------------------------------------------------" << endl;
	 cout <<  "X\t" << isumWeak << "\t" << isumStrong  << "\t =" << (isumWeak+isumStrong) << endl;

	 //Close File
	 if (pcFilename)
		 ofile->close();
}



//Uses an externally Handled Occupancy Array To Save the Average Occupancy Up to Trial t - Called By Sim signals From File
template <class T>
void reportAvgStateDistribution(vector<T*>& vCS,int piCount,const char* pcFilename,long** piOccupancy,int iTotalTrials)
{
	if (vCS.empty())
		ERREXIT(500,"Synapse Vector is empty!");

	if (iTotalTrials < 1)
	{
		cerr << "Cant Do Average Distribution Not enough Samples T:" << iTotalTrials;
		return;
	}
	assert(iTotalTrials > 0);

	ofstream *ofile=0;
	int iStateCount = vCS[0]->getCascadeSize();
	int is,js = 0;
	//2 Strength States

	long iOccupancy[2][iStateCount];
	int iThreshold[2][iStateCount];
	float fDecay[iStateCount];

	memset(iOccupancy,0,sizeof(long)*2*iStateCount); //Empty the memory Buffer
	memset(iThreshold,0,sizeof(int)*2*iStateCount); //Empty the memory Buffer
	memset(fDecay,0,sizeof(float)*iStateCount); //Empty the memory Buffer

	vector<ICascadeSynapse*>::iterator	it;

	 for (typename vector<T*>::iterator it =vCS.begin();it != vCS.end();it++)
	 {
		 is = ((int)(*it)->getStrength() == ICascadeSynapse::SYN_STRENGTH_STRONG)?1:0;
		 js = (int)(*it)->getCascadeIndex();
		 //js = 2;
			if ((*it)->getStrength() == ICascadeSynapse::SYN_STRENGTH_STRONG)
			{
				if (iThreshold[0][js] == 0)
				{
					iThreshold[0][js] 	= (*it)->getLThres();
					iThreshold[1][js] 	= (*it)->getHThres();
					fDecay[js]		 	= (*it)->getDecay();
				}
				else
				{
					//Cant Test Anymore Cause We are switching Randomly - Symmetric Assymetric
				//	assert(iThreshold[0][js] == (*it)->getFilter()->getLThres());
				}
			}

	 }

	 //Get The Avg Occupancy - Divide By the Number of Trials
	 //Report Weak Distribution
	 for(int k=0;k<iStateCount;k++)
	 {
		 cout << piOccupancy[0][k] << "->";
		 iOccupancy[0][k]+= piOccupancy[0][k]/(iTotalTrials);
		 iOccupancy[1][k]+= piOccupancy[1][k]/(iTotalTrials);

		 cout << iOccupancy[0][k] << endl;
	 }
	 cout << "AVERAGE OF:" << iTotalTrials <<endl;
	 cout << "State \t Weak \t Strong \t LTh \t HTh \t Decay" << endl;
	 if (pcFilename)
		 ofile =  new ofstream(pcFilename, ios::out ); //OPEN OUTPUT FILE

	 int isumWeak=0;
	 int isumStrong = 0;

	 //Report Weak Distribution
	 for(int k=(iStateCount-1);k>-1;k--)
	 {
		  isumWeak += iOccupancy[0][k];
		 if (pcFilename)
		  *ofile << (-(k+1)) << "\t" <<  iOccupancy[0][k] << endl;
	 }

	 for(int k=0;k<(iStateCount);k++)
	 {
			 ///Show On Screen
			 cout << (k+1) << " \t " << iOccupancy[0][k] << " \t " << iOccupancy[1][k] << "\t \t" << iThreshold[0][k] << "\t" << iThreshold[1][k] << "\t" << fDecay[k] << endl;
			 isumStrong += iOccupancy[1][k];

			 //Save to File Also
		 if (pcFilename)
		  *ofile << ((k+1)) << "\t" <<  iOccupancy[1][k] << endl;
	 } //For Each State

	 cout << "----------------------------------------------------" << endl;
	 cout <<  "X\t" << isumWeak << "\t" << isumStrong  << "\t =" << (isumWeak+isumStrong) << endl;

	 //Close File
	 if (pcFilename)
		 ofile->close();
}


//This Is called By Memory Test - Hopfield
template <class T>
void saveStateDistribution(vector<T*>& vCS,long** piOccupancy,int piCount)
{
	if (vCS.empty())
		ERREXIT(500,"Synapse Vector is empty!");

///	int iStateCount = vCS[0]->getCascadeSize();
	int is,js = 0;
	//2 Strength States

	vector<ICascadeSynapse*>::iterator	it;

	 for (typename vector<T*>::iterator it =vCS.begin();it != vCS.end();it++)
	 {
		 is = ((int)(*it)->getStrength() == ICascadeSynapse::SYN_STRENGTH_STRONG)?1:0;
		 js = (int)(*it)->getCascadeIndex();
		 //js = 2;
		 piOccupancy[is][js] +=1;

	 }

}





#endif /* UTIL_H_ */
