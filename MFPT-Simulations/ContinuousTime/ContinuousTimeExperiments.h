/*
 * ContinuousTimeExperiments.h
 *
 *  Created on: 7 Oct 2011
 *      Author: kostasl
 */

#ifndef CONTINUOUSTIMEEXPERIMENTS_H_
#define CONTINUOUSTIMEEXPERIMENTS_H_

#include "../common.h"
#include "../util.h"
#include <algorithm> //For Find
#include "PoissonSource.h"

extern int g_FilterTh;
void runAvgContinuousMemoryLifetimeSignalSimulation(int modelType,long trials,int trackedMemIndex,int CascadeSize,long synapsesPopulation,long SimTime,double dEncodingRate,string inputFile);


//Get Perceptron Signal But do not Use Synapses Pointed by track group
int testCUDAPRecallOfX(float* h_sigdata,t_inVal* W ,t_inVal** X,t_inVal* tX, uint _uiSynCount,void* d_W,void* d_X,void* d_C,t_patt_trackedtbl& vTrackedIndex,uint _uiPattsStoredSoFar);
void cleapUpCUDADeviceMem(void*& d_W,void*& d_X,void*& d_C,void*& d_odata, float*& h_odata, uint _uiSynCount);
void initCUDADeviceMem(void*& d_W,void*& d_X, void*& d_C, void*& d_odata, float*& h_odata,unsigned int _uiSynCount,uint TrackedCount);
void transferVectorsToDevice(int iNoTrackedPats,int* h_W ,t_inVal* h_X, uint _uiSynCount,void* d_W,void* d_X);

//int selectPatternToEncode(unsigned long j,uint uiNoOfPatternsStoredInTrial,bool bUseRandomPatterns,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng);
//unsigned long getNextTimestep(bool& bisRecallPeriod,bool& bisEncodingPeriod,double dTimestepSize,unsigned long currentTimestep,unsigned long& LastRecallj,unsigned long& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable);
uint getCurrentPeriodOfRecall(unsigned long j);
uint getCurrentPeriodOfRecall(unsigned long j,unsigned long peakSigTime);
double getNextTimestep(bool& bisEncodingPeriod,double dTimeStepSize,double currentTimestep,double& LastRecallj,double& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable);

int selectPatternToEncode(unsigned long j,uint uiNoOfPatternsStoredInTrial,bool bUseRandomPatterns,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng);
int makeInductionStep(float& mdRate, gsl_rng*& prng_r, double& dC);

//Get Perceptron Signal But do not Use Synapses Pointed by track group
template<class T>
int testPRecallOfX(double& _Signal,T* pCSyns, t_inVal* W ,t_inVal** X,uint _uiSynCount,uint _iStoredPatIndex)
{
	//Test Recall of _iStoredPatIndex
	//cout << "Recall index: " << _iStoredPatIndex << " Output should be :" << X[_iStoredPatIndex][_uiSynCount-1];
	double h=0.0;
	double hNTrack=0.0;

	//Do All Input vector - Assume All inputs have a required neuron output of +1
	for (uint i=0;i<(_uiSynCount);i++)
	{
			h+=W[i]*X[_iStoredPatIndex][i]; //Mu_Delta
	}

	//MEasure The normalized PostSynaptic Respose As Signal
	_Signal = h/_uiSynCount;

	//cout << X[_iStoredPatIndex][_uiSynCount-1]/(X[0][0]*X[0][0]*(_uiSynCount-1)) << endl;
	int iNeuronOut = ((h+hNTrack)>0)?1:-1; //Neuron Classifier OUtput


return iNeuronOut;
}

/* Encoding a pattern - MFPT Sim. Version
 * Handles the iteration through synapses to induce encoding While Encoding it also measures the Signal Of the 1st tracked Pattern To save on iterating again.
 *
 * t_inVal* X : The input vector
 *
 */
//Returns: 1 if Signal of 1st memory is above 0
template <class T>
uint encodeMemory(double& _Signal,vector<T*>& vpSyns,unsigned long ts,t_patt_trackedtbl& vTrackedIndex,uint uiNoOfPatternsStoredInTrial,t_inVal** X,t_inVal* W,int iPatIndex,gsl_rng*& prng_r)
{
	const uint cTrackedPatternIndex = 0; //Always return signal from 1st tracked pattern in this function
	_Signal = 0.0;
	double dC = 0.0;
	float mfRate = 1.0f;

	//uint iP_s = 0; uint iD_w = 0;

	//uint iSynCount = (uint)vpSyns.size();
	t_patt_trackedtbl::iterator itAtTrackedPatt;
	int TrackedIndex = -1;
	t_patt_trackedtbl::iterator itAtEnd = vTrackedIndex.end();
	t_patt_trackedtbl::iterator itAtStart = vTrackedIndex.begin();

	typename vector<T*>::iterator it; //Use typename Because T could be referring to a type or a member of class T

	//Search if this is a tracked pattern and obtain the index in the tracked list
	if (iPatIndex ==-1) //Not A repetition of A previous pattern? Check if this pattern is in the tracked list
	{
		int temp = 0;
		for (itAtTrackedPatt = itAtStart; itAtTrackedPatt!= itAtEnd;++itAtTrackedPatt) //Only Track 1st pattern
		{//Search if Current Pattern is tracked and map to stored pattern Index in X vector
			if (itAtTrackedPatt->first == uiNoOfPatternsStoredInTrial)
			{
				TrackedIndex = temp;
				break; //Set To the 1st Tracked Pattern
			}
			temp++;
		}
	}
///The PatIndex will be -1 when a tracked pattern is first Encountered - The trackIndex is set then to position the save of this pattern in the X vector

uint i = 0;
//Loop Through All Synapses -  //Skip This timestep if no encoding occurred bPatternArrived=false
for (it = vpSyns.begin();it!=vpSyns.end();++it)
{
	T* oCSyn = *it;


	//Whether File or Random Patterns
	//For Random Patts - The vector is here On the Fly by makeInduction Step
	int iStim;
	if (iPatIndex != -1) //If NOT random Pattern On the fly
	{
		//This is only used when selectPatternToencode Finds the repetition of a memory
		iStim = X[iPatIndex][i]; //*X[iSynCount-1]; //Assume Fixed Neuron Output of +1]
	}
	else
	{ //iPatIndex == -1 So Make pattern on the fly and store it
		//Create The Pattern One Stimulus At the time
		iStim = makeInductionStep(mfRate, prng_r, dC);
		if (TrackedIndex != -1) //But this is a tracked Pattern So Save it
		{
			X[TrackedIndex][i] = iStim; //Save The stimulus in the X vector holding The Tracked Patterns
		}
	}

	switch (iStim)
	{
		case 1:
			oCSyn->handlePOT();
		break;
		case -1:
			oCSyn->handleDEP();
		break;
		case 0:
			oCSyn->handleNOP(); //For "Sparseness"
		break;
		default:
			ERREXIT(500,"encodeMemory: Unknown Induction Stimulus")
		break;
	}

	W[i] = oCSyn->getStrength(); //Save Strength Into Weight Vector
	//Check Recall - Note: If Tracked pattern has not been stored yet then Signal Is Nonsense
	_Signal += W[i]*X[cTrackedPatternIndex][i];
	i++; //INcrement Index used for input vector
} //For Each Synapse

//Normalize Signal to Number of synapses - Note: If Tracked pattern has not been stored yet then Signal Is Nonsense
_Signal = _Signal/(double)i;

//cout << "**N:" << iSignalSynsCount << endl;
return (_Signal > 0.0);
}




//#undef USE_CUDA
/*
 * MEAN FIRST PASSAGE TIME
 * Measures Covariance of overall signal, Tracked and non-Tracked Signal.
 * Repeats a particular tracked memory at a specific interval if added to the repetitionTable.
 * Note:
 * For Random Patterns @ The Patterns Are recreated at every trial and stored in X array. The index is then used to track the signal of each pattern
 * At the Beginning of each trial there is no RESET of the synapses. It is assumed that the patterns stored will be new and stored ontop of the distribution
 * An Init Period Can be added use pattern to re-init the distribution among filter states if required. Make Sure Sufficient Init Period Exists
 * Populates a trackedIndex list with vTrackedIndex; //Key:The PattIndex - Value: 1 For Allocation Signal / 0 For no Allocation
 */
template <class T,class pFunct>
double simMFPT(double _oMFPT[],double _oMFPTVar[],pFunct pF, uint iSynCount,int iCascadeSize,uint uiInitPatterns,uint trials,double dEncodingRate,t_patt_reptbl& repetitionTable,double ts)
{
#ifdef USE_CUDA
	cout << "Using CUDA For Signal Dot Prod" << endl;
#else
	cout << "CPU For Signal" << endl;
#endif
	cout << "Mean First Passage time of 1st Tracked Pattern with repetition of tracked Signal -Net Size : " << iSynCount << " Trials:" << trials <<endl;

	const double cdLowSignalThres		= 0.0;
	uint iTrackedMemIndex		= 0;//Holds the current Tracked index when Iterating through the table
	bool bUseRandomPatterns		= true; //If no file is given Use Random Patterns
	//uint maxMemCapacity			= dEncodingRate*simTimeSeconds*2; //Max Number of Patterns to create - x2 for safety
	gsl_rng* mprng 				= g_getRandGeneratorInstance(true);

	//MEAN FIRST PASSAGE TIME VARIABLES
	const int cLowSigThresCount = 8;
	double dLowSignalMFPT[cLowSigThresCount]; //The point Where the first pattern tracked drops below noise threshold
	double dMFPTVar[cLowSigThresCount];
	int iCountThresholdsReached = 0;
	memset(dLowSignalMFPT,0.0,cLowSigThresCount*sizeof(double));
	memset(dMFPTVar,0.0,cLowSigThresCount*sizeof(double));

	char *mem_buffer 			= 0;		//This Pointer is filled by allocMem, To point to the reuseable allocated memory
	T* oCSyn;	//POinter to allocated memory
	vector<T*> vpSyns; //Vector of pointers to Synapses
	vpSyns.reserve(iSynCount);
	PoissonSource* PsMemEvent = new PoissonSource(dEncodingRate,ts,0);

	//ADD THE LIST OF TRACKED PATTERNS
	t_patt_trackedtbl vTrackedIndex; //Key:The PattIndex - Value: 1 For Allocation Signal / 0 For no Allocation
	t_patt_trackedtbl::iterator itTracked;
	t_patt_reptbl::iterator itRep;

	for (int i=0;i<=0;i+=1) //Start tracking from 1 cause 0th falls within the init patterns
		vTrackedIndex[i+uiInitPatterns] = 0; //Add Tracked Indexes to the List of tracked Patterns
	//ENDOF LIST OF TRACKED PATTERNS

	//Track All patterns That Are in the Repetition table
	unsigned long maxRepTime = 0;
	for (itRep = repetitionTable.begin(); itRep != repetitionTable.end();++itRep  )
	{
		vTrackedIndex[(*itRep).second] = (*itRep).second;  ; //Add Repetition Indexes to the List of tracked Patterns
		if (itRep->first > maxRepTime )
				maxRepTime = itRep->first;
	}


	const int ciNoOfTrackedPatterns = vTrackedIndex.size()+1;

	vector<string> slogFiles; //The list of output file names used
	//set<unsigned long> vTs; //Contains All the time point were Strength Has been Measured#
	t_inVal* 	 X[ciNoOfTrackedPatterns]; //Memory Patterns Containing The Ones Loaded from File and Random Initialization patterns
	t_inVal* 	 W = new t_inVal[iSynCount]; //Weight Vector Reflecting The state of the Synapses
	///Reserve Memory For Storing Measurements for each tracked pattern

    cout << "Number of Tracked patterns :" << (ciNoOfTrackedPatterns-1) << " Init Pat:" << uiInitPatterns << endl;
    if (ts >= 1.0 ) cout << "DISCRETE TIME"; else cout << "CONTINUOUS TIME";

	long* 	 iOccupancy[2]; //Store counters for each Strength-State Pair - Distribution

	uint t = 0; //Trial counter

	cout << " Event driven Simulation measuring signal. Synapse Init Patterns: " << uiInitPatterns  << endl;
	memset(X,0,ciNoOfTrackedPatterns*sizeof(t_inVal*)); //Setting The pointers to 0 so InitPatterns Knows to initiaze each member to a new array
	memset(W,0,iSynCount*sizeof(t_inVal)); //Setting The pointers to 0 so InitPatterns Knows to initiaze each member to a new array

	iOccupancy[0] = new long[15];
	iOccupancy[1] = new long[15];
	memset(iOccupancy[0],0,sizeof(long)*15); //Empty the memory Buffer
	memset(iOccupancy[1],0,sizeof(long)*15); //Empty the memory Buffer

	//Initialise The memory For the Tracked Patterns
	initPatternMemory(X,ciNoOfTrackedPatterns,iSynCount,iTrackedMemIndex,0.5, mprng,bUseRandomPatterns);

	/////////// LOG FILE INIT /////////////////////
	//Add the File name as the 1st entry- Used by the makeLogFileNames
	slogFiles.clear();
	string fOutName(MFPTIMES_OUTPUT_DIRECTORY);
	slogFiles.push_back(fOutName);
	makeLogFileNames<T>(slogFiles,uiInitPatterns,iCascadeSize,dEncodingRate, 0.5,trials, iSynCount,pF);
	/////////// END OF LOG FILE INIT //////////////////////////

	//ALLOCATE SYNAPSES Call Object Allocation Function - Creates Object in memory and returns pointer to 1st object
	oCSyn = (T*)(*pF)((char*)mem_buffer,iSynCount,(int)iCascadeSize,mprng,1.0);
	if (!oCSyn){//Check Failure
				ERREXIT(500,"simMFPT: Could not create synapse objects! Out Of memory?");
		}else //Save pointer
			mem_buffer = (char*)oCSyn;
	///Save pointer to a Vector only on 1st trial- Temporary Solution But Its Required
	//The pointers do not change between re-allocations
	vpSyns.clear();
	for (uint i=0;((i<iSynCount)) ;i++) //Could Do it only 1st time, POinters remain the same from then on
	{
		assert(&oCSyn[i] != NULL);
		vpSyns.push_back(&oCSyn[i]);
	}
	//Allocating Once And then Changing the Tracked Pattern on everytrial - No need to recreate Synapses - Just assume they have already been initialized by past activity
	//DistInit - Report - Only on 1st Trial
	reportStateDistribution<T>(vpSyns,iSynCount,slogFiles[0].c_str());

	uint8_t reportCycle = 0; //When Overflows we report the trial number
	///Do a Trial - First trial will allocate the memory buffer -
	t = trials;
	double dmeanSignal = 0.0;

	while (t > 0) //Start Trial Loop
	{
		//Reset the MFPT Count of  recorded a low signal conditions
		iCountThresholdsReached = 0;

		reportCycle ++;
		if (reportCycle ==0) cout << t << endl;

		uint uiNoOfPatternsStoredInTrial = 0;
		int  iPatIndex					 = 0;//The Pattern Index in focus at any particular timestep
		bool bPatternArrived 			 = false;

		double j = 0;
		double NextEncodingj = 0;
		double LastRecallj = 0;
		double lFPT = 0; //First Passage Time below 0 For this Trial
		bool bThresholdReached = false;
		double fSignalTracked = 0.0; // Signal Of 1st Pattern

		//Loop Jumping between timesteps until the signal from the last tracked memory goes below zero
		while(!bThresholdReached)
		//while(lFPT==0)
		{
			//State Machine code This Code Overrides any previous Decision on which event should occur So to accommodate the Initialisation Period
			if (uiNoOfPatternsStoredInTrial <= (uiInitPatterns))
			{
				bPatternArrived = true;
				j = NextEncodingj = 0.0; //Stuck In time Until Init Patterns and the Tracked pattern have been delivered
			} //Rest Of Timestep Control Occurs At the end of the loop

			//Encode Pattern
			if (bPatternArrived)
			{	//Select The pattern to Encode - Chooses whether pattern is to be repeated Or Ret:-1 so encodeMemory Knows to create a new Pattern on the fly
				iPatIndex = selectPatternToEncode(j,uiNoOfPatternsStoredInTrial,bUseRandomPatterns,vTrackedIndex,repetitionTable,X,ciNoOfTrackedPatterns,iSynCount,mprng);

				//Run Through Synapses Encoding Induction Signal - if iPatIndex =-1 then induction stimuli are generated on the fly
				 encodeMemory<T>(fSignalTracked,vpSyns,j,vTrackedIndex,uiNoOfPatternsStoredInTrial,X,W,iPatIndex,mprng);//Return 1 if pattern was found

				if (t==trials) //Check If it is time to Report Encoding Event
				{
					if ((vTrackedIndex.find(uiNoOfPatternsStoredInTrial)) != vTrackedIndex.end()) //Decrement One Since If Encoded It would Have increased by one
					{//After Mem Storage Distribution
						cout << endl << "T:" << t << " Tracked Memory " << uiNoOfPatternsStoredInTrial  << endl;
						reportStateDistribution<T>(vpSyns,iSynCount,slogFiles[2].c_str()); //Dist B
					}
				}

				uiNoOfPatternsStoredInTrial++;
			}//Finished Looping through all synapses - Pattern is now stored

			//MEASURE SIGNAL Obtain Distributions and Measure the Signals for all timesteps -
			for (int i=(cLowSigThresCount-iCountThresholdsReached-1);i >= 0;i--)
			{
				//Check From Highest to Lowest Threshold
				if ((fSignalTracked <= (float)(i)/100.0) && (uiNoOfPatternsStoredInTrial > uiInitPatterns)){ //If Not recorded Before THen Save time point Where signal Falls below 0
					lFPT = j;
					//dmeanSignal += fSignalTracked;
					dLowSignalMFPT[i] += lFPT; // Cummulate Mean FPT
					dMFPTVar[i]	+= lFPT*lFPT;  //Cummulate Variance
					iCountThresholdsReached++; //increment The number of thresholds crossed

					if (iCountThresholdsReached == cLowSigThresCount ) bThresholdReached = true; //Stop Loo
				//cout << (trials-t) << "-" << lFPT*ts << endl;
				}
			}

			//Go to Next Time step - The number of Patterns is incremented - The bool vars for Recall Or Encoding change depending on the type of the next event
			j = getNextTimestep(bPatternArrived,ts,j, LastRecallj,NextEncodingj, PsMemEvent,repetitionTable);

		}//Until a value is set to lFPT the time of reaching 0 -

		t--; //Decrement Trial count down to 0

		//dLowSignalMFPT += lFPT; // Cummulate Mean FPT
		//dMFPTVar	+= lFPT*lFPT; //Cummulate Variance
	}//LOOP For each trial

	cout << "---------------------------------------------------------" << endl;
///RECORD STATISTICS FROM EACH TRACKED PATTERN
//	char buffFilename[400];
//	int i = 0; //TrackedMem INdex increment
	//MFPT - Copy To output Vars
	for (int i=(cLowSigThresCount-1);i >= 0;i--)
	{
		_oMFPT[i] 	  = dLowSignalMFPT[i] = dLowSignalMFPT[i]/trials;
		_oMFPTVar[i] = dMFPTVar[i]	   = dMFPTVar[i]/(double)trials - dLowSignalMFPT[i]*dLowSignalMFPT[i]; //Calc Variance
		cout << "L.Sig :" << (float)(i/100.0) << " MFPT :" << dLowSignalMFPT[i] << " Variance :" << dMFPTVar[i] << endl;
	}
		cout << "Mean signal crossed Threshold "<< cdLowSignalThres << " with : " << dmeanSignal << endl;

//CLEAN UP///
	vpSyns.clear();
	return_temporary_buffer(mem_buffer);

	for (int i=0;i<ciNoOfTrackedPatterns;i++) //Remove All saved Patterns
		delete [] X[i];

	delete [] W;
	delete iOccupancy[0];
	delete iOccupancy[1];
	delete PsMemEvent;

	//CleaUp CUDA
#ifdef USE_CUDA
	cleapUpCUDADeviceMem(d_W,d_X,d_C,d_odata,h_odata,iSynCount);
	delete [] tX; //The Long Buffer Vector Of Tracked PAtterns
#endif
	cout << "-Fin-" << endl;

	return dLowSignalMFPT[0]; //Return the 0 Signal Cross Time
}




//#undef USE_CUDA

#endif /* CONTINUOUSTIMEEXPERIMENTS_H_ */
