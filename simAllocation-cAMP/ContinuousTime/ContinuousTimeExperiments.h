/*
 * ContinuousTimeExperiments.h
 *
 *  Created on: 7 Oct 2011
 *  Modified in Mar/2012 :
 *  This code has been modified to provide the statistics of same threshold crossings distributions under repetitious memory encoding
 *  the results of wk13-14 2012 has been obtained using this code. The sampling cut off is defined by the input variable
 *
 *      Author: kostasl
 */

#ifndef CONTINUOUSTIMEEXPERIMENTS_H_
#define CONTINUOUSTIMEEXPERIMENTS_H_

#include "../common.h"
#include "../util.h"
#include <algorithm> //For Find
#include "PoissonSource.h"
#include "../InputVectorHandling.h"

extern int g_FilterTh;
extern float g_fAllocHThres;
extern float g_fcAMPDecay;
extern float g_fcAMPMagnitude;
extern uint g_timeToSampleMetaplasticity;
extern int g_MetaplasticitySampleSize;
extern uint g_AllocRefraction;

//void runAvgContinuousMemoryLifetimeSignalSimulation(int modelType,long trials,int trackedMemIndex,int CascadeSize,long synapsesPopulation,long SimTime,double dEncodingRate,string inputFile);

//Get Perceptron Signal But do not Use Synapses Pointed by track group
int testCUDAPRecallOfX(float* h_sigdata,t_inVal* W ,t_inVal** X,t_inVal* tX, uint _uiSynCount,void* d_W,void* d_X,void* d_C,t_patt_trackedtbl& vTrackedIndex,uint _uiPattsStoredSoFar);
void cleapUpCUDADeviceMem(void*& d_W,void*& d_X,void*& d_C,void*& d_odata, float*& h_odata, uint _uiSynCount);
void initCUDADeviceMem(void*& d_W,void*& d_X, void*& d_C, void*& d_odata, float*& h_odata,unsigned int _uiSynCount,uint TrackedCount);
void transferVectorsToDevice(int iNoTrackedPats,int* h_W ,t_inVal* h_X, uint _uiSynCount,void* d_W,void* d_X);

//int selectPatternToEncode(unsigned long j,uint uiNoOfPatternsStoredInTrial,bool bUseRandomPatterns,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng);
//unsigned long getNextTimestep(bool& bisRecallPeriod,bool& bisEncodingPeriod,double dTimestepSize,unsigned long currentTimestep,unsigned long& LastRecallj,unsigned long& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable);
uint getCurrentPeriodOfRecall(unsigned long j);
uint getCurrentPeriodOfRecall(unsigned long j,unsigned long peakSigTime);

unsigned long getNextTimestep(bool& bisRecallPeriod,bool& bisEncodingPeriod,double dTimeStepSize,unsigned long currentTimestep,unsigned long& LastRecallj,unsigned long& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable);
int selectPatternToEncode(unsigned long ts,uint uiNoOfPatternsStoredInTrial,bool& bPatternIsNew,bool& bAllocatePattern,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng);
int makeInductionStep(float& mdRate, gsl_rng*& prng_r, double& dC);

//Get Perceptron Signal But do not Use Synapses Pointed by track group
//Set Allocation On If h Is above threshold
template<class T>
int testPRecallOfX(double& _SignalNTracked,double& _SignalTracked,T* pCSyns, t_inVal* W ,t_inVal** X,int* TrackedSynFlag,uint _uiSynCount,uint _iStoredPatIndex)
{
	//Test Recall of _iStoredPatIndex
	//cout << "Recall index: " << _iStoredPatIndex << " Output should be :" << X[_iStoredPatIndex][_uiSynCount-1];
	double h=0.0;
	double hNTrack=0.0;

	//Do All Input vector - Assume All inputs have a required neuron output of +1
	for (uint i=0;i<(_uiSynCount);i++)
	{
		//cout << w[i] << " " << " " << X[_iStoredPatIndex][i] << endl;
		if (TrackedSynFlag[i] == 1) // Only Measure the Post Response due to the Tracked Synapse
		{
			h+=W[i]*X[_iStoredPatIndex][i]; //Mu_Delta
		}else //Non Tracked
			hNTrack+=W[i]*X[_iStoredPatIndex][i]; //Mu_Omega
	}

	//MEasure The normalized PostSynaptic Respose As Signal
	_SignalTracked = h/_uiSynCount;
	_SignalNTracked = hNTrack/_uiSynCount;//X[_iStoredPatIndex][_uiSynCount-1]*hNTrack/(X[0][0]*X[0][0]*(_uiSynCount-1));
	//cout << X[_iStoredPatIndex][_uiSynCount-1]/(X[0][0]*X[0][0]*(_uiSynCount-1)) << endl;
	int iNeuronOut = ((h+hNTrack)>0)?1:-1; //Neuron Classifier OUtput

	//cout << " Signal: " << _SignalTracked << " - " << _SignalNTracked  << endl;

//	if (X[_iStoredPatIndex][_uiSynCount-1] == iNeuronOut)
//	{
//		Ret = 1; //Return 1 To indicate Successful classification of input
//	}

return iNeuronOut;
}




/*
 * Handles the iteration through synapses to induce encoding
 * It also detects which synapses have been modified due to a tracked memory
 * t_inVal* X : The input vector
 * The Synapses Are marked as monitored for the 1st Tracked pattern
 * Output: uiAllocSynapses -The number of allocated-frozen Synapses
 */
//Returns: The number of synapses that encoded the signal
template <class T>
uint encodeMemory(vector<T*>& vpSyns,unsigned long ts,t_patt_trackedtbl& vTrackedIndex,uint uiNoOfPatternsStoredInTrial,t_inVal** X,t_inVal* W,
				 int** iTrackedSynFlag,	int iPatIndexRelativeToX, bool& bEncodeNewPattern, bool& bAllocatePattern, gsl_rng*& prng_r,
				 uint& uiAllocSynapses, map<uint,uint>& pMDistrib, map<uint,uint>& pMDistribinSamples)
{
	double dC = 0.0;
	double h  = 0.0; //Neuron Depolarization
	float mfRate = 1.0f;
	uint iSynCount = vpSyns.size();
	uint iSignalSynsCount = 0; //Reset the Memory Synapses Count

	bool bAllocationSignal = bAllocatePattern;
	bool bresetAllocationThreshold = false;

	int TrackedIndex = iPatIndexRelativeToX;
	t_patt_trackedtbl::iterator itAtTrackedPatt;
	typename vector<T*>::iterator it; //Use typename Because T could be referring to a type or a member of class T
	///1st PERIOD - Before Storage - Start to monitor all synapse Changes

/*
	If a pattern is repeated Then it has been saved in the X vector Before Because it would have been a tracked pattern too */
	if (bEncodeNewPattern && iPatIndexRelativeToX != -1) //Pattern Is tracked and this is the first Encoding
	{
		bresetAllocationThreshold = true;
		//Reset The Time-limited Threshold Counter Distribution Before the storage of the 1st Tracked Memory
		pMDistrib.clear(); //This Distribution Is saved When the TimeLimit Is reached and Reset again at everyTrial
	}

uint uiAllocCounter = 0;
uint i = 0;
//Loop Through All Synapses -  //Skip This timestep if no encoding occurred bPatternArrived=false
for (it = vpSyns.begin();it!=vpSyns.end();++it)
{
	T* oCSyn = *it;
	//Only reset the state of non Monitored Synapses By some previous Tracked Pattern
	if (TrackedIndex >= 0) //Reached Tracked Pattern - Set State As reference
		oCSyn->startMonitoring(); //Set current state as Reference So we can detect changes at the synapse (Used to happen just before tracked index)

	//Begin new Trial ?
	if (ts == 0 && uiNoOfPatternsStoredInTrial==0){ //Start Of New Trial? Unfreeze Synapses
		//Unlock Plasticity At the beginning - But should we start from Random Point?  No Need cause Plasticity lock does not stop internal filter transitions only strength state
		oCSyn->reset(); // Unlock Plasticity,Randomize Strength/;
		oCSyn->enableMetaplasticCounting(); //Obtain Histogram Over Whole Trial
		W[i] = oCSyn->getStrength(); //Reset the weight Vector Too To match Synapses
	}

	if (bresetAllocationThreshold) //Just Before 1st Encoding of tracked Pattern
	{
		//oCSyn->resetMetaplasticCounter();
		oCSyn->enableMetaplasticCounting();

#ifdef _DEBUG
		if (i==0)
			cout << oCSyn->getMetaplasticCount() <<  " 1st Encoding of Tracked Pattern. Reset Cycle Counters before Encoding" << endl;
#endif
	}

	if (bAllocationSignal)
	{//At the allocation signal we begin measuring Stability
		oCSyn->setAllocationThreshold(g_AllocRefraction); //When same threshold is crossed X times and allocation is On-Allocate
		oCSyn->enableMetaplasticAllocation();
		//oCSyn->resetAllocationRefraction(); //Reset The internal Counter
	}

	//Whether File or Random Patterns
	int iStim;
	//Found A repeated Pattern - Use TrackedIndex TO obtain saved pattern in X
	if (!bEncodeNewPattern && iPatIndexRelativeToX != -1) //If NOT random Pattern On the fly The Read the stimulus from the X vector
	{	//For Random vectors This is only used when selectPatternToencode Finds the repetition of a memory
		iStim = X[iPatIndexRelativeToX][i]; //*X[iSynCount-1]; //Assume Fixed Neuron Output of +1
	}
	else
	{ //iPatIndex == -1 So Make pattern on the fly and store it
	 	iStim = makeInductionStep(mfRate, prng_r, dC); //Create The Pattern One Stimulus At the time

	 	if (bEncodeNewPattern) //But this is a tracked Pattern So Save it
			X[iPatIndexRelativeToX][i] = iStim; //Save The stimulus in the X vector holding The Tracked Patterns
	}

	switch (iStim)
	{ //Set What the Correct Strength State is For Each Synapse So they Can Save to the distribution
		case 1:
			//Make Sure Tracked State is called Only Once Per Trial to label the Correct threshold Crossing
			if (bresetAllocationThreshold) oCSyn->setTrackedState(ICascadeSynapse::SYN_STRENGTH_STRONG); //Set The "desired" strength State
			oCSyn->handlePOT();
		break;
		case -1:
			if (bresetAllocationThreshold) oCSyn->setTrackedState(ICascadeSynapse::SYN_STRENGTH_WEAK);
			oCSyn->handleDEP();
		break;
		case 0:
			if (bresetAllocationThreshold) oCSyn->setTrackedState(ICascadeSynapse::SYN_STRENGTH_NOTSET);
			oCSyn->handleNOP(); //For "Sparseness"
		break;
		default:
			ERREXIT(500,"encodeMemory: Unknown Induction Stimulus")
		break;
	}

	//if (i==0 && TrackedIndex >= 0)
			//cout << "i:" << oCSyn->getRunningValue() << " " << oCSyn->getMetaplasticCount() <<  " Set Tracked State S:"<< iStim << endl;

	h	+= W[i]*iStim; //Get Neural Output Due to Pattern Being Stored - Before Any Weights are updated
	W[i] = oCSyn->getStrength(); //Save new Strength  due to Encoding Into Weight Vector

	if (!oCSyn->isPlastic()) //Count Locked Synapses
		uiAllocCounter++;

	//2nd PERIOD  If this synapse is tracked then Set the flag if this synapse is a mu_delta or M_omega one -
	if (TrackedIndex >= 0)
	{
		oCSyn->stopMonitoring(); //Stop Tracking all synapses - To refine the tracked group
	// Form Memory Signal if this is the 1st run and the Synapse has changed State or Index
		if ((oCSyn->hasStrengthModified() || oCSyn->hasIndexModified()))
		{
			iTrackedSynFlag[TrackedIndex][i] = 1; //Mark As MU_DELTA Synapse - On repetition of pattern this would refine the Mu_Delta Group
			//OLD Way of a marking Synapse:oCSyn->startMonitoring(); //Record new Start Strength to monitor Change - Monitored synapses belong to the Initial Signal
			iSignalSynsCount++;//Now Save How many have changed strength to Save equilibrium point
		}//If Synapse Modified and At Tracked Memory- Add To Signal Set
		else
		{
			iTrackedSynFlag[TrackedIndex][i] = -1; //unmark Synapse for this Tracked Group
		}
	}//Reached Tracked Pattern

	//Check if Limit on the number of ThresCycles has been obtained And Save a snapshot of the distribution - Runs Only Once
	//g_MetaplasticitySampleSize*iSynCount Wait Until All Synapses Have completed their required number of Cycle Samples
	if (oCSyn->getMetaplasticDistributionSampleSize() == g_MetaplasticitySampleSize*iSynCount && pMDistribinSamples[0] == 0 && g_MetaplasticitySampleSize > 0)
	{ //Copy To Sample Limited Distribution
		pMDistribinSamples.insert(pMDistrib.begin(),pMDistrib.end());//Copy Snapshot of current Distrib IN time
		//pMDistribinSamples[0] = pMDistrib[0]; //This is Used in the Time loop To know if the Sample Has been obtained
		pMDistribinSamples[0] = ts+1; //Save the Time Of When The sample Limit Was Reached
		cout << "ThreshCycle Sample :" << g_MetaplasticitySampleSize << " Per Synapse Reached At:" << pMDistribinSamples[0] << endl;
	}

	if (pMDistribinSamples[0] != 0)
		oCSyn->disableDistributionSampleLimit(); //Once the distribution is Copied - Remove Sample Size Constraint

	i++; //INcrement Index used for input vector
} //For Each Synapse

h = h/iSynCount;

//Fixed Allocation Threshold on h. h is measured on currently encoded memory before plasticity changes due to its encoding have been expressed
if (h > g_fAllocHThres) {
	bAllocatePattern = true;
	//cout << ts << " allocation on" << endl;
}

//cout << "**N:" << iSignalSynsCount << endl;
uiAllocSynapses = uiAllocCounter;
//cout << "t:" << ts << " NAlloc:" << uiAllocCounter << endl;
return iSignalSynsCount;
}



//Copies The Time Trial Distribution To the Simulation Aggregate for the Avg Distribution Over a Fixed Time
template<class T>
void saveMetaplasticCounters(vector<T*>& vpSyns,map<uint,uint>& mpMDistribDst,map<uint,uint>& mpMDistribSrc)
{
	//Call Each Synapse To deliver Its Current Running State
	for (typename vector<T*>::iterator it = vpSyns.begin(); it != vpSyns.end(); ++it)
	{
		//uint cnt = (*it)->getMetaplasticCount();
		//cout << cnt << endl;
		(*it)->saveMetaplasticDistribution(); //Synapses Save into The mpMDistribSrc
	}

	//Save Running Distribution Into Aggregate
	for (typename map<uint,uint>::iterator it = mpMDistribSrc.begin(); it != mpMDistribSrc.end(); ++it)
	{
		//cout << it->first << endl;
		mpMDistribDst[it->first] += it->second; //At the samples of the new Trial To the previous Ones
	}

}

//#undef USE_CUDA
/*
 * Mean Signal Memory lifetime in continuous time Jumping In time Between Events
 * Measures Covariance of overall signal, Tracked and non-Tracked Signal.
 * Repeats a particular tracked memory at a specific interval if added to the repetitionTable.
 * Note:
 * For Random Patterns @ The Patterns Are recreated at every trial and stored in X array. The index is then used to track the signal of each pattern
 * At the Beginning of each trial there is no RESET of the synapses. It is assumed that the patterns stored will be new and stored ontop of the distribution
 * An Init Period Can be added use pattern to re-init the distribution among filter states if required. Make Sure Sufficient Init Period Exists
 * Populates a trackedIndex list with vTrackedIndex; //Key:The PattIndex - Value: 1 For Allocation Signal / 0 For no Allocation
 * Returns: The final SNR of the last tracked Pattern
 */
template <class T>
pair<double,double> simRepetitionAllocation(T* oCSyn, uint iSynCount,int iCascadeSize,uint uiInitPatterns, char* pinputFile,uint trials,uint simTimeSeconds,double dEncodingRate,t_patt_reptbl& repetitionTable,double ts,vector<string>& slogFiles)
{
#ifdef USE_CUDA
	cout << "Using CUDA For Signal Dot Prod" << endl;
#else
	cout << "CPU For Signal" << endl;
#endif
	cout << "C.T Signal Lifetime With repetition of tracked Signal -Net Size : " << iSynCount << " Trials:" << trials <<endl;
	//const double ts				= 1.0;//0.001;//dEncodingRate/2.0; //TimeStep is sampling at twice the frequency
	gsl_rng* mprng 				= g_getRandGeneratorInstance(false);
	//const unsigned long peakSig  =g_FilterTh*g_FilterTh*tsPerSec*0.375;
	const unsigned long lTotalTimesteps	= (double)simTimeSeconds/ts; //The Total Number of TimeSteps in the simulation- (Number of Samples )
	uint iSignalSynsCount 		= 0; //Signal Variable
	uint uiAllocSynapses		= 0; //Stores the count of Allocated Synapses returned by encodeMemory
	uint iTrackedMemIndex		= 0;//Holds the current Tracked index when Iterating through the table
	bool bShuffleLoadedVectors 	= true;
	bool bUseRandomPatterns		= (pinputFile == 0) || (pinputFile[0]=='\n'); //If no file is given Use Random Patterns
	uint maxMemCapacity			= dEncodingRate*simTimeSeconds*2; //Max Number of Patterns to create - x2 for safety
	uint uiPatCount 			= uiInitPatterns + maxMemCapacity; //Max Patterns loaded from File-Valid Only File- And Thus Max Storage Capacity is NeuronCount/5

	set<unsigned long> vTs; //Contains All the time point were Strength Has been Measured#

	// DOPAMINE SIGNAL //
	double dDALevel 					= 0.0; //The current level of Dopamine-initiated signal
	double dDAInjected 					= 0.0; //The current level of Dopamine-initiated signal
	const int iNoOfDAInjections			= 7;

	unsigned long ulLastDASignalTime[iNoOfDAInjections];//''	= 0.0; //Injection Times of DA Signal

	const double dDAStep		= g_fcAMPMagnitude;
	const double dDAThreshold 	= 1.0; //Point after which dDALevel

	// MEAN FIRST PASSAGE TIME VARIABLES //
	double dLowSignalMFPT 		= -1.0; //The point Where the first pattern tracked drops below noise threshold
	double dLastSignalValue 	= -1.0; //The signal Value at the end of simulation is assumed to be the allocated signal
	double dLastVarValue 		= -1.0; //The signal Value at the end of simulation is assumed to be the allocated signal

	PoissonSource* PsMemEvent = new PoissonSource(dEncodingRate,ts,0.0,mprng);

	// ADD THE LIST OF TRACKED PATTERNS //
	t_patt_trackedtbl vTrackedIndex; //Key:The PattIndex(No Of Patterns Stored Up To Tracked one) - Value: 1 For Allocation Signal / 0 For no Allocation
	t_patt_trackedtbl::iterator itTracked;
	t_patt_reptbl::iterator itRep;

	// HISTOGRAM of ThresholdCycles //
	map<uint,uint> mpMDistrib; //The distribution Of Metaplastic Transitions Used by All Synapses -Live
	map<uint,uint> mpMDistribinSamples; //The Avg distribution Of Metaplastic Transitions among Synapses for a fixed number of cycle samples - Sample Limited
	map<uint,uint> mpMDistribinTime; // Time Limited Average over All trials - Contains the distribution of Synapses At the cut-off point too

	mpMDistribinTime.clear();
	mpMDistrib.clear();
	mpMDistribinSamples.clear();
	mpMDistribinSamples[0] = 0;

	//ADD REPEATED PATTERNS TO TRACKED LIST AND MARK AS ALLOCATED - Track All patterns That Are in the Repetition table
	unsigned long maxRepTime = 0;
	//	vTrackedIndex[i+uiInitPatterns] = 0;  ; //Setting to 0 means Disable Allocation / 1 Means Allocate
	for (itRep = repetitionTable.begin(); itRep != repetitionTable.end();++itRep  )
	{
		///Tracked pattern Is allocated - Give Allocation Signal
		vTrackedIndex[(*itRep).second] = 1;  //Add Repeated Mem Index to the List of tracked Patterns
		if (itRep->first > maxRepTime)
			maxRepTime = itRep->first; //Save the Last Repetition Time into MaxRep
	}

	//vTrackedIndex[13] = 0; //Just Testing
	unsigned long cSampleMetaplasticCounters 	= g_timeToSampleMetaplasticity; //Sample the metaplastic Counters At Fixed Time
	cout << "Sampling Metaplastic Counters at :" << cSampleMetaplasticCounters << endl;

	const int ciNoOfTrackedPatterns = vTrackedIndex.size();
	if (bUseRandomPatterns) //For On The Fly Vectors
		uiPatCount = ciNoOfTrackedPatterns+1; //No Need to have Large X vectors - Only Tracked patterns Saved

	assert(simTimeSeconds*(1/ts) > maxRepTime); //Check If the Simulation Time Allows for the number of repetitions

	t_inVal* 	 X[uiPatCount]; //Memory Patterns Containing The Ones Loaded from File and Random Initialization patterns
	t_inVal* 	 W = new t_inVal[iSynCount]; //Weight Vector Reflecting The state of the Synapses
	///Reserve Memory For Storing Measurements for each tracked pattern
	double**	 dVar 	= new double*[ciNoOfTrackedPatterns]; //Memory For Variance Vector Of Each Tracked Pattern
	double**	 dPtSn 	= new double*[ciNoOfTrackedPatterns]; //THe Perceptron Signal of the tracked pattern at each memory storage step
	double**  	 dPntSn = new double*[ciNoOfTrackedPatterns]; //THe Perceptron Signal of Non-tracked pattern measured on the non tracked synapses set
	double** 	 dDASn	= new double*[ciNoOfTrackedPatterns]; //The Timecourse of the DA signal
	int**		 iTrackedSynFlag = new int*[ciNoOfTrackedPatterns]; //Holds An array of flags indicating if a synapse is part of the DElta or the Omega group in a tracked pattern
	uint** 		 iAllocSynapses	 = new uint*[ciNoOfTrackedPatterns];

    cout << "Number of Tracked patterns :" << ciNoOfTrackedPatterns << endl;
	//Init Time Slots For each Tracked Pattern
	for (int i=0;i<ciNoOfTrackedPatterns;i++)
	{
		dVar[i] 	= new double[lTotalTimesteps]; //Memory For Variance Vector for this Particular tracked Pattern
		dPtSn[i]	= new double[lTotalTimesteps]; //Tracked Signal
		dPntSn[i]	= new double[lTotalTimesteps]; //Non Tracked Signal
		dDASn[i]	= new double[lTotalTimesteps]; //The DA signal

		iAllocSynapses[i]	= new uint[lTotalTimesteps]; //The Number of allocated synapses at a particular time step
		iTrackedSynFlag[i] 	= new int[iSynCount]; //Array the size of Number of Synapses each element holding 1 for Mdelta or -1 for Omega

		memset(dPtSn[i],0,lTotalTimesteps*sizeof(double));
		memset(dPntSn[i],0,lTotalTimesteps*sizeof(double));
		memset(dDASn[i],0,lTotalTimesteps*sizeof(double));
		memset(dVar[i],0,lTotalTimesteps*sizeof(double));
		memset(iAllocSynapses[i],0,lTotalTimesteps*sizeof(uint));
		memset(iTrackedSynFlag[i],0,iSynCount*sizeof(int));
	}

	long  iOccupancy[2][15]; //Store counters for each Strength-State Pair - Distribution
	uint t = 0; //Trial counter

	cout << "Event driven Continuous Time Simulation measuring signal after "<< uiInitPatterns << "memories using Perceptron measures " << endl;
	memset(X,0,uiPatCount*sizeof(t_inVal*)); //Setting The pointers to 0 so InitPatterns Knows to initiaze each member to a new array
	memset(W,0,iSynCount*sizeof(t_inVal)); //Setting The pointers to 0 so InitPatterns Knows to initiaze each member to a new array

	//iOccupancy[0] = new long[15];
	//iOccupancy[1] = new long[15];
	memset(iOccupancy[0],0,sizeof(long)*15); //Empty the memory Buffer
	memset(iOccupancy[1],0,sizeof(long)*15); //Empty the memory Buffer

	//Initialise The memory For Patterns
	initPatternMemory(X,uiPatCount,iSynCount,iTrackedMemIndex,0.5, mprng,bUseRandomPatterns);

	//ALLOCATE SYNAPSES Call Object Allocation Function - Creates Object in memory and returns pointer to 1st object
	if (!oCSyn) //Check Failure
		ERREXIT(500,"simMemRepetitionAllocation: Could not create synapse objects! Out Of memory?");

	///Save pointer to a Vector only on 1st trial- Temporary Solution But Its Required
	//The pointers do not change between re-allocations
	vector<T*> vpSyns; //Vector of pointers to Synapses
	vpSyns.reserve(iSynCount);
	vpSyns.clear();
	for (uint i=0;((i<iSynCount)) ;i++) //Could Do it only 1st time, POinters remain the same from then on
	{
		assert(&oCSyn[i] != NULL);
		oCSyn[i].setMetaplasticDistribution(&mpMDistrib); //Pass Histogram pointer
		vpSyns.push_back(&oCSyn[i]);
	}

	//Allocating Once And then Changing the Tracked Pattern on everytrial - No need to recreate Synapses - Just assume they have already been initialized by past activity
	//DistInit - Report - Only on 1st Trial
	reportStateDistribution<T>(vpSyns,iSynCount,slogFiles[0].c_str());

	uint8_t reportCycle = 0; //When Overflows we report the trial number
	///Do a Trial - First trial will allocate the memory buffer -
	t = trials;

	while (t > 0)
	{
		reportCycle ++; //Report Trial Progress
		if (reportCycle ==0) cout << t << endl;

		//INIT PATTERNS Read in The Vectors from the File to X - Less are available then uiPatCount Is reduced
		//Doing it This way or just random On each synapse in the synapse loop appears to have no impact on speed.
		if (!bUseRandomPatterns)
			::readTestVectorsFromFile(pinputFile,X,iSynCount,uiPatCount,0,bShuffleLoadedVectors); //Load At 0 All Required PAtterns

		uint uiNoOfPatternsStoredInTrial = 0;
		int  iPatIndex					 = 0;//The Pattern Index in focus at any particular timestep
		bool bPatternArrived 			 = false;
		bool bTimeToTestRecall 			 = false;
		///PATTERN Loop For alluiPatCount Patterns - Induce Stimuli to synapses and Track Fusi And Perceptron Signal
		//for (long j=0;j<lTotalTimesteps;j++) //Do each TimeSlot
		unsigned long j 				= 0;
		unsigned long NextEncodingj 	= 0;
		unsigned long LastRecallj 		= 0;
		unsigned long lFPT 				= lTotalTimesteps; //First Passage Time below 0 For this Trial

		dDALevel 						= 0.0; //Reset DA Level
		memset(ulLastDASignalTime,0,sizeof(unsigned long)*iNoOfDAInjections); //Reset Time of Last DA Signal
		int iDAStimuliCount				= 0; //Counter Index used to track the number of DA injections - Reset on every Trial

		//Loop Jumping between timesteps until A Memory Arrives after the end of the recording period -- With the addition of a 1st cycle:
		//Obtain a Metacycle distribution sample If required - Once Obtained Each trial ends with the Loop Trial Until Timestep Condition
		bool bAllocatePattern = false; //This Flag will be set On-by encodeMemory when h>0.5
		//First Wait enough cycles so the Sampled Distribution is obtained (If set )
		while ((j < lTotalTimesteps) || (mpMDistribinSamples[0] == 0 && g_MetaplasticitySampleSize > 0))
		{
			bool bDASignal = false;
			bool bThresCycleSampleObtained = (mpMDistribinSamples[0] != 0 || g_MetaplasticitySampleSize == -1); //Fixed Sample SIze Distribution Is enabled And has been obtained
			//Check if time to obtain sample of metaplastic Counters - But only after Fixed Sample Size has been Obtained (If Enabled g_MetaplasticitySampleSize > 0)
			if ((j == cSampleMetaplasticCounters && cSampleMetaplasticCounters > 0) && bThresCycleSampleObtained ) //Save after this step at j Completes and The Fixed Sample Distribution Has been Obtained
			{	//Save the Current Running Cycles Into Distribution And Copy Distribution To Aggregate mpMDistribinTime
				saveMetaplasticCounters(vpSyns,mpMDistribinTime,mpMDistrib);
				//cout << j << " Save Metaplastic Counters " << mpMDistrib[0] << endl;
			}

			//State Machine code This Code Overrides any previous Decision on which event should occur So to accommodate the Initialisation Period
			if (uiNoOfPatternsStoredInTrial <= (uiInitPatterns))
			{
				bPatternArrived = true;
				bTimeToTestRecall = false; //Havent Finished INit So Dont Waste Time Measuring Signal
				j = NextEncodingj = 0; //Stuck In time Until Init Patterns and the Tracked pattern have been delivered
				if ((uiInitPatterns) == uiNoOfPatternsStoredInTrial) //Once InitPeriod Is over - Start the 1st recall period
				{
					bTimeToTestRecall = true;
					LastRecallj = 0;
				}
			} //Rest Of Timestep Control Occurs At the end of the loop
			//else
				//bTimeToTestRecall =true;

			//Encode Pattern
			if (bPatternArrived)
			{
				//Select The pattern to Encode -
				bool bEncodeNewPattern = false; //SelectPatternTOEncode will set this true if this pattern should be encoded for the first time and saved in X
				//Handle DA Signal Timecourse
				dDALevel = 0.0;
				for (int k=0;k<iDAStimuliCount;k++)
				{ //Sum all contributions
					dDALevel  				+= dDAStep*(double)(j-ulLastDASignalTime[k])*exp(-g_fcAMPDecay*(double)(j-ulLastDASignalTime[k]));
					//cout << dDALevel;
				}
				if (bUseRandomPatterns){
					iPatIndex = selectPatternToEncode(j,uiNoOfPatternsStoredInTrial,bEncodeNewPattern,bDASignal,vTrackedIndex,repetitionTable,X,uiPatCount,iSynCount,mprng);
					// Handle DA Signalling //
					if (bDASignal){
						//if (bEncodeNewPattern) //Need to reflect injection at time of 1st encoding
						//	dDAInjected += dDAStep;
						//else
						//dDAInjected			= dDALevel + dDAStep; //INject more DA->cAMP
						dDAInjected			=  dDAStep; //INject more DA->cAMP
						ulLastDASignalTime[iDAStimuliCount] 	= j;
						iDAStimuliCount++;

					}
					assert(iDAStimuliCount < iNoOfDAInjections);
					//cout << dDALevel << endl;
				}
				else//Load the next Index from File
					iPatIndex = uiNoOfPatternsStoredInTrial; //Just move to the next pattern in the file
//#ifdef _DEBUG
				//cout << j << " " << iPatIndex << " Smpl:" << mpMDistrib[0] << " Sampls:" << mpMDistribinSamples[0] << endl;
//#endif
				//Run Through Synapses Encoding Induction Signal - if iPatIndex =-1 then induction stimuli are generated on the fly
				//Sets bAllocatePattern=true if h>0.5 - Alloc Is reset to false for each trial-Synapses plasticity is reset on every trials
				if (dDALevel < dDAThreshold) //Override the allocFlag if it has been switched on
					bAllocatePattern = false; //Switch off allocation when DA below threshold - Alloc Could have been switched on in encodeM function by h > thresh.

				iSignalSynsCount = encodeMemory<T>(vpSyns,j,vTrackedIndex,uiNoOfPatternsStoredInTrial,X,W,iTrackedSynFlag,iPatIndex,bEncodeNewPattern,bAllocatePattern,mprng,uiAllocSynapses, mpMDistrib, mpMDistribinSamples);

				if (t==trials) //Check If it is time to Report Encoding Event
				{
					if ((vTrackedIndex.find(uiNoOfPatternsStoredInTrial)) != vTrackedIndex.end()) //Decrement One Since If Encoded It would Have increased by one
					{//After Mem Storage Distribution
						cout << endl << "T:" << t << " Tracked Memory " << uiNoOfPatternsStoredInTrial << "* Stored in :" << iSignalSynsCount << endl;
						reportStateDistribution<T>(vpSyns,iSynCount,slogFiles[2].c_str()); //Dist B
					}
				}

				uiNoOfPatternsStoredInTrial++;
			}//Finished Looping through all synapses - Pattern is now stored

			//Check for Error COndition
			if (uiPatCount == uiNoOfPatternsStoredInTrial && !bUseRandomPatterns)
			{
				cerr << "Patterns Used: " << uiPatCount;
				ERREXIT(500,"Insufficient number of patterns loaded for simulation - Increase uiPatCount:");
			}

			//MEASURE SIGNAL Obtain Distributions and Measure the Signals for all timesteps -
			double fSignalTracked = 0.0; //Tracked population Signal
			double fSignalNTracked = 0.0; //Not Tracked Population Signal

			//Measure Signal on monitored Synapses if Init Period Is over and we have not exceeded the total array size for storing the signal
			if (bTimeToTestRecall && (j < lTotalTimesteps)) //3rd Period - (j < lTotalTimesteps) is required because array bounds may be exceeded when on the period of measuring ThresCycle Histograms
			{
				//cout << "R@"<< j << endl; //Report Recall Times
				LastRecallj = j; //Save Time of Last Recording - Measurements are taken on every tsPeriodOfRecall
				//Measure Perceptron Sig For Each Tracked Pattern
				int i=0;
				lFPT = lTotalTimesteps; //Reset The Point in time Where 1st Tracked pattern Goes 0
				for (itTracked = vTrackedIndex.begin();itTracked!=vTrackedIndex.end();++itTracked)
				{
					fSignalTracked = fSignalNTracked = 0.0;

					//iTrackedMemIndex = itTracked->first; //The second value Contains The index Relative to the storage of the X vector
					testPRecallOfX<T>(fSignalNTracked,fSignalTracked,oCSyn,W,X,iTrackedSynFlag[i],iSynCount,i);
					//Save signal At this memory storage Time Step
					double tSignal = fSignalTracked + fSignalNTracked;
					//if (j==0 && iTrackedMemIndex == 0)
					//	cout << "TSn:" << fSignalTracked << " NTSn:" << fSignalNTracked << endl;
					//cout << i << "-" << LastRecallj << ":" << tSignal << " Sz:"<< (int)vTrackedIndex.size() << endl;
					vTs.insert(LastRecallj); //Add to List Of Measurement Times
					dPtSn[i][LastRecallj] 	+= fSignalTracked;
					dPntSn[i][LastRecallj]	+= fSignalNTracked;
					dDASn[i][LastRecallj]	+= dDALevel;
					dVar[i][LastRecallj]	+= (tSignal)*(tSignal); //Store h^2 to obtain E[h^2]
					iAllocSynapses[i][LastRecallj]	+=uiAllocSynapses; //Count The number of Allocated Synapses
					//dCoVar[LastRecallj] += dSingleCoVar; //This E[X_i*X_j]

					i++; //Increment Tracked Index
				} //And Of Loop through Tracked Patterns

			}//If Time To test Recall

			//Go to Next Time step - The number of Patterns is incremented - The bool vars for Recall Or Encoding change depending on the type of the next event
			if (uiNoOfPatternsStoredInTrial >= uiInitPatterns)
				j = getNextTimestep(bTimeToTestRecall,bPatternArrived,ts,j, LastRecallj,NextEncodingj, PsMemEvent,repetitionTable);

			if (uiNoOfPatternsStoredInTrial > (uiPatCount-1) && !bUseRandomPatterns)
				ERREXIT(500,"We run out of patterns! Increase The capacity");

		}//For each timestep

/*
		if (!(t%100))
		{
			saveStateDistribution<T>(vpSyns,(long**)iOccupancy,iSynCount); //Augment Distribution
			iDistSamples++;
			//reportAvgStateDistribution<T>(vpSyns,iSynCount,slogFiles[3].c_str(),(int**)iOccupancy,iDistSamples); //Report Avg Up to now
		}
*/
		t--; //Decrement Trial count down to 0
	}//LOOP For each trial

	//Save Avg Distribution To File
	//reportAvgStateDistribution<T>(vpSyns,iSynCount,slogFiles[3].c_str(),(long**)iOccupancy,iDistSamples); //Report Avg Up to now

/* WRITE OUTPUT TO FILES */
	//OPEN OUTPUT FILES

	char buffObjName[250];
	T::getTypeName(buffObjName);
	//////LOG File Opened////


	double dCovar = 0.0;
	double dsqE = 0.0;
///RECORD STATISTICS FROM EACH TRACKED PATTERN
	char buffFilename[400];
	int i = 0; //TrackedMem INdex increment

	//MFPT
	for (itTracked = vTrackedIndex.begin();itTracked != vTrackedIndex.end();++itTracked)
	{
		//strcpy(buffFilename,slogFiles[4].c_str());
		sprintf(buffFilename,(const char*)slogFiles[4].c_str(), itTracked->first-uiInitPatterns);
		cout << "Signal Output Files: " <<  buffFilename << endl; //Tell User Which Output file we are using
		ofstream ofile(buffFilename, ios::out ); //Open Data File
		if (!ofile.is_open())
			ERREXIT(100,"Could Not Open output files. Check directories");
		//Write Header
		ofile << "#" << buffObjName << " Event Driven Memory Lifetime simulation Ts:" << ts << " Total samples: " << lTotalTimesteps << " Signal Sampling every :"<< "On everyEncoding" <<endl;
		ofile << "#t\tE[h]\tPerceptronSigNonEncodingSynapses\tPerceptronSigEncodingSynapses\tVariance\tCoVariance\tE[h^2]\tAllocFraction\tDALevel" << endl;

		double simTime = 0.0; //Start from ts so GnuPlot Can plot the 1st timepoint on a logscale that does not start from 0.
		double dEhsquared = 0.0;
		double dAllocSignal = 0.0;
		double dSignal = 0.0;
		//tsPeriodOfRecall = 1;

		//Iterate Through All time Points of Measurement
		for (set<unsigned long>::iterator it = vTs.begin();it!=vTs.end();++it)
		{
			//tsPeriodOfRecall = getCurrentPeriodOfRecall(j);
			unsigned long j = *it;
			dSignal 	= (dPntSn[i][j]+dPtSn[i][j])/trials;
			dPntSn[i][j] 	= (dPntSn[i][j]/trials); //NonTrack Perceptron - Avg
			dPtSn[i][j] 	= (dPtSn[i][j]/trials);
			dDASn[i][j]		= (dDASn[i][j]/trials);//Avg DA signal level
			dsqE 		= (dSignal)*(dSignal); //E[X]^2
			dEhsquared	= (dVar[i][j]/trials);
			dVar[i][j] 	= dEhsquared-dsqE; // Var = E[X^2]-E[X]^2
			dCovar 		= dVar[i][j]-(1+dsqE)/(iSynCount);
			dAllocSignal = ((double)iAllocSynapses[i][j]/trials)/iSynCount;
			//Calc Mean Sig From Each Synapse At this timepoint

			simTime=ts*(double)j;

			//if (i==0) 				cout << "S:" << (dSignal) << " N:" <<  sqrt(dVar[i][j]) << endl;
			//While Above Water Update the Lifetime of the 1st tracked Memory
			if (i==0 && (dSignal > sqrt(dVar[i][j])))
			{
				dLowSignalMFPT = simTime; //Save Until last point when Avg Memory Signal is above Noise
				//cout << dSignal/sqrt(dVar[i][j]) << endl;
			}

			ofile << (simTime) << "\t" << (float)(dSignal) << "\t" << (float)(dPntSn[i][j]) << "\t" << (float)(dPtSn[i][j]) << "\t" <<(float)(dVar[i][j]) << "\t" << dCovar << "\t" << dEhsquared << "\t" << dAllocSignal << "\t" << (float)dDASn[i][j] << endl;
			dLastSignalValue = dSignal;
			dLastVarValue = (dVar[i][j]);
		}

		i++; //Next Tracked Mem Index
		ofile << "#EOF" << endl;
		ofile.close();
	}
	cout << "*Mean Signal Lifetime : " <<  dLowSignalMFPT << endl;

/*
	//Now Save the Distribution Of Metaplastic Counters sampled at a single point
	cout << "Single point Meta Distribution Output File: " <<  slogFiles[3] << endl; //Tell User Which Output file we are using
	ofstream ofile(slogFiles[3].c_str(), ios::out ); //Open Data File
	ofile << "#SampleTime :" << cSampleMetaplasticCounters << endl;
	ofile << "#MetaplasticTransitions\tNumberOfSyns" << endl;
	for (map<uint,uint>::iterator it = mpMDistrib.begin();it!= mpMDistrib.end();++it)
	{
		ofile << it->first << "\t"<< it->second/(double)trials << endl;
	}
	ofile.close();
*/

	saveCycleHistogramToFile(mpMDistribinTime, slogFiles[9],trials,cSampleMetaplasticCounters);

	saveCycleHistogramToFile(mpMDistribinSamples, slogFiles[10],trials,mpMDistribinSamples[0]);

	//CLEAR MEMORY
	vpSyns.clear();
	for (int i=0;i<ciNoOfTrackedPatterns;i++)
	{

		delete [] dVar[i]; //Memory For Variance Vector for this Particular tracked Pattern
		delete [] dPtSn[i]; //Tracked Signal
		delete [] dPntSn[i]; //Non Tracked Signal
		delete [] dDASn[i];
		delete [] iAllocSynapses[i];
		delete [] iTrackedSynFlag[i];
	}

	for (uint i=0;i<uiPatCount;i++) //Remove All saved Patterns
		delete [] X[i];

	delete [] W;
	delete [] dVar;
	delete [] dPtSn;
	delete [] dPntSn;
	delete [] dDASn;
	delete [] iTrackedSynFlag;
	delete [] iAllocSynapses;
	//delete [] iOccupancy[0];
	//delete [] iOccupancy[1];
	delete PsMemEvent;

	//CleaUp CUDA
#ifdef USE_CUDA
	cleapUpCUDADeviceMem(d_W,d_X,d_C,d_odata,h_odata,iSynCount);
	delete [] tX; //The Long Buffer Vector Of Tracked PAtterns
#endif
	cout << "-Fin-" << endl;

	//return dLowSignalMFPT;
	pair<double,double> Ret;
	Ret.first = dLastSignalValue;
	Ret.second = dLastVarValue;

	return Ret;
}





//#undef USE_CUDA

#endif /* CONTINUOUSTIMEEXPERIMENTS_H_ */
