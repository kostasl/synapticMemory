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
extern uint g_timeToSampleMetaplasticity;
extern uint g_MetaplasticitySampleSize;
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

unsigned long getNextTimestep(bool& bisRecallPeriod,bool& bisEncodingPeriod,double dTimeStepSize,unsigned long currentTimestep,unsigned long& LastRecallj,unsigned long& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable);
int selectPatternToEncode(unsigned long j,uint uiNoOfPatternsStoredInTrial,bool bUseRandomPatterns,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng);
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
uint encodeMemory(vector<T*>& vpSyns,unsigned long ts,t_patt_trackedtbl& vTrackedIndex,uint uiNoOfPatternsStoredInTrial,t_inVal** X,t_inVal* W,int** iTrackedSynFlag,int iPatIndex,gsl_rng*& prng_r,uint& uiAllocSynapses,map<uint,uint>& pMDistribinTime, map<uint,uint>& pMDistribinSamples)
{
	double dC = 0.0;
	double h = 0.0; //Neuron Depolarization
	long cMetaplasticCountersSampleSize = g_MetaplasticitySampleSize; //Sample the metaplastic Counters At Fixed Time
	float mfRate = 1.0f;
	uint iSignalSynsCount = 0; //Reset the Memory Synapses Count
	bool bAllocationSignal = false;
	bool bresetAllocationThreshold = false;
	t_patt_trackedtbl::iterator itAtTrackedPatt;
	int TrackedIndex = -1;
	t_patt_trackedtbl::iterator itAtEnd = vTrackedIndex.end();
	t_patt_trackedtbl::iterator itAtStart = vTrackedIndex.begin();
	typename vector<T*>::iterator it; //Use typename Because T could be referring to a type or a member of class T
	///1st PERIOD - Before Storage - Start to monitor all synapse Changes

/*
	If a pattern is repeated The it has been saved in the X vector Before Because it would have been a tracked pattern too
	In that case the iPatIndex will not be -1 but there is no need to Set TrackedIndex as the pattern has been stored before.
	Search if this is a tracked pattern and obtain the index in the tracked list
	//itAtTrackedPatt = vTrackedIndex.find(uiNoOfPatternsStoredInTrial);// find(vTrackedIndex.begin(),vTrackedIndex.end(),);
 */
	uint uiPatternToCheck; //Which Pattern Index Should be Checked if it has been Allocated

	if (iPatIndex ==-1) //Not A repetition of A previous pattern? Check if current pattern is in the tracked list
	{
		uiPatternToCheck = uiNoOfPatternsStoredInTrial; //Search TrackedList And Get Index to store Vector in X
	}
	else
		uiPatternToCheck = iPatIndex; //This is a Repetition of Tracked pattern So check the Given Index To see if the Pattern is Has an Allocation Signal

	////////FIND IF THIS PATTERN IS TRACKED And Allocated- Obtain Tracked list Index and Allocation Flag////////////
	int temp = 0;
	for (itAtTrackedPatt = itAtStart; itAtTrackedPatt!= itAtEnd;++itAtTrackedPatt) //Only Track 1st pattern
	{
		if (itAtTrackedPatt->first == uiPatternToCheck) //Pattern being Encoded In the tracked list?
		{//Check If it should be allocated
			TrackedIndex = temp; //TrackedIndex tells Patt generating code below to store this pattern
			//State Machine here to change behaviour between repetitions
			if (itAtTrackedPatt->second == 2) //If Flag Set and this is the repetition and not the 1st Encoding -Allocate
			{
				itAtTrackedPatt->second = 1; //Return To 1st state
			}
			///1st Repetition : itAtTrackedPatt->second == 1 <-Allocation Signal For this memory is On
			///iPatIndex == -1 1st Encoding of this memory
			if (itAtTrackedPatt->second == 1) //If Flag Set and this is the repetition and not the 1st Encoding -Allocate
			{   //Reset the Refraction counter
				bAllocationSignal = true;
				//cout << "Reset Refraction Counter" << endl;
				if (iPatIndex == -1) //Ist Encoding Of this Pattern
					bresetAllocationThreshold = true; //Reset Counters On 2nd Repetition
				//itAtTrackedPatt->second = 2;
			}
			//cout << iPatIndex << " Trck:" << itAtTrackedPatt->first <<" Found At X[" << TrackedIndex << endl;
			break; //Set To the 1st Tracked Pattern
		}
		temp++; //Used As Relative Index In the X[] Pattern Vector
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

	if (ts == 0 && uiNoOfPatternsStoredInTrial==0){ //Start Of New Trial? Unfreeze Synapses
		//Unlock Plasticity At the beginning - But should we start from Random Point?  No Need cause Plasticity lock does not stop internal filter transitions only strength state
		oCSyn->reset(); //Begin new Trial -> Unlock Plasticity, Reset Counters
	}

	if (bresetAllocationThreshold)
		oCSyn->resetMetaplasticCounter();

	if (bAllocationSignal)
	{//At the allocation signal we begin measuring Stability
		oCSyn->enableMetaplasticAllocation();
		//oCSyn->resetAllocationRefraction(); //Reset The internal Counter
	}

	if (!oCSyn->isPlastic()) //Count Locked Synapses
		uiAllocCounter++;

	//Whether File or Random Patterns
	int iStim;
	//Found A repeated Pattern - Use TrackedIndex TO obtain saved pattern in X
	if (iPatIndex != -1) //If NOT random Pattern On the fly The Read the stimulus from the X vector
	{
		//For Random vectors This is only used when selectPatternToencode Finds the repetition of a memory
		iStim = X[TrackedIndex][i]; //*X[iSynCount-1]; //Assume Fixed Neuron Output of +1
	}
	else
	{ //iPatIndex == -1 So Make pattern on the fly and store it
	  //Create The Pattern One Stimulus At the time
		iStim = makeInductionStep(mfRate, prng_r, dC);
		if (TrackedIndex != -1) //But this is a tracked Pattern So Save it
			X[TrackedIndex][i] = iStim; //Save The stimulus in the X vector holding The Tracked Patterns
	}

	switch (iStim)
	{
		case 1:
			if (TrackedIndex >= 0) oCSyn->setTrackedState(ICascadeSynapse::SYN_STRENGTH_STRONG); //Set The "desired" strength State
			oCSyn->handlePOT();
		break;
		case -1:
			if (TrackedIndex >= 0) oCSyn->setTrackedState(ICascadeSynapse::SYN_STRENGTH_WEAK);
			oCSyn->handleDEP();
		break;
		case 0:
			if (TrackedIndex >= 0) oCSyn->setTrackedState(ICascadeSynapse::SYN_STRENGTH_NOTSET);
			oCSyn->handleNOP(); //For "Sparseness"
		break;
		default:
			ERREXIT(500,"encodeMemory: Unknown Induction Stimulus")
		break;
	}

	W[i] = oCSyn->getStrength(); //Save Strength Into Weight Vector
	h	+= W[i]*iStim; //Get Neural Output Due to Pattern Being Stored

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


	//Check if Limit on the number of MetaCycles has been obtained And Save a copy of the distribution
	if (oCSyn->getMetaplasticDistributionSampleSize() == cMetaplasticCountersSampleSize)
	 pMDistribinSamples.insert(pMDistribinTime.begin(),pMDistribinTime.end());//Copy Snapshot of current Distrib IN time

	i++; //INcrement Index used for input vector
} //For Each Synapse

h = h/vpSyns.size();

//Suppose Avg Signal is 0.1 -- This could be changed to a BCM Like rule.
//if (h > 0.1) //Set STATIC Variable
//	bAllocationSignal = true;
//else
//	bAllocationSignal = false;
//cout << "**N:" << iSignalSynsCount << endl;
uiAllocSynapses = uiAllocCounter;
//cout << "t:" << ts << " NAlloc:" << uiAllocCounter << endl;
return iSignalSynsCount;
}


//Function that Incremends the number of synapses found with a specific metaplastic transition count -
//Constructs distribution Vector
template<class T>
void saveMetaplasticCounters(vector<T*>& vpSyns,map<uint,uint>& mpMDistrib)
{

	for (typename vector<T*>::iterator it = vpSyns.begin(); it != vpSyns.end(); ++it)
	{
		//uint cnt = (*it)->getMaxMetaplasticCount();
		uint cnt = (*it)->getMetaplasticCount();

		//If Key Does not Exist Then Reset This Distribution Bin to Zero
		if (mpMDistrib.find(cnt) == mpMDistrib.end())
			mpMDistrib[cnt] = 0;

		mpMDistrib[cnt]++; //Increment The number of synapses Found in this Counter State
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
 */
template <class T,class pFunct>
double simRepetitionAllocation(pFunct pF, uint iSynCount,int iCascadeSize,uint uiInitPatterns, char* pinputFile,uint trials,uint simTimeSeconds,double dEncodingRate,t_patt_reptbl& repetitionTable,double ts,vector<string>& slogFiles)
{
#ifdef USE_CUDA
	cout << "Using CUDA For Signal Dot Prod" << endl;
#else
	cout << "CPU For Signal" << endl;
#endif
	cout << "C.T Signal Lifetime With repetition of tracked Signal -Net Size : " << iSynCount << " Trials:" << trials <<endl;
	//const double ts				= 1.0;//0.001;//dEncodingRate/2.0; //TimeStep is sampling at twice the frequency
	gsl_rng* mprng 				= g_getRandGeneratorInstance(true);
	const long tsPerSec			= (1.0/ts);
	//const unsigned long peakSig  =g_FilterTh*g_FilterTh*tsPerSec*0.375;
	const unsigned long lTotalTimesteps	= (double)simTimeSeconds/ts; //The Total Number of TimeSteps in the simulation- (Number of Samples )
	uint iSignalSynsCount 		= 0; //Signal Variable
	uint uiAllocSynapses		= 0; //Stores the count of Allocated Synapses returned by encodeMemory
	uint iDistSamples			= 0;
	uint tsPeriodOfRecall 		= 1; //(double)lTotalSamples*dEncodingRate/(simTimeSeconds*10); //The period in Timesteps of recall measurements
	uint iTrackedMemIndex		= 0;//Holds the current Tracked index when Iterating through the table
	bool bShuffleLoadedVectors 	= true;
	bool bUseRandomPatterns		= (pinputFile == 0) || (pinputFile[0]=='\n'); //If no file is given Use Random Patterns
	uint maxMemCapacity			= dEncodingRate*simTimeSeconds*2; //Max Number of Patterns to create - x2 for safety
	uint uiPatCount 			= uiInitPatterns + maxMemCapacity; //Max Patterns loaded from File-Valid Only File- And Thus Max Storage Capacity is NeuronCount/5

	set<unsigned long> vTs; //Contains All the time point were Strength Has been Measured#

	//MEAN FIRST PASSAGE TIME VARIABLES
	double dLowSignalMFPT 		= -1.0; //The point Where the first pattern tracked drops below noise threshold
	char *mem_buffer 			= 0;		//This Pointer is filled by allocMem, To point to the reuseable allocated memory
	T* oCSyn;	//POinter to allocated memory
	vector<T*> vpSyns; //Vector of pointers to Synapses
	vpSyns.reserve(iSynCount);
	PoissonSource* PsMemEvent = new PoissonSource(dEncodingRate,ts,0);

	//ADD THE LIST OF TRACKED PATTERNS
	t_patt_trackedtbl vTrackedIndex; //Key:The PattIndex - Value: 1 For Allocation Signal / 0 For no Allocation
	t_patt_trackedtbl::iterator itTracked;
	t_patt_reptbl::iterator itRep;

	map<uint,uint> mpMDistrib; //The distribution Of Metaplastic Transitions Among Synapses At a particular point in time
	map<uint,uint> mpMDistribinSamples; //The Avg distribution Of Metaplastic Transitions among Synapses for a fixed number of cycle samples
	map<uint,uint> mpMDistribinTime; //The Avg distribution Of Metaplastic Transitions among Synapses sampled over a fixed amount of time


	mpMDistribinTime.clear();
	mpMDistrib.clear();

	mpMDistribinSamples.clear();

	for (int i=0;i<=0;i++) //Start tracking from 1 cause 0th falls within the init patterns
		vTrackedIndex[i+uiInitPatterns] = 0;  ; //Setting to 0 means Disable Allocation / 1 Means Allocate

	//ADD REPEATED PATTERNS TO TRACKED LIST AND MARK AS ALLOCATED - Track All patterns That Are in the Repetition table
	unsigned long maxRepTime = 0;

	for (itRep = repetitionTable.begin(); itRep != repetitionTable.end();++itRep  )
	{
		///Tracked pattern Is allocated - Give Allocation Signal
		vTrackedIndex[(*itRep).second] = 1;  //Add Repeated Mem Index to the List of tracked Patterns

		if (itRep->first > maxRepTime)
			maxRepTime = itRep->first; //Save the Last Repetition Time into MaxRep
	}
	long cSampleMetaplasticCounters 	= g_timeToSampleMetaplasticity; //Sample the metaplastic Counters At Fixed Time


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
	int**		 iTrackedSynFlag = new int*[ciNoOfTrackedPatterns]; //Holds An array of flags indicating if a synapse is part of the DElta or the Omega group in a tracked pattern
	uint** 		 iAllocSynapses	 = new uint*[ciNoOfTrackedPatterns];

    cout << "Number of Tracked patterns :" << ciNoOfTrackedPatterns << endl;
	//Init Time Slots For each Tracked Pattern
	for (int i=0;i<ciNoOfTrackedPatterns;i++)
	{
		dVar[i] 	= new double[lTotalTimesteps]; //Memory For Variance Vector for this Particular tracked Pattern
		dPtSn[i]	= new double[lTotalTimesteps]; //Tracked Signal
		dPntSn[i]	= new double[lTotalTimesteps]; //Non Tracked Signal
		iAllocSynapses[i]	= new uint[lTotalTimesteps]; //The Number of allocated synapses at a particular time step
		iTrackedSynFlag[i] 	= new int[iSynCount]; //Array the size of Number of Synapses each element holding 1 for Mdelta or -1 for Omega

		memset(dPtSn[i],0,lTotalTimesteps*sizeof(double));
		memset(dPntSn[i],0,lTotalTimesteps*sizeof(double));
		memset(dVar[i],0,lTotalTimesteps*sizeof(double));
		memset(iTrackedSynFlag[i],0,iSynCount*sizeof(int));
	}

	long* 	 iOccupancy[2]; //Store counters for each Strength-State Pair - Distribution
	uint t = 0; //Trial counter

	cout << "Event driven Continuous Time Simulation measuring signal after "<< uiInitPatterns << "memories using Perceptron measures " << endl;
	memset(X,0,uiPatCount*sizeof(t_inVal*)); //Setting The pointers to 0 so InitPatterns Knows to initiaze each member to a new array
	memset(W,0,iSynCount*sizeof(t_inVal)); //Setting The pointers to 0 so InitPatterns Knows to initiaze each member to a new array

	iOccupancy[0] = new long[15];
	iOccupancy[1] = new long[15];
	memset(iOccupancy[0],0,sizeof(long)*15); //Empty the memory Buffer
	memset(iOccupancy[1],0,sizeof(long)*15); //Empty the memory Buffer

	//Initialise The memory For Patterns
	initPatternMemory(X,uiPatCount,iSynCount,iTrackedMemIndex,0.5, mprng,bUseRandomPatterns);

	//ALLOCATE SYNAPSES Call Object Allocation Function - Creates Object in memory and returns pointer to 1st object
	oCSyn = (T*)(*pF)((char*)mem_buffer,iSynCount,(int)iCascadeSize,mprng,1.0);
	if (!oCSyn) {//Check Failure
			ERREXIT(500,"simMemRepetitionAllocation: Could not create synapse objects! Out Of memory?");
		}else //Save pointer
			mem_buffer = (char*)oCSyn;
	///Save pointer to a Vector only on 1st trial- Temporary Solution But Its Required
	//The pointers do not change between re-allocations
	vpSyns.clear();
	for (uint i=0;((i<iSynCount)) ;i++) //Could Do it only 1st time, POinters remain the same from then on
	{
		assert(&oCSyn[i] != NULL);
		oCSyn[i].setMetaplasticDistribution(&mpMDistribinTime); //Pass Histogram pointer
		vpSyns.push_back(&oCSyn[i]);
	}
	//Allocating Once And then Changing the Tracked Pattern on everytrial - No need to recreate Synapses - Just assume they have already been initialized by past activity
	//DistInit - Report - Only on 1st Trial
	reportStateDistribution<T>(vpSyns,iSynCount,slogFiles[0].c_str());

	uint8_t reportCycle = 0; //When Overflows we report the trial number
	///Do a Trial - First trial will allocate the memory buffer -
	t = trials;


#ifdef USE_CUDA
	///ALLOCATE CUDA VARIABLES
		 void* d_W,*d_X,*d_C,*d_odata; //Pointers To device Memory
		 float* h_odata; ////Buffer used To output Data from the reduction
		 t_inVal* 	tX = new t_inVal[ciNoOfTrackedPatterns*iSynCount]; //The Joined Array of all Tracked Patterns To pass to CUDA
		 memset(tX,0,ciNoOfTrackedPatterns*iSynCount*sizeof(t_inVal));
		 initCUDADeviceMem(d_W,d_X,d_C,d_odata,h_odata,iSynCount,ciNoOfTrackedPatterns);
#endif
	while (t > 0)
	{
		reportCycle ++;
		if (reportCycle ==0)
			cout << t << endl;

		//INIT PATTERNS Read in The Vectors from the File to X - Less are available then uiPatCount Is reduced
		//Doing it This way or just random On each synapse in the synapse loop appears to have no impact on speed.
		if (!bUseRandomPatterns)
			readTestVectorsFromFile(pinputFile,X,iSynCount,uiPatCount,0,bShuffleLoadedVectors); //Load At 0 All Required PAtterns

		uint uiNoOfPatternsStoredInTrial = 0;
		int  iPatIndex					 = 0;//The Pattern Index in focus at any particular timestep
		bool bPatternArrived 			 = false;
		bool bTimeToTestRecall 			 = false;
		///PATTERN Loop For alluiPatCount Patterns - Induce Stimuli to synapses and Track Fusi And Perceptron Signal
		//for (long j=0;j<lTotalTimesteps;j++) //Do each TimeSlot
		unsigned long j = 0;
		unsigned long NextEncodingj = 0;
		unsigned long LastRecallj = 0;
		unsigned long lFPT = lTotalTimesteps; //First Passage Time below 0 For this Trial

		//Loop Jumping between timesteps until A Memory Arrives after the end of the recording period
		while(j < lTotalTimesteps)
		//while (mpMDistribinSamples[0] < g_MetaplasticitySampleSize) //Check If There are enough distribution cycle-samples to complete this trial
		{
			//CHECK IF TIME TO OBTAIN Instanteneous SAMPLE of METAPLASTIC COUNTERS -
			if (j == cSampleMetaplasticCounters) //Save after this step at j Completes
			{
				//This will also trigger the EndCycle on next threshold Event
				//If the simulation time runs out, the cycles that did not reach threshold are not counted
				saveMetaplasticCounters(vpSyns,mpMDistrib);
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
				//Select The pattern to Encode - Creates the next random Pattern if required Or Chooses a pattern to be repeated
				iPatIndex = selectPatternToEncode(j,uiNoOfPatternsStoredInTrial,bUseRandomPatterns,vTrackedIndex,repetitionTable,X,uiPatCount,iSynCount,mprng);

				//Run Through Synapses Encoding Induction Signal - if iPatIndex =-1 then induction stimuli are generated on the fly
				iSignalSynsCount = encodeMemory<T>(vpSyns,j,vTrackedIndex,uiNoOfPatternsStoredInTrial,X,W,iTrackedSynFlag,iPatIndex,mprng,uiAllocSynapses, mpMDistribinTime, mpMDistribinSamples);

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
			if (bTimeToTestRecall && (j < lTotalTimesteps)) //3rd Period
			{
				LastRecallj = j; //Save Time of Last Recording - Measurements are taken on every tsPeriodOfRecall
				//Measure Perceptron Sig For Each Tracked Pattern
#ifdef USE_CUDA
				//Check Again the Vectors Sent To CUDA Because the Vectors in X now contain only the tracked patterns
				int iTrackedCount = testCUDAPRecallOfX(h_odata,W ,X,tX, iSynCount,d_W,d_X,d_C,vTrackedIndex, uiNoOfPatternsStoredInTrial);
#endif
				int i=0;
				lFPT = lTotalTimesteps; //Reset The Point in time Where 1st Tracked pattern Goes 0
				for (itTracked = vTrackedIndex.begin();itTracked!=vTrackedIndex.end();++itTracked)
				{
					fSignalTracked = fSignalNTracked = 0.0;

#ifndef USE_CUDA ///Use i instead of TrackedIndex To preserve the order Tracked Patterns Were stored
					//iTrackedMemIndex = itTracked->first; //The second value Contains The index Relative to the storage of the X vector
					testPRecallOfX<T>(fSignalNTracked,fSignalTracked,oCSyn,W,X,iTrackedSynFlag[i],iSynCount,i);
#else //NOCUDA
					if (i < iTrackedCount) //Do not exceed the available tracked patterns
						fSignalTracked = h_odata[i]/iSynCount; //Save Each Result from the CUDA Call
#endif
					//Save signal At this memory storage Time Step
					double tSignal = fSignalTracked + fSignalNTracked;
					//if (j==0 && iTrackedMemIndex == 0)
					//	cout << "TSn:" << fSignalTracked << " NTSn:" << fSignalNTracked << endl;
					//cout << i << "-" << LastRecallj << ":" << tSignal << " Sz:"<< (int)vTrackedIndex.size() << endl;
					vTs.insert(LastRecallj); //Add to List Of Measurement Times
					dPtSn[i][LastRecallj] 	+= fSignalTracked;
					dPntSn[i][LastRecallj]	+= fSignalNTracked;
					dVar[i][LastRecallj]	+= (tSignal)*(tSignal); //Store h^2 to obtain E[h^2]
					iAllocSynapses[i][LastRecallj]	+=uiAllocSynapses;
					//dCoVar[LastRecallj] += dSingleCoVar; //This E[X_i*X_j]
					#ifdef VERBOSE
						cout << LastRecallj << "," << iSn << "," << fSignal << endl;
					#endif

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
	//string strDir(MFPTIMESSNR_OUTPUT_DIRECTORY);
	//strDir.append(buffObjName);

	char buffFilename[400];
	//////LOG File Opened////


	double dCovar,dsqE;
///RECORD STATISTICS FROM EACH TRACKED PATTERN

	int i = 0; //TrackedMem INdex increment

	//MFPT
	for (itTracked = vTrackedIndex.begin();itTracked != vTrackedIndex.end();++itTracked)
	{
		//strcpy(buffFilename,slogFiles[4].c_str());
		sprintf(buffFilename,(const char*)slogFiles[4].c_str(), itTracked->first-uiInitPatterns);
		cout << "Signal Output Files: " <<  buffFilename << endl; //Tell User Which Output file we are using

		ofstream* ofile = openfile(buffFilename,ios::out);

		if (!ofile->is_open())
			ERREXIT(100,"Could Not Open output files. Check directories");
		//Write Header
		*ofile << "#" << buffObjName << " Event Driven Memory Lifetime simulation Ts:" << ts << " Total samples: " << lTotalTimesteps << " Signal Sampling every :"<< "On everyEncoding" <<endl;
		*ofile << "#t\tE[h]\tPerceptronSigNonEncodingSynapses\tPerceptronSigEncodingSynapses\tVariance\tCoVariance\tE[h^2]\tAllocFraction" << endl;

		double simTime = 0.0; //Start from ts so GnuPlot Can plot the 1st timepoint on a logscale that does not start from 0.
		double dEhsquared = 0.0;
		//tsPeriodOfRecall = 1;

		//Iterate Through All time Points of Measurement
		for (set<unsigned long>::iterator it = vTs.begin();it!=vTs.end();++it)
		{
			//tsPeriodOfRecall = getCurrentPeriodOfRecall(j);
			unsigned long j = *it;
			double dSignal 	= (dPntSn[i][j]+dPtSn[i][j])/trials;
			dPntSn[i][j] 	= (dPntSn[i][j]/trials); //NonTrack Perceptron - Avg
			dPtSn[i][j] 	= (dPtSn[i][j]/trials);
			dsqE 		= (dSignal)*(dSignal); //E[X]^2
			dEhsquared	= (dVar[i][j]/trials);
			dVar[i][j] 	= dEhsquared-dsqE; // Var = E[X^2]-E[X]^2
			dCovar 		= dVar[i][j]-(1+dsqE)/(iSynCount);
			double dAllocSignal = ((double)iAllocSynapses[i][j]/trials)/iSynCount;
			//Calc Mean Sig From Each Synapse At this timepoint

			simTime=ts*(double)j;

			//if (i==0) 				cout << "S:" << (dSignal) << " N:" <<  sqrt(dVar[i][j]) << endl;
			//While mean Above Water Update the Lifetime of the 1st tracked Memory -
			//Beware that in Low Averaging Fluctuations will mess up this calculation
			if (i==0 && (dSignal > sqrt(dVar[i][j]) ))
			{
				dLowSignalMFPT = simTime; //Save Until last point when Avg Memory Signal is above Noise
				//cout << dSignal/sqrt(dVar[i][j]) << endl;
			}

			*ofile << (simTime) << "\t" << (float)(dSignal) << "\t" << (float)(dPntSn[i][j]) << "\t" << (float)(dPtSn[i][j]) << "\t" <<(float)(dVar[i][j]) << "\t" << dCovar << "\t" << dEhsquared << "\t" << dAllocSignal << endl;
		}

		i++; //Next Tracked Mem Index
		*ofile << "#EOF" << endl;
		ofile->close();
	}
	cout << "*Mean Signal Lifetime : " <<  dLowSignalMFPT << endl;

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


	//Now Save the Distribution Of Metaplastic Counters Sampled Over a Time interval
	cout << "In Time Meta Distribution Output File: " <<  slogFiles[9] << endl; //Tell User Which Output file we are using
	ofstream ofile2(slogFiles[9].c_str(), ios::out ); //Open Data File For Wrong Cycles
	ostringstream oss;
	oss << "#SampleDuration :" << cSampleMetaplasticCounters << " Trials:" << trials << std::endl;
	oss	<< "#Cycle-Size\tOverallFrequency\tCorrectStateFq\tWrongStateFq\tTotalSamples" << std::endl;
	ofile2 << oss.str();

	//Split Into Two distributions to ease saving into one File
	unsigned long lDistribSum = 0;
	int maxOccupancy =0; //Find the Maximum State Occupied
	int iOccupancyIndex;
	map<uint,uint> mDistribCorrect;
	map<uint,uint> mDistribWrong;
	//Set-Up The Distributions before Printing
	lDistribSum = mpMDistribinTime[0]; //The Total Number of Samples Is saved in 0
	for (map<uint,uint>::iterator it = mpMDistribinTime.begin();it!= mpMDistribinTime.end();++it)
	{
		if (it->first == 0)
			continue; //Skip the 0 Occupancy. It is used to hold the sum of samples used in the distribution

		//lDistribSum += it->second; //Sum the Total Number of samples
		if (it->first > 999) //Is this the Fq of correct one or Wrong one?
		{
			iOccupancyIndex = it->first-1000;
			mDistribCorrect[iOccupancyIndex] = it->second;
		}
		else
		{
			iOccupancyIndex = it->first;
			mDistribWrong[iOccupancyIndex] = it->second;
		}

		if  (maxOccupancy < iOccupancyIndex) maxOccupancy = iOccupancyIndex;
		//Replace With Sum Of Both
		mpMDistribinTime[iOccupancyIndex] = mDistribWrong[iOccupancyIndex] + mDistribCorrect[iOccupancyIndex]; //Fix To total Value
	}
	//Output To File
	for (int i=0;i<=maxOccupancy;i++)
	{
		ofile2 << i << "\t"<< mpMDistribinTime[i] << "\t" << mDistribCorrect[i] << "\t"<<  mDistribWrong[i] << "\t" << lDistribSum << endl;
	}
	ofile2.close();

	//if (mpMDistribinSamples.size() == 0)
	//	cerr << "MetaCycle Distribution in Time finished Before Distribution in Sample Size could complete! Increase the Sample Time" << endl;
	mDistribCorrect.clear();
	mDistribWrong.clear();
	lDistribSum	=	0;
	for (map<uint,uint>::iterator it = mpMDistribinSamples.begin();it!= mpMDistribinSamples.end();++it)
	{
		if (it->first == 0)
			continue; //Skip the 0 Occupancy. It is used to hold the sum of samples used in the distribution

		lDistribSum += it->second; //Sum the Total Number of samples
		if (it->first > 999) //Is this the Fq of correct one or Wrong one?
		{
			iOccupancyIndex = it->first-1000;
			mDistribCorrect[iOccupancyIndex] = it->second;
		}
		else
		{
			iOccupancyIndex = it->first;
			mDistribWrong[iOccupancyIndex] = it->second;
		}

		if  (maxOccupancy < iOccupancyIndex) maxOccupancy = iOccupancyIndex;
		//Replace With Sum Of Both
		mpMDistribinSamples[iOccupancyIndex] = mDistribWrong[iOccupancyIndex] + mDistribCorrect[iOccupancyIndex]; //Fix To total Value
	}
		//Now Save the Distribution Of Metaplastic Counters Over a fixed Sampled Size
		cout << " Meta Distribution Output File: " <<  slogFiles[10] << endl; //Tell User Which Output file we are using
		ofstream ofile3(slogFiles[10].c_str(), ios::out ); //Open Data File For Wrong Cycles
		ostringstream oss2;
		oss2 << "#SampleSize :" << lDistribSum << " Trials:" << trials << std::endl;
		oss2	<< "#Cycle-Size\tOverallFrequency\tCorrectStateFq\tWrongStateFq\tSamplesPerTrial" << std::endl;
		ofile3 << oss2.str();

		//Output To File
		for (int i=0;i<=maxOccupancy;i++)
		{
			ofile3 << i << "\t"<< (float)mpMDistribinSamples[i]/trials << "\t" << (float)mDistribCorrect[i]/trials << "\t"<<  (float)mDistribWrong[i]/trials<< "\t" << (float)lDistribSum/trials << endl;
		}
		ofile3.close();

	//CLEAR MEMORY
	vpSyns.clear();
	return_temporary_buffer(mem_buffer);
	for (int i=0;i<ciNoOfTrackedPatterns;i++)
	{

		delete [] dVar[i]; //Memory For Variance Vector for this Particular tracked Pattern
		delete [] dPtSn[i]; //Tracked Signal
		delete [] dPntSn[i]; //Non Tracked Signal
		delete [] iAllocSynapses[i];
	}

	for (uint i=0;i<uiPatCount;i++) //Remove All saved Patterns
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

	return dLowSignalMFPT;
}



//#undef USE_CUDA

#endif /* CONTINUOUSTIMEEXPERIMENTS_H_ */
