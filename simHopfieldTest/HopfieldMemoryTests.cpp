//@Author Kostas Lagogiannis 17/12/10
//Test Code to conduct experiments on Memory Lifetimes Using allocation in Autoassociative neural networks

#include "common.h"
#include "util.h"

#include "synapseAllocators.h"
#include "InputVectorHandling.h"
#include "HopfieldMemoryTests.h"

#include <iostream>



//Global Vars
uint g_time = 0;
bool bRecallInProgress = false;
extern float g_fAllocHThres;
extern int g_FilterTh; //Used for Single Filter Experiments
extern double g_FilterDecay; //0.0916986;
extern float g_fAllocHThres; //Default Post Synaptic depol. Threshold to switch on Allocation
extern float g_fcAMPDecay; //The timeconstant for the cAmp alpha process (With 0.5 it takes approx 10 tsteps for a complete wave)
extern float g_fcAMPMagnitude;
extern double g_dcAMPMax; // A globally set  saturation value of cAMP.
extern float g_fPKAAllocThres;
extern float g_fAllowedRecallError;

extern gsl_rng* g_rng_r;

//Simply Leanrning Rule - No Synapse Model Class
//Does A dotProd to produce the correlation matrix for a number of training Patterns patCount
float** makeWeightMatrix(int NetSize, float** Xin, int patCount) {
	float** w = new float*[NetSize];
	float norm = 1;

	for (int i = 0; i < NetSize; i++) {
		w[i] = new float[NetSize];

		for (int k = 0; k < NetSize; k++) {
			w[i][k] = 0.0;
			for (int j = 0; j < patCount; j++) //For Each Pattern
				w[i][k] += Xin[j][i] * Xin[j][k];

			if (i == 0 && k == 0)
				norm = w[0][0]; //Normalize
			w[i][k] = 100 * w[i][k] / norm; //TrappenBerg
		}
	}

	return w;
}

///Uses externally Managed Array store New pattern Ontop Of the previous Weights
//Superposition Classic Hopfield Learning
float** makeWeightMatrix(int NetSize, int* Xin, float** W) {

//	float norm = 1;

	for (int i = 0; i < NetSize; i++) {

		for (int k = i; k < NetSize; k++) {
			//if (i !=k)
			W[k][i] = W[i][k] += Xin[i] * Xin[k]; //Symmetric
			//W[i][k] += Xin[i]*Xin[k];
			//else
//				W[i][k] = W[k][i] = 0;

			//if (i==0 && k==0) norm = W[0][0];//Normalize ???What is This??
//				W[k][i] = W[i][k] = NetSize*W[i][k]/norm; //TrappenBerg

			if (i == k)
				W[i][k] = 0;
		}
	}
	return W;
}

template<class T>
inline void deleteMemoryBuffer(uint ArrSize, T**& buffer) {
	assert(buffer != 0);

	for (uint i = 0; i < ArrSize; i++) {
		return_temporary_buffer(buffer[i]);
		//delete[] buffer[i];
		//delete the memory for each char array in array first
	}
	//then delete the memory to the array of pointers
	delete[] buffer;
	buffer = 0;

}

//For Weight Matrix
template<>
inline void deleteMemoryBuffer<float>(uint ArrSize, float**& buffer) {
	assert(buffer != 0);

	for (uint i = 0; i < ArrSize; i++) {
		delete[] buffer[i];
		//delete the memory for each char array in array first
	}
	//then delete the memory to the array of pointers
	delete[] buffer;
	buffer = 0;

}
//For Weight Matrix
template<>
inline void deleteMemoryBuffer<t_inVal>(uint ArrSize, t_inVal**& buffer) {
	assert(buffer != 0);

	for (uint i = 0; i < ArrSize; i++) {
		delete[] buffer[i];
		//delete the memory for each char array in array first
	}
	//then delete the memory to the array of pointers
	delete[] buffer;
	buffer = 0;

}
//void deleteWeightMatrixCasc(float** w,int NetSize,char** buffer)
//{
//
//	for(int i = 0; i < NetSize; i++)
//	{
//		delete [] w[i];
//	//	delete [] buffer[i];
//	//delete the memory for each char array in array first
//	}
//	//then delete the memory to the array of pointers
//	delete [] w;
//	//delete [] buffer;
//
//}

/*
   USES LEGACY CASCADE SYNAPSE

 A Binary Learning Rule For Associative learning as Described by Barret & Rossum 2008
 Without stochastic updates the forgetting is too fast and only the last pattern is remembered
 Attempting a recall of any other pattern Fails - Setting q low 0.1 then the most recent pattern can be recalled accuratelly
 But it is also possible to recall past ones for pats 10 in a Net of 500 neur.
 Can be Given A memory pointer which it initialises to Hold cascadeSynapse Objects

* Called During each Trial -
* Makes A symmetric Weight Matrix As Required by a Hopfield Memory recurrent network
1 if bReconstructObjects =true The Objects Are reconstructed everytime - But the memory space is retained -To optimize speed
2 The Vector X is randomized EveryTime function is called
3 Then synObjects Are stimulated
4 Weight Matrix Constructed from synObjects Strength

 * There is the option of allocating the 1st memory to all available synapses (Changes the uniform Distribution)
 * and then testing lifetime
 * @buffer pointer to where objects of CascadeSynapse Are stored - The memory is allocated the 1st time the function runs and then it is re-used. This improves speed
 * @Weight Array Pointer // The memory is allocated the 1st time the function runs and then it is re-used. This improves speed
 * @function To Allocate Array of Synapses - Can be any of Filter Allocation ones Found in CascadeModelSim.cpp
 * The Synapse Objects may be retained in memory and re-used if bReconstructObjects = false.
 * This Option Should be used when we are sure the distribution among cascade states is stationary after memory storage
 */
template<class T>
float** makeWeightMatrix(int NetSize, t_inVal** Xin, uint patCount,
		uint iTrackedIndex, int _CascadeSize, T**& buffer, float**& w,
		vector<T*>& vSyn, vector<string>& slogfiles) {

	bool bReconstructObjects = false; //Destruct and Re-construct objects in memory buffer every time function is called

	//char *buffer; //Temp Memory Allocation pointer
	//If CascadeObjects Have been Previously Inited then pointer Will not be zero and thus no allocation is required
	bool NeedToInitMemory = (buffer == 0); //The Memory Is initialized with the synapses 1once for a whole experiment
	bool bAllocateWVector = (w == 0); //Memory reallocation Switch - Every Time testHopfieldNet is called, the Weight vector is deleted

	T** s;

	T* pseg = 0;
	//float StartStrength = 0.0f;
	size_t sizeBlock = sizeof(T);
	gsl_rng* mprng = g_getRandGeneratorInstance(false);

	//initialize Memory If This is the 1st Call to the function

	if (bAllocateWVector) {
		w = new float*[NetSize];
		assert( w != NULL);
	}

	//Check if the Memory for Vectors and the Buffer for the Objects has not been Initialised
	if (NeedToInitMemory) {
		vSyn.clear(); //Clear Vector of Pointers to Synapses
		vSyn.reserve(0.55*NetSize * NetSize); //Size Vector To save Time
		//Buffer Inited In Alloc Function - Here Make The Array of pointers to buffers
		//buffer = new char*[NetSize]; //Make Pointer Array to Pointers of buffers where the each array of synapses is to stored
		//pair<T*,ptrdiff_t> pmem;
		//pmem = get_temporary_buffer<T*>(NetSize); //new t_inVal[_uiNeuronCount];
		buffer = new T*[NetSize];

		assert( buffer != NULL);
		memset(buffer,0,NetSize*sizeof(T*));
		cout << "Memory for synapses in bytes :" << NetSize*sizeBlock <<  endl;

	} else {//Memory Has been Inited So delete Old Objects to create New over allocated Memory  //Call Destructors Manually

		if (bReconstructObjects)
				for (typename vector<T*>::iterator it = vSyn.begin();it != vSyn.end(); ++it)(*it)->~T();
		}


	//<ALLOCATION> new Synapses And Store pointers In vector
	uint iSynCount = 0;
	for (int i = 0; i < (NetSize); i++) //For EACH NEURON
		{

		if (bAllocateWVector) {
			w[i] = new float[NetSize]; //The Float Strength Matrix - New One On every Cycle
			assert( w[i] != NULL);
		}

		if (NeedToInitMemory) {
			buffer[i] = 0; //Mem Will be re-init in Alloc Function
			//buffer[i] = pseg = new char[NetSize * sizeBlock]; //Alloc Memory Block to Array Of Objects

		} else {
			memset(w[i], 0, NetSize*sizeof(float)); //Reset Contents of weight matrix to 0 for New Round
			pseg = (T*) buffer[i]; //memory Had Been Initialised on Previous Run so just reposition pointer
		}

		s = (T**) (buffer); ///Cast Memory Location Pointer for Objects

		if (bReconstructObjects || NeedToInitMemory) {
			///Initialise Synapses Using Allocation Function - Only 1/2 of NetSize^2 - As there are symmetric connections
			s[i] = (T*)allocSynapseArray<T>(vSyn,(T*)buffer[i], i,(int) _CascadeSize, mprng, 1.0);
			buffer[i] = s[i]; //Buffer Is not by ref so set address here
			iSynCount += i;

		} //If Init Memory Or Objects Required

		if (!s[i] && i>0) //Check Failure
			ERREXIT(500,"makeWeightMatrix<T>: Could not create synapse objects! Out Of memory?");

	} // FOR EACH NEURON // END OF ALLOCATION

	//Number of Synapses Should Be Half a square Matrix - Diagonal which is 0 - No self Connection
	assert((uint)vSyn.size() == 0.5*NetSize*NetSize-0.5*NetSize);

	///--<LEARNING>--///  - LEarn Patterns using synapse Objects And Return Float Weight matrix -
	//1st Time Around Init the synapses with x number of patterns found before the Tracked Memory So the Distribution moves to Equilibrium
	uint istartPoint = 0;
	if (NeedToInitMemory || bReconstructObjects) //Reconstructing Object Requires To run INit Patterns.
	{
		cout << NeedToInitMemory << " Learning " << patCount << " Patterns from " << istartPoint << endl;
	}
	else { //If Objects are not being reconstructed on Every Run, Then Store again only the tracked memory and onwards - The Distribution of Synapses Should Remain Unaltered
		istartPoint = iTrackedIndex;
		//cout << "Learning Pats from :" << istartPoint << endl;
	}

	for (uint j = istartPoint; j < patCount; j++) //For Each Pattern
	{
		///Report Distribution First Time this is Run
		if (NeedToInitMemory && j == iTrackedIndex) {
			cout << "Inited Synapses. Distribution Before Memory Storage is:"	<< endl;
			const char * fname = slogfiles[0].c_str(); //  slogFiles[0].c_str();
			reportStateDistribution(vSyn, NetSize * NetSize, fname);
		}

		double r; //Randomize Vector
		for (int i = 0; i < NetSize; i++)
		{	r	= gsl_rng_uniform(mprng);
			Xin[j][i] = (r < 0.5) ? 1 : -1;
		}//Make Test Vector

		for (int i = 0; i < NetSize; i++) {
			//Randomize DESIGN - Throw in Noise?
//			double r = gsl_rng_uniform(mprng );
//			if (r>0.5 && j > istartPoint)
//				Xin[j][i] = -Xin[j][i];

			for (int k = 0; k <= i; k++) {
				//w[i][k] = 0.0;
				//W_ik From i->k weight is the output of neuron i at training pattern j
				//Binary learning Rule according to barret&Rossum Cannot Work For AutoAssociative Nets
				//w[i][k] = (Xin[j][i] > 0)?1:-1;
				//Symmetric - Refer to Same Objects
				if (k < i) {
					if (Xin[j][i] * Xin[j][k] > 0) {
						s[i][k].handlePOT();
					} else
						s[i][k].handleDEP();

					w[k][i] = w[i][k] = s[i][k].getStrength(); //Copy To Float matrix
				} else
					w[i][k] = 0; //Connection to Self is 0
			} //END OF For Each Afferent

		} //END OF For Each Neuron
	} //For Each Pattern

	return w;
}

//A Binary Learning Rule For Associative learning as Described by Barret&Rossum 2008
//Without stochastic updates the forgetting is too fast and only the last pattern is remembered
//Attempting a recall of any other pattern Fails - Setting q low 0.1 then the most recent pattern can be recalled accuratelly
//But it is also possible to recall past ones for pats 10 in a Net of 500 neur.
float** makeWeightMatrixBin(int NetSize, float** Xin, int patCount) {
	float q[patCount];
	///Probability of Plasticity for each memory Pattern - mem 0 is Allocated so q[0] = 1
	q[0] = 1.0;
	for (int i = 1; i < patCount; i++)
		q[i] = 0.125;

	float** w = new float*[NetSize];
	gsl_rng* mprng = g_getRandGeneratorInstance(true);
	double r;

	//Allocate new Float
	for (int i = 0; i < NetSize; i++) {
		w[i] = new float[NetSize];
		for (int k = 0; k < NetSize; k++) //INITIALIZE RANDOM wEIGHT
				{
			r = gsl_rng_uniform(mprng);
			w[i][k] = (r > 0.5) ? 1 : -1;
		}
	}

	//Learn Patterns
	for (int j = 0; j < patCount; j++) //For Each Pattern
			{
		for (int i = 0; i < NetSize; i++) {
			for (int k = 0; k < NetSize; k++) {
				//w[i][k] = 0.0;
				//W_ik From i->k weight is the output of neuron i at training pattern j
				//Binary learning Rule according to barret&Rossum Cannot Work For AutoAssociative Nets
				//w[i][k] = (Xin[j][i] > 0)?1:-1;
				//Do Stochastic Update
				r = gsl_rng_uniform(mprng);
				if (r < q[j]) {
					w[i][k] = (Xin[j][i] * Xin[j][k] > 0) ? 1 : -1;
				}
				if (i == k)
					w[i][k] = 0; //No Self Connection
				//cout << w[i][k] << " ";
			}
			//cout << endl;
		}
	}
	return w;
}

inline uint calcHammingDistance(t_inVal* X1, t_inVal* X2, uint NetSize) {
	uint HammingDistance = 0;
	for (uint i = 0; i < NetSize; i++) {
		if (X1[i] != X2[i])
			HammingDistance++;
	}

	return HammingDistance;
}

/*
 ////////10/12/12 - Transfered from  LEGACY CASCADESYNAPSE  /////////

 //Training And Testing A Binary Hopfield Net Using Binary Synapses And Binary Neurons As Amit & Fusi 1994
 // Return number of times stored pattern was recalled successfully
 //@iPatCount:			 Total number of patterns to be stored in the weight matrix
 //@NeuronCount : 		 The number of neurons in the recurrent net
 //@ProbeNoiseLevel :	 The probability of bit inversion in the vector used as a memory recall cue.
 //@StoredPatIndex:	 	 The pattern on which to attempt recall after storage.
 //@pAllocationFunct : 	 Options Are: allocCascadeSynArray, allocCascadeSynDoubleThresFiltArray, allocCascadeSynDualFiltArray
 //@vector<cascadeSynapse*> vSyn; //Vector Of Pointers To CascadeSynapses - Filled After Weigth Matrix is Created
 //For Binary Neurons and Synapses Capacity Grows as 6.75826677× log (_uiNeuronCount)

  *
  * For Each Trial
  *  makeWeightMatrix
  *  "Inject Probe Vector with Noise:"
  *  Run n= NetSize * uiNetUpdateCycles Ouput Neuron Updates
  *  Check if Output vector is matches stored pattern within Acceptable Error
  * Next Trial
 */
template<class T>
int testHopfieldBinarySyns(int _iPatCount, uint _uiNeuronCount, t_inVal** X,
		t_inVal* tX, float _fProbeNoiseLevel, int _iStoredPatIndex,
		int _iCascadeSize, float& AvgSignal, uint trials,
		vector<string>& slogFiles, vector<T*>& vSyn, T**& mem_buffer) {
	///Simulation Statistics//
	const uint uiNetUpdateCycles = 20;
	const uint ciInterruptCondition = trials * 0.10; // If the first n trials give 0% or 100% The simulation result is taken as min or Max and stops.

	AvgSignal = 0.0; //Reset
	int AvgSignalSamples = 0;

	gsl_rng* mprng = g_getRandGeneratorInstance(false);

	uint NetSize = _uiNeuronCount;
	int tPatCount = _iPatCount;
	float ProbeNoise = _fProbeNoiseLevel; //The percentage of Inverted bits to the trained pattern to construct probe pattern
	const float fAcceptedError = NetSize * g_fAllowedRecallError; //NetSize*0.05; //Accepted Deviation From Stored pattern - Or when Examining the stability of output pattern
	int StoredPatIndex = _iStoredPatIndex; //The pattern index that is being tracked
	int iCascadeSize = _iCascadeSize;
	//BinaryNeuron aBN[NetSize];

	int HammingDistance 	= 0;
	int HammingDistance2 	= 0;

	uint t = trials; //Timer
	uint rcallHits = 0;
	clock_t start, finish; //TIME ESTIMATION
	double duration = 0;
	double Totalduration = 0;

	float** W = 0; //The Weight Matrix - Initialized in makeWeight function

#ifdef MEM_TEST_VERBOSE
	//Run Network
	cout << "- Test Binary Synapse Learning Hopfield Net - Simple Associative Learning Rule" << endl;
	cout << "----Pat Count: " << _iPatCount << " Cascade Size:" << _iCascadeSize << "----- Net Size:" << NetSize << endl;
	cout.unsetf(ios_base::floatfield);
#endif

	///DO RECALL TRIALS TO OBTAIN STATISTICS
	while (t > 0)
	{
		start = clock();

		if (t%(trials/10))
		{cout << endl << endl << " Cascade Size:" << _iCascadeSize
				<< " Pat Count: " << _iPatCount << " Tracking:"
				<< (StoredPatIndex + 1) << " Trial:" << (trials - t)
				<< " Recall Hits Up to Now:" << rcallHits << endl;}


		if (duration > 0)
			cout << " ETL:" << (t * (Totalduration / (trials - t))) / 60
					<< "mins" << endl;

		///WEIGHT MATRIX INIT
		//cout << (trials-t) << " Making Weight Matrix..." << endl;
		//Construct Weight Matrix
		//Update Allocated Weight MAtrix - 1st Call Allocates Memory - subsequent reuses the memory
		makeWeightMatrix<T>(NetSize, X, tPatCount, _iStoredPatIndex,
				iCascadeSize, mem_buffer, W, vSyn, slogFiles);


		///LOG DISTRIBUTION
		// Record Distribution After Learning around the Initial Trials t=10
		if (t == 10) {
			const char * fname = slogFiles[2].c_str();
			reportStateDistribution<T>(vSyn, NetSize * NetSize, fname);
		}

		//INJECT Probe VECTOR X[0] - Recall First Pattern Stored
		//cout << "Inject Probe Vector with Noise:" << ProbeNoise 	<< " Update Cycles:" << uiNetUpdateCycles << endl;
		//Copy Tracked PAttern To Test Output Vector
		for (uint j = 0; j < NetSize; j++)
		{
			tX[j] = X[StoredPatIndex][j];
			//cout << ((tX[j]>0)?"+":"-");
			if (ProbeNoise > 0) //If the Probe Is noisy
			{
				double r = gsl_rng_uniform(mprng);
				if (r > ProbeNoise)
					tX[j] = -tX[j];
			}
		} //End Loop Creating The imposed(probe) Output Vector
		HammingDistance 	= calcHammingDistance(tX, X[StoredPatIndex], NetSize);
		//END OF PROBE INJECTION

		//// ---  LET IT FREE RUN --//
		//cout << "Free Search of Stable Point" << endl;
		//cout << "Showing Change in Distance to stored pattern E:" << endl;
		uint t_stable = 0;
		uint iSearchTime = 1;

		HammingDistance2 = 0;
		HammingDistance = 0;
		//Update N*uiNetUpdateCycles Randomly picked OutputNeurons State
		for (uint rep = 0; rep < NetSize * uiNetUpdateCycles; rep++) //Maximum Number of Network Updates
		{
			//Randomly Update Network
			double r = gsl_rng_uniform(mprng);
			r = (r == 0.0) ? 0.001 : r; //Check 0 Boundary case
			int i = ceil(r * NetSize) - 1;
			tX[i] = 0;
			//UPDATE NEURON i
			for (uint j = 0; j < NetSize; j++)
				tX[i] += tX[j] * W[j][i]; //Sum All Inputs To this Neuron

			tX[i] = (tX[i] > 0) ? 1.0 : -1.0; //Save Neuron's output - Activation Function
			//HammingDistance += (tX[i] != X[StoredPatIndex][i])?1:-1; //Update in StepWise Manner Starting From 0

			//Once a cycle Is complete Start Measuring the time(updates) during which the output does not change (Accounting the Error Too)
			if (rep > NetSize * iSearchTime) {
				iSearchTime++;
				HammingDistance = calcHammingDistance(tX, X[StoredPatIndex],
						NetSize);
				if (abs(HammingDistance - HammingDistance2) > fAcceptedError) {
					//cout << "E:" << HammingDistance << endl;
					HammingDistance2 = HammingDistance;
					t_stable = 0; //Reset Time We have been Stable
				} else
					t_stable++; //Count Time output has been stable

			}

			if (t_stable > uiNetUpdateCycles)
				break; //Stop Running If we have been stable for 30 Cycles
		}///END Loop Updating the OUtput Neurons

		bool FoundStored = false;
		float Signal = 0.0;

		//Read Output
		//cout << "Calculate Stored pattern's Distance To Output" << endl;
#ifdef MEM_TEST_VERBOSE
		cout << "PatNo\tDist.\tSignal" << endl;
		//MEasure Distance To Stored PAttern#
		//Measure Distance of Output to Stored Patterns
		for (int k=0;k<tPatCount;k++)
		{

			HammingDistance = calcHammingDistance(tX,X[k],NetSize);
			//According To Kempter Leibold 2010 GAmma Is Signal
			Signal = (float)(NetSize-HammingDistance)/NetSize - (float)HammingDistance/NetSize;
			//Check If Pattern Of Interest -
			if (k == StoredPatIndex)
			{
				AvgSignal+= Signal;
				if (HammingDistance < fAcceptedError) //Recall Distance to Tracked Pattern
				{
					FoundStored = true;
					rcallHits++;
				}
			}
			cout << k << "\t" << HammingDistance << "\t" << Signal << endl;
		}
#else //Speed Up Check Only Against Tracked Pattern
		HammingDistance = calcHammingDistance(tX, X[StoredPatIndex], NetSize);
		//According To Kempter Leibold 2010 GAmma Is Signal - This Appears to Match Well Against Recall Probability
		Signal = (float) (NetSize - HammingDistance) / (float)NetSize;	//- (float) HammingDistance / (float) NetSize;
		AvgSignal += Signal;
		AvgSignalSamples++;
		if (HammingDistance < fAcceptedError) //Recall Distance to Tracked Pattern
		{
			FoundStored = true;
			rcallHits++;
		}
#endif

		if (FoundStored)
			cout << "Stored Pattern Found!" << " Search Iter.: " << iSearchTime
					<< endl;
		else
			cout << "Pat.Index:" << StoredPatIndex << " Not Found."
					<< " Search Iter.: " << iSearchTime << endl;
		t--;
		finish = clock();
		duration = (double) (finish - start) / CLOCKS_PER_SEC;
		printf("\n Runtime was %2.1f seconds\n", duration);
		Totalduration += duration;

		//CHECK INTERRUPT CONDITIONs - Out Of Recorded Behaviour - So Stop Now.
		if ((trials - t) == rcallHits && (trials - t) > ciInterruptCondition) { //100% up to now so
			rcallHits = trials; //ShortCut
			cout << "***Interrupting - Looks like a 100%" << endl;
			break;

		}
		if (((trials - t) > ciInterruptCondition) && rcallHits == 0) {
			cout << "***Interrupting - Looks like a 0%" << endl;
			break;
		}

	} //END OF MAIN TRIAL LOOP

	cout << "--------------------------" << endl;
	cout << "Stored :" << (tPatCount - 1)
			<< " Patterns over traced memory. Recall Total Hits : " << rcallHits
			<< "/" << trials << endl;

	cout << " Delete Weight Matrix..." << endl;

	deleteMemoryBuffer(NetSize, W);

	AvgSignal = AvgSignal / (float)AvgSignalSamples;

	return rcallHits;
}

/*
 ///Converted from Legacy CascadeSynapse Is Used
 //HopField Memory Test
  *	Init Pattern Memory for Max 100 patterns
  * for Each Pattern (uint i = 1 + iTrackedMemIndex; i < PatCount; i++)
  * 	testHopfieldBinarySyns
  * 	RecallDuration - The number of patterns during which the tracked pattern was recalled
 */
template<class T>
int SearchForNetCapacity(uint _uiNeuronCount, uint iTrackedMemIndex,
		float _fProbeNoiseLevel, int _iCascadeSize, uint trials,float& AvgRecallSignal, int& RecallSuccessAtCapacityLimit,
		int& RecallDuration,//Counter of consecutive number of Successfule patterns recalled - Filter Rising Signal could start recalling at middlepoint in sequence of patterns
		int NoRecallPatternsCountThreshold, int iMemoryReps = 0)
{
	//T* oCSyn;
	const int iRecallSuccessHitsThreshold 	= (float) trials / 2;


	double dRepIntervalsecs = 1;
	const bool bUseRandomPatterns 			= true;
	uint maxMemCapacity 					= 100+NoRecallPatternsCountThreshold;//0.2*sqrt(_uiNeuronCount); //Number of Patterns to Expect After tracked pattern
	uint PatCount 							= 1 + iTrackedMemIndex + maxMemCapacity; //Max Patterns Created And Thus Max Storage Capacity is NeuronCount/5

	const float Fp 							= 0.5; //Probability of POT Signal
	int NoRecallPatternsCount 				= 0; //Counter of consecutive number of failed patterns to recall

	//const float corrPercent = 0.5;
	int recallHits;

	int iCapacity = 0;
	float AvgSignal;
	double r; //Random Var.

	//This Function Will Handle the Memory Of the Whole Experiment- Improves speed and adds monitoring Control of synapse Population
	T** mem_buffer = 0; //Pointer To Memory Allocated for CascadeSynapses
	vector<T*> vSyn; //Vector Of Pointers To CascadeSynapses - Filled After Weight Matrix is Created
	t_inVal** X = new t_inVal*[PatCount]; //Memory PAtterns Containing The Ones Loaded from File and Random Initialization patterns
	t_inVal tX[_uiNeuronCount]; //Aux Vector

	gsl_rng* mprng = g_getRandGeneratorInstance(false);

	//Init Memory For Patterns -- Vectors Re-Randomized At MakeWeight MAtrix
	for (uint i = 0; i < PatCount; i++) {
		X[i] = new t_inVal[_uiNeuronCount];
		//Make Random Patterns to Fill the Gap to the Ones loaded from the File - Or All Random Patterns
		if (i < iTrackedMemIndex || bUseRandomPatterns) {
			for (uint j = 0; j < _uiNeuronCount; j++) {
				r = gsl_rng_uniform(mprng);
				X[i][j] = (r < Fp) ? 1 : -1; //Make Test Vector
			}
		} else
			memset(X[i], 0, sizeof(float) * _uiNeuronCount); //Just Init to zero to detect creepy Bugs
	}

	char fname[200];
	//sprintf(fname,"randXVectorN_%dPat_%d_Corr_%2.2f.bin",_uiNeuronCount,maxMemCapacity,corrPercent);
	//Load Patterns From File
	if (!bUseRandomPatterns) {
		sprintf(fname, "HD%d-N300-Set67.dat", _uiNeuronCount);
		int iInsertionIndex = 0; //iTrackedMemIndex;
		readTestVectorsFromFile(fname, X, _uiNeuronCount, PatCount,
				iInsertionIndex, false); //Load At 0 All Required PAtterns

		//	//Fix Vector to Scale 0.5 For Hopfield Because this is what we used in previous Hopfield Nets
		for (uint i = iInsertionIndex; i < maxMemCapacity; i++) {
			for (uint k = 0; k < _uiNeuronCount; k++)
				X[i][k] = 0.5 * X[i][k];
		}
	}

	vector<string> slogFiles;
	//Add the File name as the 1st entry- Used by the makeLogFileNames
	slogFiles.push_back(HOPFIELD_OUTPUT_DIRECTORY);

	makeLogFileNames<T>(slogFiles, iTrackedMemIndex, _iCascadeSize, iMemoryReps,
			dRepIntervalsecs, 0.5, trials, _uiNeuronCount * _uiNeuronCount);
	//makeLogFileNames(slogFiles,trials,_iCascadeSize,iTrackedMemIndex,0.5,_uiNeuronCount*_uiNeuronCount,trials,pF);

	const char * foutname = slogFiles[4].c_str();
	if (!foutname)
		ERREXIT(500, "Could not open Output file-Empty filename");
	ofstream ofile(foutname, ios::out); //OPEN OUTPUT FILE

	if (!ofile.is_open()) {
		cerr << "Could not open Output file-Directory Missing? " << foutname
				<< endl;
		ERREXIT(500, "Could not open Output file-Directory Missing?")
	}
	ofile << "PatStored\t R.Hits\t AvgSignal" << endl;
#ifdef MEM_TEST_VERBOSE
	cout << "No.St.Patt\t R.Hits\t AvgSignal" << endl;
#endif

	//Increment Stored PAtterns And MEasure Recall Hits
	AvgRecallSignal 				= 0;
	RecallSuccessAtCapacityLimit 	= 0;
	RecallDuration					= 0;
	for (uint i = 1 + iTrackedMemIndex; i < PatCount; i++)
	{	//Always Recall Pattern zero And Incrementally Increase the stored number of overlayed patterns
		//An Empty Memory buffer pointer and Vector for the Synapses is Passed - Which is filled by the called functions
		//Return the number of times Recall was within the error margin in T Trials
		recallHits = testHopfieldBinarySyns<T>(i, _uiNeuronCount, X, tX, 0.0,
												iTrackedMemIndex, _iCascadeSize, AvgSignal, trials, slogFiles,vSyn, mem_buffer);
		// deleteMemoryBuffer(_uiNeuronCount,mem_buffer);
		//Write To Output File

		ofile << i << "\t " << recallHits << "\t " << AvgSignal << endl;
#ifdef MEM_TEST_VERBOSE
		cout << i << "\t " << recallHits << "\t " << AvgSignal << endl;
#endif
		//Check if Recall Succesσfull for i patterns stored
		if (recallHits < iRecallSuccessHitsThreshold)
		{
			NoRecallPatternsCount++; //Increment Counter of non recall - If No memory is stored for 3  then Stop Trying
		}
		else { //Reset the Norecall possible for X patterns - and store more patterns on top

			iCapacity = i - iTrackedMemIndex; //Store Capacity Up to last successful recall
			AvgRecallSignal 			 = AvgSignal; //The last one saved will be when The PatCount = Capacity
			RecallSuccessAtCapacityLimit = recallHits;

			if (NoRecallPatternsCount == 0) //Is this Success following previous ones?
					RecallDuration++; //Yes:Increment Recall streak Count
			else
					RecallDuration = 0; //No, its new:ReStart counting the streak of consecutive recalls
			//Now Reset The Number of No Recall patterns
			NoRecallPatternsCount = 0;
		}

		//Stop If Recall has been impossible for X number of patterns now
		if (NoRecallPatternsCount >= NoRecallPatternsCountThreshold)
			break; //Normally We expect Capacity to decrease With the Number of Patterns Stored - But with filters this aint so

	}
	//CLEAN UP
	cout << "Cleaning Up Memory..." << endl;
	ofile.close();

	cout << " Call Destructors..." << endl;
	//Call Destructors Manually
	for (typename vector<T*>::iterator it = vSyn.begin(); it != vSyn.end();	++it)
		(*it)->~T();
	vSyn.clear();

	if (mem_buffer != 0)
		deleteMemoryBuffer(_uiNeuronCount, mem_buffer);

	//Delete Pattern Memory
	deleteMemoryBuffer(PatCount, X);


	//gsl_rng_free(mprng); -> Now done in static instance
	g_rng_r =  g_getRandGeneratorInstance(false, true); //Free INstance in static Var

	char buff[150];
	T::getTypeName(buff);
	cout << "-:Fin:- Max Patterns stored: " << iCapacity << " Recall Width:" << RecallDuration << " with Recall " << RecallSuccessAtCapacityLimit << "/" << trials
			<< "  Object Size :" << _iCascadeSize << " Object Type: " << buff << " AllowedRecallError:" << g_fAllowedRecallError
			<< " PatternsStored after last:" << NoRecallPatternsCountThreshold 	<< endl;

	return iCapacity;
}

//Measure Capacity of 1000 Synapses Using Hadamard Vectors as input
void doHopfieldCapacityTest(int modelType, string modelName, uint iNeuronCount,
		uint trials, uint initPatterns, int maxCascSize, int startIndex) {
	const uint iSynCount = iNeuronCount * iNeuronCount;
	const uint itrials = trials; //1000;
	const uint iNoOfInitPatterns = initPatterns;

	int C[maxCascSize];
//	 int  (*Tfunct)(ICascadeSynapse* oCSyn,uint _uiNeuronCount,uint iTrackedMemIndex,float _fProbeNoiseLevel,int _iCascadeSize, uint itrials,int iMemoryReps);

	char fname[200];
	strcpy(fname, modelName.c_str());

	cout << fname << " Hopfield Capacity Report Per Cascade Size:" << endl;

	char params[150];
	sprintf(params, "-HOPFCapacity-T%d-RNDV_N%d_E%1.2f.N%d-%ddat", itrials,
			iNeuronCount, g_fAllowedRecallError, startIndex, maxCascSize);
	strcat(fname, params);

	string fOutName(HOPFIELD_OUTPUT_DIRECTORY);
	fOutName.append(fname);
	// ofstream* ofile = new ofstream(fOutName.c_str(), ios::out ); //OPEN OUTPUT FILE
	string strfname(fname);
	string strDir(HOPFIELD_OUTPUT_DIRECTORY);
	ofstream* ofile = openfile(strDir, fname, ios::app);

	float fRecallSignal;
	int   iRecallDuration = 0;

	if (!(ofile->is_open())) {
		cerr << "Could not open Output file-Directory Missing?"
				<< fOutName.c_str() << endl;
		ERREXIT(500, "Could not open Output file-Directory Missing?")
	}
	(*ofile) << "#" << modelName
			<< "  Synapse Model HOPFIELD Capacity Report Trials:" << itrials
			<< " Init Patts:" << iNoOfInitPatterns << " Neurons:"
			<< iNeuronCount << " Allowed Recall Error:" << g_fAllowedRecallError << endl;

	(*ofile) << "#CascSize\tCapacity\tAvgSignal\tRecallSuccessInTrials\tRecallPatDuration" << endl;
	int rHits; //Recall Counts out of T trials At capacity Limi
	for (int i = startIndex; i <= maxCascSize; i++) {
		g_FilterTh = i;
		g_UpdaterQ = 1.0 / (g_FilterTh * g_FilterTh);
		g_fAllocHThres = getCAthres(i, 0, modelType);

		cout << "SynSz:" << g_FilterTh << " Decay Fc:" << g_fcAMPDecay
				<< " cAMPInj:" << g_fcAMPMagnitude << " h_thres:"
				<< g_fAllocHThres << endl;

		int NoRecallPatternsCountThreshold = 10;
		switch (modelType) {
		case 1: //synapseCascade
			NoRecallPatternsCountThreshold = 10;
			C[i - 1] = SearchForNetCapacity<synapseCascade>( iNeuronCount, initPatterns, 0.0f, i, trials, fRecallSignal, rHits, iRecallDuration, NoRecallPatternsCountThreshold);
			break;
		case 2: //Cascade Filter
			NoRecallPatternsCountThreshold = 0.375*(float)i*i+10;
			C[i - 1] = SearchForNetCapacity<synapseCascadeFilterUnified>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;
		case 3: //Cascade Filter With Decay
			NoRecallPatternsCountThreshold = 0.375*(float)i*i+10;
			C[i - 1] = SearchForNetCapacity<synapseCascadeFilterUnifiedWithDecay>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;
		case 5: // //synapseCascadeFilterDual DUAL Filter
			NoRecallPatternsCountThreshold = 2*i+10; //Peak Is mix OF All DualFilters But largest one is 4
			C[i - 1] = SearchForNetCapacity<synapseFilterDual>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;
		case 7: //Single DUAL Filter
			NoRecallPatternsCountThreshold = 2*i+10; //2 Theta - 1
			C[i - 1] = SearchForNetCapacity<synapseSingleFilterDual>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;
		case 8: //A Single Filter Synapse
			NoRecallPatternsCountThreshold = 0.375*(float)i*i+10;
			C[i - 1] = SearchForNetCapacity<synapseSingleFilterUnifiedWithDecay>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;
		case 9: //A Stochastic Updater Synapse
			NoRecallPatternsCountThreshold = 3;
			C[i - 1] = SearchForNetCapacity<synapseSingleUpdater>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;
		case 11:
			 NoRecallPatternsCountThreshold = 0.375*(float)i*i + 10;
			C[i - 1] = SearchForNetCapacity<synapseSingleFilterUnifiedWithDecayReflecting>( iNeuronCount, initPatterns, 0.0f, i, trials,fRecallSignal,rHits,iRecallDuration,NoRecallPatternsCountThreshold);
			break;

		default:
			cerr << modelName << endl;
			ERREXIT(100, "Unhandled Model Type");
			break;
		};
		///Write To output File
		(*ofile) << i << "\t" << C[i - 1] << "\t" << fRecallSignal <<"\t" << (float)rHits/trials << "\t" << iRecallDuration  << endl;
	}

	for (int i = startIndex; i <= maxCascSize; i++) {
		cout << i << "\t" << C[i - 1] << endl;
	}

	cout << fname << endl;

	ofile->close();
	delete ofile;

}

//Training And Testing A Binary Hopfield Net Using Binary Synapses And Binary Neurons As Amit & Fusi 1994
// Return number of times stored pattern was recalled successfully
//@iPatCount:			 Total number of patterns to be stored in the weight matrix
//Takes tX copies to X and then Updates X by trial*Netsize times - Random selection
int recallHopfieldNetPattern(uint _uiNeuronCount, uint StartNeuron, t_inVal* tX,
		t_inVal* X, float** W, float& AvgSignal, uint trials) {
	///Simulation Statistics//
	const uint uiNetUpdateCycles = trials;

	AvgSignal = 0.0; //Reset
	//pAllocationFunct = &allocCascadeSynDoubleThresFiltArray; //Assign Value To Function Pointer

	//IF THE NET IS INITIALIZED WITH SYNS OF 0 STRENGTH THEN Binary Synapses 0-1 can be used
	gsl_rng* mprng = g_getRandGeneratorInstance(false,false);

	uint NetSize = _uiNeuronCount;
	float fAccError = _uiNeuronCount * 0.001;
	//BinaryNeuron aBN[NetSize];
	int HammingDistance = 0;
	int HammingDistance2 = 0;

	uint rcallHits = 0;
	clock_t start, finish; //TIME ESTIMATION
	double duration = 0;
	double Totalduration = 0;

	start = clock();

	//INJECT Probe VECTOR X[0] - Recall First Pattern Stored
	cout << "Inject Probe Vector  Update Cycles:" << uiNetUpdateCycles << endl;
	memcpy(X, tX, NetSize * sizeof(int));
	//END OF PROBE INJECTION

//It Appears that although the above is instructive it is very hard to accurately recall any pattern  other than the last one!
	//LET IT FREE RUN
	cout << "Free Search of Stable Point" << endl;
	//cout << "Showing Change in Distance to stored pattern E:" << endl;

	HammingDistance2 = 0;
	HammingDistance = 0;
	long ldiff = 0; //Number of Differences of update cycle to previous state
	long timeStable = 0;
	for (uint rep = 0; rep < NetSize * uiNetUpdateCycles; rep++) //Maximum Number of Network Updates
			{
		//Randomly Update Network
		//uint i = round(gsl_rng_uniform(mprng)*(NetSize-1)); //Select A Random Neuron
		uint i = gsl_ran_flat(mprng, StartNeuron, StartNeuron + NetSize - 1);

		if (!bRecallInProgress)
			exit(0); //pthread_exit normally

		int tmp = 0;
		//X[i] = 0; ///Set Output to Zero
		//UPDATE NEURON i
		for (uint j = 0; j < NetSize; j++)
			tmp += X[j] * W[j][i]; //Sum All Inputs To this Neuron

		tmp = (tmp > 0) ? 1 : -1; //Clip output

		if (tmp != X[i])
			ldiff++;

		timeStable++;

		//Now Update Output Neuron
		X[i] = tmp; //Save Neuron's output - Activation Function
		//HammingDistance += (tX[i] != X[StoredPatIndex][i])?1:-1; //Update in StepWise Manner Starting From 0
		if ((rep % NetSize) == 0) //Reset Integration Time and Sum of Error
				{
			timeStable = 0;
			ldiff = 0;
		}

		if ((ldiff < fAccError) && (timeStable > NetSize * 0.5)) {
			cout << "Stable.. Stop." << endl;
			break; //Stop Search
		}

	} //Update for a number of maximum network cycles - Loop Can be interrupted

	cout << "Diff:" << ldiff;
	//Read Output
	finish = clock();
	duration = (double) (finish - start) / CLOCKS_PER_SEC;
	printf("\n Runtime was %2.1f seconds\n", duration);
	Totalduration += duration;

	return rcallHits;
}
