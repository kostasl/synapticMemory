/*
 * main.cpp
 *  Author: kostasl
 *  Created on: Dec 15, 2011
 *
 *  Modified in Mar/2012 :
 *  This code has been modified to provide the statistics of same threshold crossings distributions under repetitious memory encoding
 *  the results of wk13-14 2012 has been obtained using this code. The sampling cut-off is defined by the input variable metaSampleSize
 *  which defines the number of threshold cycle samples required from each synapse before the simulation ends. After this limit is reached
 *  each synapse stops writing to the global instance distribution map.
 *
 * PKA Branch 14/9/12:
 * Code Modified to implement Fused h(t) and DA(t) signalling with cAMP saturation for PKA production
 * Here cAMP production is:
 *  u'(t) = G (u_max- u(t)) DA(t)h(t) - Fc u(t)
 *
 * u_max =1 setting the saturation upper bound (Increasing Bound only scales the final PKA)
 *
 *
 */
#include <boost/lexical_cast.hpp>
#include "common.h"
#include "ContinuousTime/ContinuousTimeExperiments.h"
//#include "../synapseModels/common.h"

#include "InputVectorHandling.h"
#include "synapseAllocators.h"

#include <boost/program_options.hpp> //Located in usr/include --Using /usr/local/boost_1_42_0/stage/lib/libboost_program_options.a statically linked lib

namespace po = boost::program_options;
using namespace std;

//GLOBAL VARS
int g_FilterTh 			= 7; //Used for Single Filter Experiments
double g_FilterDecay 	= 0.0; //0.0916986;

uint g_timeToSampleMetaplasticity	= 0; //Used by Sim code as the time to sample the number of metaplastic transitions
int g_MetaplasticitySampleSize		= 0;//Sim Code Stops saving to the distribution of same threshold crossings Once this number of samples has been gathered
double g_UpdaterQ 		= 1.0/(g_FilterTh*g_FilterTh); //The single Updater Transitions - Make sure its in double format
float g_fAllocHThres	= 0.0; //Default Post Synaptic depol. Threshold to switch on Allocation
float g_fcAMPDecay		= 0.01; //The timeconstant for the cAmp alpha process (With 0.5 it takes approx 10 tsteps for a complete wave)
float g_fcAMPMagnitude	= 0.0;
double g_dcAMPMax		= 1.0;// A globally set  saturation value of cAMP.
float g_fPKAAllocThres	= -1000; //Threshold beyond which the integrating PKA signal switches allocation ON - Set -Ve means no allocation
float g_fInjectionGain	= 1.0; //The GAIN of the cAMP production Process
int g_iHillOrder		= 4; //The threshold function hill order
uint g_AllocRefraction	= 1;//The Same Threshold Counter Limit required to allocate a synapse - Set to -1 For No Allocation
string g_outputTag;

// Signal h_thresholds {Theta,repetitions,MaxSignal<-Used as Threshold}
const float g_fCAthresUFilter[][4][3] = {
					{{2, 1, 0.335173}, {2, 2, 0.62145}, {2, 4, 0.771769}, {2, 8,0.793688}},
					{{3, 1, 0.243198}, {3, 2, 0.437438}, {3, 4,0.598619}, {3, 8, 0.651432}},
					{{4, 1, 0.187282}, {4, 2, 0.334954}, {4, 4, 0.473701}, {4, 8, 0.531825}},
					{{5, 1, 0.151413}, {5, 2, 0.270386}, {5, 4, 0.388254}, {5, 8, 0.443233}},
					{{6, 1, 0.126832}, {6, 2, 0.226367}, {6, 4, 0.327752}, {6, 8, 0.377805}},
					{{7, 1, 0.109031}, {7, 2, 0.194552}, {7, 4, 0.2831}, {7, 8, 0.328315}},
					{{8, 1, 0.0955744}, {8, 2, 0.170522}, {8, 4, 0.248942}, {8, 8,0.289863}},
					{{9, 1, 0.0850569}, {9, 2, 0.151749}, {9, 4, 0.22203}, {9, 8, 0.259249}},
					{{10, 1, 0.0766153}, {10, 2,0.136684}, {10, 4, 0.200308}, {10, 8, 0.234358}},
					{{11, 1, 0.0696926}, {11, 2, 0.124332}, {11, 4, 0.182421}, {11, 8, 0.213751}},
					{{12, 1, 0.0639141}, {12, 2, 0.114022}, {12, 4, 0.167445}, {12, 8, 0.196428}},
					{{13, 1, 0.0590183}, {13, 2, 0.105287}, {13, 4, 0.154726}, {13, 8, 0.181671}}
					};
const float g_fCAthresSUSpaced[][4][3] = {
		{  {2, 1,0.220624}, {2, 2, 0.334349}, {2, 4, 0.423187}, {2, 8,0.453064}},
		   {{3, 1, 0.0853399}, {3, 2, 0.137476}, {3, 4, 0.188786}, {3, 8, 0.215084}},
		   {{4, 1, 0.045726}, {4, 2, 0.0751888}, {4, 4, 0.106405}, {4, 8, 0.124745}},
		   {{5, 1, 0.0286135}, {5, 2, 0.0474927}, {5, 4, 0.0681678}, {5, 8, 0.0810866}},
		   {{6, 1, 0.0196291}, {6, 2, 0.0327453}, {6, 4, 0.0473656}, {6, 8, 0.0568081}},
		   {{7, 1, 0.0143155}, {7, 2, 0.0239536}, {7, 4, 0.0348114}, {7, 8, 0.0419639}},
		   {{8, 1, 0.010908}, {8, 2, 0.0182878}, {8, 4, 0.0266585}, {8, 8, 0.0322437}},
		   {{9, 1, 0.00859046}, {9, 2, 0.0144217}, {9, 4, 0.0210669}, {9, 8, 0.0255396}},
		   {{10, 1, 0.00694197}, {10, 2, 0.0116654}, {10, 4, 0.0170661}, {10, 8, 0.020724}},
		   {{11, 1, 0.00572721}, {11, 2, 0.00963094}, {11, 4, 0.0141054}, {11, 8, 0.01715}},
		   {{12, 1, 0.0048061}, {12, 2, 0.00808635}, {12, 4, 0.0118532}, {12, 8, 0.0144253}},
		   {{13, 1, 0.00409094}, {13, 2, 0.00688595}, {13, 4, 0.0101003}, {13, 8, 0.0123011}},
		   {{14, 1, 0.00352451}, {14, 2, 0.00593452}, {14, 4, 0.00870925}, {14, 8,  0.0106132}},
		   {{15, 1, 0.00306823}, {15, 2, 0.00516761}, {15, 4, 0.00758697}, {15, 8, 0.00924995}}
};


double runContinuousMemoryRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, vector<double> vpReptimes, int CascadeSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile);
void runAllocSignalVsRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, int iRepMemoryCount, int FilterSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile);


//Calcium threshold Lookup table Mapping
double getCAthres(int theta, int reps,int modelType)
{
if (modelType == 8) //synapseSingleFilterUnifiedWithDecay
{
	switch(reps)
	{
	case 1:
		return g_fCAthresUFilter[theta-2][0][2]; //obtain CA filter
	case 2:
		return g_fCAthresUFilter[theta-2][1][2]; //obtain CA filter
	case 4:
		return g_fCAthresUFilter[theta-2][2][2]; //obtain CA filter
	case 8:
		return g_fCAthresUFilter[theta-2][3][2]; //obtain CA filter
	default:
		//Return r=4 by default
		return g_fCAthresUFilter[theta-2][2][2]; //obtain CA filter
	}
}
if (modelType == 9 || modelType == 1) //SU Synapse
{
	switch(reps)
	{
	case 1:
		return g_fCAthresSUSpaced[theta-2][0][2]; //obtain CA filter
	case 2:
		return g_fCAthresSUSpaced[theta-2][1][2]; //obtain CA filter
	case 4:
		return g_fCAthresSUSpaced[theta-2][2][2]; //obtain CA filter
	case 8:
		return g_fCAthresSUSpaced[theta-2][3][2]; //obtain CA filter
	default:
		//Return r=4 by default
		return g_fCAthresSUSpaced[theta-2][2][2]; //obtain CA filter
	}
}

return 0.0;
}

#define OUTPUT_FILENAME "_AllocSignalVsRepTime-NoSat-PKA_n"

int main(int argc, char* argv[])
{

	clock_t start, finish;

	po::options_description all("Allowed options"); //The group Of All options
	po::options_description basicSim("Simulation Averaging options ");
	po::options_description singleFilterSim("Single Filter Simulation options");
	po::options_description cascadeSim("Cascade simulation options");
	po::options_description AllocationOptions("PKA Allocation Experiments - simulation options");

	string inputFile,modelName = "synapseSingleFilterUnifiedWithDecay"; //Default
	string simulationName = "simRepetition";
	int startIndex,endIndex,simulationType,modelType,synapsesPopulation,trackedMemIndex,initPeriod;
	unsigned int trials;

	long lSimtimeSeconds = 250;
	int RepMemoryIndex = 0;
	int RepMemoryCount = 0;

	double ts = 1.000;//When Set to 1 simu. is in discrete Time
	double dEncodingRate = 1.0;
	vector<uint> pviAllocIndex; //The list of Memories to Allocate - Allocation Signal Enabled
	vector< double > vdRepTime;//Vector of repetition times
	map<string,unsigned int> mapSimType;
	map<string,int> mapSynapseAllocator; //An association of a the target object name With the allocation Function for the synapse Population


	basicSim.add_options()
	    ("help", "produce help message")
	    ("model,M", po::value<string>(&modelName), "The model to run the simulation on")
		("simulation,S", po::value<string>(&simulationName)->default_value(simulationName), "The simulation name to run")
		("trials,T", po::value<unsigned int>(&trials)->default_value(100), "Number of iteration to average over")
	    ("cSimTimeSecs", po::value<long>(&lSimtimeSeconds)->default_value(lSimtimeSeconds), "Duration of continuous time simulation in seconds-For AllocTests its the recording time after the last repetition.")
		("synapsesSize", po::value<int>(&synapsesPopulation)->default_value(10000), "The number of synapses to use - Has to match the vector file size where required")
		("inputFile,V", po::value<string>(&inputFile)->default_value("\n"), "The vector input file to use from directory MemoryInputVectors. If No file given then Random Vectors are used.")
		("startSize", po::value<int>(&startIndex)->default_value(7), "The range of model size parameter to begin testing - interpretation is model dependent")
		("endSize", po::value<int>(&endIndex)->default_value(7), "The range of model size parameter to end testing - interpretation is model dependent")
		("metaSampleTime", po::value<uint>(&g_timeToSampleMetaplasticity)->default_value(g_timeToSampleMetaplasticity), "Time to sample metaplasticity distribution")
		("metaSampleSize", po::value<int>(&g_MetaplasticitySampleSize)->default_value(g_MetaplasticitySampleSize), "The number of samples to obtain for the metaplasticity cycle distribution");

	singleFilterSim.add_options()
	   ("FilterDecay,D", po::value<double>(&g_FilterDecay)->default_value(g_FilterDecay), "The Single Filter's Decay Value");

	cascadeSim.add_options()

	   ("initPeriod,I", po::value<int>(&initPeriod)->default_value(0), "Memories used to initialise Synapses")
	   ("StochUpdQ,Q", po::value<double>(&g_UpdaterQ)->default_value(g_UpdaterQ), "Stochastic Updater Transition probability q")
	   ("Timestep,ts", po::value<double>(&ts)->default_value(ts), "Sim. Timstep in seconds. Set to 1.0 for a discrete time simulation.")
	   ("cEncodingRate,R", po::value<double>(&dEncodingRate)->default_value(dEncodingRate), "Encoding rate - Or Rate of Incoming Patterns Default: 1sec");

	AllocationOptions.add_options()
		("repPatIndex,RI", po::value<int>(&RepMemoryIndex)->default_value(RepMemoryIndex), "The index of the pattern to repeat relative to the 1st tracked pattern")
		("repPatCount,RC", po::value<int>(&RepMemoryCount)->default_value(RepMemoryCount), "For PKA vs Rep. Interval experiments it sets the number of repetitions")
		("repTimes,RT", po::value< vector<double> >(&vdRepTime)->multitoken(), "The relevant time intervals a pattern will be repeated after initial encoding")
		("AllocDepolThres,RT", po::value< float >(&g_fAllocHThres)->default_value(g_fAllocHThres), "SignalThreshold For Allocation-Set Automatically for simulation: AllocSignalVsRepetitionTime")
		("AllocRefrac,RP", po::value<uint>(&g_AllocRefraction)->default_value(g_AllocRefraction), "The period a synapse needs to be stable before it is allocated-Threshold Counter Tagging-")
		("PKAAllocThres,PK", po::value<float>(&g_fPKAAllocThres)->default_value(g_fPKAAllocThres), "The PKA level above which global allocation is switched on.")
		("cAMPDecay,Fc", po::value<float>(&g_fcAMPDecay)->default_value(g_fcAMPDecay), "cAMP decay F_c rate. Std Vals : 0.5,0.05 or 0.01");

	if (g_MetaplasticitySampleSize == 0)
		g_MetaplasticitySampleSize = -1;///Disable This limit

	all.add(basicSim).add(singleFilterSim).add(cascadeSim).add(AllocationOptions);
	po::variables_map vm;
	po::store(po::parse_command_line(argc,  argv, all), vm);
	po::notify(vm);

	///Add List Of Simulation Types
	mapSimType["simMemSignalinTime"] = 1;
	mapSimType["simMemSignalsFromFile"] = 2;
	mapSimType["PerceptronTest"] = 3;
	mapSimType["HopfieldTest"] = 4;
	mapSimType["simMemSignalinContinuousTime"] = 5;
	mapSimType["simEscTime"] = 6;
	mapSimType["MeanMemoryLifetime"] = 7;
	mapSimType["simRepetition"] = 8;
	mapSimType["ThresholdCycleFq"] = 9;
	mapSimType["AllocSignalVsRepetitionTime"] = 10;


	mapSynapseAllocator["synapseCascade"] = 1;// (pAllocationFunct)allocSynapseArrayCascade<synapseCascade>;
	mapSynapseAllocator["synapseCascadeFilterUnified"] = 2;
	mapSynapseAllocator["synapseCascadeFilterUnifiedWithDecay"] = 3;// (pAllocationFunct)allocSynapseArrayCascade<synapseCascadeFilterUnifiedWithDecay>;
	mapSynapseAllocator["synapseSingleFilterDual"] = 4;//(pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterDual>;
	//mapSynapseAllocator["cascadeDelayed"]
	//mapSynapseAllocator["CascadeSamplingFilter"]
	//mapSynapseAllocator["synapseSingleFilterDual"]
	mapSynapseAllocator["synapseSingleFilterUnifiedWithDecay"] = 8;//(pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecay>;
	mapSynapseAllocator["synapseSingleUpdater"] = 9;//(pAllocationFunct)allocSynapseArraySingleQ<synapseSingleUpdater>;
	//mapSynapseAllocator["synapseSerialCascade"]
	mapSynapseAllocator["synapseSingleFilterUnifiedWithDecayRefl"] = 11;

	trackedMemIndex = initPeriod;

	if (vm.count("simulation"))
		simulationType = mapSimType[simulationName];
	else
	{
		simulationType = mapSimType["AllocSignalVsRepetitionTime"];
		cerr << "No simulation argument Specified, Using default: AllocSignalVsRepetitionTime" << endl;
	}

	if (simulationType==0)
		ERREXIT(2,"SimulationType not recognized");

	if (vm.count("model") == 0)
		ERREXIT(1,"No model argument Specified");

	if (mapSynapseAllocator.find(modelName) == mapSynapseAllocator.end()  )
	{
		cout << "Model Name Can be on of:" << endl;
		for (map<string,int>::iterator it = mapSynapseAllocator.begin(); it!=mapSynapseAllocator.end();++it)
			cout << it->first << endl;
		ERREXIT(100,"Model name not recognized");
	}
	else
	{
		modelType = mapSynapseAllocator[modelName];
	}

	start = clock();
	float minRequiredEncodings = 3.0f;
	double cPeakcAMPToTheta = 0.138384;
	if (simulationType == 1) //Simulate Signal In Time MLT
	{

		////OPEN OUTPUT FILES To save The point When MEAN signal Drops below SNR=1
		string buffFilename(ALLOCCTEVNT_OUTPUT_DIRECTORY);
		buffFilename.append(modelName);
		char buff[100];
		std::sprintf(buff,"_FPT-N%d_%d-%d_T%d_ts%.2f.dat",synapsesPopulation,startIndex,endIndex,trials,ts);
		buffFilename.append(buff);

		cout << "@ Simulation " << simulationName << " Output File:" << buffFilename.c_str() << endl;
		ofstream ofile(buffFilename.c_str(), ios::app ); //Open Data File for Appending So you dont Overwrite Previous Results

		if (!ofile.is_open())
			ERREXIT(errno,"Could not Open output file");

		//////LOG File Opened////
		ofile << "#First Passage Time Is where SNR=1 - That is the point where on Avg Signal crosses 0" << endl;
		ofile << "#Size\tMSFPT" << endl;


		cout << "****MEMORY LIFETIME SIMULATION******" << endl;
		double dMSFPT;
		//For Cascade Indexes/Symmetric Filter Sizes
		for (int i=startIndex;i<=endIndex;i+=2)
		{
			 g_FilterTh 		= i;
			 g_UpdaterQ 		= 1.0/(g_FilterTh*g_FilterTh);
			 g_fAllocHThres 	=  getCAthres(i,(int)vdRepTime.size(),modelType);
			 /*Set the decay so Max Amplitude is at peak h(T) time */
			 /* Magnitude of cAMP */
			 g_fcAMPMagnitude	= 1.0;  //
			 cout << "SynSz:"<< g_FilterTh << " Decay Fc:" << g_fcAMPDecay << " cAMPInj:" << g_fcAMPMagnitude << " h_thres:" << g_fAllocHThres << endl;
			 dMSFPT = runContinuousMemoryRepetition(modelType,ts,trials,trackedMemIndex,RepMemoryIndex,vdRepTime,i,synapsesPopulation,lSimtimeSeconds,dEncodingRate,inputFile);
			 cout << i << "\t Lifetime of Mean signal :" << dMSFPT << endl;

			 ofile << i << "\t" << dMSFPT << endl;
		 }//Loop For Each Cascade Index

		ofile.close();


	}//END IF SIMULATION TYPE (1) MLT

	if (simulationType == 10) //Run simulation To obtain Alloc Signal Size Vs Reptime for each Theta
	{
		cout << "****MEMORY ALLOCATION PER REP.INTERVAL SIMULATION******" << endl;
		//For Cascade Indexes/Symmetric Filter Sizes
		for (int i=startIndex;i<=endIndex;i+=2)
		{
			 g_FilterTh			= i;
			 g_UpdaterQ 		= 1.0/(g_FilterTh*g_FilterTh);
			 g_fAllocHThres 	=  getCAthres(g_FilterTh,4,modelType);
			 //g_fcAMPDecay		= 0.01; //The timeconstant for the cAMP Exp Decay
			 g_fcAMPMagnitude	= 1.0;  //The magnitude for the cAmp alpha process beta=1/0.138384Theta^2r where r is the number of repetitions desired to reach PKA thres
			 cout << "Theta:"<< g_FilterTh << " Decay Fc:" << g_fcAMPDecay << " cAMPInj:" << g_fcAMPMagnitude << " h_thres (Fixed To 4 reps) :" << g_fAllocHThres << endl;

			 runAllocSignalVsRepetition(modelType,ts,trials,trackedMemIndex,RepMemoryIndex,RepMemoryCount ,i,synapsesPopulation,lSimtimeSeconds,dEncodingRate,inputFile);
		 }//Loop For Each Cascade Index
	}

  ///Measure Duration
  finish = clock();

  ///Print Duration of Run -
  double duration = (double)(finish - start) / CLOCKS_PER_SEC;//Eclipse Reports Problem with  CLOCKS_PER_SEC But it compiles normally -Eclipse Bug
  std::printf( "\n Runtime was %2.1f seconds\n", duration );
 // std::exit(0);
return 0;
}



/*
 * This Is the Main Function that sets Up and Calls the Allocation Continuous Time Experiments
 * A repetition table is setup with the time of repetition and the memory index required.
 * Each Pattern added for repetition is automatically added to the list of tracked Pattern
 * RepMemoryIndex : Ignore The InitPatterns Just assume 0 is the first tracked Memory
 * Selects the appropriate Continuous Time Simulation To Run from the various objects and allocation Functions Available
 *
 * Returns the First Passage Time Of the Mean Signal under the noise
 */

double runContinuousMemoryRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, vector<double> vpReptimes, int CascadeSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile)
{
	t_simRet simRet; //The Time of Mean signal Crosses the noise

	double dRepIntervalsecs = 0;
	pair<double,double> dAllocSignal;
	int iRepCount = vpReptimes.size();

	cout << "MLT of Repeated memory pattern:" << endl;
	pAllocationFunct pF;
	t_patt_reptbl repetitionTable;

	//Create A map of repetition time (after init period) and pattern Number in the total patterns stored
	//repetitionTable[Timepoint of repetition] = RepMemoryIndex+InitPeriod (NoOfPatternsStored)
	///Rep Times Need to Account For the init period
	for ( vector<double>::iterator it = vpReptimes.begin();it !=vpReptimes.end(); ++it)
	{
		dRepIntervalsecs += *it;
		//Add the INdex of the memory with  A key Being the Time When it should be repeated and the value -> The Pattern Number cInitPeriod+repIndex
		repetitionTable[dRepIntervalsecs +(RepMemoryIndex)*(1.0/ts)*dEncodingRate] = trackedMemIndex+RepMemoryIndex; //Use the Absolute Pattern Number
		cout << "ts:" << ts << " Repetition of Memory " << RepMemoryIndex << " at t:" << (dRepIntervalsecs+RepMemoryIndex) << endl;
	}
	dRepIntervalsecs = vpReptimes[0];
	cout << "#########" << endl;
	cout <<  " h_thres >" << g_fAllocHThres << endl;
/////////// LOG FILE INIT /////////////////////
	//Add the File name as the 1st entry- Used by the makeLogFileNames
	vector<string> slogFiles; //The list of output file names used
	slogFiles.clear();
	string fOutName(ALLOCCTEVNT_OUTPUT_DIRECTORY);
	slogFiles.push_back(fOutName);
	/////////// END OF LOG FILE INIT //////////////////////////

	char* mem_buffer 			= 0;		//This Pointer is filled by allocMem, To point to the reuseable allocated memory
	gsl_rng* mprng 				= g_getRandGeneratorInstance(true);



switch (modelType)
{
case 1: //synapseCascade
{
	 pF =  (pAllocationFunct)allocSynapseArray<synapseCascade>;
	 synapseCascade* oCSyn; //Local To this block
	 oCSyn = (synapseCascade*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"simMemRepetitionAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 makeLogFileNames<synapseCascade>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 simRet =
	 simRepetitionAllocation<synapseCascade>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 2: //Cascade Filter
{
	 pF =  (pAllocationFunct)allocSynapseArray<synapseCascadeFilterUnified>;
	 synapseCascadeFilterUnified* oCSyn;
	 oCSyn = (synapseCascadeFilterUnified*)(*pF)(mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);

	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

	 makeLogFileNames<synapseCascadeFilterUnified>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 simRet =
	 simRepetitionAllocation<synapseCascadeFilterUnified>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);

}
break;

case 3: //Cascade Filter With Decay
{
	 pF =  (pAllocationFunct)allocSynapseArray<synapseCascadeFilterUnifiedWithDecay>;
	 synapseCascadeFilterUnifiedWithDecay* oCSyn;
	 oCSyn = (synapseCascadeFilterUnifiedWithDecay*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

	 makeLogFileNames<synapseCascadeFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 simRet =
	 simRepetitionAllocation<synapseCascadeFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 7:
{ //Single DUAL Filter
	 pF =  (pAllocationFunct)allocSynapseArray<synapseSingleFilterDual>;
	 synapseSingleFilterDual* oCSyn;
	 oCSyn = (synapseSingleFilterDual*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

	 makeLogFileNames<synapseSingleFilterDual>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 simRet =
	 simRepetitionAllocation<synapseSingleFilterDual>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 8: //A Single Filter Synapse
{
	 synapseSingleFilterUnifiedWithDecay* oCSyn;
	slogFiles.clear();
	slogFiles.push_back(fOutName);
	makeLogFileNames<synapseSingleFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs,0.5,trials, synapsesPopulation);

	oCSyn = (synapseSingleFilterUnifiedWithDecay*)allocSynapseArray<synapseSingleFilterUnifiedWithDecay>((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space
	 simRet =
	simRepetitionAllocation<synapseSingleFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}

break;

case 9: //A Stochastic Updater Synapse
{
	 pF =  (pAllocationFunct)allocSynapseArray<synapseSingleUpdater>;
	 synapseSingleUpdater* oCSyn;
	 oCSyn = (synapseSingleUpdater*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

	 makeLogFileNames<synapseSingleUpdater>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //simMemSignalinContinuousTime<synapseSingleUpdater,pAllocationFunct>(pF, synapsesPopulation,i,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate);
	 simRet =
	 simRepetitionAllocation<synapseSingleUpdater>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 11: //U Filter Reflecting Boundary
{
	 pF =  (pAllocationFunct)allocSynapseArray<synapseSingleFilterUnifiedWithDecayReflecting>;
	 synapseSingleFilterUnifiedWithDecayReflecting* oCSyn;
	 oCSyn = (synapseSingleFilterUnifiedWithDecayReflecting*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

	 makeLogFileNames<synapseSingleFilterUnifiedWithDecayReflecting>(slogFiles,trackedMemIndex,CascadeSize,iRepCount,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 simRet =
	 simRepetitionAllocation<synapseSingleFilterUnifiedWithDecayReflecting>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

default:
	cerr << "Don't Have a repetition simulation associated with this type "<< modelType << " of object." << endl;
	ERREXIT(500,"Object not recognized for repetition simulation");
	break;
}

//Clear Object Memory
return_temporary_buffer(mem_buffer);


return simRet.dMeanSignalLifetime; //Return Time When SNR <= 1
}


/*
 * Simulation Scans different single Repetition Times To obtain the final allocated signal Size.
 * These are stored to a file so we can plot AllocSignalSize Vs Repetition Time For each Theta assuming Alloc Threshold Is set to peak signal
 */
void runAllocSignalVsRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, int iRepMemoryCount, int FilterSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile)
{
	t_simRet  AllocSignal; //The Time of Mean signal Crosses the noise
	double dRepIntervalsecs = 0;
	const int iMemoryReps = iRepMemoryCount; //The number of applied repetitions
	g_fAllocHThres 	=  getCAthres(FilterSize,4,modelType);

	cout << "FUSED h(t) and DA(t) and Saturation for PKA production." << endl;
	cout << "Allocation SNR Against Repetition Interval " << endl;

	pAllocationFunct pF;
	t_patt_reptbl repetitionTable;

	float PeakTime = 0.375*(float)FilterSize*(float)FilterSize;

	//Create A map of repetition time (after init period) and pattern Number in the total patterns stored
	//repetitionTable[Timepoint of repetition] = RepMemoryIndex+InitPeriod (NoOfPatternsStored)
	///Rep Times Need to Account For the init period

	///////// LOG FILE INIT /////////////////////
	//Add the File name as the 1st entry- Used by the makeLogFileNames
	vector<string> slogFiles; //The list of output file names used
	slogFiles.clear();
	string fOutName(ALLOCCTEVNT_OUTPUT_DIRECTORY);
	slogFiles.push_back(fOutName);
	/////////// END OF LOG FILE INIT //////////////////////////

	char buff[300];
	sprintf(buff,"%d%s%d_N%d_T%d_Fc%.2f_r%d.dat",modelType,OUTPUT_FILENAME,FilterSize,synapsesPopulation,trials,g_fcAMPDecay,iMemoryReps);
	string sAggregateFile(buff);
//	sAggregateFile += boost::lexical_cast<std::string>(modelType);
//	sAggregateFile.append(OUTPUT_FILENAME);  sAggregateFile += boost::lexical_cast<std::string>(FilterSize);
//	sAggregateFile.append("_N"); sAggregateFile += boost::lexical_cast<std::string>(synapsesPopulation);
//	sAggregateFile.append("_T"); sAggregateFile += boost::lexical_cast<std::string>(trials);
//	sAggregateFile.append("_Fc");sAggregateFile += boost::lexical_cast<std::string>(g_fcAMPDecay);
//	sAggregateFile.append("_r"); sAggregateFile += boost::lexical_cast<std::string>(iMemoryReps);
//	sAggregateFile.append(".dat");

	cout << "Signal Output Files: " <<  sAggregateFile << endl; //Tell User Which Output file we are using
	ofstream* pfile = openfile(fOutName,sAggregateFile,ios::out);

	//ofstream* ofile(sAggregateFile.c_str(), ios::out ); //Open Data File
	if (!pfile->is_open())
		ERREXIT(101,"Could Not Open output files. Check directories");
	//Write Header
	(*pfile) << "#RepTime\tAllocSNR\tAllocSignal\tAllocVariance\tAllocThreshold\tPKALevel\tPKAVariance\tSNR_FPT" << endl;

	const float MaxRepTime = 100+PeakTime;
	int iRepIntervalStep = 5; //Is in Numerical Model Mathematica Results
	dRepIntervalsecs = 1;

	while (dRepIntervalsecs <= MaxRepTime)
	{
		cout << "#########" << endl;
		 clock_t start2 = clock();

		repetitionTable.clear();
		//Add the Index of the memory with  A key Being the Time When it should be repeated and the value -> The Pattern Number cInitPeriod+repIndex
		//Save RepTime in timesteps not seconds.
		int  iabsRepTime = (dRepIntervalsecs +((double)RepMemoryIndex)*dEncodingRate)*(1.0/ts);
		for (int i=1;i<=iMemoryReps;i++)
		{
			repetitionTable[iabsRepTime*i] = trackedMemIndex+RepMemoryIndex; //Use the Absolute Pattern Number
			cout << "ts:" << ts << " Repetition of Memory " << RepMemoryIndex << " at t:" << (dRepIntervalsecs+RepMemoryIndex) << endl;
		}

		//Fix Simulation Time To be after the last rep time and when cAMP has decayed by 99%
		//lSimtimeSeconds = iMemoryReps*iabsRepTime +  + 3.1/g_fcAMPDecay;
		//No Fix to test Allocated signal a fixed time after last repetition
		long lAllocSignalMeasureTime =  iMemoryReps*iabsRepTime +  lSimtimeSeconds; //Add Sim Time Parameter to last rep

		cout <<"TOTAL Sim.Time " << lAllocSignalMeasureTime << " ts:" << ts << " Rep. Memory " << RepMemoryIndex << " " << iMemoryReps << " times with interval: " << (dRepIntervalsecs) << "secs" << " h_thres >" << g_fAllocHThres << endl;
		cout << "Stability Threshold : " << g_AllocRefraction << endl;

		char* mem_buffer 			= 0;		//This Pointer is filled by allocMem, To point to the reuseable allocated memory
		gsl_rng* mprng 				= g_getRandGeneratorInstance(true);

		switch (modelType)
		{
		case 1: //synapseCascade
		{
			synapseCascade* oCSyn;
			slogFiles.clear();
			slogFiles.push_back(fOutName);
			makeLogFileNames<synapseCascade>(slogFiles,trackedMemIndex,FilterSize,iMemoryReps,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
			oCSyn = (synapseCascade*)allocSynapseArray<synapseSingleFilterUnifiedWithDecay>((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

			AllocSignal = simRepetitionAllocation<synapseCascade>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
		}
		break;
		case 8: //A Single Filter Synapse
		{
			synapseSingleFilterUnifiedWithDecay* oCSyn;
			slogFiles.clear();
			slogFiles.push_back(fOutName);
			makeLogFileNames<synapseSingleFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,FilterSize,iMemoryReps,dRepIntervalsecs,0.5,trials, synapsesPopulation);

			oCSyn = (synapseSingleFilterUnifiedWithDecay*)allocSynapseArray<synapseSingleFilterUnifiedWithDecay>((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space
			AllocSignal = simRepetitionAllocation<synapseSingleFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lAllocSignalMeasureTime,dEncodingRate,repetitionTable,ts,slogFiles);
		}
		break;
		case 11: //U Filter Reflecting Boundary
		{
			synapseSingleFilterUnifiedWithDecayReflecting* oCSyn;
			slogFiles.clear();
			slogFiles.push_back(fOutName);
			makeLogFileNames<synapseSingleFilterUnifiedWithDecayReflecting>(slogFiles,trackedMemIndex,FilterSize,iMemoryReps,dRepIntervalsecs,0.5,trials, synapsesPopulation);

			oCSyn = (synapseSingleFilterUnifiedWithDecayReflecting*)allocSynapseArray<synapseSingleFilterUnifiedWithDecayReflecting>((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

			makeLogFileNames<synapseSingleFilterUnifiedWithDecayReflecting>(slogFiles,trackedMemIndex,FilterSize,iMemoryReps,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
			AllocSignal = simRepetitionAllocation<synapseSingleFilterUnifiedWithDecayReflecting>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lAllocSignalMeasureTime,dEncodingRate,repetitionTable,ts,slogFiles);
		}
		break;
		case 9: // "synapseSingleUpdater"
		{

			synapseSingleUpdater* oCSyn;
			slogFiles.clear();
			slogFiles.push_back(fOutName);
			makeLogFileNames<synapseSingleUpdater>(slogFiles,trackedMemIndex,FilterSize,iMemoryReps,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
			oCSyn = (synapseSingleUpdater*)allocSynapseArray<synapseSingleUpdater>((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space
			AllocSignal = simRepetitionAllocation<synapseSingleUpdater>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lAllocSignalMeasureTime,dEncodingRate,repetitionTable,ts,slogFiles);
		}
		break;

		case 3: //Cascade Filter
		{
			 pF =  (pAllocationFunct)allocSynapseArray<synapseCascadeFilterUnifiedWithDecay>;
			 synapseCascadeFilterUnifiedWithDecay* oCSyn;
			 oCSyn = (synapseCascadeFilterUnifiedWithDecay*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space

			 makeLogFileNames<synapseCascadeFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,FilterSize,iMemoryReps,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
			 simRepetitionAllocation<synapseCascadeFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
		}

		break;
		default:
			cout << "Don't Have a repetition simulation associated with this type "<< modelType << " of object." << endl;
			cerr << "Don't Have a repetition simulation associated with this type "<< modelType << " of object." << endl;
			ERREXIT(100,"Object not recognized for repetition simulation");
			break;
		}

		(*pfile) << iabsRepTime << "\t"<< AllocSignal.pairAllocSignalVal.first /sqrt(AllocSignal.pairAllocSignalVal.second) <<
				"\t" << AllocSignal.pairAllocSignalVal.first << "\t"
				<< AllocSignal.pairAllocSignalVal.second << "\t"
				<< g_fAllocHThres << "\t"
				<< AllocSignal.pairPKAVal.first << "\t"
				<< AllocSignal.pairPKAVal.second << "\t"
				<< AllocSignal.dMeanSignalLifetime << endl;

		//Clear Object Memory
		cout << "R.I:" << iabsRepTime << " PKA:" << AllocSignal.pairPKAVal.first << endl;

		return_temporary_buffer(mem_buffer);

		std::clock_t finish2 = clock();
		float duration = (double)(finish2 - start2) / CLOCKS_PER_SEC;//Eclipse Reports Problem with  CLOCKS_PER_SEC But it compiles normally -Eclipse Bug
		std::printf( "\n SubSim Runtime was %2.1f secs. Est.Time left: %2.1f mins \n\n", duration, (duration*(MaxRepTime - dRepIntervalsecs)/iRepIntervalStep)/60);

		dRepIntervalsecs += iRepIntervalStep; //increment the repetition Time
	}

pfile->close();

}
