/*
 * main.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: kostasl
 *
 *   This is the Mean First Passage Time Project - Adapted from the Lifetime SNR Project
 *   The simulation Currently supports the measuring MFPT of the 1st tracked Memory. Simulation Runs until Signal Drop Below 0
 */

#include "common.h"
#include "ContinuousTime/ContinuousTimeExperiments.h"
//#include "../synapseModels/common.h"

#include "InputVectorHandling.h"
#include "synapseAllocators.h"

#include <boost/program_options.hpp> //Located in usr/include --Using /usr/local/boost_1_42_0/stage/lib/libboost_program_options.a statically linked lib

namespace po = boost::program_options;
using namespace std;

//GLOBAL VARS
int g_FilterTh = 6; //Used for Single Filter Experiments
double g_FilterDecay = 0.0; //0.0916986;
double g_UpdaterQ = 1.0/(g_FilterTh*g_FilterTh); //The single Updater Transitions - Make sure its in double format
string g_outputTag;



void runMFPTWithRepetition(double _oMFPT[],double _oMFPTVar[],int modelType,double ts,long trials,int trackedMemIndex,int RepMemoryIndex,int RepCount,double dRepIntervalsecs ,int CascadeSize,long synapsesPopulation,double dEncodingRate,int iRepetitions);

int main(int argc, char* argv[])
{

	clock_t start, finish;
	start = clock();
	map<string,unsigned int> mapSimType;
	map<string,int> mapSynapseAllocator; //An association of a the target object name With the allocation Function for the synapse Population

	po::options_description all("Allowed options"); //The group Of All options
	po::options_description basicSim("Simulation Averaging Options ");
	po::options_description singleFilterSim("Single Filter Simulation Options");
	po::options_description cascadeSim("Cascade simulation option");
	po::options_description AllocationOptions("Allocation Experiment simulation options");

	char file[200] = "\n";

	string inputFile,simulationName,modelName;
	simulationName = "MFPT"; //Default Value
	int startIndex,endIndex,simulationType,modelType,synapsesPopulation,trackedMemIndex,initPeriod,iRepetitions;
	initPeriod = 0; //Fix -Otherwise Tracked Signal not measured
	unsigned int trials;
	//bool bUseCascadeParadigm;

	long lSimtimeSeconds = 200;
	int RepMemoryIndex = 10;
	double dRepTime = g_FilterTh*g_FilterTh*0.375;//Default Val is Peak Signal
	double ts = 1.000;//When Set to 1 simu. is in discrete Time
	double dEncodingRate = 1.0;
	vector<uint> pviTrackedIndex; //The List of Memories to Track - No Allocation Signal
	vector<uint> pviAllocIndex; //The list of Memories to Allocate - Allocation Signal Enabled

	file[0] = 0; //Set to Null So we can detect if input has been given

	basicSim.add_options()
	    ("help", "produce help message")
	    ("model,M", po::value<string>(&modelName), "The model to run the simulation on")
		("simulation,S", po::value<string>(&simulationName)->default_value(simulationName), "The simulation name to run")
	    ("startSize", po::value<int>(&startIndex)->default_value(g_FilterTh), "A parameter Representing the Model Size to start testing. For cascade this is the number of states in each Cascade strength")
	    ("endSize", po::value<int>(&endIndex)->default_value(g_FilterTh), "A parameter representing the Model Size to end testing. For cascade this is the number of states in each Cascade strength")
		("trials,T", po::value<unsigned int>(&trials)->default_value(10000), "Number of iteration to average over")
	    ("cSimTimeSecs,secs", po::value<long>(&lSimtimeSeconds)->default_value(lSimtimeSeconds), "Duration of continuous time simulation in seconds")
		("synapsesSize", po::value<int>(&synapsesPopulation)->default_value(10000), "The number of synapses to use - Has to match the vector file size where required")
		("inputFile,V", po::value<string>(&inputFile)->default_value("\n"), "The vector input file to use from directory MemoryInputVectors. If No file given then Random Vectors are used.")
		("UnAllocTID", po::value< vector<uint> >(&pviTrackedIndex)->multitoken(), "Pattern indexes to Track that will not be allocated. These patterns will signal No allocation")
		("AllocateTID", po::value< vector<uint> >(&pviAllocIndex)->multitoken(), "Pattern Index on which to enable the global Allocation Signal");

	singleFilterSim.add_options()
	   ("FilterTh", po::value<int>(&g_FilterTh)->default_value(6), "The Single Filter's upper&Lower threshold Value")
	   ("FilterDecay,D", po::value<double>(&g_FilterDecay)->default_value(g_FilterDecay), "The Single Filter's Decay Value");

	cascadeSim.add_options()
	    ("initPeriod,I", po::value<int>(&initPeriod)->default_value(0), "Memories used to initialise Synapses")
	    ("StochUpdQ,Q", po::value<double>(&g_UpdaterQ)->default_value(g_UpdaterQ), "A small string to attach to the outputFiles so they The experiment is identified")
	    ("Timestep,ts", po::value<double>(&ts)->default_value(ts), "Sim. Timstep in seconds. Set to 1.0 for a discrete time simulation.")
	    ("cEncodingRate,R", po::value<double>(&dEncodingRate)->default_value(dEncodingRate), "Encoding rate - Or Rate of Incoming Patterns Default: 1sec");

	AllocationOptions.add_options()
		("repPatIndex,RI", po::value<int>(&RepMemoryIndex)->default_value(RepMemoryIndex), "The index of the pattern to repeat relative to the 1st tracked pattern")
		("repCount,RC", po::value<int>(&iRepetitions)->default_value(0), "The number of times a pattern should be repeatedly encoded at RepTimeIntervals")
		("repInterval,RT", po::value<double>(&dRepTime)->default_value(dRepTime), "The number of times a pattern should be repeatedly encoded at RepTimeIntervals");


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
	mapSimType["MFPT"] = 9;

	mapSynapseAllocator["synapseCascade"] = 1;// (pAllocationFunct)allocSynapseArrayCascade<synapseCascade>;
	//mapSynapseAllocator["synapseFilterUnified"]
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
		ERREXIT(1,"No simulation argument Specified");
	}

	if (simulationType==0)
		ERREXIT(2,"SimulationType not recognized");

	if (vm.count("model") == 0)
		ERREXIT(1,"No model argument Specified");

	if (mapSynapseAllocator.find(modelName) == mapSynapseAllocator.end()  )
	{	ERREXIT(2,"Model name not recognized");}
	else
	{
		modelType = mapSynapseAllocator[modelName];
		cout << modelName << endl;

		if (modelType == 9)
			cout << "Stoch. Updater Q:" << g_UpdaterQ << endl;
		if (modelType == 8)
			cout << "Single Filter. st. Theta:" << startIndex  << " Decay" << g_FilterDecay << endl;
	}


	//simulationType=9 MFPT Save Results To file
	//OPEN OUTPUT FILES
	string buffTargetDir(MFPTIMES_OUTPUT_DIRECTORY);
	string buffFilename;
	buffTargetDir.append("_").append(modelName).append("/");
	char buff[200];

	ofstream* ofile[8];

	for (int j=0;j<8;j++)
	{
		sprintf(buff,"MFPT-N%d_%d-%d_T%d_lst%.2f_ts%.2f.dat",synapsesPopulation,startIndex,endIndex,trials,(float)(j/100.0),ts);
		buffFilename.assign(buffTargetDir);
		buffFilename.append(buff);

		ofile[j] = new ofstream(buffFilename.c_str(), ios::app ); //Open Data File for Appending So you dont Overwrite Previous Results
		if (!ofile[j]->is_open()) ERREXIT(errno,"Could not Open output file");
		//////LOG File Opened////
		*ofile[j] << "#Size\tMFTP\tSTDVar\tlSigThres" << endl;
	}

	cout << "Running Simulation: " << simulationName << " Output File:" << buffFilename.c_str() << endl;


	//For Cascade Indexes
	double dMFPT[8],dMFPTVar[8];
	for (int i=startIndex;i<=endIndex;i++)
	{
	 g_FilterTh =i;
	 //g_UpdaterQ = 1.0/(g_FilterTh*g_FilterTh);

	 runMFPTWithRepetition(dMFPT,dMFPTVar,modelType,ts,trials,trackedMemIndex,RepMemoryIndex,iRepetitions,dRepTime,i,synapsesPopulation,dEncodingRate,iRepetitions);

	 for (int j=0;j<8;j++)
		 *ofile[j] << i << "\t" << dMFPT[j] << "\t" << sqrt(dMFPTVar[j]) << "\t" << (float)j/100 << endl;
	}//Loop For Each Cascade Index

	 for (int j=0;j<8;j++)
		 ofile[j]->close();

	 ///Measure Duration
  finish = clock();
  ///Print Duration of Run - //TODO: This gives the wrong Time When using Threads!

  double duration = (double)(finish - start) / (CLOCKS_PER_SEC*60); //CLOCKS_PER_SEC;//Eclipse Reports Problem with  CLOCKS_PER_SEC But it compiles normally -Eclipse Bug
  printf( "\n Runtime was %2.1f minutes\n", duration );
 // std::exit(0);
return 0;
}



/*
 * This Is the Main Function that sets Up and Calls the Allocation Continuous Time Experiments
 * A repetition table is setup with the time of repetition and the memory index required.
 * Each Pattern added for repetition is automatically added to the list of tracked Pattern
 *
 * Selects the appropriate Continuous Time Simulation To Run from the various objects and allocation Functions Available
 */

void runMFPTWithRepetition(double _oMFPT[],double _oMFPTVar[],int modelType,double ts,long trials,int trackedMemIndex,int RepMemoryIndex,int RepCount,double dRepIntervalsecs ,int CascadeSize,long synapsesPopulation,double dEncodingRate,int iRepetitions)
{
cout << "MEAN FIRST PASSAGE TIME Continuous/Discrete time With a Repeated memory pattern" << endl;
pAllocationFunct pF;
t_patt_reptbl repetitionTable;

unsigned long lrepInterval = dRepIntervalsecs*(1.0/ts); //Normalize from Seconds to Timesteps

///Rep Times Need to Account For the init period
for (int i=1;i<=RepCount; i++)
{
	repetitionTable[i*lrepInterval+(RepMemoryIndex)*(1.0/ts)] = RepMemoryIndex+trackedMemIndex; //Set Repetition Time At the peak of the tracked pattern
	cout << "ts:" << ts << " Repetition of " << RepMemoryIndex << " at t:" << (i*lrepInterval+RepMemoryIndex) << endl;
}
cout << "####" << endl;
switch (modelType)
{
case 1: //synapseCascade
	 pF =  (pAllocationFunct)allocSynapseArrayCascade<synapseCascade>;
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 simMFPT<synapseCascade,pAllocationFunct>(_oMFPT,_oMFPTVar,pF, synapsesPopulation,CascadeSize,trackedMemIndex, trials,dEncodingRate,repetitionTable,ts);
break;
case 3: //Cascade
	 pF =  (pAllocationFunct)allocSynapseArrayCascade<synapseCascadeFilterUnifiedWithDecay>;
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 simMFPT<synapseCascadeFilterUnifiedWithDecay,pAllocationFunct>(_oMFPT,_oMFPTVar,pF, synapsesPopulation,CascadeSize,trackedMemIndex, trials,dEncodingRate,repetitionTable,ts);
break;
case 7:
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterDual>;
	 simMFPT<synapseSingleFilterDual,pAllocationFunct>(_oMFPT,_oMFPTVar,pF, synapsesPopulation,CascadeSize,trackedMemIndex, trials,dEncodingRate,repetitionTable,ts);
break;
case 8: //A Single Filter Synapse
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecay>;
	 simMFPT<synapseSingleFilterUnifiedWithDecay,pAllocationFunct>(_oMFPT,_oMFPTVar,pF, synapsesPopulation,CascadeSize,trackedMemIndex, trials,dEncodingRate,repetitionTable,ts);
break;
case 9: //A Stochastic Updater Synapse
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleUpdater>;
	 //simMemSignalinContinuousTime<synapseSingleUpdater,pAllocationFunct>(pF, synapsesPopulation,i,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate);
	 simMFPT<synapseSingleUpdater,pAllocationFunct>(_oMFPT,_oMFPTVar,pF, synapsesPopulation,CascadeSize,trackedMemIndex, trials,dEncodingRate,repetitionTable,ts);
break;
case 11: //U Filter Reflecting Boundary
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecayReflecting>;
	 simMFPT<synapseSingleFilterUnifiedWithDecayReflecting,pAllocationFunct>(_oMFPT,_oMFPTVar,pF, synapsesPopulation,CascadeSize,trackedMemIndex, trials,dEncodingRate,repetitionTable,ts);
break;

default:
	cerr << "Don't Have a repetition simulation associated with this type "<< modelType << " of object." << endl;
	ERREXIT(500,"Object not recognised for repetition simulation");
	break;
}


}




