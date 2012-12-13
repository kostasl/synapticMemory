/*
 * main.cpp
 *  Author: kostasl
 *  Created on: Dec 13, 2011
 *
 * Runs the Hopfield Memory Test
 *
 *
 */
#include <boost/lexical_cast.hpp>
#include "common.h"
#include "synapseAllocators.h"
#include "HopfieldMemoryTests.h"

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
float g_fPKAAllocThres	= 1000; //Threshold beyong which the integrating PKA signal switches allocation ON
float g_fInjectionGain	= 1.0; //The GAIN of the cAMP production Process
int g_iHillOrder		= 4; //The threshold function hill order
uint g_AllocRefraction	= 0;//The Same Threshold Counter Limit required to allocate a synapse --0.375*g_FilterTh*g_FilterTh;

float g_fAllowedRecallError = 0.05;
string g_outputTag;


double runContinuousMemoryRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, vector<double> vpReptimes, int CascadeSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile);
void runAllocSignalVsRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, int iRepMemoryCount, int FilterSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile);


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
	int startIndex,endIndex,simulationType,modelType,iNeuronPopulation,trackedMemIndex,initPeriod;
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
	    ("cSimTimeSecs", po::value<long>(&lSimtimeSeconds)->default_value(lSimtimeSeconds), "Duration of continuous time simulation in seconds")
		("NetSize", po::value<int>(&iNeuronPopulation)->default_value(100), "The number of synapses to use - Has to match the vector file size where required")
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
		("AllocDepolThres,RT", po::value< float >(&g_fAllocHThres)->default_value(g_fAllocHThres), "The relevant time intervals a pattern will be repeated after initial encoding.Set Automatically for simulation: AllocSignalVsRepetitionTime")
		("AllocRefrac,RP", po::value<uint>(&g_AllocRefraction)->default_value(g_AllocRefraction), "The period a synapse needs to be stable before it is allocated-Threshold Counter Tagging")
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
	if (simulationType == 4) //Simulate Signal In Time MLT
	{

		////OPEN OUTPUT FILES To save The point When MEAN signal Drops below SNR=1
//		string buffFilename(ALLOCCTEVNT_OUTPUT_DIRECTORY);
//		buffFilename.append(modelName);
//		char buff[100];
//		std::sprintf(buff,"_FPT-N%d_%d-%d_T%d_ts%.2f.dat",synapsesPopulation,startIndex,endIndex,trials,ts);
//		buffFilename.append(buff);
//
//		cout << "@ Simulation " << simulationName << " Output File:" << buffFilename.c_str() << endl;
//		ofstream ofile(buffFilename.c_str(), ios::app ); //Open Data File for Appending So you dont Overwrite Previous Results
//
//		if (!ofile.is_open())
//			ERREXIT(errno,"Could not Open output file");
//
//		//////LOG File Opened////
//		ofile << "#First Passage Time Is where SNR=1 - That is the point where on Avg Signal crosses 0" << endl;
//		ofile << "#Size\tMSFPT" << endl;


		cout << "****HOPFIELD MEMORY LIFETIME TEST******" << endl;
		double dMSFPT;
		//For Cascade Indexes/Symmetric Filter Sizes
			 g_fcAMPMagnitude	= 1.0;  //
			 doHopfieldCapacityTest(modelType,modelName, iNeuronPopulation,trials, trackedMemIndex,endIndex,startIndex);

	}//END IF SIMULATION TYPE (1) MLT

  ///Measure Duration
  finish = clock();

  ///Print Duration of Run -
  double duration = (double)(finish - start) / CLOCKS_PER_SEC;//Eclipse Reports Problem with  CLOCKS_PER_SEC But it compiles normally -Eclipse Bug
  std::printf( "\n Runtime was %2.1f seconds\n", duration );
 // std::exit(0);
return 0;
}


