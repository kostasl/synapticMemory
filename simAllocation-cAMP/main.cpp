/*
 * main.cpp
 *
 *  Created on: Dec 15, 2011
 *
 *  Modified in Mar/2012 :
 *  This code has been modified to provide the statistics of same threshold crossings distributions under repetitious memory encoding
 *  the results of wk13-14 2012 has been obtained using this code. The sampling cut off is defined by the input variable metaSampleSize
 *  which defines the number of threshold cycle samples required from each synapse before the simulation ends. After this limit is reached
 *  each synapse stops writing to the global instance distribution map.
 *
 *      Author: kostasl
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
int g_FilterTh 			= 6; //Used for Single Filter Experiments
double g_FilterDecay 	= 0.0; //0.0916986;

uint g_timeToSampleMetaplasticity	= 0; //Used by Sim code as the time to sample the number of metaplastic transitions
int g_MetaplasticitySampleSize		= 1;//Sim Code Stops saving to the distribution of same threshold crossings Once this number of samples has been gathered
double g_UpdaterQ 		= 1.0/(g_FilterTh*g_FilterTh); //The single Updater Transitions - Make sure its in double format
float g_fAllocHThres	= 0.0; //Default Post Synaptic depol. Threshold to switch on Allocation
float g_fcAMPDecay		= 0.0; //The timeconstant for the cAmp alpha process (With 0.5 it takes approx 10 tsteps for a complete wave)
float g_fcAMPMagnitude	= 0.0; //The magnitude for the cAmp alpha process b = (g_cAMPDecay^2)/r where r is the number of repetitions desired to reach PKA thres

uint g_AllocRefraction	= 0;//The Same Threshold Counter Limit required to allocate a synapse --0.375*g_FilterTh*g_FilterTh;
string g_outputTag;

double runContinuousMemoryRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, vector<double> vpReptimes, int CascadeSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile);
void runAllocSignalVsRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, vector<double> vpReptimes, int FilterSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile);

int main(int argc, char* argv[])
{

	clock_t start, finish;

	po::options_description all("Allowed options"); //The group Of All options
	po::options_description basicSim("Simulation Averaging Options ");
	po::options_description singleFilterSim("Single Filter Simulation Options");
	po::options_description cascadeSim("Cascade simulation option");
	po::options_description AllocationOptions("Allocation Experiment simulation options");

	string inputFile,modelName;
	string simulationName = "simRepetition";
	int startIndex,endIndex,simulationType,modelType,synapsesPopulation,trackedMemIndex,initPeriod;
	unsigned int trials;

	long lSimtimeSeconds = 250;
	int RepMemoryIndex =0;

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
		("trials,T", po::value<unsigned int>(&trials)->default_value(10000), "Number of iteration to average over")
	    ("cSimTimeSecs", po::value<long>(&lSimtimeSeconds)->default_value(lSimtimeSeconds), "Duration of continuous time simulation in seconds")
		("synapsesSize", po::value<int>(&synapsesPopulation)->default_value(10000), "The number of synapses to use - Has to match the vector file size where required")
		("inputFile,V", po::value<string>(&inputFile)->default_value("\n"), "The vector input file to use from directory MemoryInputVectors. If No file given then Random Vectors are used.")
		("startSize", po::value<int>(&startIndex)->default_value(1), "The range of model size parameter to begin testing - interpretation is model dependent")
		("endSize", po::value<int>(&endIndex)->default_value(1), "The range of model size parameter to end testing - interpretation is model dependent")
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
		("repTimes,RT", po::value< vector<double> >(&vdRepTime)->multitoken(), "The relevant time intervals a pattern will be repeated after initial encoding")
		("AllocDepolThres,RT", po::value< float >(&g_fAllocHThres)->default_value(g_fAllocHThres), "The relevant time intervals a pattern will be repeated after initial encoding.Set Automatically for simulation: AllocSignalVsRepetitionTime")
		("AllocRefrac,RP", po::value<uint>(&g_AllocRefraction)->default_value(g_AllocRefraction), "The period a synapse needs to be stable before it is allocated-Threshold Counter Tagging");

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
	mapSynapseAllocator["synapseFilterUnified"] = 2;
	mapSynapseAllocator["synapseFilterUnifiedWithDecay"] = 3;// (pAllocationFunct)allocSynapseArrayCascade<synapseCascadeFilterUnifiedWithDecay>;
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
	}

	start = clock();
	float minRequiredEncodings = 3.0f;
	double cPeakcAMPToTheta = 0.138384;
	if (simulationType == 1) //Switch the Simulation Type
	{
		double dMSFPT;
		//For Cascade Indexes/Symmetric Filter Sizes
		for (int i=startIndex;i<=endIndex;i++)
		{
			 g_FilterTh =i;
			 g_UpdaterQ = 1.0/(g_FilterTh*g_FilterTh);
			 g_fAllocHThres = 	2*0.766/(float)i; //According to beta/Theta law
			 /*Set the decay so Max Amplitude is at peak h(T) time */
			 g_fcAMPDecay		= 1.0/(0.375*(pow(i,2))); //The timeconstant for the cAmp alpha process (With 0.5 it takes approx 10 tsteps for a complete wave)
			 /* Magnitude of cAMP set so 3 properly spaced repetitions are required to reach PKA threshold*/
			 g_fcAMPMagnitude	= 1.0/(minRequiredEncodings*cPeakcAMPToTheta*(pow(i,2)));  //The magnitude for the cAmp alpha process beta=1/0.138384Theta^2r where r is the number of repetitions desired to reach PKA thres

			 dMSFPT = runContinuousMemoryRepetition(modelType,ts,trials,trackedMemIndex,RepMemoryIndex,vdRepTime,i,synapsesPopulation,lSimtimeSeconds,dEncodingRate,inputFile);
			 cout << i << "\t MFPT :" << dMSFPT << endl;
		 }//Loop For Each Cascade Index
	}

	if (simulationType == 10) //Run simulation To obtain Alloc Signal Size Vs Reptime for each Theta
	{
		//For Cascade Indexes/Symmetric Filter Sizes
		for (int i=startIndex;i<=endIndex;i++)
		{
			 g_FilterTh = i;
			 g_UpdaterQ = 1.0/(g_FilterTh*g_FilterTh);
			 g_fAllocHThres = 	2*0.766/(float)i; //According to beta/Theta law

			 /*Set the decay so Max Amplitude is at peak h(T) time - Needs a fix cause r=3 does not allocated for  */
			 g_fcAMPDecay		= 1.0/(0.375*(pow(i,2))); //The timeconstant for the cAmp alpha process (With 0.5 it takes approx 10 tsteps for a complete wave)
			 /* Magnitude of cAMP set so 3 properly spaced repetitions are required to reach PKA threshold*/
			 g_fcAMPMagnitude	= 1.0/(4*0.138384*(pow(i,2)));  //The magnitude for the cAmp alpha process beta=1/0.138384Theta^2r where r is the number of repetitions desired to reach PKA thres
			 cout << "PKA alpha:" << g_fcAMPDecay << " beta:" << g_fcAMPMagnitude << endl;

			 runAllocSignalVsRepetition(modelType,ts,trials,trackedMemIndex,RepMemoryIndex,vdRepTime,i,synapsesPopulation,lSimtimeSeconds,dEncodingRate,inputFile);
		 }//Loop For Each Cascade Index
	}

  ///Measure Duration
  finish = clock();

  ///Print Duration of Run -
  double duration = (double)(finish - start) / CLOCKS_PER_SEC;//Eclipse Reports Problem with  CLOCKS_PER_SEC But it compiles normally -Eclipse Bug
  printf( "\n Runtime was %2.1f seconds\n", duration );
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
	double dMFPT; //The Time of Mean signal Crosses the noise
	double dRepIntervalsecs = 0;
	pair<double,double> dAllocSignal;

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
	cout << "#########" << endl;
	cout << "Allocation at h>" << g_fAllocHThres << endl;
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
	 pF =  (pAllocationFunct)allocSynapseArrayCascade<synapseCascade>;
	 synapseCascade* oCSyn; //Local To this block
	 oCSyn = (synapseCascade*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"simMemRepetitionAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 makeLogFileNames<synapseCascade>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 //dMFPT =
	 simRepetitionAllocation<synapseCascade>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 2: //Cascade Filter
{
	 pF =  (pAllocationFunct)allocSynapseArrayCascade<synapseCascadeFilterUnified>;
	 synapseCascadeFilterUnified* oCSyn;
	 oCSyn = (synapseCascadeFilterUnified*)(*pF)(mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);

	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 //mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 makeLogFileNames<synapseCascadeFilterUnified>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 //dMFPT =
	 simRepetitionAllocation<synapseCascadeFilterUnified>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);

}
break;

case 3: //Cascade Filter
{
	 pF =  (pAllocationFunct)allocSynapseArrayCascade<synapseCascadeFilterUnifiedWithDecay>;
	 synapseCascadeFilterUnifiedWithDecay* oCSyn;
	 oCSyn = (synapseCascadeFilterUnifiedWithDecay*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 makeLogFileNames<synapseCascadeFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //(pFunct pF, int iSynCount,int iCascadeSize,unsigned int iSimTime,unsigned int ciInitPeriod,double mdRate, double dFp=0.5)
	 //Also Available : simMemSignalinContinuousTime
	 //dMFPT =
	 simRepetitionAllocation<synapseCascadeFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 7:
{ //Single DUAL Filter
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterDual>;
	 synapseSingleFilterDual* oCSyn;
	 oCSyn = (synapseSingleFilterDual*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 makeLogFileNames<synapseSingleFilterDual>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //dMFPT =
	 simRepetitionAllocation<synapseSingleFilterDual>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 8: //A Single Filter Synapse
{
	 synapseSingleFilterUnifiedWithDecay* oCSyn;
	slogFiles.clear();
	slogFiles.push_back(fOutName);
	makeLogFileNames<synapseSingleFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs,0.5,trials, synapsesPopulation);

	oCSyn = (synapseSingleFilterUnifiedWithDecay*)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecay>((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally
	//dMFPT =
	simRepetitionAllocation<synapseSingleFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}

break;

case 9: //A Stochastic Updater Synapse
{
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleUpdater>;
	 synapseSingleUpdater* oCSyn;
	 oCSyn = (synapseSingleUpdater*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 ::makeLogFileNames<synapseSingleUpdater>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //simMemSignalinContinuousTime<synapseSingleUpdater,pAllocationFunct>(pF, synapsesPopulation,i,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate);
	 //dMFPT =
	 ::simRepetitionAllocation<synapseSingleUpdater>(oCSyn, synapsesPopulation,CascadeSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
}
break;

case 11: //U Filter Reflecting Boundary
{
	 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecayReflecting>;
	 synapseSingleFilterUnifiedWithDecayReflecting* oCSyn;
	 oCSyn = (synapseSingleFilterUnifiedWithDecayReflecting*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)CascadeSize,mprng,1.0);
	 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
	 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

	 makeLogFileNames<synapseSingleFilterUnifiedWithDecayReflecting>(slogFiles,trackedMemIndex,CascadeSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
	 //dMFPT =
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


return dMFPT;
}


/*
 * Simulation Scans different single Repetition Times To obtain the final allocated signal Size.
 * These are stored to a file so we can plot AllocSignalSize Vs Repetition Time For each Theta assuming Alloc Threshold Is set to peak signal
 */

void runAllocSignalVsRepetition(int modelType,double ts, long trials, int trackedMemIndex, int RepMemoryIndex, vector<double> vpReptimes, int FilterSize, long synapsesPopulation, long lSimtimeSeconds, double dEncodingRate, string inputFile)
{
	pair<double,double> AllocSignal; //The Time of Mean signal Crosses the noise
	double dRepIntervalsecs = 0;


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

	string sAggregateFile = fOutName;
	sAggregateFile.append("AllocSignalVsRepTime-cAMP_n");
	sAggregateFile += boost::lexical_cast<std::string>(FilterSize);
	sAggregateFile.append("_T");
	sAggregateFile += boost::lexical_cast<std::string>(trials);
	sAggregateFile.append(".dat");

	cout << "Signal Output Files: " <<  sAggregateFile << endl; //Tell User Which Output file we are using
	ofstream ofile(sAggregateFile.c_str(), ios::out ); //Open Data File
	if (!ofile.is_open())
		ERREXIT(101,"Could Not Open output files. Check directories");
	//Write Header
	ofile << "#RepTime\tAllocSNR\tAllocSignal\tAllocVariance\tAllocThreshold" << endl;

	const float MaxRepTime = PeakTime*2;
	int iRepIntervalStep = (MaxRepTime/20.0);
	if (iRepIntervalStep == 0) iRepIntervalStep = 1; //Fix For Small Filters
	dRepIntervalsecs = iRepIntervalStep; //Start Rep

	while (dRepIntervalsecs <= MaxRepTime)
	{
		cout << "#########" << endl;
		 clock_t start2 = clock();

		repetitionTable.clear();
		//Add the INdex of the memory with  A key Being the Time When it should be repeated and the value -> The Pattern Number cInitPeriod+repIndex

		int  iabsRepTime = dRepIntervalsecs +(RepMemoryIndex)*(1.0/ts)*dEncodingRate;
		for (int i=1;i<5;i++)
			repetitionTable[iabsRepTime*i] = trackedMemIndex+RepMemoryIndex; //Use the Absolute Pattern Number

		cout << "ts:" << ts << " Repetition of Memory " << RepMemoryIndex << " 4 times with interval :" << (dRepIntervalsecs) << endl;

		cout << "Allocation at h >" << g_fAllocHThres << endl;

		char* mem_buffer 			= 0;		//This Pointer is filled by allocMem, To point to the reuseable allocated memory
		gsl_rng* mprng 				= g_getRandGeneratorInstance(true);

		switch (modelType)
		{
		case 8: //A Single Filter Synapse
		{
			 synapseSingleFilterUnifiedWithDecay* oCSyn;
			slogFiles.clear();
			slogFiles.push_back(fOutName);
			makeLogFileNames<synapseSingleFilterUnifiedWithDecay>(slogFiles,trackedMemIndex,FilterSize,dRepIntervalsecs,0.5,trials, synapsesPopulation);

			oCSyn = (synapseSingleFilterUnifiedWithDecay*)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecay>((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally
			AllocSignal = simRepetitionAllocation<synapseSingleFilterUnifiedWithDecay>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
		}
		break;
		case 11: //U Filter Reflecting Boundary
		{
			 pF =  (pAllocationFunct)allocSynapseArraySingleQ<synapseSingleFilterUnifiedWithDecayReflecting>;
			 synapseSingleFilterUnifiedWithDecayReflecting* oCSyn;
			 oCSyn = (synapseSingleFilterUnifiedWithDecayReflecting*)(*pF)((char*)mem_buffer,synapsesPopulation,(int)FilterSize,mprng,1.0);
			 if (!oCSyn) ERREXIT(500,"MemAllocation: Could not create synapse objects! Out Of memory?");
			 mem_buffer = (char*)oCSyn; //If New Allocation Then Update The Mem Buff pointer to newly Allocated Space \\TODO:Could Happen Internally

			 makeLogFileNames<synapseSingleFilterUnifiedWithDecayReflecting>(slogFiles,trackedMemIndex,FilterSize,dRepIntervalsecs, 0.5,trials, synapsesPopulation);
			 AllocSignal = simRepetitionAllocation<synapseSingleFilterUnifiedWithDecayReflecting>(oCSyn, synapsesPopulation,FilterSize,trackedMemIndex,(char*)inputFile.c_str(), trials,lSimtimeSeconds,dEncodingRate,repetitionTable,ts,slogFiles);
		}
		break;
		default:
			cerr << "Don't Have a repetition simulation associated with this type "<< modelType << " of object." << endl;
			ERREXIT(500,"Object not recognized for repetition simulation");
			break;
		}

		ofile << iabsRepTime << "\t"<< AllocSignal.first/sqrt(AllocSignal.second) << "\t" << AllocSignal.first << "\t" << AllocSignal.second << "\t" << g_fAllocHThres << endl;
		//Clear Object Memory
		return_temporary_buffer(mem_buffer);

		clock_t finish2 = clock();
		float duration = (double)(finish2 - start2) / CLOCKS_PER_SEC;//Eclipse Reports Problem with  CLOCKS_PER_SEC But it compiles normally -Eclipse Bug
		printf( "\n SubSim Runtime was %2.1f secs. Est.Time left: %2.1f mins \n\n", duration, (duration*(MaxRepTime - dRepIntervalsecs)/iRepIntervalStep)/60);

		dRepIntervalsecs += iRepIntervalStep; //increment the repetition Time
	}

ofile.close();

}
