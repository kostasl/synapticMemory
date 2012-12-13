// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#ifndef _STD_INC
#define _STD_INC

#ifdef __linux__
	#define _REENTRANT
	#define _POSIX_SOURCE
#endif

#define _USE_MATH_DEFINES //IN order to use Math COnstants
//#define USE_CUDA //If Set GPU is used for measuring the Signal in Continuous Time Allocation

#include <assert.h> //For Debug
#include <stdio.h>
#include <stdint.h> //For Small DataTypes
#include <iomanip> //For Set Precision
#include <string.h> //strcat
#include <stdlib.h> // for srand ( ) and rand ( ) and _itoa
#include <time.h> // for time ( ) and time_t
#include <iostream>
#include <fstream> //File Streams
//#include <direct.h> // for getcwd
#include <math.h>     // for exp(), log(), and log10()
#include <stdlib.h> //To have abort()
#include <new> //for parameter new allocation
#include <memory>

#include <vector>
#include <set>
#include <map>
#include <ctime>


// ###GSL Note: For the library to work in MSVC, I had to change to the Multithreaded version WinGsl_md.lib
// Also under Properties->C/C++->Code GEneration->Run Time Library Change to Multithreaded Debug
//#include <WinGsl.h >
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


typedef int t_inVal; //Input Memory Vector value types (-1,1) used so int

///Program Parameters

#define MAX_AFFERENTS 1202
#define MAX_SPIKES	2500 //Number of simultanuous injected spikes that can be monitored by the neuron
#define G_MAX		0.0151  //0.015 //The Song Conductance Max VAlue used in IFNeuron And SynapseEnsemble
#define G_INH		0.05 //0.05 //The Song INH Conductance Fixed VAlue used in IFNeuron And SynapseEnsemble
#define IFFIRERATE_PERIOD 0.5 ///The Period which the IFNeuron calculates the Fire rate and Refreshes the SpikeList
#define SPIKE_EXPIRATION_TIMECONSTS 2.2 //Sets When A Spike is considered Expired
//#define DEBUG_LOG

#define USE_SONG_CONDUCTANCE //When Defined the simpler Song method is used to calculate gex
///Need to change this and make a separate Song Synapse
//#define USE_SONG_LEARNING	// Synaptic modification implemented as the double exponential rule and switch rule is ignored.
//#define VERBOSE
//#define MEM_TEST_VERBOSE //Hopfield Mem Test will Write recall hits to output

#define _MAX_PATH 550
#ifndef IRIDIS
	#define HOPFIELD_OUTPUT_DIRECTORY "..//simResults/HopFieldMemRND//"
	#define PERCEPTRON_OUTPUT_DIRECTORY "..//simResults//PerceptronMemResultsV5RND//"
	#define MTL_OUTPUT_DIRECTORY "//..//simResults//MemoryLifetimeResultsV2//"
	#define MTLC_OUTPUT_DIRECTORY "//..//simResults//MemoryLifetimeContinuousResults//"
	#define MTLCEVNT_OUTPUT_DIRECTORY "//..//simResults//MemoryLtEventDrivenContinuousResults//"
	#define ALLOCCTEVNT_OUTPUT_DIRECTORY "..//simResults//MemAllocationResults//"
	#define PCAPACITY_CONTINUOUS_OUTPUT_DIRECTORY "..//simResults//PerceptronCapacityContinuousResults//"

	#define MTLMEANF_OUTPUT_DIRECTORY "//..//simResults//MemoryLifetimeResultsMeanField//"
	#define SIGNALS_OUTPUT_DIRECTORY "//..//simResults//SignalsComparison3//"
	#define INPUT_VECTORS_DIRECTORY "..//simResults//MemoryInputVectors//"
	#define ESCAPETIMES_OUTPUT_DIRECTORY "..//simResults//EscapeTimesResults//"
	#define MFPTIMESSNR_OUTPUT_DIRECTORY "..//simResults//MFPT-SNR//"
#else
	#define HOPFIELD_OUTPUT_DIRECTORY "results//HopFieldMem//"
	#define PERCEPTRON_OUTPUT_DIRECTORY "results//PerceptronMem//"
	#define MTL_OUTPUT_DIRECTORY "results//MemoryLifetime//"
	#define MTLMEANF_OUTPUT_DIRECTORY "results//MemoryLifetimeResultsMeanField//"
	#define SIGNALS_OUTPUT_DIRECTORY "results//SignalsComparison//"
	#define INPUT_VECTORS_DIRECTORY "results//MemoryInputVectors//"
	#define ALLOCCTEVNT_OUTPUT_DIRECTORY "results//MemAllocationResults//"
	#define MFPTIMESSNR_OUTPUT_DIRECTORY "results//MFPT-SNR//"
#endif

	#define NTHREADS 30 //THREADING ?? MAXIMUM THREADS
//#define USE_NATIVE_RAND
//#define TH_FILT_FREEZE_THRESHOLD //Freezes the growth of the thresholds
#define USE_CASCADE_DEFAULT_PROBABILITIES //Cascade Synapses Initialize default geometric progression P and Q transition probabilities
#define UNIFILTER_RESET_ZERO


#define ERREXIT(code, str) errexit((int)(code), __LINE__ ,__FILE__, str);


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


static void errexit(int code,uint lineno ,const char* srcFile,const char* str)
{
	//fprintf(stderr,"line:%d",__LINE__);
	fprintf(stderr,"%s line %d : %s: %s\n",srcFile,lineno,(str),strerror(code));
	exit(1);
}


using namespace std;


static std::ofstream* openfile(string strDir,string strFile,ios::openmode omode)
{
	string strbuff(strDir);
	strbuff.append(strFile);

	cout << strbuff << endl;
	std::ofstream* file = new ofstream(strbuff.c_str(), omode ); //Open Data File for Appending So you dont Overwrite Previous Results
		if (!file->is_open())
		{
			cerr << strDir;
			string cmd = "mkdir ";
			cmd.append(strDir);
			cout <<  "Create Output directory" << endl;
			int ret = system(cmd.c_str());
			cout << "Ret:" << ret << endl;
			if (!ret == 0)
				ERREXIT(ret,"Missing path to output directory-could not create Model output directory");

			file = new ofstream(strFile.c_str(), ios::app ); //Open Data File for Appending So you dont Overwrite Previous Results
			if (!file->is_open())
				ERREXIT(errno,"Could not Open output file");
		}

return file;
}




///Simulation Global Variables
//gsl_rng * rng_r; //Used by GSL Rand Num Generator
extern char FilePath[_MAX_PATH]; // _MAX_PATH represents the longest possible path on this OS

//extern gsl_rng * rng_r; //Used by GSL Rand Num Generator



///GLOBAL Instance and FunctioN
static gsl_rng* g_rng_r = 0;

//class cascadeSynapse; //Need this Prototype Here For Template Allocator Functions - Implementation of the class is however now Legacy
class ICascadeSynapse;
//Synapse Model Allocation Function Templated to Fit Old And New Synapse Classes
//Wrapper Struct required as no Template function pointers are allowed
template <typename T>
struct synapseAllocator {
    typedef T* (*pFunct)(char*,int,int,gsl_rng*,float);
    typedef std::vector<T*>* pvSyns;
    typedef pvSyns (*vpFunct)(char*,int,int,gsl_rng*,float); //The vector Equivalent
    pFunct pF; //Accessory Member Declared
};

//Declare for Legacy reasons
//typedef synapseAllocator<cascadeSynapse>::pFunct pAllocationFunct; //Old Naming style for the two allocator functions
typedef synapseAllocator<ICascadeSynapse>::pFunct pAllocationFunct;
typedef synapseAllocator<ICascadeSynapse>::vpFunct pVectorAllocationFunct2;
using namespace std;


typedef map<unsigned long,unsigned int> t_patt_reptbl; //Holds pairs of time, pattern index - For repeating input patterns
typedef map<unsigned int,unsigned int>  t_patt_trackedtbl; //Holds the list of index numbers of Tracked Patterns

//Found in mltExperiments


#endif // _STD_INC
