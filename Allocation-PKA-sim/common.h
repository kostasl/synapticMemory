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
	#define HOPFIELD_OUTPUT_DIRECTORY "..//simResults/HopFieldMemResultsV5RND//"
	#define PERCEPTRON_OUTPUT_DIRECTORY "..//simResults//PerceptronMemResultsV5RND//"
	#define MTL_OUTPUT_DIRECTORY "//..//simResults//MemoryLifetimeResultsV2//"
	#define MTLC_OUTPUT_DIRECTORY "//..//simResults//MemoryLifetimeContinuousResults//"
	#define MTLCEVNT_OUTPUT_DIRECTORY "//..//simResults//MemoryLtEventDrivenContinuousResults//"
	#define ALLOCCTEVNT_OUTPUT_DIRECTORY "..//simResults//MemAllocationPKAResults//"
	#define PCAPACITY_CONTINUOUS_OUTPUT_DIRECTORY "..//simResults//PerceptronCapacityContinuousResults//"

	#define MTLMEANF_OUTPUT_DIRECTORY "//..//simResults//MemoryLifetimeResultsMeanField//"
	#define SIGNALS_OUTPUT_DIRECTORY "//..//simResults//SignalsComparison3//"
	#define INPUT_VECTORS_DIRECTORY "..//simResults//MemoryInputVectors//"
	#define ESCAPETIMES_OUTPUT_DIRECTORY "..//simResults//EscapeTimesResults//"
	#define MFPTIMESSNR_OUTPUT_DIRECTORY "..//simResults//MFPT-SNR//"
	#define THRESCYCLE_OUTPUT_DIRECTORY "..//simResults//ThresCycles//"
#else
	#define HOPFIELD_OUTPUT_DIRECTORY "results//HopFieldMem//"
	#define PERCEPTRON_OUTPUT_DIRECTORY "results//PerceptronMem//"
	#define MTL_OUTPUT_DIRECTORY "results//MemoryLifetime//"
	#define MTLMEANF_OUTPUT_DIRECTORY "results//MemoryLifetimeResultsMeanField//"
	#define SIGNALS_OUTPUT_DIRECTORY "results//SignalsComparison//"
	#define INPUT_VECTORS_DIRECTORY "results//MemoryInputVectors//"
	#define ALLOCCTEVNT_OUTPUT_DIRECTORY "results//MemAllocationPKAResults//"
#endif

	#define NTHREADS 30 //THREADING ?? MAXIMUM THREADS
//#define USE_NATIVE_RAND
//#define TH_FILT_FREEZE_THRESHOLD //Freezes the growth of the thresholds
#define USE_CASCADE_DEFAULT_PROBABILITIES //Cascade Synapses Initialize default geometric progression P and Q transition probabilities
#define UNIFILTER_RESET_ZERO


#define ERREXIT(code, str) errexit((int)(code), __LINE__ ,__FILE__, str);

static void errexit(int code,uint lineno ,const char* srcFile,const char* str)
{
	//fprintf(stderr,"line:%d",__LINE__);
	fprintf(stderr,"%s line %d : %s: %s\n",srcFile,lineno,(str),strerror(code));
	exit(1);
}


///Simulation Global Variables
//gsl_rng * rng_r; //Used by GSL Rand Num Generator
extern char FilePath[_MAX_PATH]; // _MAX_PATH represents the longest possible path on this OS


///GLOBAL Instance and FunctioN
static gsl_rng* g_rng_r = 0;

class cascadeSynapse; //Need this Prototype Here For Template Allocator Functions - Implementation of the class is however now Legacy
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
using namespace std;

typedef int t_inVal; //Input Memory Vector value types (-1,1) used so int

typedef pair< pair<double,double>, pair<double,double> > t_simRet; //The Return type for simRepetitionAllocation
//Declare for Legacy reasons
//typedef synapseAllocator<cascadeSynapse>::pFunct pAllocationFunct; //Old Naming style for the two allocator functions
typedef synapseAllocator<ICascadeSynapse>::pFunct pAllocationFunct;
typedef synapseAllocator<ICascadeSynapse>::vpFunct pVectorAllocationFunct2;



typedef map<unsigned long,unsigned int> t_patt_reptbl; //Holds pairs of time, pattern index - For repeating input patterns
typedef map<unsigned int,unsigned int>  t_patt_trackedtbl; //Holds the list of index numbers of Tracked Patterns

//Found in mltExperiments


#endif // _STD_INC
