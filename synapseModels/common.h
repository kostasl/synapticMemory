
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
#include <iostream>
#include <fstream> //File Streams
//#include <direct.h> // for getcwd
#include <cmath>     // for exp(), log(), and log10()
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


//#define USE_NATIVE_RAND
//#define TH_FILT_FREEZE_THRESHOLD //Freezes the growth of the thresholds
#define USE_CASCADE_DEFAULT_PROBABILITIES //Cascade Synapses Initialize default geometric progression P and Q transition probabilities
#define UNIFILTER_RESET_ZERO



#define liberrexit(code,str)                          \
   {fprintf(stderr,"%s: %s\n",(str),strerror(code)); \
    exit(1);}



///GLOBAL Instance and FunctioN - Called by some Default Constructors
//The calling code must implemend the Rand Instance function
extern gsl_rng* g_getRandGeneratorInstance(bool newInstance);


using namespace std;



#endif // _COMMON_INC
