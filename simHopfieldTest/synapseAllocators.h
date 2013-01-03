/*
 * SynapseAllocators.h
 * File Contains the Functions Called by simulations to initialize the arrays synapse objects
 * //TODO:A Better implementation of Allocators would be to take the a sample object as parameter and use Copy constructors to initialise the population
 * This would remove the required Global Variables
 *  Created on: Dec 19, 2011
 *      Author: kostasl
 */

#ifndef SYNAPSEALLOCATORS_H_
#define SYNAPSEALLOCATORS_H_


#include "../synapseModels/ICascadeSynapse.h"
#include "../synapseModels/synapseCascade.h"
#include "../synapseModels/synapseSingleUpdater.h"
#include "../synapseModels/synapseSingleFilterUnifiedWithDecay.h"
#include "../synapseModels/synapseSingleFilterDual.h"
#include "../synapseModels/synapseSingleFilterUnifiedWithDecayReflecting.h"


extern double g_UpdaterQ;
extern int g_FilterTh;
extern double g_FilterDecay;
extern uint g_AllocRefraction;
extern int g_MetaplasticitySampleSize;



//A Generic Allocation Function For all ICascadeSynapse Type Objects
//NOTES: This Can be converted to return a T*, But.. Issues with older non template functions arise
template <class T>
ICascadeSynapse* allocSynapseArray(vector<T*> &vSyn,T* buffer, int iSynCount,int iParamSize,gsl_rng* prng_r,float StimRate)
{
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<T*,ptrdiff_t> pmem;

	  if (buffer == 0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<T>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (T*)buffer;
		  pmem.second = iSynCount; //Assume Count is correct from Previous init
	  }

   if (!pmem.first && iSynCount > 0)
	   ERREXIT(100,"Could not allocate memory");

   if (pmem.second < iSynCount)
  	 ERREXIT(100,"Could not allocate all the required memory");

//	   long p =  (long)pmem.first;
//   cout << p << " Next:" << (long)(pmem.first+1) <<  " diff:" << (long)(pmem.first+1)-p << endl;
//   cout << "Bytes Alloc:" << ((long)(pmem.first+1)-p)*pmem.second << endl;
//   cout << "Bytes Required:" << sizeof(synapseCascade)*iSynCount << endl;
   	   T* pObj;
	for (int i =0;i<iSynCount;i++)
	{
		pseg = (char*)(pmem.first+i);

		//Check If Fixed Starting Point Of strengths is used
		if (bFixedStartState)
		{
			if (i < iSynCount/2)
				startStrength = ICascadeSynapse::SYN_STRENGTH_STRONG;
			else
				startStrength = ICascadeSynapse::SYN_STRENGTH_STRONG;

			pObj = new(pseg) T(iParamSize,startStrength, prng_r);
		}
		else//Starting Strength is Random
			pObj = new(pseg) T(iParamSize,prng_r); //Do not Fix the Start index - Uniform Distrib

		vSyn.push_back(pObj);//Add to Target Vector
	}


	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}

//////ALLOCATORS OF SINGLE FILTERS ///
//Allocating Stochastic Updaters At a particular Cascade Index

/*
 * Refactoring Made this redundant
 template <class T>
ICascadeSynapse* allocSynapseArray(char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const int iCascadeSize = 1;
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<T*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<T>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (T*)buffer;
		  pmem.second = iSynCount; //Assume Count is correct from Previous init
	  }

   if (!pmem.first) ERREXIT(100,"Could not allocate memory");
   if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");

//	 long p = (long)pmem.first;
//   cout << p << " Next:" << (long)(pmem.first+1) <<  " diff:" << (long)(pmem.first+1)-p << endl;
//   cout << "Bytes Alloc:" << ((long)(pmem.first+1)-p)*pmem.second << endl;
//   cout << "Bytes Required:" << sizeof(synapseCascade)*iSynCount << endl;

	for (int i =0;i<iSynCount;i++)
	{
		pseg = (char*)(pmem.first+i);

		//Check If Fixed Starting Point Of strengths is used
		if (bFixedStartState)
		{
			if (i < iSynCount/2)
				startStrength = ICascadeSynapse::SYN_STRENGTH_STRONG;
			else
				startStrength = ICascadeSynapse::SYN_STRENGTH_WEAK;
			//CascadeSize
			new(pseg) T(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting Strength is Random
			new(pseg) T(iCascadeSize,IndexOfTransitionProb,prng_r); //Do not Fix the Start index - Uniform Distrib
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}
*/


#endif /* SYNAPSEALLOCATORS_H_ */
