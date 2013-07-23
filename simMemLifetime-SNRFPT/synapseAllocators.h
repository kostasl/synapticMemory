/*
 * SynapseAllocators.h
 * File Contains the Functions Called by simulations to initialize the arrays synapse objects
 * Init Details Are found in the Allocation Functions which call the constructors
 *
 * //TODO:A Better implementation of Allocators would be to take the a sample object as parameter and use Copy constructors to initialise the population
 * This would remove the required Global Variables
 *  Created on: Dec 19, 2011
 *      Author: kostasl
 */

#ifndef SYNAPSEALLOCATORS_H_
#define SYNAPSEALLOCATORS_H_


#include "ICascadeSynapse.h"
#include "synapseCascade.h"
#include "synapseSingleUpdater.h"
#include "synapseSingleFilterUnifiedWithDecay.h"
#include "synapseSingleFilterDual.h"
#include "synapseFilterDual.h" //This is a cascade of Dual Filters
#include "synapseSingleFilterUnifiedWithDecayReflecting.h"

extern double g_UpdaterQ;
extern int g_FilterTh;
extern double g_FilterDecay;
extern uint g_AllocRefraction;

//A Generic Allocation Function For all ICascadeSynapse Type Objects
//NOTES: This Can be converted to return a T*, But.. Issues with older non template functions arise
template <class T>
ICascadeSynapse* allocSynapseArray(char*buffer,int iSynCount,int iCascadeSize,gsl_rng* prng_r,float StimRate)
{
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
   if (pmem.second < iSynCount)
  	 ERREXIT(100,"Could not allocate all the required memory");


//	   long p =  (long)pmem.first;
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
				startStrength = ICascadeSynapse::SYN_STRENGTH_STRONG;

			new(pseg) T(iCascadeSize,startStrength, prng_r);
		}
		else//Starting Strength is Random
			new(pseg) T(iCascadeSize,prng_r); //Do not Fix the Start index - Uniform Distrib
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}



//////ALLOCATORS OF SINGLE FILTERS ///
//Allocating Stochastic Updaters At a particular Cascade Index
//GENERIC
/*
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
/*

//Allocating Stochastic Updaters
template <>
ICascadeSynapse* allocSynapseArray<synapseCascade>(char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const int iCascadeSize = 1;
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseCascade*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {  //Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseCascade>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseCascade*)buffer;
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
			new(pseg) synapseCascade(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting Strength is Random
			new(pseg) synapseCascade(g_UpdaterQ,prng_r);
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}
*/

//Allocating Stochastic Updaters
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleUpdater>(char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const int iCascadeSize = 1;
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleUpdater*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleUpdater>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseSingleUpdater*)buffer;
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
			new(pseg) synapseSingleUpdater(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting Strength is Random
			new(pseg) synapseSingleUpdater(g_UpdaterQ,prng_r,g_MetaplasticitySampleSize);
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}




//Allocating Single U Filters: Sets Upper Lower Threshold And the allocation Threshold

//TODO:A Better implementation of Allocators would be to take the a sample object as parameter and use Copy constructors to initialiaze the population
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleFilterUnifiedWithDecay>(char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleFilterUnifiedWithDecay*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialise cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleFilterUnifiedWithDecay>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseSingleFilterUnifiedWithDecay*)buffer;
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
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting State is Random not linked to Cascade States
		{
			new(pseg) synapseSingleFilterUnifiedWithDecay(-g_FilterTh,g_FilterTh,g_FilterDecay,g_MetaplasticitySampleSize, prng_r); //Constructor sets Running Value - Filter State
			//cout << ((synapseSingleFilterUnifiedWithDecay*)(pseg))->getMetaplasticCount() << endl;
		}
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);

	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release

}

////ALLOCATOR OF Single Filters With Fixed Filter StartState Used to Obtain Threshold Cycle Distribution
ICascadeSynapse* allocSynapseArraySingleUFilter(char*& buffer,int iSynCount,int iFilterSize,double dFilterDecay,int iFilterState,int iRequiredCycleSamples, gsl_rng* prng_r, float StimRate)
{
	  //const int iCascadeSize = 1; //Used So File Name Reflects the Filter Size
	 bool bFixedStartState = false;
	 ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleFilterUnifiedWithDecay*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {	//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleFilterUnifiedWithDecay>(iSynCount);
		 buffer = (char*)pmem.first;
		 cout << iSynCount << " MemAllocated :" << iSynCount*sizeof(synapseSingleFilterUnifiedWithDecay) << endl;
	  }
	  else
	  {
		  pmem.first = (synapseSingleFilterUnifiedWithDecay*)buffer;
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
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting State is Random not linked to Cascade States
			new(pseg) synapseSingleFilterUnifiedWithDecay(-iFilterSize,iFilterSize,dFilterDecay,iFilterState,iRequiredCycleSamples, prng_r); //Constructor sets Running Value - Filter State
			//new(pseg) ICascadeSynapse();
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);

	}
	cout << "Stability Threshold : " << g_AllocRefraction << endl;
	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release

}



//Allocating Stochastic Updaters
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleFilterUnifiedWithDecayReflecting>(char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  //const int iCascadeSize = 1; //Used So File Name Reflects the Filter Size
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleFilterUnifiedWithDecayReflecting*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleFilterUnifiedWithDecayReflecting>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseSingleFilterUnifiedWithDecayReflecting*)buffer;
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

		new(pseg) synapseSingleFilterUnifiedWithDecayReflecting(-g_FilterTh,g_FilterTh,g_FilterDecay);
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}




//Allocating Single DUAL FILTER
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleFilterDual>(char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const int iCascadeSize = 1;
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleFilterDual*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleFilterDual>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseSingleFilterDual*)buffer;
		  pmem.second = iSynCount; //Assume Count is correct from Previous init
	  }

   if (!pmem.first) ERREXIT(100,"Could not allocate memory");
   if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");

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
			new(pseg) synapseSingleFilterDual(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting Strength is Random
			//new(pseg) synapseSingleFilterDual(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);
			new(pseg) synapseSingleFilterDual(g_FilterTh,g_FilterTh,0.0); //Using Global VAriable(Bad Practice but needed)
		//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}

/*
 *Allocate Single Object - According to the Required Constuctor
 */



#endif /* SYNAPSEALLOCATORS_H_ */
