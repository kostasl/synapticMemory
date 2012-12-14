/*
 * synapseAllocators.cpp
 *
 *  Created on: 13 Dec 2012
 *      Author: kostasl
 */

#include "common.h"
#include "synapseAllocators.h"
#include "ICascadeSynapse.h"
#include "synapseCascade.h"
#include "synapseSingleFilterUnifiedWithDecay.h"


//Allocating Fusi Cascade
template <>
ICascadeSynapse* allocSynapseArray<synapseCascade>(vector<synapseCascade*> &vSyn,char* buffer,int iSynCount,int iCascadeSize, gsl_rng* prng_r, float StimRate)
{
	 // const int iCascadeSize = 1;
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseCascade*,ptrdiff_t> pmem;

	  if (buffer == 0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseCascade>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseCascade*)buffer;
		  pmem.second = iSynCount; //Assume Count is correct from Previous init
	  }

	  if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");
	  if (!pmem.first && iSynCount > 0) ERREXIT(100,"Could not allocate memory");


//	 long p = (long)pmem.first;
//   cout << p << " Next:" << (long)(pmem.first+1) <<  " diff:" << (long)(pmem.first+1)-p << endl;
//   cout << "Bytes Alloc:" << ((long)(pmem.first+1)-p)*pmem.second << endl;
//   cout << "Bytes Required:" << sizeof(synapseCascade)*iSynCount << endl;
   synapseCascade* pObj;
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
			pObj = new(pseg) synapseCascade(iCascadeSize,1, prng_r);
		}
		else//Starting Strength is Random
			pObj = new(pseg) synapseCascade(iCascadeSize,prng_r);

		vSyn.push_back(pObj);//Add to Target Vector

	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}

template <>
ICascadeSynapse* allocSynapseArray<synapseSingleFilterUnifiedWithDecay>(vector<synapseSingleFilterUnifiedWithDecay*> &vSyn,char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	 // const int iCascadeSize = 1; //Used So File Name Reflects the Filter Size
	 ICascadeSynapse::SYN_STRENGTH_STATE startStrength;

	  const bool bFixedStartState = false;

	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleFilterUnifiedWithDecay*,ptrdiff_t> pmem;

	  if (buffer ==0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleFilterUnifiedWithDecay>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseSingleFilterUnifiedWithDecay*)buffer;
		  pmem.second = iSynCount; //Assume Count is correct from Previous init
	  }

   if (!pmem.first && iSynCount > 0) ERREXIT(100,"Could not allocate memory");
   if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");

//	 long p = (long)pmem.first;
//   cout << p << " Next:" << (long)(pmem.first+1) <<  " diff:" << (long)(pmem.first+1)-p << endl;
//   cout << "Bytes Alloc:" << ((long)(pmem.first+1)-p)*pmem.second << endl;
//   cout << "Bytes Required:" << sizeof(synapseCascade)*iSynCount << endl;
   synapseSingleFilterUnifiedWithDecay* pObj;
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

			pObj = new(pseg) synapseSingleFilterUnifiedWithDecay(-g_FilterTh, g_FilterTh, g_FilterDecay, startStrength, g_MetaplasticitySampleSize, prng_r);
		}
		else//Starting State is Random not linked to Cascade States
			pObj = new(pseg) synapseSingleFilterUnifiedWithDecay(-g_FilterTh,g_FilterTh,g_FilterDecay,g_AllocRefraction,0.0);
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);

		vSyn.push_back(pObj);//Add to Target Vector
	}
	cout << "Stability Threshold : " << g_AllocRefraction << endl;
	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release

}


//Allocating Stochastic Updaters
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleUpdater>(vector<synapseSingleUpdater*> &vSyn,char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const int iCascadeSize = 1;
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
	  char* pseg; //Generic Pointer to allocated memory
	  pair<synapseSingleUpdater*,ptrdiff_t> pmem;

	  if (buffer == 0) //If null Then Allocate, Otherwise Create Objects over Previously Allocated Memory
	  {
		//Allocate memory Buffer to initialize cascadeSynapse objects
		 pmem = get_temporary_buffer<synapseSingleUpdater>(iSynCount);
	  }
	  else
	  {
		  pmem.first = (synapseSingleUpdater*)buffer;
		  pmem.second = iSynCount; //Assume Count is correct from Previous init
	  }

   if (!pmem.first  && iSynCount > 0) ERREXIT(100,"Could not allocate memory");
   if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");

//	 long p = (long)pmem.first;
//   cout << p << " Next:" << (long)(pmem.first+1) <<  " diff:" << (long)(pmem.first+1)-p << endl;
//   cout << "Bytes Alloc:" << ((long)(pmem.first+1)-p)*pmem.second << endl;
//   cout << "Bytes Required:" << sizeof(synapseCascade)*iSynCount << endl;
    synapseSingleUpdater* pObj;
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
			pObj = new(pseg) synapseSingleUpdater(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting Strength is Random
			pObj = new(pseg) synapseSingleUpdater(g_UpdaterQ,prng_r,g_MetaplasticitySampleSize);

			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);
		vSyn.push_back(pObj);//Add to Target Vector
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}


//Allocating Single U Filters: Sets Upper Lower Threshold And the allocation Threshold
//TODO:A Better implementation of Allocators would be to take the a sample object as parameter and use Copy constructors to initialise the population

//Allocating Stochastic Updaters
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleFilterUnifiedWithDecayReflecting>(vector<synapseSingleFilterUnifiedWithDecayReflecting*> &vSyn,char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
{
	  const int iCascadeSize = 1; //Used So File Name Reflects the Filter Size
	  const bool bFixedStartState = false;
	  ICascadeSynapse::SYN_STRENGTH_STATE startStrength;
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

   if (!pmem.first && iSynCount > 0) ERREXIT(100,"Could not allocate memory");
   if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");

//	 long p = (long)pmem.first;
//   cout << p << " Next:" << (long)(pmem.first+1) <<  " diff:" << (long)(pmem.first+1)-p << endl;
//   cout << "Bytes Alloc:" << ((long)(pmem.first+1)-p)*pmem.second << endl;
//   cout << "Bytes Required:" << sizeof(synapseCascade)*iSynCount << endl;
   synapseSingleFilterUnifiedWithDecayReflecting* pObj;
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
			pObj = new(pseg) synapseSingleFilterUnifiedWithDecayReflecting(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting State is Random not linked to Cascade States
			pObj = new(pseg) synapseSingleFilterUnifiedWithDecayReflecting(-g_FilterTh,g_FilterTh,g_FilterDecay);
			//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);

		vSyn.push_back(pObj);//Add to Target Vector
	}


	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}




//Allocating Stochastic Updaters DUAL FILTER
template <>
ICascadeSynapse* allocSynapseArray<synapseSingleFilterDual>(vector<synapseSingleFilterDual*> &vSyn,char*buffer,int iSynCount,int IndexOfTransitionProb, gsl_rng* prng_r, float StimRate)
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

   if (!pmem.first && iSynCount > 0) ERREXIT(100,"Could not allocate memory");
   if (pmem.second < iSynCount) ERREXIT(100,"Could not allocate all the required memory");

   synapseSingleFilterDual* pObj;
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
			pObj = new(pseg) synapseSingleFilterDual(iCascadeSize,IndexOfTransitionProb,startStrength, prng_r);
		}
		else//Starting Strength is Random
			//new(pseg) synapseSingleFilterDual(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);
			pObj = new(pseg) synapseSingleFilterDual(g_FilterTh,g_FilterTh,0.0); //Using Global VAriable(Bad Practice but needed)
		//new(pseg) synapseSingleFilterUnifiedWithDecay(iCascadeSize,IndexOfTransitionProb,ICascadeSynapse::SYN_STRENGTH_NOTSET, prng_r);

		vSyn.push_back(pObj);//Add to Target Vector
	}

	return (ICascadeSynapse*)pmem.first;
	//Use return_temporary_buffer(buffer) To release
}

/*
 *Allocate Single Object - According to the Required Constuctor
 */

