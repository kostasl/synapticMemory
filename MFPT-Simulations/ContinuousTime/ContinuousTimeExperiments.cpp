

#include "ContinuousTimeExperiments.h"


/*
 * This Function Generates a random +/-1 response using information from the last call.
 * It is used to generate the random stimuli induction train.
 *
 * dC is correlation //Prob[I2 = ±I1] = (1 ± C2)/2,
 * When using a correlation This function makes use of a static variable to remember the last induction step.
 *
 * This way it may generate driftless induction stimuli by Setting the correlation term to dC = -1 . so -1, 1,-1,1 sequences can be generated
 *	NOT MultiThread Safe - Uses static Variable
 * dC : Correlation term 0 gives random sequence, -1 gives  toggle of induction stimuli
 * Returns : 1 for a POT step -1 for DEP
 *
 * NOTE: If the function is inline then AllocationExperiment.cpp cannot find it!
 */
int makeInductionStep(float& mdRate, gsl_rng*& prng_r, double& dC)
{
	static int lastInductionStep = 1;
	int newInductionStep = 0;
	double fp,r;//Random Num for rates of f_+, f_- POT and DEP rates
	int iSn 	= 0; //Signal Variable

	#ifdef USE_NATIVE_RAND
		fp = rand()/(double)RAND_MAX;
		r = rand()/(double)RAND_MAX;
	#else
		fp 		= gsl_rng_uniform(prng_r);
		r = (mdRate < 1.0)?gsl_rng_uniform(prng_r):0.0;
	#endif

	///3 Possibilities:
	//Either Both DEP & POT Occur and we do DEP THEN POT
	//or both occur and we call in the other order probabilistically using fp < fm
	///Or One of the two occurs and then call sequence does not matter
	if (r < mdRate)
	{
			if (fp < (1.0 + dC)/2.0)
				newInductionStep = lastInductionStep;
			else
				newInductionStep = -1*lastInductionStep;

			if (newInductionStep == 1 )
			{
				iSn = 1;
			}
			else
			{
				iSn = -1;
			}

			lastInductionStep = newInductionStep;
			return iSn;
	}
	else //No stimulus induction Call NOP in case decay is required
	{
		iSn = 0;
	}

return iSn;
}
/*
 * Returns -1 If the next Pattern is to A random one i.e Not a repeated Pattern
 * Checks if it is time to repeat a pattern and sets the index to the pattern to be re-encoded
 * Returns Index of repeated pattern since all past patterns are stored In X and are overwritten on every trial
 *
 * Tracked patterns are regenerated on every trial and stored in X array so recall can be tested
 * bisRecallPeriod = true; //FOR MFPT ALWAYS RECALL
 */
int selectPatternToEncode(unsigned long j,uint uiNoOfPatternsStoredInTrial,bool bUseRandomPatterns,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng)
{
	t_patt_trackedtbl::iterator		itt;
	int Ret = -1; //Return -1 if random Pattern and not a tracked one.

	//Search Repetition Table 1st
	t_patt_reptbl::iterator it = repetitionTable.find(j);

	//If this is a time to repeat a pattern then return its index from the TrackedList So we can locate it in the previously encoded patterns
	if (it != repetitionTable.end())
	{
		int i=0;
		//Find the index of the repeated pattern in the Tracked List
		for (itt = vTrackedIndex.begin();itt != vTrackedIndex.end(); ++itt)
		{
			if (repetitionTable[j] == itt->first) //Locate Index Of this Pattern relative to the Tracked List
			{
				Ret = i;//Found a pattern that needs to be repeated at time J Return the relative index of the tracked List
				break;
			}
		}
		i++;
	}

return Ret;
}


double getNextTimestep(bool& bisEncodingPeriod,double dTimeStepSize,double currentTimestep,double& LastRecallj,double& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable)
{
	const bool bDiscreteTime = (dTimeStepSize >= 1.00);
	double Nextj,peakSigTime; //The value Returned

	/* Change this to make Simulation in Discrete Time */
	if (bDiscreteTime)
	{
		Nextj = currentTimestep + 1.0;
	}
	else
		Nextj = PsMemEvent->getTimeUntilNextEvent() + currentTimestep; //This will be the next time step. If next timestep exceeds period of recording then record in this cycle at the regular interval

	//Search If a pattern is to be repeated before the next Encoding Event is to occur
	//Assumes the rep.time values are sorted in ascending Order
	t_patt_reptbl::iterator it;
	for (it = repetitionTable.begin();it!=repetitionTable.end();++it)
	{
		peakSigTime = it->first;
		//If rep.time is after this point in time AND before the next regular encoding event
		//Change Next encoding time for the repetition time.
		if ((currentTimestep < peakSigTime) && (Nextj >= peakSigTime))
		{
			Nextj = peakSigTime; //Change to encode on repetition
			break; //Stop Searching
		}
	}
	//Every Step is an Encoding/Recall Period for MFPT
	bisEncodingPeriod  =  true;

	return Nextj;
}



/*
 * Calculates the frequency of recall periods in number of timesteps until next recall using the point we are in the simulation now -
 * Currently assumes a log scale frequency, because results are plotted in logscale then only those points on the plot matter.
 */
uint getCurrentPeriodOfRecall(unsigned long j)
{
	uint period =  pow(10,floor(log10(j)));
	period = (period == 0)?1:period;

	period = (period  > 1)?period/2:period; //Double the frequency if we are at times above 10 timesteps - Smoother curves )We cant do this below 10 timesteps offcourse

	uint n = j/period; ///Align to Period Intervals
	period = (n+1)*period - j;


	return period;
}



/*
 * Same as above but makes sure a recall occurs at the peak of a filters signal
 */

uint getCurrentPeriodOfRecall(unsigned long j,unsigned long peakSigTime)
{
	uint period =  pow(10,floor(log10(j+1)));
	period = (period  > 1)?period /2:period; //Double the frequency if we are at times above 10 timesteps - Smoother curves )We cant do this below 10 timesteps offcourse

	if ((j < peakSigTime) && (j+period > peakSigTime)) //Check If Signal Peak Is missed And Correct
		period = peakSigTime-j;

	return period;
}

