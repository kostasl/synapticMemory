

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
 * Returns Index of repeated pattern - Non Repeated or tracked patterns are Saved on the same position in the X vector and are overwritten on every trial
 *
 * Tracked patterns are regenerated on every trial and stored in X array so recall can be tested
 * bisRecallPeriod = true; //FOR MFPT ALWAYS RECALL
 */
int selectPatternToEncode(unsigned long j,uint uiNoOfPatternsStoredInTrial,bool bUseRandomPatterns,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng)
{
	t_patt_trackedtbl::iterator		itt;
	int Ret = -1; //Return -1 if random Pattern and not a tracked one.

	if (!bUseRandomPatterns) //Load the next Index from File
		Ret = uiNoOfPatternsStoredInTrial;

	//Search Repetition Table 1st to see if it is time to repeat a pattern
	t_patt_reptbl::iterator it = repetitionTable.find(j);

	//If this is a time to repeat a pattern then return its index from the TrackedList So we can locate it in the previously encoded patterns
	//If Rep is not at 1st repetition (0) then return an index - Otherwise the pattern has not being stored before
	if (it != repetitionTable.end() && it->second < uiNoOfPatternsStoredInTrial)
	{
		int i=0;
		//Find the index of the repeated pattern in the Tracked List
		Ret = it->second;// repetitionTable[j];//Found a pattern that needs to be repeated at time J Return the relative index of the tracked List - THis Maps onto the Index In the X vector
		//cout << " Ret:" << Ret << endl;
	}

return Ret;
}

//Next Timestep Predicts that we need a measurement just before pattern repetition so the sudden signal change can be seen.
unsigned long getNextTimestep(bool& bisRecallPeriod,bool& bisEncodingPeriod,double dTimeStepSize,unsigned long currentTimestep,unsigned long& LastRecallj,unsigned long& NextEncodingj, PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable)
{
	const bool bDiscreteTime = (dTimeStepSize >= 1.0);
	unsigned long Nextj, peakSigTime; //The value Returned
	uint tsPeriodOfRecall;
	bool bNextEncodingIsSignalRepetition = false;

	if (bisEncodingPeriod) //If this was a recall period obtain the time of the nextEncoding Event
	{
		/* Change this to make Simulation in Discrete Time */
		if (bDiscreteTime)
			NextEncodingj = currentTimestep + dTimeStepSize;
		else
			NextEncodingj = PsMemEvent->getTimestepsUntilNextEvent()+currentTimestep; //This will be the next time step. If next timestep exceeds period of recording then record in this cycle at the regular interval


		//Search If a pattern is to be repeated before the next Encoding Event is to occur
		//Assumes the rep.time values are sorted in ascending Order
		t_patt_reptbl::iterator it;
		for (it = repetitionTable.begin();it!=repetitionTable.end();++it)
		{
			peakSigTime = it->first;
			//If rep.time is after this point in time AND before the next regular encoding event
			//Change Next encoding time for the repetition time.
			if ((currentTimestep < (peakSigTime-1)) && (NextEncodingj >= (peakSigTime-1))) //Obtain A measurement Just Before Repetition
			{
				NextEncodingj = (peakSigTime-1); //Change to encode on repetition
				bNextEncodingIsSignalRepetition = false;
				break; //Stop Searching
			}

			if ((currentTimestep < peakSigTime) && (NextEncodingj >= peakSigTime))
			{
				NextEncodingj = peakSigTime; //Change to encode on repetition
				bNextEncodingIsSignalRepetition = true;
				break; //Stop Searching
			}
		}
	}

	//Next Period Is recording only?-Or recording too
	bisEncodingPeriod = bisRecallPeriod =  false; //Reset Flags for next Timestep Cycle
	tsPeriodOfRecall = 1;//getCurrentPeriodOfRecall(currentTimestep);

//Now check if a measurement Period occurs before the next encoding one
unsigned long nextRecallt = currentTimestep + tsPeriodOfRecall;//LastRecallj+tsPeriodOfRecall;

	if ((nextRecallt) > NextEncodingj)
	{
		bisEncodingPeriod = true;
		Nextj = NextEncodingj;

		// Next encoding is a signal repetition event or Just Before A rep Event Do Measurement
		if (bNextEncodingIsSignalRepetition || ((peakSigTime-1) == NextEncodingj))
			bisRecallPeriod = true;


	}else{
		Nextj = nextRecallt;
		bisRecallPeriod = true; //Do measurement And save it to index j=LastRecallj+tsPeriodOfRecall
		bisEncodingPeriod = false;
	}

	//Do Encoding And Measurement if Encoding period coincides with Recall Period
	if ((nextRecallt == NextEncodingj))
	{
		bisRecallPeriod = true;
		bisEncodingPeriod = true;
	}

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



