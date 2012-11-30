

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


/* Single U Filter cAMP Calculator used in SignalRep Function -
 * Template Specialization for calculating the cAMP signal as required to match
 * The numerical integration results in MATHEMATICA. Here I take the solution to u'(t)=dIn(u_m-u(t)) to
 * Inject the correct amount of cAMP during 1 timestep*/
/*
Remove Because it does not work for high decay rates or high repetitions r
template<>
double getcAMP<synapseSingleFilterUnifiedWithDecay>(double ts,double dcAMPLevel,double dCA, double dDA,double h_thres)
{
const int iPwr = g_iHillOrder;
double h_hill = pow(dCA,iPwr)/( pow(dCA,iPwr)+pow(h_thres,iPwr) );

assert(!isnan(dCA));

double dIn				= g_fInjectionGain*dDA*h_hill;

//Problem with Fast decays.. Better split decay step in two
dcAMPLevel				-= 0.5*(double)g_fcAMPDecay*dcAMPLevel; //Do Half Decay before;
//Here Used to save time not to calculate without inpu
if (dIn > 0.0) //The solution to u'[t] = s(um-u[t]) with u[0]=H where H is the value just before repetition
	dcAMPLevel					= 1.0 + exp(-dIn*ts)*(dcAMPLevel - 1.0); //This Works For Filter

dcAMPLevel				-= 0.5*(double)g_fcAMPDecay*dcAMPLevel; //Do Half Decay After;

return dcAMPLevel;
}
*/

/*
 * Returns -1 if the pattern has not been stored in the X vector Before
 * 			Otherwise returns an index relevant to the X vector So the pattern can be repeated on stored
 * 			Sets bAllocatePattern=true if pattern in tracked LIst is allocated
 * 			Sets bPatternIsNew - If pattern is tracked But this is the 1st encoding and thus Does not Exist in the the X vector
 *
 *	If a tracked Pattern say 4 is given at a timepoint when repetition occurs 4x1, then repetition takes precedence
 *
 * 1st: Check if time to repeat a pattern
 * 2nd: Find sequence relevant index in trackedlist to use as index in X[].
 *
 *
 *bPatternIsNew is set to true When a pattern needs to be encoded for the 1st time- If pattern index is also returned then this pattern needs to be saved in the given index.
 * Tracked patterns are regenerated on every trial and stored in X array so recall can be tested
 * bAllocatePattern - Is the DopAmine signal
 */
int selectPatternToEncode(unsigned long ts,uint uiNoOfPatternsStoredInTrial,bool& bPatternIsNew,bool& bAllocatePattern,t_patt_trackedtbl& vTrackedIndex,t_patt_reptbl& repetitionTable,t_inVal** X,uint uiPatCount,uint iSynCount,gsl_rng*& prng)
{
	t_patt_trackedtbl::iterator		itt;
	int Ret = -1; //Return -1 if random Pattern and not a tracked one.
	uint uiPatternToCheck; //Which Pattern Index Should be Checked if it has been Allocated

	t_patt_trackedtbl::iterator itAtTrackedPatt;
	t_patt_trackedtbl::iterator itAtEnd = vTrackedIndex.end();
	t_patt_trackedtbl::iterator itAtStart = vTrackedIndex.begin();

	//Search Repetition Table 1st to see if it is time to repeat a pattern
	t_patt_reptbl::iterator it = repetitionTable.find(ts);
	//Check That this Pattern is in the past and then return the absolute Pattern Index
	if (it != repetitionTable.end() && it->second < uiNoOfPatternsStoredInTrial)
	{
		Ret = it->second;// repetitionTable[j];//Found a pattern that needs to be repeated at time J Return the relative index of the tracked List -
		uiPatternToCheck = Ret; //This is a Repetition of Tracked pattern So check the Given Index To see if the Pattern is Has an Allocation Signal
	}else
			//if (Ret ==-1) //Not A repetition of A previous pattern? Check if current pattern is in the tracked list anyway
	{
		uiPatternToCheck = uiNoOfPatternsStoredInTrial; //Search TrackedList And Get Index to store Vector in X
	}

	////////FIND IF A PATTERN IS TRACKED And Allocated- Obtain Tracked list Index and Allocation Flag////////////
	int tempIndex = 0; //The sequence found in the tracked list is also the sequence that patterns are stored in the X vector -/Used As Relative Index In the X[] Pattern Vector
	for (itAtTrackedPatt = itAtStart; itAtTrackedPatt!= itAtEnd; ++itAtTrackedPatt) //Only Track 1st pattern
	{
		//cout <<  itAtTrackedPatt->first << " Check P:" << uiPatternToCheck << endl;
		if (itAtTrackedPatt->first == uiPatternToCheck) //Is the Pattern being Encoded In the tracked list?
		{//Check If it should be allocated

			///1st Repetition : itAtTrackedPatt->second == 1 <-Allocation Signal For this memory is On (DA)
			//Do not Set Flag On 1st repetition because h signal is zero on average - Thus any +ve contributions set cAMP
			//- But -ve do not decrease  && uiNoOfPatternsStoredInTrial > 0
			if (itAtTrackedPatt->second == 1)
			{
				//DA Signal Here - This Flag will be set to false again if the DA threshold is not set
			 	bAllocatePattern = true; //Allow post-Syn Depol To Se Alloc Signal
			}
			//Note: You can set a State Machine here (ex. (itAtTrackedPatt->second = 2) to change behaviour between repetitions

			if(Ret == -1) //If pattern in Track list but not repeated then Flag As a new Pattern to Be stored
			{
				bPatternIsNew = true; //Pattern is Tracked And It is the first time it will be encoded - Encode will store it
			}

			Ret = tempIndex; //Set to index to use in the X vector of Patterns
			break; //Set To the 1st Tracked Pattern found
		}
		tempIndex++; //Used As Relative Index In the X[] Pattern Vector
	}

return Ret;
}

//Next Timestep Predicts that we need a measurement just before pattern repetition so the sudden signal change can be seen.
unsigned long getNextTimestep(bool& bisRecallPeriod,bool& bisEncodingPeriod,double dTimeStepSize,unsigned long currentTimestep,unsigned long& LastRecallj,
								unsigned long& NextEncodingj, unsigned long lTotalTimesteps,
								PoissonSource*& PsMemEvent,t_patt_reptbl&  repetitionTable)
{
	const bool bDiscreteTime = (dTimeStepSize >= 1.0);
	bool bLastTimeStep = currentTimestep == (lTotalTimesteps-1); //
	unsigned long Nextj,peakSigTime; //The value Returned
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
			peakSigTime = it->first; //PeakSigTime Is stored in secs*ts while Encoding j is in timesteps by default
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
	tsPeriodOfRecall = 1.0;// getCurrentPeriodOfRecall(currentTimestep);

//Now check if a measurement Period occurs before the next encoding one
unsigned long nextRecallt = currentTimestep + tsPeriodOfRecall;//LastRecallj+tsPeriodOfRecall;

	if ((nextRecallt) > NextEncodingj) //If An encoding period precedes the next period of recall
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

	//Check we are crossing over Last timestep and make sure we recall last timepoint signal
	//If not at the currently lasttimestep already but we are exceeding it then move back and recall the last time point.
	if ((Nextj >= (lTotalTimesteps-1)) && (currentTimestep < (lTotalTimesteps-1)) && !bLastTimeStep)
	{
			Nextj = lTotalTimesteps -1; //Go to Last timestep and test recall
			bLastTimeStep = true;
	}

	if (bLastTimeStep)
		bisRecallPeriod = true;

	//uncomment following line to record signal At each timestep
	//bisRecallPeriod = true;

	return Nextj;
}




/*
 * Calculates the frequency of recall periods in number of timesteps until next recall using the point we are in the simulation now -
 * Currently assumes a log scale frequency, since results are plotted in logscale then only those points on the plot matter.
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



