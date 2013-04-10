#include  "../stdafx.h"
#include "../InputVectorHandling.h"
#include "PerceptronUnifiedWDecayFilter.h"


/*
 * Spaggeti Version Of A Double threshold Filter with Decay - Running A Perceptron Test
 * This code loads the vectors once and then randomizes over the starting point
 */

void doNoddyPerceptronTest(float _fLowSignalThreshold)
{

	const bool bRandomizeOverSequenceOfPatterns = false; //Determines if Patterns Should Be shuffled - Or just presented at the file Sequence
	const unsigned int ciTrials					= 100000;
	const unsigned int ciNumberOfInitPatterns 	= 1;
	unsigned int iTrackedMemIndex 				= 0; //Will Be set to the Index of a random Pattern that occured After  ciNumberOfInitPatterns
	unsigned int ciPatCount 					= 1000;
	const unsigned ciInputSize			 		= 1000; //The size of the input vector = 1+Afferent Synapses
	//const int iRandomVectorSet					= 1000;//The Set Of Hadamard Vectors Used
	const int ciStartIndex						=  0; //The cascade index to start capacity tests from
	const int ciMaxCscIndex						= 14; // Biggest cascade To be tested
	const float cfLowSignalThreshold			= _fLowSignalThreshold; //The point to stop adding More Memories
	///Patterns Buffer
	t_inVal* X[ciPatCount];
	float* W; //Weight Vector
	double C[ciMaxCscIndex+1];//The Capacity Found Per Cascade Size

	char fname[250];
	sprintf(fname,"HD%d.dat",ciInputSize);

///Synapse state Variables
	int iFRun[ciInputSize]; //Running Sum for each Synapse
	int iCascadeIndex[ciInputSize]; //CascadeIndex of Each Synapse
	int iFLHThres[ciInputSize][2]; //The current L and H Threshold For each Synapse
	int iStrength[ciInputSize];
	double dDecayRate[ciInputSize];///Decay rate Of Each Synapse
	gsl_rng* mprng = getRandGeneratorInstance(true);

	//Statistics
	double dlastTrialSignal,dAvgSignal;
	unsigned int iRecallHits = 0;

	///Init Memory
	W = new float[ciInputSize];
	memset((void*)W,0,sizeof(float)*ciInputSize);

	for (uint i=0;i<ciPatCount;i++)
	{
		X[i] = new t_inVal[ciInputSize];
		//Zero Out To detect any initialization errors
		memset(X[i],0,sizeof(t_inVal)*ciInputSize);

	}
	//Read Test Vectors

	//sprintf(fname,"HD%d-N300-Set%d.dat",ciInputSize,iRandomVectorSet);
	//sprintf(fname,"randXVectorN_%dPat_300_Corr_0.50.bin",ciInputSize);

	readTestVectorsFromFile(fname,(t_inVal**)&X,ciInputSize,ciPatCount,0); //Load At 0 All Required PAtterns

	//Open Output file
	ofstream* ofile;
	char fparams[150];
	strcpy(fname,PERCEPTRON_OUTPUT_DIRECTORY);
	sprintf(fparams,"DblFilt_RevNODDYTrial%dRST0HV%d_I%d_SThres%1.2f.N%d-%ddat",ciTrials,ciInputSize,ciNumberOfInitPatterns,cfLowSignalThreshold,ciStartIndex+1,ciMaxCscIndex+1);
	strcat(fname,fparams);
	ofile = new ofstream(fname, ios::out ); //OPEN OUTPUT FILE
	cout << "OutFile: " << fname << endl;
	*ofile << "#unified double threshold filter with decay - Noddy code Revision" << endl;

	//Loop For all CascadeSizes
for(int itestIndex = ciStartIndex;itestIndex<=ciMaxCscIndex;itestIndex++)
{

	//Loop For all trials
	C[itestIndex] = 0.0;
	dAvgSignal = 0.0;
	for (unsigned int t=0;t<ciTrials;t++) ///TRIAL LOOP
	{
		if (t%1000 ==0) //Time showing
			cout << "*";
		cout.flush();

		iRecallHits = 0;
		//<INIT>Initialise Each Synapse CascadeIndex & Low and High threshold
		for(unsigned int i=0;i<ciInputSize;i++)
		{

			double r = gsl_rng_uniform(mprng);
			iStrength[i] = (r < 0.5)?-1:1;//Set Random Strength
			//iCascadeIndex[i] = round(itestIndex*gsl_rng_uniform(mprng));//Set Random Cascade INdex
			iCascadeIndex[i] = floor(itestIndex*gsl_rng_uniform(mprng)*0.99);//Set Random Cascade INdex
			//Set Thresholds for the Index and strength set
			set_thresholds(iCascadeIndex[i],iStrength[i],iFLHThres[i],dDecayRate[i],itestIndex,mprng);

			//Initialize Filter State According to PDF
			iFRun[i] = set_injection_filter_state(iFRun[i], iCascadeIndex[i],iStrength[i],iFLHThres[i],dDecayRate[i],itestIndex,mprng);

		}//END OF Init Each Synapse

		//<LEARNING> Store Each Pattern Until Signal Drops below low threshold
		for(unsigned int k = 0;k<ciPatCount;k++)
		{
			int PatIndex;
			if (bRandomizeOverSequenceOfPatterns)
				PatIndex = floor(gsl_rng_uniform(mprng)*0.99*ciPatCount); ///Randomize sequence of patterns
			else
				PatIndex = k;

			for (uint i=0;i<(ciInputSize-1);i++) //Go over every synapse
			{ //Assume Desired output is  X[j][_uiSynCount-1] last vector value
				//Do POT AND DEP
				int iSig = 0;
				if (X[PatIndex][i]*X[PatIndex][ciInputSize-1] > 0) //Correlate WIth desired OUtput to obtain induction Signal
					iSig = +1;//POT Signal
				else
					iSig = -1;
				//Pass To function That modifies and checks Filter (//Add To running Sum.
				//Check if Threshold Reached //Check if terminal State 	// if not Change Index To new State)
				addSampleToFilter(iFRun[i],iSig, iCascadeIndex[i], iStrength[i],iFLHThres[i],dDecayRate[i],itestIndex,mprng);

				W[i] = iStrength[i]; //Save new Strength To weight Matrix
			}//END OF For Each Afferent
			//<TEST RECALL> Saves Signal Into dlastTrialSignal variable
			if (k >= ciNumberOfInitPatterns) //Check if required number of init patterns have been stored
			{
				if (k==ciNumberOfInitPatterns) //Once Reached Store The Tracked Pattern Index
					iTrackedMemIndex = PatIndex; //Save this pattern as the tracked index -

				int ret = testPerceptronRecallOfStoredPattern( dlastTrialSignal, W,(t_inVal**)X, ciInputSize, iTrackedMemIndex);
				assert (ret <2);
				iRecallHits +=  ret;
				//cout << "t:" << t << " Pat:" << (k-ciTrackedMemIndex)  << " S:"<< dlastTrialSignal << endl;
				if (dlastTrialSignal < cfLowSignalThreshold)
					{
						break; 	//If signal Is below Threshold Break
						cout << "LSIG @" << k << endl;
					}
			}
				//iRecallHits = k;

//			float h=0.0;
//			int Ret = 0;
//			for (uint i=0;i<(ciInputSize-1);i++)
//			{
//				//cout << w[i] << " " << " " << X[_iStoredPatIndex][i] << endl;
//				h+=W[i]*X[ciTrackedMemIndex][i];
//			}
//
//			//MEasure The normalized PostSynaptic Respose As Signal
//			dlastTrialSignal = X[ciTrackedMemIndex][ciInputSize-1]*h/(ciInputSize-1);
//
//			float iNeuronOut = (h>0)?1:-1; //Neuron Classifier OUtput
//			//cout << " Signal: " << _Signal << endl;
//			if (X[ciTrackedMemIndex][ciInputSize-1] == iNeuronOut)
//			{
//				Ret = 1; //Return 1 To indicate Successful classification of input
//			}
			//Check if Max Capacity Reached - Low Signal Condition
		}//END OF For each Pattern

		C[itestIndex] += iRecallHits;//Save the Number of Patterns Stored On the lAst Trial
		dAvgSignal +=dlastTrialSignal;

	}//END OF TRIALS LOOP
	C[itestIndex] = C[itestIndex]/ciTrials; //Calc Avg Capacity
	//Output Capacity for this cascade Size
	cout << "Size: " << (itestIndex+1) << " Capacity:" << C[itestIndex] << endl;
	*ofile << (itestIndex+1) << "\t" << C[itestIndex] << endl;

}	//Loop For all CascadeSizes

//<Clean Up>
ofile->close();
delete ofile;

delete [] W;
for (uint i=0;i<ciPatCount;i++)
	{
		delete [] X[i];


	}

gsl_rng_free(mprng);

}


/*
 * @_Signal Reference to variable to Return the signal for tracked pattern
 * @_iStoredPatIndex The index of the pattern we are tracking
 * @Returns : 1 For successful recall , 0 for failure
 */
inline int testPerceptronRecallOfStoredPattern(double& _Signal,float* W,int** X,uint _uiSynCount,uint _iStoredPatIndex)
{
	//Test Recall of _iStoredPatIndex
	//cout << "Recall index: " << _iStoredPatIndex << " Output should be :" << X[_iStoredPatIndex][_uiSynCount-1];
	float h=0.0;
	int Ret = 0;
	for (unsigned int i=0;i<(_uiSynCount-1);i++)
	{
		//cout << w[i] << " " << " " << X[_iStoredPatIndex][i] << endl;
		h+=W[i]*X[_iStoredPatIndex][i];
	}

	//MEasure The normalized PostSynaptic Respose As Signal
	_Signal = X[_iStoredPatIndex][_uiSynCount-1]*h/(X[0][0]*X[0][0]*(_uiSynCount-1));
	float iNeuronOut = (h>0)?1:-1; //Neuron Classifier OUtput
	//cout << " Signal: " << _Signal << endl;
	if (X[_iStoredPatIndex][_uiSynCount-1] == iNeuronOut)
	{
		Ret = 1; //Return 1 To indicate Successful classification of input
	}

	return Ret;
}


/*
 * Sets the thresholds by Clearly listing each case
 * It then Adds the Random switching of thresholds which cancels
 * The previous Ordering of p&q thresholds BUT it is still explicitly
 * coded to clearly demonstrate the ideas that have been superimposed
 */
inline void set_thresholds(int _iCascadeIndex,int _iStrength,int* _piFLHThres,double& _dDecayRate,int _iTerminalIndex,gsl_rng* _rng)
{
	if (_iStrength == 1) //STRONG SYNAPSE
	{	//Check if Synapse is in terminal State
		if (_iCascadeIndex < _iTerminalIndex)
		{
			_piFLHThres[0] 	= -miThreshold_r100[_iCascadeIndex][1];//Low Thres q transitions
			_piFLHThres[1]	= miThreshold_r100[_iCascadeIndex][0];//High Thres p transitions
			_dDecayRate 		= mdDecay_r100[_iCascadeIndex][0];
		}
		else //Terminal States
		{
			_piFLHThres[0] 	= -miThresholdTerminal_r100[_iCascadeIndex][1];//Low Thres q transitions
			_piFLHThres[1]	=  miThresholdTerminal_r100[_iCascadeIndex][0];//High Thres p transitions
			_dDecayRate 	=  mdDecayTerminal_r100[_iCascadeIndex][0];
		}
//WEAK SYNAPSE
	}else{
		if (_iCascadeIndex < _iTerminalIndex)
		{
			_piFLHThres[0] 	= -miThreshold_r100[_iCascadeIndex][0];//Low Thres q transitions
			_piFLHThres[1]	= miThreshold_r100[_iCascadeIndex][1];//High Thres p transitions
			_dDecayRate 	= mdDecay_r100[_iCascadeIndex][0];
		}
		else //Terminal States
		{
			_piFLHThres[0] 	= -miThresholdTerminal_r100[_iCascadeIndex][0];//Low Thres q transitions
			_piFLHThres[1]	= miThresholdTerminal_r100[_iCascadeIndex][1];//High Thres p transitions
			_dDecayRate 	= mdDecayTerminal_r100[_iCascadeIndex][0];
		}

	} //ENDOF IF STRONG


	//Implement Symmetric Assymetric Idea By Randomly Switching Between thres
	if (_iCascadeIndex != _iTerminalIndex)
	{
		//g_rng_r = getRandGeneratorInstance(false);
		double r = gsl_rng_uniform(_rng);
		if (r < 0.5)
		{	//SWAP
			int dummy4 = _piFLHThres[0];
			_piFLHThres[0] = -_piFLHThres[1];
			_piFLHThres[1] = -dummy4;
		}
	}

}


inline int set_injection_filter_state(int& _iFilterState,int _iCascadeIndex,
									  int _iStrength,int* _piFLHThres,
									  double& _dDecayRate,int _iTerminalIndex,
									  gsl_rng* _rng)
{
	const int ciMaxInternalStates = 7; //The Sum of Filter states -> add Low + High threshold +1
	double r = gsl_rng_uniform(_rng);
	double p = 0.0; //Aux Var. CDF
	int i; //Aux Var. Iterator
	int countUpOrDown = 0; //Aux Variable to assist in counting Up/Down the Filter State
	const double (*dPDFUsed)[15][ciMaxInternalStates]; //Pointer to PDF. If at terminal State then Use The reflecting Boundary PDF

	_iFilterState = 0;
	//return;
	if (((_iTerminalIndex) < 1)) return 0; //Not For 1,2 thresholds filter

	if (_iCascadeIndex < _iTerminalIndex )
		dPDFUsed = (&mdPDF_r100);
	else
		dPDFUsed = (&mdPDFTerminal_r100);


	if ((_iCascadeIndex) == 1) //Small Thresholds Use Uniform (Assymetric Threshold)
	{
		//Assign Random Start State --Add +1 To Include State 0 bUT LIMIT TO BELOW THRESHOLD
		_iFilterState = round(r*(-(_piFLHThres[0]+1) + (_piFLHThres[1]-1))) +(_piFLHThres[0]+1);
	}
	else //Use init Distribution
	{

		for ( i = 0; i< ciMaxInternalStates;i++)
		{
			if (_iStrength == 1) //STRONG SYNAPSE
				countUpOrDown = i; // If STRONG Synapse - The PDF is directed Correctly - Low Bound Is Reflecting
			else
				countUpOrDown = ciMaxInternalStates - i-1; //-1 cause We require an index

			p += (*dPDFUsed)[_iCascadeIndex][i]; //Accumulate the Pdf
			if (p > r){
				//Found the spot since r was just exceeded
				_iFilterState = countUpOrDown-3; //Remove Offset so i=0 becomes state -3 floor(ciMaxInternalStates/2)
				break;
			}
		}
	}

	//Note That miRVal can be = to Threshold If the Threshold Is a Holding and Not Absorbing
if ((_iFilterState < _piFLHThres[0]) || (_iFilterState > _piFLHThres[1]))
{
	cout << "RVal:" << _iFilterState <<" LTh:" << _piFLHThres[0] << " HTh:" << _piFLHThres[1] << endl;
	assert((_iFilterState > _piFLHThres[0]) && ((_iFilterState < _piFLHThres[1])));
}
return _iFilterState;
}

inline void switchResetSynapse(int& _iFilterState,int& _iCascadeIndex,
							  int& _iStrength,int* _piFLHThres,
							  double& _dDecayRate,int _iTerminalIndex,
							  gsl_rng* _rng)
{
	if (_iStrength == 1){
		_iStrength = -1;
	}else{
		_iStrength = 1;
		}

		//Set To Top Of opposite Cascade
		_iCascadeIndex = 0;

		 set_thresholds(_iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
		 _iFilterState = 0; //Re-Inject to Zero In new Cascade State
}

/*
 * Function Handling the change of a filters internal sum
 * Calls the decay, Check of thresholds and Change CascadeState if required
 * Function Has been designed to work for Rate = 1.0 So Time since last induction is always 1
 */
int addSampleToFilter(int& _iFilterState,int _InductionSignal,int& _iCascadeIndex,
		  int& _iStrength,int* _piFLHThres,
		  double& _dDecayRate,int _iTerminalIndex,
		  gsl_rng* _rng)
{
	int Ret = 0; //Aux Var. Return Value of Function
	int iTimeSinceLastIndunction = 1;
	//Do Decay
	_iFilterState += drawDecaySteps(_iFilterState, _dDecayRate, iTimeSinceLastIndunction, _rng);
	//Add New Signal
	_iFilterState +=_InductionSignal;
	//Check Thresholds  Condition
	if (_iFilterState <= _piFLHThres[0]) //LTHres Reached?
	{
		//If Weak Synapse Then p Transition  -
		if (_iStrength == -1)
		{
			if (_iCascadeIndex < _iTerminalIndex){
				_iCascadeIndex++; //Do A p Transitions
				 Ret = 1; //Signal METAplastic Transition Occurred
				 //Set New Thresholds
				 set_thresholds(_iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
				 _iFilterState = 0;
				 //RE-INJECTION
				 //set_injection_filter_state(_iFilterState, _iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
			}else{//Terminal State //No p Transitions
				Ret = 0;//No Transition
				_iFilterState = _piFLHThres[0]; //Holding Barrier
			}
		}else{///Strong Synapse + Low Thresh Condition = Switch Cascade Strength
			//Change to STRONG And Reset Cascade Index - Reset Runnning Sum
			switchResetSynapse(_iFilterState,_iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
			Ret = -1; //Plastic Transition Occured
			}
	}else //Otherwise High Threshold Has been reached
		if (_iFilterState >= _piFLHThres[1])
		{ //p THRESHOLD REACHED?
			if (_iStrength == 1)
			{///Strong Synapse  So P transition
				if (_iCascadeIndex < _iTerminalIndex){
					_iCascadeIndex++; //Do A p Transitions
					 set_thresholds(_iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
					 _iFilterState = 0;
					 //set_injection_filter_state(_iFilterState, _iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
					 Ret = 1; //METAplastic Transition Occurred
				}else{//Terminal State //No p Transitions
					Ret = 0;//No Transition
					_iFilterState = _piFLHThres[1]; //Holding Barrier
				}
			}else{ //Q Threshold Reached - WEAK SYNAPSE -> SWITCH TO STRONG
				//Change to weak And Reset Cascade Index Reset Runnning Sum
				switchResetSynapse(_iFilterState,_iCascadeIndex,_iStrength,_piFLHThres,_dDecayRate,_iTerminalIndex,_rng);
				Ret = -1; //Plastic Transition Occurred
				}
		}

	return Ret;
}

/*
 * Calculate The number of Decay Steps that have occured
 * Since the last induction signal
 */
int drawDecaySteps(int _miRFilterValue,double _dDecayRate,int _iTimeSinceLastInduction,gsl_rng* _rng)
{

	double p;
	unsigned int rDecaySteps = 0;


	if (_miRFilterValue == 0) return 0;

	if (_miRFilterValue > 0)
	{
		p = 1.0-exp(-_dDecayRate*_iTimeSinceLastInduction);
		rDecaySteps -= gsl_ran_binomial(_rng,p,_miRFilterValue);
		//miRFilterValue -= rDecaySteps; //Subtract to move to zero
	}
	else //Increment Sum - -Ve side is decaying back to zero
	{
		p = 1.0-exp(-_dDecayRate*_iTimeSinceLastInduction);
		rDecaySteps = gsl_ran_binomial(_rng,p,-_miRFilterValue);
		//miRFilterValue += rDecaySteps; //Add To move to zero
	}


	return rDecaySteps;
}
