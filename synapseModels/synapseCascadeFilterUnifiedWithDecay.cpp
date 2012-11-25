/*
 * synapseFilterUnifiedWithDecay.cpp
 *
 *  Created on: 8 Mar 2011
 *      Author: kostasl
 *
 */

#include "synapseCascadeFilterUnifiedWithDecay.h"

///RE-INjection At ZERO - PDFs
const double synapseCascadeFilterUnifiedWithDecay::mdPDF_r100[20][7] = {
		 {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, //n=1
		 {0.0, 0.0 ,0.5, 0.5, 0.0, 0.0, 0.0}, //n=2 - Code ignores these values - Uniform Distrib used
		 {0.0, 0.0 ,0.256193, 0.480363, 0.256193, 0.0, 0.0}, //n=3
		 {0.0, 0.0, 0.368089, 0.262159, 0.368089, 0.0, 0.0}, //n=4 i0
	     {0.06250000000000000000, 0.12500000000000000000, 0.18750000000000000000, 0.25000000000000000000, 0.187500000000000000, 0.125000000000000000, 0.06250000000000000}, //n=5 i0
	     {0.05173184305966837970, 0.12056985594910431200, 0.20712844504212299500, 0.24113971189820862500, 0.207128445042122995, 0.120569855949104312, 0.05173184305966838},//n=6 0i
	     {0.04178090320605485768, 0.11991970205845851808, 0.22669658910714549203, 0.22320561125668226443, 0.226696589107145492, 0.119919702058458518, 0.041780903206054857}, //n=7
	     {0.03033980324165782516, 0.11375865233380612532, 0.24617387244966487268, 0.21945534394974235366, 0.246173872449664872, 0.113758652333806125, 0.030339803241657825},//n=8
	     {0.02207650931777412286, 0.10759264781989243877, 0.26477172580743867559, 0.21111823410978952554, 0.264771725807438676, 0.107592647819892438, 0.022076509317774123},//n=9
	     {0.01551202508390696529, 0.10039869969193056834, 0.28466441487178715927, 0.19884972070475061420, 0.284664414871787159, 0.100398699691930568, 0.015512025083906965},//n=10
	     {0.01107689622611658533, 0.09328862481463681768, 0.30283461593474448079, 0.18559972604900423240, 0.302834615934744481, 0.093288624814636817, 0.011076896226116585},//n=11
	     {0.00794007027933235005, 0.08620047067753524242, 0.31990807418748821973, 0.17190276971128837562, 0.319908074187488220, 0.086200470677535242, 0.007940070279332350},//n=12
	     {0.00561104924088536600, 0.07892128935747792290, 0.33667028676943866612, 0.15759474926439608994, 0.336670286769438666, 0.078921289357477923, 0.005611049240885366}, //n=13
	     {0.00390308244102358892, 0.07158541208445645002, 0.35298599021679465694, 0.14305103051545060824, 0.352985990216794657, 0.071585412084456450, 0.00390308244102359}, //n=14
	     {0.00279223634637394967, 0.06516152558458873724, 0.36691538387921754998, 0.13026170837963952624, 0.366915383879217550, 0.065161525584588737, 0.00279223634637395}, //n=15
	     {0.001968,					0.058866,				0.380300,				0.117732,				0.380300,			0.058866,				0.001968},
         {0.001390,0.053069,0.392473,0.106137,0.392473,0.053069,0.001390},
	    {0.000982,0.047731,0.403555,0.095463,0.403555,0.047731,0.000982},
        {0.000694,0.042851,0.413604,0.085702,0.413604,0.042851,0.000694}
 };



//With the View that the lower theshold is the Reflecting one
//		//Reflecting								//Absorbing
const double synapseCascadeFilterUnifiedWithDecay::mdPDFTerminal_r100[][7] = {
						{0.0000, 0.0000, 0.0000,  1.0000, 	0.0000, 	0.0000,	0.0000}, //n=1 - But never Terminal
						{0.0000, 0.0000 ,0.0000,  1.0000,   0.0000,		0.0000, 0.0000}, //n=2 - Sim i0
						{0.0000, 0.0000 ,0.5000,  0.5000,	0.0000,		0.0000,	0.0000}, //n=3 //Thres 2,1
						{0.250000,0.250000,0.250000,0.250000,0.000000,0.000000, 0.0000}, //n=4 i0 //Thres 4,1

						{0.000000,0.090158,0.353180,0.209320,0.347342,0.000000,0.000000}, //n=5 i0 Thres 3,2
						{0.000000,0.036118,0.433009,0.098241,0.432632,0.000000,0.000000},//n=6 0
						{0.000000,0.016705,0.467707,0.047918,0.467670,0.000000,0.000000}, //n=7 Thres 3,2
						//Thres 4,3
						{0.010107,0.091144,0.313765,0.185852,0.311854,0.087278,0.000000},//n=8
						{0.004560,0.074516,0.348584,0.151336,0.348010,0.072994,0.000000},//n=9
						{0.002148,0.060309,0.378044,0.121938,0.377867,0.059695,0.000000},//n=10 Simulation: 0		0.03592		0.4216		0.1172		0.3365		0.08507		0.003781
						{0.001036,0.048478,0.402318,0.097679,0.402263,0.048228,0.000000},//n=11
						{0.000507,0.038794,0.422026,0.077973,0.422008,0.038692,0.000000},//n=12
						{0.000250,0.030954,0.437888,0.062111,0.437883,0.030913,0.000000},//n=13
						{0.000124,0.024653,0.450589,0.049411,0.450588,0.024636,0.000000},//n=14
						{0.000062,0.019610,0.460726,0.039274,0.460725,0.019603,0.000000},//n=15
						{0.000031,0.015586,0.468800,0.031201,0.468799,0.015584,0.000000},
						{0.000015,0.012382,0.475222,0.024778,0.475222,0.012381,0.000000},
						{0.000008,0.009833,0.480326,0.019674,0.480326,0.009833,0.000000},
						{0.000004,0.007808,0.484381,0.015619,0.484381,0.007807,0.000000},
						{0.000002,0.006198,0.487602,0.012398,0.487602,0.006198,0.000000}
						/* Some of my old injection values - Changed them to match T.E
						{0.00464808937484574594, 0.0748993804197276717, 0.347784022265603529, 	  0.152133069779678630,   0.347192837188120,   0.073342600972023,	0.0},//n=9
						{0.0023538540965529977,  0.061932317432758907,   0.37469405428023903, 	  0.125284389548946356,   0.374489542593762,   0.0612458420477402, 	0.0},//n=10
						{0.00141733266077160678, 0.053316276325584176,   0.39241757230475938, 	  0.107574555123429597,   0.39232617721386379, 0.05294808637159145,	0.0},//n=11
						{0.00051651803985112326, 0.039029305280931016,   0.42154776460888664, 	  0.078451179945527004,   0.42152979889170465, 0.03892543323309956,	0.0},//n=12
						{0.00025185174365122620, 0.031032426527322594,   0.43773065851853435,	  0.062269089548515900,   0.43772508641419692, 0.03099088724777900,	0.0}, //n=13
						{0.00012365502526760665, 0.024638552656549110,   0.45061763483271835,    0.049382304279847994,   0.45061589936851583, 0.02462195383710111,	0.0}, //n=14
						{0.00006935118818222132, 0.020388589671107832,   0.45916198988909251,    0.040837990932378369,   0.45916132011373434, 0.0203807582055047,	0.0} //n=15*/
				 }; //Changed them on 31/5 to Match T.E header File


//alpha = 0.5
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
 const int8_t synapseCascadeFilterUnifiedWithDecay::miThreshold_r100[][2] ={
		 	 	 	 	 	 	 	 	 	 	{1,1},{1,2},
												{2,2},{2,2},
												{4,4},{4,4}, //n=6
												{4,4},{4,4},
												{4,4},{4,4}, //n=10
												{4,4},{4,4}, //n=12
												{4,4},{4,4}, //n=14
												{4,4},{4,4},
												{4,4},{4,4},
												{4,4},{4,4}}; ////n=20

 //TODO: UNIFIED FILTER WITH DECAY: NEED TO MAKE THE REFLECTING VALUES FOR R=0.1 and r=0.01
 //First Val is Reflecting 2nd is Absorbing -- The escape time of the reflecting boundary is the same as the n-1 Cascade state
 const int8_t synapseCascadeFilterUnifiedWithDecay::miThresholdTerminal_r100[][2] ={
		 	 	 	 	 	 	 	 	 	 	{1,1},{1,1},
 												{2,1},{4,1}, //n=3
 												{3,2},{3,2}, //n=5
 												{3,2},{4,3}, //n=7
 												{4,3},{4,3}, //n=9
 												{4,3},{4,3}, //n=11
 												{4,3},{4,3}, //n=13
 												{4,3},{4,3},
 												{4,3},{4,3},
 												{4,3},{4,3} }; ////n=20


 const double synapseCascadeFilterUnifiedWithDecay::mdDecayTerminal_r100[][2] ={ {0.0,0.0},{0.0,0.0},//n=2
 												{0.0,0.0},{0.0,0.0},//n=4
 												{1.02199,1.02199},{1.9347,1.9347}, //n=6
 												{2.705,2.705},{0.86011,0.86011},//n=8
 												{1.1173263161,1.1173263161},{1.36333,1.36333}, //n=10
 												{1.603255,1.603255},{1.839678,1.839678},//n=12
 												{2.074023,2.074023},{2.3071099136,2.3071099136},//n=14
 												{2.539424,2.539424},{2.7712634713,2.7712634713},
 												{3.0028062354,3.0028062354},{3.2341645842,3.2341645842},
 												{3.4654076648,3.4654076648},{3.6965785572,3.6965785572} }; //nReflecting Boundary Values


 const double synapseCascadeFilterUnifiedWithDecay::mdDecay_r100[][2] =
 { {0.0000000000,0.0000000000},{0.0000000000,0.0000000000},//n=2
   {0.0000000000,0.0000000000},{1.0986122887,1.0986122887},//n=4
   {0.0000000000,0.0000000000},{0.1675660078,0.1675660078}, //n=6
   {0.3121871366,0.3121871366},{0.4443631040,0.4443631040},//n=8
   {0.5694798840,0.5694798840},{0.6904734513,0.6904734513}, //n=10
   {0.8090210052,0.8090210052},{0.9261117805,0.9261117805},//n=12
   {1.0423410960,1.0423410960},{1.1580711338,1.1580711338},//n=14
   {1.2735229248,1.2735229248},{1.3888307621,1.3888307621},//n=15 - Its always reflecting ANYWAY so this wont be used
   {1.5040751772,1.5040751772},{1.6193032749,1.6193032749},
   {1.7345414227,1.7345414227},{0.0,0.0} };


//alpha - 0.05
const double synapseCascadeFilterUnifiedWithDecay::mdDecay_r010[][2] ={ {0.0,0.0},{0.0,0.0},//n=2
												{0.0,0.0},{0.182322,0.182322},//n=4
												{0.0,0.0},{0.0192914,0.0192914},//n=6
												{0.0412691,0.0412691},{0.0670579,0.0670579},//n=8
												{0.0974952,0.0974952},{0.133265,0.133265},//n=10
												{0.174915,0.174915},{0.222847,0.222847}, //n=12
												{0.277292,0.277292}, {0.338298,0.338298}, //n=14
												{0.405733,0.405733},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}};

//alpha = 0.05
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
const int8_t synapseCascadeFilterUnifiedWithDecay::miThreshold_r010[][2] ={ {1,1},{2,1},//n=2
																	{2,2},{2,2},//n=4
																	{4,4},{4,4}, //n=6
																	{4,4},{4,4},//n=8
																	{4,4},{4,4}, //n=10
																	{4,4},{4,4}, //n=12
																	{4,4},{4,4}, //n=14
																	{4,4},////n=15
																	{4,4},{4,4},{4,4},{4,4},{4,4}};

 //TODO: UNIFIED FILTER WITH DECAY: NEED TO MAKE THE REFLECTING VALUES FOR R=0.1 and r=0.01
const int8_t synapseCascadeFilterUnifiedWithDecay::miThresholdTerminal_r010[][2] ={ {0,1},{0,1},
 												{1,1},{3,1},
 												{2,2},{2,2}, //n=6
 												{2,2},{3,3},
 												{3,3},{3,3}, //n=10
 												{3,3},{3,3}, //n=12
 												{3,3},{3,3}, //n=14
 												{3,3},
 												{3,3},{3,3},{3,3},{3,3},{3,3}	}; ////n=15


const double synapseCascadeFilterUnifiedWithDecay::mdDecayTerminal_r010[][2] ={ {0.0,0.0},{0.0,0.0},//n=2
 												{0.0,0.0},{0.0,0.0},//n=4
 												{1.022,1.022},{1.9347,1.9347}, //n=6
 												{2.7,2.7},{0.8601,0.8601},//n=8
 												{1.11733,1.11733},{1.36333,1.36333}, //n=10
 												{1.60326,1.60326},{1.83968,1.83968},//n=12
 												{2.07402,2.07402},{2.3120,2.3120},//n=14
 												{2.53942,2.53942},
 												{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; //Reflecting Boundary Values


//alpha - 0.005
const double synapseCascadeFilterUnifiedWithDecay::mdDecay_r001[][2] ={ {0.0,0.0},{0.0,0.0},//n=2
												{0.0,0.0},{0.0198026,0.0198026},//n=4
												{0.0582689,0.0582689},{0.131028,0.131028},//n=6
												{0.262364,0.262364},{0.482426,0.482426},//n=8
												{0.0449093,0.0449093},{0.0664888,0.0664888},//n=10
												{0.0958382,0.0958382},{0.135234,0.135234}, //n=12
												{0.0358837,0.0358837}, {0.0463981,0.0463981}, //n=14
												{0.0593447,0.0593447},//n=15
												{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0}}; //I have not Calc these values for this rate

//alpha = 0.005
//Array of Upper and Lower Threshold values used by strong and weak cascade, 1st column is p escape and 2nd is q
///Matched Escape times to 1/q+p
//N=15 Double Decay Filter
const int8_t synapseCascadeFilterUnifiedWithDecay::miThreshold_r001[20][2] ={  {1,1},{2,1}, //n=2
															{2,2},{2,2}, //n=4
															{2,2},{2,2}, //n=6
															{2,2},{2,2}, //n=8
															{3,3},{4,4}, //n=10
															{3,3},{3,3}, //n=12
															{4,4},{4,4}, //n=14
															{4,4},//n=15
															{4,4},{4,4},{4,4},{4,4}};


 //TODO: MISSING VALUES UNIFIED FILTER WITH DECAY: NEED TO MAKE THE REFLECTING VALUES FOR R=0.1 and r=0.01
 //The set for A Strong Synapse is {p,q} So Low threshold Escape is the 2nd Value in the set
const int8_t synapseCascadeFilterUnifiedWithDecay::miThresholdTerminal_r001[][2] ={ {0,1},{0,1},
 												{1,1},{3,1}, //Low thres is 1 and High 3 For a Strong Synapse
 												{2,2},{2,2}, //n=6
 												{2,2},{3,3},
 												{3,3},{3,3}, //n=10
 												{3,3},{3,3}, //n=12
 												{3,3},{3,3}, //n=14
 												{3,3},////n=15
 												{3,3},{3,3},{3,3},{3,3},{3,3}};


const double synapseCascadeFilterUnifiedWithDecay::mdDecayTerminal_r001[][2] ={ {0.0,0.0},{0.0,0.0},//n=2
 												{0.0,0.0},{0.0,0.0},//n=4
 												{1.022,1.022},{1.9347,1.9347}, //n=6
 												{2.7,2.7},{0.8601,0.8601},//n=8
 												{1.11733,1.11733},{1.36333,1.36333}, //n=10
 												{1.60326,1.60326},{1.83968,1.83968},//n=12
 												{2.07402,2.07402},{2.3120,2.3120},//n=14
 												{2.53942,2.53942},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0},{0.0,0.0} }; //Reflecting Boundary Values



//Set Up the Array with the sets indexed by rate r
//This is a pointer to an array - Not An array of pointers
const PtrThresholdSetSmall synapseCascadeFilterUnifiedWithDecay::miThresholdsSet[] = {&synapseCascadeFilterUnifiedWithDecay::miThreshold_r100,
																 &synapseCascadeFilterUnifiedWithDecay::miThreshold_r010,
																 &synapseCascadeFilterUnifiedWithDecay::miThreshold_r001};

//This is a pointer to an array of thresholds- Not An array of pointers
const PtrThresholdSetSmall synapseCascadeFilterUnifiedWithDecay::miTerminalThresholdsSet[] = {&synapseCascadeFilterUnifiedWithDecay::miThresholdTerminal_r100,
																		 &synapseCascadeFilterUnifiedWithDecay::miThresholdTerminal_r010,
																		 &synapseCascadeFilterUnifiedWithDecay::miThresholdTerminal_r001};


//Set Up the Array with the sets indexed by rate r For Decay and Terminal States Decay
const PtrDecaySet synapseCascadeFilterUnifiedWithDecay::mdDecaysSet[] = {&synapseCascadeFilterUnifiedWithDecay::mdDecay_r100,
																 &synapseCascadeFilterUnifiedWithDecay::mdDecay_r010,
																 &synapseCascadeFilterUnifiedWithDecay::mdDecay_r001};
//
////This is a pointer to an array of thresholds- Not An array of pointers
const PtrDecaySet synapseCascadeFilterUnifiedWithDecay::mdTermDecaySet[] = {&synapseCascadeFilterUnifiedWithDecay::mdDecayTerminal_r100,
																	 &synapseCascadeFilterUnifiedWithDecay::mdDecayTerminal_r010,
																	 &synapseCascadeFilterUnifiedWithDecay::mdDecayTerminal_r001};


synapseCascadeFilterUnifiedWithDecay::synapseCascadeFilterUnifiedWithDecay():super() {
	//
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	mprng = g_getRandGeneratorInstance(false);
	miCascadeIndex = 0;
	miTerminalIndex = 0;
	miParameterSetUsed = 0;
	setFilterThresholds();
	initialiseFilterState();
	mbIsMonitored = false;

	//errexit(400,"synapseFilterUnifiedWithDecay: Do not Use This Constructor"); //Base class does not initialize
}
synapseCascadeFilterUnifiedWithDecay::synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,gsl_rng * rng_r)
{
	//:super(piCascadeSize,rng_r)

	mprng = rng_r;
	miTerminalIndex = piCascadeSize - 1;
	//RAND INDEX
	double r = gsl_rng_uniform(mprng)*0.999;
	miCascadeIndex = miStartIndex = floor(piCascadeSize*r);//Take Cascade Index as concat double to give 0-Count

	super::init( piCascadeSize,miCascadeIndex,SYN_STRENGTH_NOTSET,0);

	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	//setFilterThresholds();
	//initialiseFilterState();
}

synapseCascadeFilterUnifiedWithDecay::synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,
															ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
															gsl_rng * rng_r):super(piCascadeSize,penumStartStrength,rng_r)
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	setFilterThresholds();
	initialiseFilterState();
}
//Used to Start From A fixed point and test how distribution evolves
synapseCascadeFilterUnifiedWithDecay::synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,gsl_rng * rng_r):super(piCascadeSize,piStartIndex,rng_r)
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	setFilterThresholds();
	initialiseFilterState();

}

//Default Value is rate=1.0
synapseCascadeFilterUnifiedWithDecay::synapseCascadeFilterUnifiedWithDecay(int piCascadeSize,int piStartIndex,
															ICascadeSynapse::SYN_STRENGTH_STATE penumStartStrength,
															gsl_rng *  prng_r,
															int iRateDependentParameterSet):super(piCascadeSize,piStartIndex,penumStartStrength,prng_r,iRateDependentParameterSet)
{
	//Need to Explicitly call these in the constructor - In other cases They are handled by the virtual pointer
	setFilterThresholds();
	initialiseFilterState();
}


int synapseCascadeFilterUnifiedWithDecay::addSample(int iValue)
{
	miTimeSinceLastInduction++; //Increment time since last induction
	if (iValue != 0) //Some Induction step - So Calculate Decay Steps in between
	{
		doStochasticDecay();
		miTimeSinceLastInduction = 0; //Reset time since last event
	}
	return super::addSample(iValue);
}

////////DECAY FUNCTIONS///////////
/// There is a single decay rate depending on the sign of the running sum
//Decrement running sum based on decay probability eta * filter state (running sum)
//As we approach the threshold the decay increases proportionally
int synapseCascadeFilterUnifiedWithDecay::doStochasticDecay()
{

	double p;
	unsigned int rDecaySteps;

	//g_rng_r = getRandGeneratorInstance();
	if (miRFilterValue == 0) return 0;

	if (miRFilterValue > 0) //p side decay
	{
		p = 1.0-exp(-mdDecayRate*miTimeSinceLastInduction);
		rDecaySteps = gsl_ran_binomial(mprng,p,miRFilterValue);
		miRFilterValue -= rDecaySteps; //Subtract to move to zero
	}
	else //Increment Sum - q side is decaying back to zero
	{
		p = 1.0-exp(-mdDecayRate*miTimeSinceLastInduction);
		rDecaySteps = gsl_ran_binomial(mprng,p,-miRFilterValue);
		miRFilterValue += rDecaySteps; //Add To move to zero
	}


	return rDecaySteps;
}


//Sets where the running sum should be initialized after every change to Cascade state -
void synapseCascadeFilterUnifiedWithDecay::reInjectFilterStateToCascadeState()
{

#ifdef	UNIFILTER_RESET_ZERO
	miRFilterValue = 0;
#else
	initialiseFilterState();
#endif
}

////Sets where the running sum should be initialise When the Object is first created
void synapseCascadeFilterUnifiedWithDecay::initialiseFilterState()
{

	double r = gsl_rng_uniform(mprng);
	double p = 0.0;
	int i; //Iterator
	int countUpOrDown = 0;
	const double (*dPDFUsed)[20][ciMaxInternalStates]; //Pointer to PDF. If at terminal State then Use The reflecting Boundary PDF

	mbIsMonitored = false;
	miRFilterValue = 0;
	//return;
	if (((miTerminalIndex) < 1)) return; //Not For 1,2 thresholds filter

	//TODO: Add PDF Sets For All Rates -
	if (miCascadeIndex < miTerminalIndex )
		dPDFUsed = (&mdPDF_r100);
	else
		dPDFUsed = (&mdPDFTerminal_r100);


	if ((miCascadeIndex) == 1) //Small Thresholds Use Uniform
	{
		//Assign Random Start State --Add +1 To Include State 0 bUT LIMIT TO BELOW THRESHOLD
		miRFilterValue = round(r*(-(miLThres+1) + (miHThres-1))) +(miLThres+1);
	}
	else //Use init Distribution
	{

		for ( i = 0; i< ciMaxInternalStates;i++)
		{
			if (penumStrength == ICascadeSynapse::SYN_STRENGTH_WEAK)
				countUpOrDown = i; // miStartPStepsDone > miStartNStepsDone ???If Weak Synapse - The PDF is directed Correctly - Low Bound Is Reflecting
			else
				countUpOrDown = ciMaxInternalStates - i-1; //-1 cause We require an index

			p += (*dPDFUsed)[miCascadeIndex][i]; //Accumulate the Pdf
			if (p > r){
				//Found the spot since r was just exceeded
				miRFilterValue = countUpOrDown-3; //Remove Offset so i=0 becomes state -3 floor(ciMaxInternalStates/2)
				break;
			}
		}
	}
	//miRValue = 0;
	//Note That miRVal can be = to Threshold If the Threshold Is a Holding and Not Absorbing
if(( (miRFilterValue < miLThres) || (miRFilterValue > miHThres) ))
{
	assert((miRFilterValue > miLThres) && ((miRFilterValue < miHThres)));

}




}


//Use Index / strength and the Lookup arrays to set new L/H thresholds
//Called after changes to cascade index
//TODO : Add the decay
void  synapseCascadeFilterUnifiedWithDecay::setFilterThresholds()
{
	if (penumStrength == SYN_STRENGTH_STRONG)
	{
		if (miCascadeIndex < miTerminalIndex)
		{
			miLThres = -(*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1]; //Q
			miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0]; //P
			mdDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
		}
		else //Terminal States
		{
			miLThres = -(*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
			miHThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
			mdDecayRate = (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][0];
		}

	}else{//WEAK SYNAPSE
		if (miCascadeIndex < miTerminalIndex) //Non Terminal
		{
			miLThres = -(*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][0]; //Q
			miHThres = (*miThresholdsSet[miParameterSetUsed])[miCascadeIndex][1]; //P
			mdDecayRate = (*mdDecaysSet[miParameterSetUsed])[miCascadeIndex][0];
		}
		else //Terminal States
		{
			miLThres = -(*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][0];
			miHThres = (*miTerminalThresholdsSet[miParameterSetUsed])[miCascadeIndex][1];
			mdDecayRate = (*mdTermDecaySet[miParameterSetUsed])[miCascadeIndex][0];
		}

	} //ENDOF IF STRONG

	//implement Symmetric Assymetric Idea By Randomly Switching Between thres
	//Should We Randomly Switch Between Thresholds ?
	if (miCascadeIndex != miTerminalIndex)
	{
		//g_rng_r = getRandGeneratorInstance(false);
		double r = gsl_rng_uniform(mprng);
		if (r < 0.5)
		{
			dummy4 = miLThres;
			miLThres = -miHThres;
			miHThres = -dummy4;
		}
	}
}

double synapseCascadeFilterUnifiedWithDecay::getHDecay() const
{
	return mdDecayRate;
}

double synapseCascadeFilterUnifiedWithDecay::getLDecay() const
{
	return mdDecayRate;
}

double synapseCascadeFilterUnifiedWithDecay::getDecay() const
{
	return mdDecayRate;
}



int8_t synapseCascadeFilterUnifiedWithDecay::handlePOT()
{
	return super::handlePOT();
}

int8_t synapseCascadeFilterUnifiedWithDecay::handleDEP()
{
	return super::handleDEP();
}

void synapseCascadeFilterUnifiedWithDecay::handleNOP()
{
	super::handleNOP();
}

void synapseCascadeFilterUnifiedWithDecay::reset()
{
	super::reset(); //If we Re-init from Non zero position Then The escape time through any boundary Reduces
	 //The reset is only called by the testEscapetime Function - And thus it affects the result of this function
}

void synapseCascadeFilterUnifiedWithDecay::getTypeAsString(char* buff)
{
	getTypeName(buff);
}

void synapseCascadeFilterUnifiedWithDecay::getTypeName(char* buff)
{
	strcpy(buff,"_synapseCascadeFilterUnifiedWDecay");
}
synapseCascadeFilterUnifiedWithDecay::~synapseCascadeFilterUnifiedWithDecay() {
	//  Auto-generated destructor stub
}
