/*
 * util.cpp
 *
 *  Created on: Dec 19, 2011
 *      Author: kostasl
 */
#include "common.h"
#include "util.h"


///GLOBAL Instance and FunctioN - Called by some Default Constructors
//static gsl_rng* g_rng_r = 0;
//Default value is False
gsl_rng* g_getRandGeneratorInstance(bool newInstance)
{
		 ///Setup the GSL Random number Generator
		  // select random number generator
		//Return Pointer to global instance if it has been initialized
		if (g_rng_r && !newInstance)
			  return g_rng_r;

		  time_t t; //Used for random num generation
		  unsigned int seed = unsigned(time(&t)) + clock() + rand();
		  //unsigned int seed = 1;
		  gsl_rng* mrng_r = gsl_rng_alloc (gsl_rng_mt19937);

		  if (!mrng_r)
		  {
			  ERREXIT(1,"GSL RNG Init Failed. Out Of Memory?");
		  }

		  gsl_rng_set(mrng_r,seed);
		  //Seed random number generator
		  srand(seed);
		  //cout << "New Generator Seeded:" << seed;
	//END OF GSL SETUP

		  //SET TO GLOBAL INSTANCE
		 //if (!newInstance)
			  g_rng_r = mrng_r;

		  return mrng_r;
}

/* Called by makelogfilenames
   Uses Simulation Parameters to create a vector of strings holding file names to be used for loging distributions
   or data. vLogFiles[0] needs to contain the base file name part to use to construct the filename
   Ex. vLogFiles[0] = "mltCasc"; //THE CORPUS
*/
void MakeListOfFiles(vector<string>& vLogFiles,unsigned int ciInitPeriod,int iCascadeSize,double mdRate,double dFp, int iSynCount)
{
	char buffDFile[_MAX_PATH];
	memset(buffDFile, 0, _MAX_PATH);
	//Check if Base File Name Given
	if (vLogFiles.size() == 0 )
	{
		cout << "Base Filename not passed to makeLogFileNames" << endl;
		cerr << "Base Filename not passed to makeLogFileNames" << endl;
		ERREXIT(110,"makeLogFileNames- Corpus at [0] is Missing");
	}

	//File Path is Global So have to be carefull when using threads
	if (!getcwd(buffDFile, _MAX_PATH))// reads the current working directory into the array FilePath
		ERREXIT(50,"makeLogFileNames:Could not read working dir");


	strcat(buffDFile,"//");
	strcat(buffDFile,vLogFiles[0].c_str());


	string sbuffDistInit(buffDFile);
	sbuffDistInit.append(".distInit");//0
	string sbuffDistA(buffDFile); //Buffer Filename Strings - After Init Period
	sbuffDistA.append(".distA");//1
	string sbuffDistB(buffDFile); //After Memory Storage
	sbuffDistB.append(".distB");//2
	string sbuffDistC(buffDFile); //End Of simulation Distribution
	sbuffDistC.append(".distC");//3
	string sbuffDFile(buffDFile);
	sbuffDFile.append(".dat");//4
	string sbuffEFile(buffDFile);

	sbuffEFile.append(".distMemStor");//5

	string sbuffFDist(buffDFile);
	sbuffFDist.append(".distFilter");//6 //For Internal Filter State Reports

	string sbuffCov(buffDFile);
	sbuffCov.append(".Covariance");//6 //For Internal Filter State Reports

	string sbuffSig(buffDFile);
	sbuffSig.append(".MeanSignal");//6 //For Internal Filter State Reports

	string sbuffCovMatrix(buffDFile);
	sbuffCovMatrix.append(".CovarianceMatrix");//6 //For Internal Filter State Reports

	string sbuffMetaFile(buffDFile);
	sbuffMetaFile.append(".distMetaT");//9

	string sbuffCorrectMetaFile(buffDFile); //File to Save the metaplastic Cycle Frequency of synapses that are in the desired state for Memory
	sbuffCorrectMetaFile.append(".distMetaC");//10


	///OPEN LOG FILE - Make File Names///
	///////Get File Name
	////Make Distribution Data File Names
	vLogFiles.clear();
	vLogFiles.push_back(sbuffDistInit); //0
	vLogFiles.push_back(sbuffDistA);
	vLogFiles.push_back(sbuffDistB);
	vLogFiles.push_back(sbuffDistC);
	vLogFiles.push_back(sbuffDFile); //4

	vLogFiles.push_back(sbuffFDist);
	vLogFiles.push_back(sbuffCov);
	vLogFiles.push_back(sbuffSig);
	vLogFiles.push_back(sbuffCovMatrix); //8
	vLogFiles.push_back(sbuffMetaFile);
	vLogFiles.push_back(sbuffCorrectMetaFile); //10

}

//////////OUTPUT THRES-CYCLES OVER FIXED SAMPLE
void createCycleHistogramFile(map<uint,uint>& mDistrib, string outputFilename,unsigned long trials,unsigned long sampleTime)
{
	//Now Save the Distribution Of Metaplastic Counters Over a fixed Sampled Size
	cout << " Meta Distribution Output File: " <<  outputFilename << endl; //Tell User Which Output file we are using
	ofstream ofile3(outputFilename.c_str(), ios::out ); //Open Data File For Wrong Cycles
	ostringstream oss2;
	oss2 << "#Distribution of consecutive threshold crossings" << endl;
	oss2 << "#Trials:" << trials << " SampleLimit Reached At t:"<< sampleTime << endl;
	oss2 << "#Cycle-Size\tOverallFrequency\tCorrectStateFq\tWrongStateFq\tTotalSamples\tTrialNo\tSampleTime" << endl;
	ofile3 << oss2.str();

//USE APPEND To CYCLE FILE To write output
//		//Output To File
//		for (int i=1;i<=maxOccupancy;i++)
//		{
//			ofile3 << i << "\t"<< mDistrib[i] << "\t" << mDistribCorrect[i] << "\t"<<  mDistribWrong[i] << "\t" << lDistribSum << endl;
//		}
//		ofile3 << "#-\t-\t" << lDistribSumCorrect << "\t"<< lDistribSumWrong << endl;
//		ofile3.close();
}

///APPEND Threshold Cycles To file
void appendCycleHistogramToFile(map<uint,uint>& mDistrib, string outputFilename,unsigned long totalTrials,unsigned long sampleTime,unsigned long trialNo)
{
	//Save the Distribution Of Metaplastic Counters Sampled Over a Time interval
	//Split Into Two distributions to ease saving into one File
	unsigned long lDistribSum = 0; //The Total Number of Samples Is saved in 0
	unsigned long lDistribSumCorrect = 0; //The Total Number of Samples Is saved in 0
	unsigned long lDistribSumWrong = 0; //The Total Number of Samples Is saved in 0

	int maxOccupancy =0; //Find the Maximum State Occupied
	int iOccupancyIndex;
	map<uint,uint> mDistribCorrect;
	map<uint,uint> mDistribWrong;

	mDistribCorrect.clear();
	mDistribWrong.clear();

	for (map<uint,uint>::iterator it = mDistrib.begin();it!= mDistrib.end();++it)
	{
		if (it->first == 0)
			continue; //Skip the 0 Occupancy. It is used to hold the sum of samples used in the distribution

		lDistribSum += it->second; //Sum the Total Number of samples
		if (it->first > 999) //Correct Ones are given INdexes with +1000
		{
			iOccupancyIndex = it->first-1000;
			lDistribSumCorrect += mDistribCorrect[iOccupancyIndex] = it->second;
		}
		else
		{
			iOccupancyIndex = it->first;
			lDistribSumWrong+= mDistribWrong[iOccupancyIndex] = it->second;
		}

		if  (maxOccupancy < iOccupancyIndex) maxOccupancy = iOccupancyIndex;
		//Replace With Sum Of Both
		mDistrib[iOccupancyIndex] = mDistribWrong[iOccupancyIndex] + mDistribCorrect[iOccupancyIndex]; //Fix To total Value
	}
		//Now Save the Distribution Of Metaplastic Counters Over a fixed Sampled Size
		ofstream ofile3(outputFilename.c_str(), ios::app ); //Open Data File For Wrong Cycles

		//Output To File
		for (int i=1;i<=maxOccupancy;i++)
		{
			ofile3 << i << "\t"<< mDistrib[i] << "\t" << mDistribCorrect[i] << "\t"<<  mDistribWrong[i]
			            << "\t" << lDistribSum << "\t" << trialNo <<  "\t" << sampleTime << endl;
		}
		//ofile3 << "#-\t-\t" << lDistribSumCorrect << "\t"<< lDistribSumWrong << endl;
		ofile3.close();
}
