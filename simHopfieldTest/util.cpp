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
			  liberrexit(1,"GSL RNG Init Failed. Out Of Memory?");
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
	vLogFiles.push_back(sbuffDistInit);
	vLogFiles.push_back(sbuffDistA);
	vLogFiles.push_back(sbuffDistB);
	vLogFiles.push_back(sbuffDistC);
	vLogFiles.push_back(sbuffDFile);

	vLogFiles.push_back(sbuffFDist);
	vLogFiles.push_back(sbuffCov);
	vLogFiles.push_back(sbuffSig);
	vLogFiles.push_back(sbuffCovMatrix);
	vLogFiles.push_back(sbuffMetaFile);
	vLogFiles.push_back(sbuffCorrectMetaFile);

}


double getCAthres(int theta, int reps,int modelType)
{
if (modelType == 8) //synapseSingleFilterUnifiedWithDecay
{
	switch(reps)
	{
	case 1:
		return g_fCAthresUFilter[theta-2][0][2]; //obtain CA filter
	case 2:
		return g_fCAthresUFilter[theta-2][1][2]; //obtain CA filter
	case 4:
		return g_fCAthresUFilter[theta-2][2][2]; //obtain CA filter
	case 8:
		return g_fCAthresUFilter[theta-2][3][2]; //obtain CA filter
	default:
		//Return r=4 by default
		return g_fCAthresUFilter[theta-2][2][2]; //obtain CA filter
	}
}
if (modelType == 9 || modelType == 1) //SU Synapse
{
	switch(reps)
	{
	case 1:
		return g_fCAthresSUSpaced[theta-2][0][2]; //obtain CA filter
	case 2:
		return g_fCAthresSUSpaced[theta-2][1][2]; //obtain CA filter
	case 4:
		return g_fCAthresSUSpaced[theta-2][2][2]; //obtain CA filter
	case 8:
		return g_fCAthresSUSpaced[theta-2][3][2]; //obtain CA filter
	default:
		//Return r=4 by default
		return g_fCAthresSUSpaced[theta-2][2][2]; //obtain CA filter
	}
}

return 0.0;
}

