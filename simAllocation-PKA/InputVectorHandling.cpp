/*
 * InputVectorHandling.cpp
 *
 *  Created on: 10 Mar 2011
 *      Author: kostasl
 */


#include "common.h"
#include "util.h"
//#include "../synapseModels/common.h"
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <stdio.h>

#include "InputVectorHandling.h"
using namespace std;
/*

Takes the default Hadamard Matrix File - Shuffles a required number of vectors and writes
the shuffled vectors to a new file given a named random set number

*/

inline double dotprod(t_inVal* X1,t_inVal* X2,uint VecSize)
{
	double cdot = 0.0;
	for (uint j=0;j<VecSize;j++)
	{
		cdot += X1[j]*X2[j];
	}
return cdot;
}
void makeHadamardShuffledSet(uint PatCount,uint VecSize)
{
	 //INPUT PATTERN
	t_inVal* X[PatCount]; //Train

	//Construct Memory for VECTOR X
	for (uint i=0;i<PatCount;i++)
	{
		X[i] = 0;

		while (X[i] == 0) //Retry In case of error
		{
		  X[i] = new (nothrow) t_inVal[VecSize];
		  if (X[i] == 0)
		  {
		    cerr << "Error: memory could not be allocated for pattern:" << i;
		  }
		}
	}

		char fname[100];
		sprintf(fname,"HD%d.dat",VecSize);
		readTestVectorsFromFile(fname,X,VecSize,PatCount,0,true);
		g_getRandGeneratorInstance(false);
		sprintf(fname,"HD%d-N%d-Set%d.dat",VecSize,PatCount,(int)(100*gsl_rng_uniform(g_rng_r)));
		writeTestVectorsToFile(fname,X,VecSize,PatCount);


}


/*
* Reads Vectors from a file and copies them to extrernally initialized vector X
* Shuffles the read File Vectors and Also tests loaded vectors for correlation
* -- The Number of Patterns is passed by reference so it can be changed. The calling function may be calling readVectors a few times
* and we don't want messages appearing repeatedly. Once the PatToLoad changes the warning message does not appear again.
*
* bRecycleFile : set this flag if you want the X vector filled regardless if the file contains less vectors. The file we be re-read sequencially
* from the beginning. This is used with Hadamart vectors.
*/

void readTestVectorsFromFile(const char* fname,t_inVal** X,uint VecSize, uint& iPatsToLoad,uint indexToInsertVectors, bool bShuffleVectors)
{
	const bool bRecycleFile = false; //Once we run out of patterns in the file - The file is read from the begining again
	const bool bTestAfterLoad = false;
	//const bool bShuffleVectors = false; //Default Param Value is False
	assert(indexToInsertVectors >= 0 && VecSize > 0 && X!=0);
	gsl_rng* prng = g_getRandGeneratorInstance(true);
	//Used for Testing Data
	float cdotTotal = 0;
	float maxcorrelation;
	float cdot = 0;
	const float cMaxCorrPercent = 0.51;

	map<uint,uint> vecIndex;
	string fInName(INPUT_VECTORS_DIRECTORY);
	fInName.append(fname);

	ifstream ifile(fInName.c_str(), ios::binary | ios::in ); //OPEN OUTPUT FILE

	if (!ifile.is_open())
	{
		cerr <<  "Could not open vector file : " << fInName.c_str();
		ERREXIT(500, "Could not open vector file for reading" );
	}
	//Sample
	if (bTestAfterLoad)
		cout << "Sample."<< endl;

	uint pos=0;
	ifile.seekg(0, ios::end); //Get File Length
	uint size =  ifile.tellg()/(sizeof(t_inVal)*VecSize);
	//cout << "Reading Vector File:" << fname << endl;
	//
	if (size < iPatsToLoad)
	{
		if (size ==0)
		{
			cerr << " No vectors found in file or file does not exist" << endl;
			ERREXIT(150,"Failed to load vectors from file.");
		}
		else
			cout << "File" << fname <<  " contains:"<< size << " vectors Requested to load:" << iPatsToLoad << ((bShuffleVectors)?" Shuffled":" in-order ") << endl;
			if (bRecycleFile)
			{
				//cout << "I Will be repeating the patterns in the file." << endl;
				//iPatsToLoad = size;
				bShuffleVectors = false;
			}
			else
				iPatsToLoad = size; //Limit Loaded patterns to the Available one in the file
	}

	const uint VecCount =  iPatsToLoad + indexToInsertVectors;

	ifile.seekg(0, ios::beg);
	for (uint i=indexToInsertVectors;i<VecCount;i++)
	{
		if (bShuffleVectors) //If Shuffling
		{
			do{ //Look for Unique Random vector not Picked Previously
				pos = std::ceil(gsl_rng_uniform(prng)*size)-1;
			}while (vecIndex.find(pos) != vecIndex.end());

			ifile.seekg(sizeof(t_inVal)*VecSize*pos, ios::beg);
			vecIndex[pos] = pos;

		} //Otherwise Read Sequentially

		ifile.read((char*)X[i],sizeof(t_inVal)*VecSize);

		///Readout A sample if required
		for (int k=0;k<20 && (i<8)&& bTestAfterLoad;k++)
			cout << ((X[i][k]>0)?"+":"-");
			if ((i<8) && bTestAfterLoad) cout << endl;

		uint bytesread = ifile.gcount();
		if (ifile.eof())
		{
			if (bRecycleFile)
				ifile.seekg(0, ios::beg);
			else
			{
				ifile.close();
				cerr << " Found " << (i-indexToInsertVectors) << " Instead of " << iPatsToLoad << endl;
				iPatsToLoad = (i-indexToInsertVectors);
				cerr << "Lowered Max Capacity to:" <<  iPatsToLoad << endl;
				break;
			}
			//ERREXIT(1,"Vectors in file are less than required"); //Changed to Non Fatal Error

		}else if (bytesread != sizeof(t_inVal)*VecSize )
		{
			cerr << i << " Invalid Vector Size Read - File : " << (int)(bytesread)/(int)sizeof(t_inVal) << " Expected:" << VecSize << endl;
			ERREXIT(200,"Invalid Vector Size Read - File does not match vector size");
		}
	}

	ifile.close();
	fInName.clear();
	////////TEST LOADED DATA//////////

	maxcorrelation = VecSize * X[0][0]*X[0][0] * cMaxCorrPercent;
	 for (uint i=indexToInsertVectors;(i<VecCount) && bTestAfterLoad;i++)
	 {
		 cdot = 0;
		for (uint k=indexToInsertVectors;k<i;k++)
		{
			cdotTotal +=cdot = dotprod(X[k],X[i],VecSize);

			//cout << i << "Vs" << k << " cdot:" << cdot << endl;

			if  (abs(cdot) > maxcorrelation)
			{
				cerr << "Max Correlation " << maxcorrelation << " Exceeded between patterns " << i  << " and " << k << endl;
				ERREXIT(1,"Max Correlation Exceeded ");
			}
		}

	 }


	if (bTestAfterLoad)
	{
		cout << "Loaded " << (VecCount-indexToInsertVectors);
			cout << " Vectors with All-to-All correlation :" << cdotTotal << endl;
	}
}

/*
 * Function Created for TE, to convert my binary file to txt
 */
void convertVectorFileToTxt(uint uiPatCount,uint iSynCount,char* pInputFile)
{
	t_inVal* X[uiPatCount]; //Memory PAtterns Containing The Ones Loaded from File and Random Initialization patterns

	gsl_rng* mprng = g_getRandGeneratorInstance(true);
	cout << "Making Text Version of Hadamard Vector File" << pInputFile << endl;
	//Initialise The memory For Patterns
	initPatternMemory(X,uiPatCount,iSynCount,0,0.5, mprng,false);

	//Read in The Vectors from the File to X
	readTestVectorsFromFile(pInputFile,X,iSynCount,uiPatCount,0,false); //Load At 0 All Required PAtterns

	char filename[250];
	strcpy(filename,pInputFile);
	strcat(filename,".txt");
	cout << uiPatCount;
	writeTestVectorsToFileTxt(filename,X, iSynCount, uiPatCount);
}


/*
///Allocates memory for each pattern if required if pointer to Array is Null
 Otherwise it assumes memory has been initialized and it just creates a new population of random vectors
*/
void initPatternMemory(t_inVal **X,uint PatCount,uint _uiNeuronCount,uint iTrackedMemIndex,float Fp,gsl_rng* mprng,bool bUseRandomPatterns)
{
	double r;
	//Init Memory For Patterns
	for (uint i=0;i<PatCount;i++)
	{
		if (X[i] == 0) //Memory May Have been initialized at a previous run.
			X[i] = new t_inVal[_uiNeuronCount];
		//Random Patterns if Required to fill the gap before the tracked pattern  - Or if Tesing with Random - Not Hadamard
		if (i <= iTrackedMemIndex || bUseRandomPatterns)
		{
			 for (uint j=0;j<_uiNeuronCount; j++)
			 {
				 r = gsl_rng_uniform(mprng);
				 X[i][j] = (r < Fp)?1:-1; //Make Test Vector
			 }
		}
		else
			memset(X[i],0,sizeof(t_inVal)*_uiNeuronCount); //Just Init to zero to detect creepy Bugs
	}
}

void writeTestVectorsToFileTxt(char* cfOutName,t_inVal**X,uint VecSize,uint PatCount)
{
	string fOutName(INPUT_VECTORS_DIRECTORY);
	fOutName.append(cfOutName);

	ofstream* ofile = new ofstream(fOutName.c_str(),  ios::out ); //OPEN OUTPUT FILE

	 if (!ofile->is_open())
	 {
		 cerr << "Could not open Output file-Directory Missing?" << cfOutName << endl;
		 ERREXIT(500,"writeTestVectorsToFileTxt:Could not open Output file-Directory Missing?");
	 }


	 for(uint i = 0; i < PatCount;i++)
	 {
		 int sum = 0;
		 for (uint j=0;j<VecSize;j++)
		 {
			 (*ofile) << X[i][j] << " ";
			 //cout <<  X[i][j] << " ";
			 sum += X[i][j];
		 }
		 (*ofile) << endl;
		 cout << sum << endl; //Sum should be 0 - Except one line
	 }

	//Test DOT PRODUCT
	ofile->close();

	cout << "Wrote " << PatCount << " to file :" << cfOutName << endl;
}

void writeTestVectorsToFile(char* cfOutName,t_inVal** X,uint VecSize,uint PatCount)
{
	string fOutName(INPUT_VECTORS_DIRECTORY);
	fOutName.append(cfOutName);

	ofstream* ofile = new ofstream(fOutName.c_str(), ios::binary | ios::out ); //OPEN OUTPUT FILE

	 if (!ofile->is_open())
	 {
		 cerr << "Could not open Output file-Directory Missing?" << cfOutName << endl;
		 ERREXIT(500,"Could not open Output file-Directory Missing?");
	 }


	 for(uint i = 0; i < PatCount;i++)
	 {
		 ofile->write((char*)X[i],sizeof(t_inVal)*VecSize);


	 }

	//Test DOT PRODUCT
	ofile->close();

	cout << "Wrote " << PatCount << " to file :" << cfOutName << endl;
}


/*
 * FUNCTIONS THAT SEARCH THE INPUT VECTOR FIELD
 */



//SINGLE POINT MUTATIONS OF A POPULATION OF VECTORS - Reverse mutation if fitness not improved
//RANDOM SEARCH FOR CONSTRUCTING CLOSE TO ORTHOGONAL VECTORS within 1% error
void makeGN2TestVectorsInFile(int PatCount,uint VecSize,float fbitBalance)
{
	//const uint sTimeout = 500000000;
	float corrThres = VecSize*0.25*0.01;
	//float corrThres = 0.25;
	gsl_rng* mprng =  g_getRandGeneratorInstance(true);

	uint p = 0;
	//string fOutName("..//HopFieldMemResultsV4//");
	char fname[100];
	sprintf(fname,"randGN2XVectorN_%dPat_%d.bin",(int)(VecSize),PatCount);

	string fOutName(INPUT_VECTORS_DIRECTORY);
	fOutName.append(fname);
	ofstream* ofile = new ofstream(fOutName.c_str(), ios::binary | ios::out ); //OPEN OUTPUT FILE

	 if (!ofile->is_open())
	 {
		 cerr << "Could not open Output file-Directory Missing?" << fOutName.c_str() << endl;
		 ERREXIT(500,"Could not open Output file-Directory Missing?")
	 }
	 //INPUT PATTERN
	 t_inVal* X[PatCount]; //Train
	 t_inVal  T[VecSize]; //Train
	float cdot = 0;
	float cdotprev = 0;


	uint SearchTime = 0;
	//Init Population
	for (int i=0;i<PatCount;i++)
	{
		X[i] = new t_inVal[VecSize];
		//Random Init
		for (uint j=0;j<VecSize;j++)
		{
				double r = gsl_rng_uniform(mprng );
				T[j] = X[i][j] = (r < fbitBalance)?0.5:-0.5; //Make Test Vector
			//cout << X[i][j];
		}
	}

	//INitial Fitness
	int c;
	cdot =0;
	for (int i=0;i<PatCount;i++)
	{
		cdot += abs(dotprod(X[0],X[i],VecSize));
	}
//Start Search
	while(cdot > corrThres)
	{
		SearchTime++;
		do
		{
		c = ceil((gsl_rng_uniform(mprng)+0.001)*(PatCount))-1;
		}while (c==PatCount);

		do
		{
			p = ceil((gsl_rng_uniform(mprng)+0.001)*(VecSize))-1;
		}while (p==VecSize);

		X[c][p] = -X[c][p]; //Mutation
		cdotprev = cdot;
		cdot = 0;
		//Get Population Fitness
		for (int y=0;y<PatCount;y++)
		{
			for (int i=0;i<PatCount;i++)
			{
				cdot += abs(dotprod(X[y],X[i],VecSize));
			}
		}
		//Reverse Change
		if (cdotprev < cdot)
		{
			X[c][p] = -X[c][p];
		}

		if (SearchTime%100000)
			cout << cdot<< endl;
	}

	//ofile->write((char*)X[i],sizeof(float)*VecSize);
	//Test DOT PRODUCT
	cout << "fin: cdot" << cdot;
	ofile->close();

	gsl_rng_free(mprng);
}



//RANDOM SEARCH FOR CONSTRUCTING CLOSE TO ORTHOGONAL VECTORS within 1% error
void makeGNTestVectorsInFile(int PatCount,uint VecSize,float fbitBalance)
{
	const uint sTimeout = VecSize*2000000;
	float corrThres = VecSize*0.25*0.055;
	//float corrThres = 0.25;
	gsl_rng* mprng =  g_getRandGeneratorInstance(true);
	char fname[100];
	int p = 0;
	//string fOutName("..//HopFieldMemResultsV4//");
	 //INPUT PATTERN
	t_inVal* X[PatCount]; //Train
	t_inVal  T[VecSize]; //Train
	float cdot = 0.0f;
	float cdotprev;
	int foundCount = 0;
	//Init Memory For Patterns
//	for (int i=0;i<PatCount;i++)
//	{
////		tX[i] = new float[VecSize]; //Not Required But.. Leave For now
//	}

	bool Search = true;

	uint SearchTime = 0;
	//Construct TEST VECTOR T
	for (int i=0;i<PatCount;i++)
	{
		X[i] = 0;

		while (X[i] == 0)
		{
		  X[i] = new (nothrow) t_inVal[VecSize];
		  if (X[i] == 0)
		  {
		    cerr << "Error: memory could not be allocated for pattern:" << i;
		  }
		}

		//Random Init
		for (uint j=0;j<VecSize;j++)
		{
			double r = gsl_rng_uniform(mprng );
//			if (i>0)
//			{
				T[j] = X[i][j] = (r < fbitBalance)?1:-1; //Make Test Vector
//			}
//			else
//			{
//				T[j] = X[i][j] = (j < VecSize/2)?0.5:-0.5;
//			}
			//cout << X[i][j];
		}

		int lastUpdate = 0;
		SearchTime = 0;
		Search = true;
		//Genetic Algorithm Style Random Mut.
		cdotprev = 0.0f;
		while (Search) //Recreate Pattern Until It is Unique
		{
			SearchTime++;
			lastUpdate++;
			Search = false;
			//Make A rand Vector
			//cout << endl;/home/kostasl/CodeProjects/memoryframework_ver1.01/Release/SynapticFilterMemory

			//Fitness - Search Dot Prod from Other Patterns
			cdot =0;
			for (int k=0;k<i;k++)
			{
				cdot += abs(dotprod(X[k],T,VecSize));
				//cout << SearchTime << " X_"<< k << ".X_"<< i << ": " << cdot << endl;
				//cout << i << " " << corrThres << " < " << cdot << endl;
			}
			//Save correlation of Start Vector
			if (lastUpdate == 1) cdotprev = cdot;

			Search = (cdot > corrThres);

			//Has Time Past and We are worse from where we started? Start Again with a small change
			if (Search && (lastUpdate > 60000)) //Reverse the change
			{	//cout << i << "Rst." << endl;
				///Mutate Starting Point
				p = ceil((gsl_rng_uniform(mprng)+0.001)*VecSize)-1;
				X[i][p] = -X[i][p]; //Make Test Vector

				for (uint j=0;j<VecSize;j++) //Reset Search Vector
					T[j] = X[i][j];

				lastUpdate=0;
			}


			if (!Search)///Found
				{
					///Copy to OutPut
					foundCount++;
					for (uint j=0;j<VecSize;j++)
					{
						cout << ((T[j]>0)?"+":"-") << "";
						X[i][j] = T[j];

					}

					//cout << endl << i << " " << corrThres << " < " << cdot << endl;
					break;
					lastUpdate=0;

				}

			//Do Mutation
			p = ceil((gsl_rng_uniform(mprng)+0.001)*VecSize)-1;
			T[p] = -T[p]; //Make Test Vector

			 //If Similarity not Detected Then Move to Next one
			if(SearchTime > sTimeout)//Do not Search For Ever
			{
				cout << "Tired of searching... exiting" << endl;
				Search =false;

				break;

			}
		}//Searching While

		if ((foundCount < (i+1))) break; //Pattern Not Found Stop looking
		//Write To File
		cout << i+1 << " Ts:" << SearchTime << endl;
	} //For Each Pattern

	sprintf(fname,"randGNXVectorN_%dPat_%d_Corr_%4.2f.bin",(int)(VecSize),foundCount,corrThres/0.25);

	writeTestVectorsToFile(fname,X,VecSize,foundCount);

	for(int i = 0; i < foundCount;i++)
		delete [] X[i];

	gsl_rng_free(mprng);
}


/*
 * Gram Schmitt Orthogonilization
 * Return The number of Patterns Found
 * Max Patterns if Vector Size is a multiple of 2 even numbers! powers of 2
 *
 */

uint makeGSTestVectorsInFile(uint PatCount,uint VecSize,float fbitBalance)
{
	gsl_rng* mprng =  g_getRandGeneratorInstance(true);
	char fname[200];
	 //INPUT PATTERN
	t_inVal* X[PatCount]; //Random Init Pats
	t_inVal* W[PatCount]; //Orthogonal Pats
	double c[PatCount]; //The Projection Coefficients
	float cdotTotal[PatCount];
	float cdot = 0.0;

	bool Search = true;
	float corrThres = VecSize*1.0*0.009;

	uint SearchTime = 0;

	//Construct TRAIN VECTOR //Construct Training Patterns
	//BEGIN GS LOOP
	uint StableTime = 0;
	int prevSuccess;
	uint bestMatch = 0;
	uint countSuccess = 0;

	StableTime = 0;
	unsigned short reportInterval = 0; //Once Overflow - Report Search Time

	while(Search && StableTime < 1000000)
	{
		SearchTime++;
		reportInterval++;
		//MAKE PATTERNS

		bool randomPat = true;
		for (uint i=0;i<PatCount;i++)
		{
			X[i] = new t_inVal[VecSize];
			W[i] = new t_inVal[VecSize];
			bool up =true;
			int jcount = 0;
			float split;
				//Make A rand Vector
				for (uint j=0;j<VecSize;j++)
				{

					split  = ((float)VecSize/(pow(2,(i+1))));
					if (jcount >= std::ceil(split))
					{
						up = !up;
						jcount =0;
					}
					jcount++;

					if (randomPat)
					{
						double r = gsl_rng_uniform(mprng );
						W[i][j]  = (r < fbitBalance)?1.0:-1.0; //Make Test Vector
					}
					else
					{
						 W[i][j] = (up)?1.0:-1.0;
					}
//					cout << (( W[i][j] > 0)?"+":"-") << "";

				}

//				cout << floor(split) << " - "<< randomPat <<  endl;

				//Switch To random Patterns After Structured Ones
				//randomPat = (floor(split)==0)?true:false;

		}//END OF MAKE RAND PATTERNS


//return 0;
	prevSuccess = countSuccess; //Save
	countSuccess = 0;

	//DO GRAM SCHMIT - SKIP 1st Pattern
	for (uint i=1;i<PatCount;i++)
	{
			//Do Gram-Schmit
			for (uint l=0;l<i;l++)
			{
				c[i] = dotprod(W[l],W[i],VecSize)/dotprod(W[l],W[l],VecSize);
				//cout << c[i]<< endl;
				//Multiply Coefficient to base to obtain projection
				for (uint k=0;k<VecSize;k++)
				{
					W[i][k]= (W[i][k] - c[i]*W[l][k]);
				}
			}

	}//For Each Pattern

	//cout<< SearchTime << " Fin... Now Reconstruct & measure"<< endl;

	countSuccess = 0;
	for (uint i=0;i<PatCount;i++)
	{
		cdotTotal[i] = 0;
		//cout << i << " ";
		//Reconstruct
		for (uint k=0;k<VecSize;k++)
		{
			W[i][k] = (W[i][k]>0)?1.0:-1.0;

		}
		//cout << endl;

		//cout << endl;

		for (uint l=0;l<i;l++)
		{
			cdot = abs(dotprod(W[i],W[l],VecSize));

			cdotTotal[i]  +=cdot;
//			cout << "W_"<<i<<"." << "W_"<<l<<"="<< cdot << endl;
		}
		//Count Vectors whose total correlation to other vectors is less than Threshold
		countSuccess += (cdotTotal[i] <= corrThres)?1:0;
	}//For Each Pattern

	//return 0;
	//Save if this Trial Is better
	if (countSuccess > bestMatch)
	{
		cout << "Better :" << countSuccess << endl;
		StableTime = 0;
		bestMatch = countSuccess;
		int g = 0;
		for (uint i=0;i<PatCount;i++)
		{
			if (cdotTotal[i] <= corrThres) //Save only The Orthogonal Ones
			{

				for (uint k=0;k<VecSize;k++)
				{
					X[g][k] = W[i][k];
					cout << ((X[g][k]  > 0)?"+":"-") << "";
				}
				g++;
				cout << endl;
			}

		}
		//Copy to Best matching Vectors
	}
	else
		StableTime++;

	Search = ((bestMatch < PatCount));
	if (reportInterval == 0 || reportInterval == 32768 ||  reportInterval == 16384 || reportInterval == 8192 || reportInterval == 4090)
		cout <<  SearchTime << " Solution Stabletime: "<< StableTime << " Best:" << bestMatch << " Unc.Pats Now:"<< countSuccess << endl;

}//WHile Search


//Save Results To FIle
	sprintf(fname,"randGSOrthoWVectorsFxStartN_%dPat_%d.bin",(int)(VecSize),bestMatch);

	string fOutName(INPUT_VECTORS_DIRECTORY);
	fOutName.append(fname);

	ofstream* ofile = new ofstream(fOutName.c_str(), ios::binary | ios::out ); //OPEN OUTPUT FILE

	 if (!ofile->is_open())
	 {
		 cerr << "Could not open Output file-Directory Missing?" << fOutName.c_str() << endl;
		 ERREXIT(500,"Could not open Output file-Directory Missing?");
	 }

	for (uint i=0;i<bestMatch;i++)
		ofile->write((char*)X[i],sizeof(float)*VecSize);

	//Test DOT PRODUCT
	ofile->close();

	gsl_rng_free(mprng);
	cout << "******** Vector Size:" << VecSize << "Fin**********"  << endl;
	cout <<  "t:" << SearchTime << " Solution Stabletime : "<< StableTime << " No. Of Uncorr. Pats : " << bestMatch <<endl;
	//cout << "Search:" << Search <<" Finished Stabletime: " << StableTime << endl;
	return bestMatch;
}


//RANDOM SEARCH FOR CONSTRUCTING CLOSE TO ORTHOGONAL VECTORS within 1% error
//Each Vector Randomly constructed until is orthogonal with Previous ones withing % error
void makeTestVectorsInFile(uint PatCount,uint VecSize,float fbitBalance)
{
	const uint sTimeout = 500000000;
	const float corrPercent = 1.00;
	float corrThres = VecSize*corrPercent;
	//float corrThres = 0.25;
	gsl_rng* mprng =  g_getRandGeneratorInstance(true);
	char fname[100];

	sprintf(fname,"randXVectorN_%dPat_%d_Corr_%2.2f.bin",(VecSize),PatCount,corrPercent);

	string fOutName(INPUT_VECTORS_DIRECTORY);
	fOutName.append(fname);

	ofstream* ofile = new ofstream(fOutName.c_str(), ios::binary | ios::out ); //OPEN OUTPUT FILE

	 if (!ofile->is_open())
	 {
		 cerr << "Could not open Output file-Directory Missing?" << fOutName.c_str() << endl;
		 ERREXIT(500,"Could not open Output file-Directory Missing?")
	 }
	 //INPUT PATTERN
	 t_inVal* X[PatCount]; //Train
	float cdot = 0;

	//Init Memory For Patterns
//	for (int i=0;i<PatCount;i++)
//	{
////		tX[i] = new float[VecSize]; //Not Required But.. Leave For now
//	}

	bool Search = true;

	uint SearchTime = 0;
	//Construct TRAIN VECTOR //Construct Training Patterns
	for (uint i=0;i<PatCount;i++)
	{
		X[i] = new t_inVal[VecSize];

		SearchTime = 0;
		Search = true;
		while (Search) //Recreate Pattern Until It is Unique
		{
			SearchTime++;
			Search = false;
			//Make A rand Vector
			for (uint j=0;j<VecSize;j++)
			{
				double r = gsl_rng_uniform(mprng );
				X[i][j] = (r < fbitBalance)?1:-1; //Make Test Vector
				//cout << X[i][j];
			}
			//cout << endl;

			//Search Dot Prod from Other Patterns
			cdot =0;
			for (uint k=0;k<i;k++)
			{
				cdot = dotprod(X[k],X[i],VecSize);

				//cout << SearchTime << " X_"<< k << ".X_"<< i << ": " << cdot << endl;
				Search = ((cdot > corrThres) || (cdot < -corrThres));
				//cout << i << " " << corrThres << " < " << cdot << endl;
				if (Search) break;
			}
			 //If Similarity not Detected Then Move to Next one
			if(SearchTime > sTimeout)//Do not Search For Ever
			{
				cout << "Tired of searching... exiting" << endl;
				Search =false;
				ofile->close();
				return;

			}
		}//Searching While
		//Write To File
		cout << i << " Ts:" << SearchTime << endl;

		ofile->write((char*)X[i],sizeof(t_inVal)*VecSize);
	}
	//Test DOT PRODUCT
	ofile->close();

	gsl_rng_free(mprng);
}



//Creates A driftless stimulation By making a random pattern then alternating it for the required number of vectors
/// This effect could be made much simpler just by alternating the last bit of the vector which denotes the output of the Neuron
void makeDriftlessTestVectorsInFile(uint PatCount,uint VecSize,float fbitBalance)
{
	gsl_rng* mprng =  g_getRandGeneratorInstance(true);
	char fname[150];

	sprintf(fname,"AlternatingXVectorN_%dPat_%d.bin",(VecSize),PatCount);

	string fOutName(INPUT_VECTORS_DIRECTORY);
	fOutName.append(fname);

	ofstream* ofile = new ofstream(fOutName.c_str(), ios::binary | ios::out ); //OPEN OUTPUT FILE

	 if (!ofile->is_open())
	 {
		 cerr << "Could not open Output file-Directory Missing?" << fOutName.c_str() << endl;
		 ERREXIT(500,"Could not open Output file-Directory Missing?");
	 }
	 //INPUT PATTERN
	float* X[PatCount]; //Train
	//Init 1st Pattern

	X[0] = new float[VecSize]; //Init 1st pattern Memory
	for (uint j=0;j<VecSize;j++)
	{
		double r = gsl_rng_uniform(mprng );
		X[0][j] = (r < fbitBalance)?1:-1; //Make Test Vector
	//cout << X[i][j];
	}

	//Construct TRAIN VECTOR //Construct Training Patterns
	for (uint i=1;i<=PatCount;i++)
	{
		X[i] = new float[VecSize]; //init This patterns Memory

		//Make A rand Vector
		cout << endl;
		cout << i;

		for (uint j=0;j<(VecSize-1);j++) //All But the last one should Alternate - The last one is the required output
		{
			X[i][j] = -1*X[i-1][j]; //Make -Ve of previous Pattern
			//cout << "\t" << X[i][j] ;
		}
		X[i][(VecSize-1)] = X[i-1][(VecSize-1)]; //Keep This Constant



		//Write To File
		ofile->write((char*)X[i],sizeof(float)*VecSize);
	}
	//Test DOT PRODUCT
	ofile->close();

	gsl_rng_free(mprng);
}
