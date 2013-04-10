/*
 * CameraDemo.cpp
 *
 *  Created on: 3 Mar 2011
 *      Author: kostasl
 *
 *      //This Demo Will capture An image and save it onto an array which will be used as input
 *      to train a Hopfield Network
 */

#include "../stdafx.h"
#include "CameraDemo.h"
#include <pthread.h>
#include "cv.h"
#include "highgui.h"
#include "../HopfieldMemoryTests.h"

void cpImgToInputVector(int* X,IplImage* imgSrc,uint N, CvPoint& ROITopLeft, CvPoint& ROIBottomRight);
void cpImgToImg(IplImage* imgDest,IplImage* imgSrc,CvPoint& ROITopLeft, CvPoint& ROIBottomRight);
//void     copyVectorToImg(IplImage* imgDest,int* X,uint N,CvPoint& ROITopLeft, CvPoint& ROIBottomRight); //Something is Wrong
void cpVectorToImg(IplImage* imgDest,int* X,CvPoint& ROITopLeft, CvPoint& ROIBottomRight);


//Global So All threads can access it
float** W; //The Weight Matrix
int* rX; //The recalled Pattern
int* pX; //The Probe Pattern



int tids[10];
const int cConcurrentThreads = 8;

const uint iTargetWindow = 190; //One dimension of the square focus window
uint NetSize = iTargetWindow*iTargetWindow;

bool bStorageInProgress = false;
//bRecallInProgress = false; //Found in Hopfield

void* runRecallThread(void* pdata)
{
	float Signal;
	bool foundRunningThread = false;

	int nchunk = NetSize/cConcurrentThreads;
	uint startN = ((*(int*)pdata)-1)*nchunk;
	recallHopfieldNetPattern(nchunk,startN,pX,rX,W,Signal,30);


	tids[(*(int*)pdata)-1] = 0; //Set This Threads ID to 0 To flag that it has completed
	for (int i=0;i<cConcurrentThreads;i++) //Check the list of threads to See if they have all completed
		if (tids[i] !=0)
			foundRunningThread = true;

	bRecallInProgress = foundRunningThread; //Set The recall Flag to false if all threads have completed

	pthread_exit(0);
	return pdata;
}

void showCameraInput()
{
    CvCapture *capture = 0;
    IplImage  *frame_rgb = 0; //Colour Version of camera view
    IplImage *im_gray = 0; //GrayScale Version of camera view
    IplImage *im_gray2 = 0; //GrayScale Version of camera view
    IplImage *im_recallbw=0; //The Image In the focus Area Where the Recall pattern is saved-

    //float Signal = 0;


    CvPoint ROITopLeft; //((frame_rgb->width-iTargetWindow)/2, (frame_rgb->height-iTargetWindow)/2);
    CvPoint ROIBottomRight;// cvPoint((frame_rgb->width+iTargetWindow)/2, (frame_rgb->height+iTargetWindow)/2);
    CvPoint ptMessages = cvPoint(10, 130);
    CvPoint ptMessagesSuccess = cvPoint(10, 360);
    CvPoint ptInstruction1 = cvPoint(10, 10);
    CvPoint ptInstruction2 = cvPoint(10, 450);

    ROITopLeft.x = -100; //Init To fixed Value so Initialization can be detected Afterwards
    bRecallInProgress = false;
    //THREADING
	 pthread_t threads[10 ];		  /* holds thread info */

	 int errcode;
	 int interframeDelay = 100; //msWaitKey And refresh period
	//////


	///INIT MEMORY///
    W = new float*[NetSize]; //The Weight Matrix
    for (uint i=0;i<NetSize;i++)
    {
    	W[i] = new float[NetSize];
    	memset(W[i],0,NetSize*sizeof(float));
    }
    rX = new int[NetSize];
    pX	= rX;// new int[NetSize];



    char key = '\n';

    //Command Line Instrcutions
    cout << "-Associative Network Picture Memorization and Recall Demo-" << endl;
    cout << "Press 'q' to quit"<< endl;
    cout << "Press 's' to store image in red square"<< endl;
    cout << "Press 'm' to recall an image using image in red square as a cue"<< endl;
    cout << "Press 'r' to forget all memories"<< endl;
    cout << "---- Author: Konstantinos Lagogiannis 2011--" << endl;
    cout << "Using OpenCV Robotics vision library" << endl;
    cout << "Processing an Image size of " << iTargetWindow << "x" << iTargetWindow << " using " << NetSize<< endl;
    // create a window for the video
    cvNamedWindow( "Memory", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "TheEye", CV_WINDOW_AUTOSIZE );

    /* initialize font and add text */
    CvFont font;
    CvFont fontInstruction;
    cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 1.0, 1.0, 0, 2, CV_AA);
    cvInitFont(&fontInstruction, CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5, 0, 1, CV_AA);

//     initialize camera
     capture = cvCreateCameraCapture( 0 );

    // always check
    if ( !capture ) {
        errexit(100," Cannot create Capture stream 0");
    }
    //Grey Canvas


//MAIN LOOP - RECEIVING KEYSTROKES AND SHOWING IMAGES
    while( key != 'q' ) {
        // get a frame
        frame_rgb = cvQueryFrame( capture );
        // always check
        if( !frame_rgb ) break;

        //Get Dimensions And Draw Red Area
        /* draw a red box for ROI */
        if (ROITopLeft.x == -100) //Init Only Once
        {
        	ROITopLeft 	   = cvPoint((frame_rgb->width-iTargetWindow)/2 , (frame_rgb->height-iTargetWindow)/2);
        	ROIBottomRight = cvPoint(ROITopLeft.x+iTargetWindow , ROITopLeft.y+iTargetWindow);

//        	ROITopLeft 	   = cvPoint(0 , 30);
//        	ROIBottomRight = cvPoint(iTargetWindow,iTargetWindow+30 );
        	im_gray 	= cvCreateImage(cvSize(frame_rgb->width,frame_rgb->height),IPL_DEPTH_8U,1);
        	im_gray2 	= cvCreateImage(cvSize(frame_rgb->width,frame_rgb->height),IPL_DEPTH_8U,1);
        	im_recallbw = cvCreateImage(cvSize(iTargetWindow,iTargetWindow),IPL_DEPTH_8U,1);
        }


        //Fill GrayScale Image Buffer
        cvCvtColor(frame_rgb,im_gray,CV_RGB2GRAY);

        //Threshold TO B&W
        cvAdaptiveThreshold(im_gray,im_gray,255,CV_ADAPTIVE_THRESH_GAUSSIAN_C,CV_THRESH_BINARY,301,-3);
        //cvThreshold(im_gray,im_gray , 100, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);//use Otsu algorithm to Calculate Threshold

        //cvFlip(im_gray, NULL, 0);  cvFlip(im_gray, NULL, 1);
        cvCvtColor(im_gray,frame_rgb,CV_GRAY2RGB); //Convert Back to RGB Canvas So we can Draw Colour Objects


        cvRectangle(frame_rgb,                    /* the dest image */
        		ROITopLeft,        /* top left point */
        		ROIBottomRight,       /* bottom right point */
                cvScalar(0, 0, 255, 0), /* the color; blue */
							 5, 8, 0);               /* thickness, line type, shift */


       // display current frame

        if (bRecallInProgress)
		  {
        	cpVectorToImg(im_gray2,rX,ROITopLeft,ROIBottomRight); //Copy Back to Check What is happening
        	//copyVectorToImg(im_recallbw,rX,NetSize,ROITopLeft,ROIBottomRight);
        	cvShowImage( "Memory", im_gray2 );
        	cout << " Recall..." << endl;
        	bStorageInProgress = false;
		  }

        //Show Instructions
        cvPutText(frame_rgb, "Memorize images by pressing 's'.", ptInstruction1, &fontInstruction, cvScalar(10, 10, 255, 0));
        cvPutText(frame_rgb, "Recall using the current view as que with 'm'", ptInstruction2, &fontInstruction, cvScalar(10, 10, 255, 0));

        // exit if user press 'q'
        key = cvWaitKey( interframeDelay );
        interframeDelay = 100;//Reset To Default


        //Handle Key Command
        switch (key)
        {
        case 's':
        	cout << "Learning..." << endl;

        	///Save To Memory
        	cpImgToInputVector(pX,im_gray,NetSize,ROITopLeft,ROIBottomRight);
        	cpVectorToImg(im_gray2,pX,ROITopLeft,ROIBottomRight); //Copy Back to Check What is happening
        	//cpImgToInputVector(pX,im_gray2,NetSize,ROITopLeft,ROIBottomRight);
        	//cpImgToImg(im_gray2,im_gray,ROITopLeft,ROIBottomRight);
        	cvPutText(frame_rgb, "Learning...", ptMessages, &font, cvScalar(10, 10, 255, 0));
        	cvShowImage( "TheEye", frame_rgb );
        	cvShowImage( "Memory", im_gray2);
        	cvWaitKey(50);
        	bStorageInProgress = true;
            makeWeightMatrix(NetSize,pX,W);

//            for(int i = 0;i<NetSize;i++)
//                   { cout << endl;
//                   	for(int j = 0;j<NetSize;j++)
//                   		cout << W[i][j] << " ";
//                   }
            cout <<  "Memorized. " << NetSize << endl;
        	cvPutText(frame_rgb, "Picture memorized.", ptMessagesSuccess, &font, cvScalar(10, 10, 255, 0));
        	interframeDelay = 1000;
            break;
        case 'r':
        	//Reset Memory
        	 for (uint i=0;i<NetSize;i++)
        		 memset(W[i],0,NetSize*sizeof(W[i][0]));
        	 cout<< "Forgotten. " << endl;
        	 cvPutText(frame_rgb, "Forgot all memories", ptMessages, &font, cvScalar(10, 10, 255, 0));
        	 interframeDelay = 1000;
        	break;

        case 'm':


        	bRecallInProgress =true;
        	cout << "Recall from what I see..." << endl;
        	cvPutText(frame_rgb, "Trying to remember using current view as a cue...", ptMessages, &font, cvScalar(10, 10, 255, 0));
        	cpImgToInputVector(pX,im_gray,NetSize,ROITopLeft,ROIBottomRight);
        	interframeDelay = 1000;

        	for (int i=0;i<cConcurrentThreads;i++)
        	{
				tids[i] = threads[i] = i+1;
				errcode=pthread_create(&threads[i],  // thread struct
									  0,        // default thread attributes
									  runRecallThread, //start routine
									 &tids[i]); //Arg. To Routine

	    	  	  if (errcode)
	        	   {         errexit(errcode,"pthread_create");     	    }

        	}


        	break;
        case 'c':
        	cvPutText(frame_rgb, "Stopped Recall", ptMessages, &font, cvScalar(10, 10, 255, 0));
        	interframeDelay = 1000;
        	///CANCEL RECALL PROCEDURE
        	bRecallInProgress = false;
        	break;

        }

        cvShowImage( "TheEye", frame_rgb );

    } //END OF MAIN LOOP

    // free memory
    cvReleaseImage(&im_gray);
    cvReleaseImage(&im_recallbw);
    cvReleaseImage(&frame_rgb);

    cvReleaseCapture( &capture );
    cvDestroyWindow( "TheEye" );
    cvDestroyWindow( "Memory" );

    //pthread_cancel(	threads[1]);

    for (uint i=0;i<NetSize;i++)
    	delete [] W[i];


}

void copyVectorToImg(IplImage* imgDest,int* X,uint NetSize,CvPoint& ROITopLeft, CvPoint& ROIBottomRight)
{
	cout << "Vector->Img" << imgDest->height << "*" << imgDest->width << endl;
	assert(imgDest->depth == IPL_DEPTH_8U);
	uint8_t* imgPix = (uint8_t*)(imgDest->imageData);
	uint N = imgDest->width*imgDest->height;
	assert (NetSize >= N);

	for (uint i=0;i<N;i++)
	{
		//During Recall Thread The X may be out of the -1 - 1 region
//		if (X[i] != 1 && X[i] != -1)
//		{
//			cerr<< "X["<<i<< "]=" << X[i]<<endl;
//		}
		*(imgPix+i) = (uint8_t)(X[i] == 1)?254:0;
//		 cout << X[i] << " " <<  endl;
	}



}


void cpImgToImg(IplImage* imgDest,IplImage* imgSrc,CvPoint& ROITopLeft, CvPoint& ROIBottomRight)
{
	uint8_t* imgPixSrc = (uint8_t*)(imgSrc->imageData);
	uint8_t* imgPixDest = (uint8_t*)(imgDest->imageData);

	uint startRowPixel = ROITopLeft.x;
	uint endRowPixel = ROIBottomRight.x;
	uint rowLength = imgSrc->width;

	assert(imgSrc->depth == IPL_DEPTH_8U);

	uint j = 0; //Index Of Target Vector
	long startPix,endPix;
	//Copy Column
	for (uint k=ROITopLeft.y;k<ROIBottomRight.y;k++)
	{
//		cout << k;
		//Copy Row
		startPix = rowLength*k+startRowPixel;
		endPix = rowLength*k+endRowPixel;

		assert (startPix >= 0 && endPix > 0);
		assert (startPix < endPix);
		for (uint i=startPix;i<endPix;i++)
		{
			//imgPix[i] = 120;
			imgPixDest[i] = imgPixSrc[i];
			j++; //Index of target vector
		}
//		cout << endl;
		if (j > NetSize)
		{
			cerr << "Network Too Small For captured Area!";
			break;
		}
	}

}

void cpVectorToImg(IplImage* imgDest,int* X,CvPoint& ROITopLeft, CvPoint& ROIBottomRight)
{
	uint8_t* imgPixDest = (uint8_t*)(imgDest->imageData);

	uint startRowPixel = ROITopLeft.x;
	uint endRowPixel = ROIBottomRight.x;
	uint rowLength = imgDest->width;

	assert(imgDest->depth == IPL_DEPTH_8U);

	uint j = 0; //Index Of source Vector
	long startPix,endPix;
	//Copy Column
	for (uint k=ROITopLeft.y;k<ROIBottomRight.y;k++)
	{
		//cout << k;
		//Copy Row
		startPix = rowLength*k+startRowPixel;
		endPix = rowLength*k+endRowPixel;

		assert (startPix >= 0 && endPix > 0);
		assert (startPix < endPix);
		for (uint i=startPix;i<endPix;i++)
		{
			//imgPix[i] = 120;
			imgPixDest[i] = (X[j]==1)?254:0;
			j++; //Index of target vector
		}
//		cout << endl;
		if (j > NetSize)
		{
			cerr << "Network Too Small For captured Area!";
			break;
		}
	}

}

//Has to be A Gray Scale Img Pointer
void cpImgToInputVector(int* X,IplImage* imgSrc,uint NetSize,CvPoint& ROITopLeft, CvPoint& ROIBottomRight)
{
	uint8_t* imgPix = (uint8_t*)(imgSrc->imageData);
	uint startRowPixel = ROITopLeft.x;
	uint endRowPixel = ROIBottomRight.x;
	uint rowLength = imgSrc->width; //The number of pixels to skip a row at the source img

	assert(imgSrc->depth == IPL_DEPTH_8U);

	uint j = 0; //Index Of Target Vector

	//Copy Column
	for (uint k=ROITopLeft.y;k<ROIBottomRight.y;k++)
	{
//		cout << k;
		//Copy Row
		long startPix = rowLength*k+startRowPixel;
		long endPix = rowLength*k+endRowPixel;

		assert (startPix >= 0 && endPix > 0);
		assert (startPix < endPix);
		for (uint i=startPix;i<endPix;i++)
		{
			//imgPix[i] = 120;
			X[j] = (imgPix[i]>10)?1:-1;
//			cout << X[j] ;
			j++; //Index of target vector
		}
//		cout << endl;
		if (j > NetSize)
		{
			cerr << "Network Too Small For captured Area!";
			break;
		}
	}
}


