

#include  "../stdafx.h"
//#include "scalarProd_kernel.cu"

//Prototype
__global__ void scalarSigProdGPU(
    float *d_C,
    int *d_A,
    int *d_W,
    int vectorN,
    int elementN
);

// Enable this for error checking
#define CUDA_CHECK_ERROR

#define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        if ( cudaSuccess != err )
        {
            fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }
    } while ( 0 );

#pragma warning( pop )

#endif  // CUDA_CHECK_ERROR

    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        cudaError_t err = cudaGetLastError();
        if ( cudaSuccess != err )
        {
            fprintf( stderr, "cudaCheckError() failed at %s:%i : %s.\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }

        // More careful checking. However, this will affect performance.
        // Comment if not needed.
        err = cudaThreadSynchronize();
        if( cudaSuccess != err )
        {
            fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s.\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }
    } while ( 0 );

#pragma warning( pop )

#endif // CUDA_CHECK_ERROR

    return;
}

#ifndef _SIGNAL_KERNEL_H_
#define _SIGNAL_KERNEL_H_


#define MAX_BLOCK_DIM_SIZE 65535

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template<class T>
struct SharedMemory
{
    __device__ inline operator       T*()
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }

    __device__ inline operator const T*() const
    {
        extern __shared__ int __smem[];
        return (T*)__smem;
    }
};


// Device code
__global__ void VecMult(int* X,int* W, int* C, unsigned int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N)
        C[i] = X[i]*W[i];
}


template <class T>
__global__ void
reduce3(T *g_idata, T *g_odata, uint n)
{
    T *sdata = SharedMemory<T>();
   // T *sdataW = SharedMemory<T>();
    
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

    T mySum = (i < n) ? g_idata[i] : 0;
    if (i + blockDim.x < n) 
        mySum += g_idata[i+blockDim.x];  

    sdata[tid] = mySum;
    __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) 
    {
        if (tid < s) 
        {
            sdata[tid] = mySum = mySum + sdata[tid + s];
        }
        __syncthreads();
    }

    // write result for this block to global mem 
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


void initCUDADeviceMem(void*& d_W,void*& d_X,void*& d_C,void*& d_odata,float*& h_odata,unsigned int _uiSynCount,uint _uiTrackCount)
{
	 uint size = _uiSynCount*sizeof(int); //Memory Size for Input Vector
	 int threadsPerBlock = 512;
	 int N = _uiSynCount;
	 //int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	    
	 CudaSafeCall( cudaMalloc((void**)&d_W, size) ); //Weight V
	 CudaSafeCall( cudaMalloc((void**)&d_X,size*_uiTrackCount) ); //Input V
	 cout << "P:" << d_X << endl;
	 CudaSafeCall( cudaMalloc((void**)&d_C, _uiTrackCount*sizeof(float)) ); //Result Mux Vector Memory Space
	 //CudaSafeCall( cudaMalloc((void**)&d_odata, blocksPerGrid*sizeof(float)) ); //Reduction Results Per block

	 h_odata = new float[_uiTrackCount]; //The Output of each GPU Block - Copied to the Host
	 memset(h_odata,1,sizeof(float)*_uiTrackCount);
}


void cleapUpCUDADeviceMem(void*& d_W,void*& d_X,void*& d_C,void*& d_odata,float*& h_odata,uint _uiSynCount)
{

		 delete [] h_odata;
	    //Clean Up Device
		cout << "P:" << d_X  << endl;
		CudaSafeCall( cudaFree(d_X) );
		CudaSafeCall( cudaFree(d_C) );
		//CudaSafeCall( cudaFree(d_odata) );
		CudaSafeCall( cudaFree(d_W) );
	
}

void transferVectorsToDevice(int iNoTrackedPats,int* h_W ,t_inVal* h_X, uint _uiSynCount,void*& d_W,void*& d_X)
{
	uint size = _uiSynCount*sizeof(t_inVal);
    CudaSafeCall( cudaMemcpy((void*)d_W, (const void*)h_W, size, cudaMemcpyHostToDevice) ); //Weight
   	CudaSafeCall( cudaMemcpy((void*)d_X, (const void*)h_X, size*iNoTrackedPats, cudaMemcpyHostToDevice) ); //Input Patterns
	
}

//Get Perceptron Signal But do not Use Synapses Pointed by track group
// iTrackedIndex : Give the index of the currently tested Tracked pattern from the list of tracked Patts - Copy Weight Vector only the 1st time optimization
//Return The number of Tracked Patterns That can be tested now given how many patterns have been stored
int testCUDAPRecallOfX(float* h_sigdata,t_inVal* W ,t_inVal** X,t_inVal* tX, uint _uiSynCount,void* d_W,void* d_X,void* d_C,t_patt_trackedtbl& vTrackedIndex,uint _uiPattsStoredSoFar)
{
	//Test Recall of _iStoredPatIndex
	//cout << "Recall index: " << _iStoredPatIndex << " Output should be :" << X[_iStoredPatIndex][_uiSynCount-1];
	uint N = _uiSynCount;
	//uint size = N*sizeof(int);
    //int threadsPerBlock = 512;
   // int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

	int i=0;
	//assert(uiNoOfPatternsStoredInTrial > uiInitPatterns);
	t_patt_trackedtbl::iterator itTracked;
	for (itTracked = vTrackedIndex.begin();itTracked!=vTrackedIndex.end();++itTracked) //Copy Tracked Patterns that have Occured Up To now
	{
		if (itTracked->first < (_uiPattsStoredSoFar)) //Only If Tracked pattern Has been created - Note We are Comparing Index To Patt count
			memcpy(&tX[i*_uiSynCount],X[itTracked->first],sizeof(t_inVal)*_uiSynCount); //Join All tracked Patterns into one vector
		i++;
	}
	int iTrackedCount = i;
	//cout << W[0] << "X:" << X[0] << " N:" << _uiSynCount << " TP:"<<iTrackedCount << " dW"<<d_W << " dX"<<d_X << endl;
	// Allocate vectors in device memory happens once at beginning
	// Copy vectors from host memory to device memory
    //Only Copy One Weight Vector Which is used Against all Tracked Patterns
    //CudaSafeCall( cudaMemcpy((void*)d_W, (const void*)W, size, cudaMemcpyHostToDevice) ); //Weight
   //	CudaSafeCall( cudaMemcpy((void*)d_X, (const void*)X, size*iNoTrackedPats, cudaMemcpyHostToDevice) ); //Input Patterns
     transferVectorsToDevice(iTrackedCount,W,tX,_uiSynCount,d_W,d_X);
	
//	  CudaSafeCall( cudaMemcpy((void*)d_W, (const void*)W, size, cudaMemcpyHostToDevice) ); //Weight
//	  CudaSafeCall( cudaMemcpy((void*)d_X, (const void*)X, size*iNoTrackedPats, cudaMemcpyHostToDevice) ); //Input Patterns
	// Invoke kernel for Dot Prod
	 CudaSafeCall( cudaDeviceSynchronize());
	 scalarSigProdGPU<<<128, 256>>>((float*)d_C, (int*)d_X,(int*)d_W, iTrackedCount, N);
	 CudaSafeCall( cudaDeviceSynchronize());
	 CudaCheckError();
	 cudaMemcpy((void*)h_sigdata,(const void*)d_C, iTrackedCount*sizeof(float), cudaMemcpyDeviceToHost); //Get Result
	 
//	 for (int i=0;i<N;i++)
//		 h+=W[i]*X[i];
	 
	// h = h_odata[0]; //Test 1st Patt only
	  
	//MEasure The normalized PostSynaptic Respose As Signal
    //_Signal = h_odata[0]/_uiSynCount;
	// cout << _Signal << endl;
    
//	_SignalNTracked = hNTrack/_uiSynCount;//X[_iStoredPatIndex][_uiSynCount-1]*hNTrack/(X[0][0]*X[0][0]*(_uiSynCount-1));
	//cout << X[_iStoredPatIndex][_uiSynCount-1]/(X[0][0]*X[0][0]*(_uiSynCount-1)) << endl;
	//int iNeuronOut = ((_Signal)>0)?1:-1; //Neuron Classifier OUtput
	//cout << " Signal: " << _Signal << endl;

//	if (X[_iStoredPatIndex][_uiSynCount-1] == iNeuronOut)
//	{
//		Ret = 1; //Return 1 To indicate Successful classification of input
//	}

return iTrackedCount;
}


#endif // #ifndef _SIGNAL_KERNEL_H_
