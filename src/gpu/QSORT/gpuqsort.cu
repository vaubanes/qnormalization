#include "stdio.h"
#include "gpuqsort.h"
#include "Qfunc.c"

#include "simpletimer.cu"

#include <algorithm>

#define INDICES_TMP "index_tmp.bin"

// Keep tracks of the data blocks in phase one
template <typename element>
struct BlockSize
{
	unsigned int beg;
	unsigned int end;
	unsigned int orgbeg;
	unsigned int orgend;
	element		 rmaxpiv;
	element		 lmaxpiv;
	element		 rminpiv;
	element		 lminpiv;

	bool		 altered;
	bool		 flip;
	element		 pivot;
};

// Holds parameters to the kernel in phase one
template <typename element>
struct Params
{
	unsigned int from;
	unsigned int end;
	element pivot;
	unsigned int ptr;
	bool last;
};

// Used to perform a cumulative sum between blocks.
// Unnecessary for cards with atomic operations.
// Will be removed when these becomes more common
template <typename element>
struct Length
{
	element maxpiv[MAXBLOCKS];
	element minpiv[MAXBLOCKS];

	unsigned int left[MAXBLOCKS];
	unsigned int right[MAXBLOCKS];
};

// Since we have divided up the kernel in to three
// we need to remember the result of the cumulative sum
// Unnecessary for cards with atomic operations.
// Will be removed when these becomes more common
struct Hist
{
	unsigned int left[(MAXTHREADS)*MAXBLOCKS];
	unsigned int right[(MAXTHREADS)*MAXBLOCKS];
};

struct LQSortParams
{
	unsigned int beg;
	unsigned int end;
	bool flip;
	unsigned int sbsize;
};

#include "gpuqsort_kernels.cu"

#undef THREADS
#define THREADS threads

/**
* The main sort function
* @param data		Data to be sorted
* @param size		The length of the data
* @param timerValue Contains the time it took to sort the data [Optional]
* @returns 0 if successful. For non-zero values, use getErrorStr() for more information about why it failed.
*/
template <typename element>
int GPUQSort<element>::sort(  
        struct Files *fList,
        struct params *p,
 unsigned int blockscount, unsigned int threads, unsigned int sbsize, unsigned int phase)
{
 unsigned int size=p->nG;
 unsigned int nE = p->nE;
 unsigned int * data = new unsigned int[size];
        float * dataf = new float[size];
 unsigned int * pos = new unsigned int[size]; 
 float * cum = new float[size];
 double medida, timerValue[5];	
 
 int exp,i;
 
	if(!init)
		return 1;

	if(!threads||!blockscount||!sbsize)
	{
		threads   = 1<<(int)round(log(size * TK + TM)/log(2.0));
		blockscount = 1<<(int)round(log(size * MK + MM)/log(2.0));
		sbsize    = 1<<(int)round(log(size * SK + SM)/log(2.0));
	}
 

#ifdef HASATOMICS
		unsigned int* doh;
		unsigned int oh;

		cudaGetSymbolAddress((void**)&doh,"ohtotal");
		oh=0;
		cudaMemcpy(doh,&oh,4,cudaMemcpyHostToDevice);
#endif

	if(threads>MAXTHREADS)
		return 1; 

	if(blockscount>MAXBLOCKS)
		return 1;

timerValue[0]=0.0;
timerValue[1]=0.0;
timerValue[2]=0.0;
timerValue[3]=0.0;
timerValue[4]=0.0;
  

SimpleTimer st,st1,st2;
int device;
if(p->GPU==-1)
{
//MRequena: Displaying available devices visible to host
int deviceCount;
cudaGetDeviceCount (&deviceCount);
printf("Available devices:\n");
for (device=0; device<deviceCount; device++){
cudaDeviceProp deviceProp;
cudaGetDeviceProperties(&deviceProp, device);
printf("\tDevice %d: %s \n",device,deviceProp.name);
}
device= -1;
while (device<0 || device>=deviceCount){
printf("\nIntrduce un dispositivo: ");
scanf("%d", &device);
}
}
else  device = p->GPU;
printf("\nUsing GPU device %d.\n",device);
   cudaSetDevice(device);    // CUT_DEVICE_INIT();


	// Copy the data to the graphics card and create an auxiallary array
	ddata2 = 0; ddata = 0; dcum =0;
        dpos =0; dpos2=0;  //andres
	if(!errCheck(cudaMalloc((void**)&ddata2,(size)*sizeof(element))))
		return 1;
	if(!errCheck(cudaMalloc((void**)&ddata,(size)*sizeof(element))))
		return 1;

        // andres prepare indices
	if(!errCheck(cudaMalloc((void**)&dpos,(size)*sizeof(unsigned int))))
		return 1;
	if(!errCheck(cudaMalloc((void**)&dpos2,(size)*sizeof(unsigned int))))
		return 1;
	
 	// init mean vector to zero values
        if(!errCheck(cudaMalloc((void**)&dcum,(size)*sizeof(float))))
		return 1;

	if(!errCheck(cudaMemset(dcum, 0, size*sizeof(float)) ))
		return 1;
		
	// abre fichero
	FILE * fI;
	char nombref[1024];
		
	if((fI=fopen(INDICES_TMP,"wb"))==NULL) return 1;
		
	

// bucle de experimentos
for(exp=0; exp<nE; exp++)
{

// carga fichero
  	
  	// Start measuring time
		cudaThreadSynchronize();
		
		st.start();
    printf("Loading %3d/%d File: %s Genes: %d\n",exp,nE,fList[exp].fname,size);
	
    LoadFile(fList, exp, data);
    
    // Measure the time taken by loading from HD
		timerValue[1]+= st.end();
		
		st.start();
	
 	
	if(!errCheck(cudaMemcpy(ddata, data, size*sizeof(element), cudaMemcpyHostToDevice) ))
		return 1;

    // init pos vector
        unsigned int k;
        for(k=0; k<size; k++) {pos[k]=k;}  // hacer esto en paralelo en la GPU ?
    	if(!errCheck(cudaMemcpy(dpos, pos, size*sizeof(unsigned int), cudaMemcpyHostToDevice) )) return 1;

        // init_pos <<<blockscount, threads>>> (dpos, size);   // version GPU
        
        // Measure the time taken by loading into device
		timerValue[2]+= st.end();
	

	
	
			// Start measuring time
		
		st.start();
	

	// We start with a set containg only the sequence to be sorted
	// This will grow as we partition the data
	workset[0].beg = 0;
	workset[0].end = size;
	workset[0].orgbeg = 0;
	workset[0].orgend = size;
	workset[0].altered = false;
	workset[0].flip = false;


	// Get a starting pivot
	workset[0].pivot = (min(min(data[0],data[size/2]),data[size-1]) + max(max(data[0],data[size/2]),data[size-1]))/2;
	unsigned int worksize = 1;


	unsigned int blocks = blockscount/2;
	unsigned totsize = size;
	unsigned int maxlength = (size/blocks)/4;

	unsigned int iterations = 0;
	bool flip = true;


	// Partition the sequences until we have enough
	while(worksize<blocks)
	{

		unsigned int ws = totsize/blocks;
		unsigned int paramsize = 0;
		// Go through the sequences we have and divide them into sections
		// and assign thread blocks according to their size

		for(unsigned int i=0;i<worksize;i++)
		{

			if((workset[i].end-workset[i].beg)<maxlength)
				continue;

			// Larger sequences gets more thread blocks assigned to them
            unsigned int blocksassigned = max((workset[i].end-workset[i].beg)/ws,1);

			for(unsigned int q=0;q<blocksassigned;q++)
			{
				params[paramsize].from = workset[i].beg + ws*q;
				params[paramsize].end = params[paramsize].from + ws;
				params[paramsize].pivot = workset[i].pivot;
				params[paramsize].ptr = i;
				params[paramsize].last = false;
				paramsize++;
				
			}
			params[paramsize-1].last = true;
			params[paramsize-1].end = workset[i].end;

			workset[i].lmaxpiv=0;
			workset[i].lminpiv=0xffffffff;
			workset[i].rmaxpiv=0;
			workset[i].rminpiv=0xffffffff;
		}

		if(paramsize==0)
			break;
        		// Copy the block assignment to the GPU
		if(!errCheck(cudaMemcpy(dparams, params, paramsize*sizeof(Params<element>), cudaMemcpyHostToDevice) ))
			return 1;
 if(p->Verbose) printf("part1\n");

		// Do the cumulative sum
		if(flip)
			part1<<< paramsize, THREADS, (THREADS+1)*2*4+THREADS*2*4 >>>(ddata,dparams,dhists,dlength);
		else
			part1<<< paramsize, THREADS, (THREADS+1)*2*4+THREADS*2*4 >>>(ddata2,dparams,dhists,dlength);
		if(!errCheck((cudaMemcpy(length, dlength,sizeof(Length<element>) , cudaMemcpyDeviceToHost) )))
			return 1; 

		// Do the block cumulative sum. Done on the CPU since not all cards have support for
		// atomic operations yet. 
		for(unsigned int i=0;i<paramsize;i++)
		{
			unsigned int l = length->left[i];
			unsigned int r = length->right[i];
			
			length->left[i] = workset[params[i].ptr].beg;
			length->right[i] = workset[params[i].ptr].end;
			
			workset[params[i].ptr].beg+=l;
			workset[params[i].ptr].end-=r;
			workset[params[i].ptr].altered = true;
			
			workset[params[i].ptr].rmaxpiv = max(length->maxpiv[i],workset[params[i].ptr].rmaxpiv);
			workset[params[i].ptr].lminpiv = min(length->minpiv[i],workset[params[i].ptr].lminpiv);
			
			workset[params[i].ptr].lmaxpiv = min(workset[params[i].ptr].pivot,workset[params[i].ptr].rmaxpiv); 
			workset[params[i].ptr].rminpiv = max(workset[params[i].ptr].pivot,workset[params[i].ptr].lminpiv); 

			
		}

		// Copy the result of the block cumulative sum to the GPU
		if(!errCheck((cudaMemcpy(dlength, length, sizeof(Length<element>), cudaMemcpyHostToDevice) )))
			return 1;

		// Move the elements to their correct position
 if(p->Verbose) printf("part2\n");
		if(flip)
			part2<<< paramsize, THREADS >>>(dpos,dpos2,ddata,ddata2,dparams,dhists,dlength);
		else
			part2<<< paramsize, THREADS >>>(dpos2,dpos,ddata2,ddata,dparams,dhists,dlength);

		// Fill in the pivot value between the left and right blocks
		//part3<<< paramsize, THREADS >>>(ddata,dparams,dhists,dlength);

		flip = !flip;

		// Add the sequences resulting from the partitioning
		// to set
		unsigned int oldworksize = worksize;
		totsize = 0;
		for(unsigned int i=0;i<oldworksize;i++)
		{
			if(workset[i].altered)
			{
				if(workset[i].beg-workset[i].orgbeg>=maxlength)
					totsize += workset[i].beg-workset[i].orgbeg;
				if(workset[i].orgend-workset[i].end>=maxlength)
					totsize += workset[i].orgend-workset[i].end;

				workset[worksize].beg = workset[worksize].orgbeg = workset[i].orgbeg;
				workset[worksize].end = workset[worksize].orgend = workset[i].beg;
				workset[worksize].flip=flip;
				workset[worksize].altered = false;
				workset[worksize].pivot = (workset[i].lminpiv/2+workset[i].lmaxpiv/2);

				worksize++;

				workset[i].orgbeg = workset[i].beg = workset[i].end;
				workset[i].end = workset[i].orgend;
				workset[i].flip=flip;
				workset[i].pivot = (workset[i].rminpiv/2+workset[i].rmaxpiv/2);
				workset[i].altered = false;
			}
		}
		iterations++;

	}

	// Due to the poor scheduler on some graphics card
	// we need to sort the order in which the blocks
	// are sorted to avoid poor scheduling decisions
	unsigned int sortblocks[MAXBLOCKS*2];
	for(int i=0;i<worksize;i++)
		sortblocks[i]=((workset[i].end-workset[i].beg)<<(int)round(log((float)(MAXBLOCKS*4.0f))/log(2.0f))) + i;
	std::sort(&sortblocks[0],&sortblocks[worksize]);

	if(worksize!=0)
	{
		// Copy the block assignments to the GPU
		for(int i=0;i<worksize;i++)
		{
		 	unsigned int q = (worksize-1)-(sortblocks[i]&(MAXBLOCKS*4-1));

			lqparams[i].beg =  workset[q].beg;
			lqparams[i].end = workset[q].end;
			lqparams[i].flip = workset[q].flip;
			lqparams[i].sbsize = sbsize;
//                      printf("BEFS: %3d %10d %10d %10d %10d\n",i,lqparams[i].beg,lqparams[i].end,lqparams[i].flip,sbsize);
		}

		if(!errCheck((cudaMemcpy(dlqparams, lqparams, worksize*sizeof(LQSortParams), cudaMemcpyHostToDevice) )))
			return 1;
                
 if(p->Verbose) printf("lqsort\n");
		// Run the local quicksort, the one that doesn't need inter-block synchronization
		if(phase!=1)
			lqsort<<< worksize, THREADS >>>(dpos,dpos2,ddata,ddata2,dlqparams,phase);
	
	}

	cudaThreadSynchronize();
	
 if(p->Verbose) printf("promedio .......... %d\n",exp);
        promedio <<<blockscount, threads>>> (ddata, dcum, size);
        
    // Measure the Time taken by CPU+GPU
                medida= st.end();
		timerValue[0]+=medida; 
                printf("processing time: %8.2f ms\n",medida);
	
	st.start();    
    // lee posiciones y escribelas a disco    
    if(!errCheck((cudaMemcpy(pos, dpos, size*sizeof(unsigned int), cudaMemcpyDeviceToHost) )))		return 1;
    
    // Measure the time taken by loading back into CPU
		timerValue[4]+= st.end();
	
	// write to pos file
	    st.start();    
	     
          //fseek(fI, size*exp*sizeof(unsigned int), SEEK_SET);
        fwrite(pos, sizeof(unsigned int), size, fI);    
	    
    // Measure the time taken by loading back into CPU
		timerValue[3]+= st.end();
	            
        
     } // end bucle exp   
     
     fclose(fI);
       
       
 if(p->Verbose) printf("final divide .......... \n");
        divide   <<<blockscount, threads>>> (dcum, nE, size);
 
	
	

	err = cudaThreadSynchronize();
	// Free the data
	if(err!=cudaSuccess)
	{
		cudaFree(ddata);
		cudaFree(ddata2);
		cudaFree(dpos2);
		cudaFree(dpos);
		cudaFree(dcum);
		return 1;
	}
	// Copy the result back to the CPU
	//if(!errCheck((cudaMemcpy(data, ddata, size*sizeof(element), cudaMemcpyDeviceToHost) )))		return 1;
	
     st.start();	
	if(!errCheck((cudaMemcpy(cum, dcum, size*sizeof(element), cudaMemcpyDeviceToHost) )))
		return 1;
		
		// Measure the time taken by loading back into CPU
		timerValue[4]+= st.end();
		
    
    
    FILE * fO;
    
    // write to OUTPUT file
	st.start();    
	
	if((fI=fopen(INDICES_TMP,"rb"))==NULL) return 1;
	if((fO=fopen(p->fOutName,"wb"))==NULL) return 1;

	for(exp=0;exp<nE;exp++)
	{
	  fread(pos, sizeof(unsigned int), size, fI);    
	  for(i=0;i<size;i++)
	  {
	   dataf[i]=cum[pos[i]];
	  }	
	  fwrite(dataf,sizeof(float), size, fO);
if(p->Verbose)
{
		for(int i=0;i<10;i++)
		printf("data  ::: %5d = %10.2f   ,   %5d = %10.2f   ,    %5d = %10.2f\n",i,dataf[i],size/2-5+i,dataf[size/2-5+i],size-10+i,dataf[size-10+i]);
	 printf("----------\n\n");	
}
if(p->Verbose)
{
		for(int i=0;i<10;i++)
		printf("pos  ::: %5d = %10d   ,   %5d = %10d   ,    %5d = %10d\n",i,pos[i],size/2-5+i,pos[size/2-5+i],size-10+i,pos[size-10+i]);
	 printf("----------\n\n");	
}
	  
	}
	
	fclose(fI); fclose(fO);
	    
	    
    // Measure the time taken by loading back into CPU
		timerValue[3]+= st.end();
	
	
	
	//debug salida
if(0 && p->Verbose)
{
		for(int i=0;i<1000000;i++)
		//printf("cum  ::: %5d = %10.2f   ,   %5d = %10.2f   ,    %5d = %10.2f\n",i,cum[i],size/2-500+i,cum[size/2-500+i],size-1000+i,cum[size-1000+i]);
		printf("%10.2f,%10.2f,%10.2f,%10.2f,%10.2f,%10.2f\n",cum[i],cum[i]+1000000,cum[i]+2000000, cum[i]+3000000, cum[i]+4000000, cum[i]+5000000 );
		
}
		printf("\n----------  TIEMPOS -----------------------\n");
	printf("CPU+GPU  : %10.2f\n",timerValue[0]);
	printf("HD->CPU  : %10.2f\n",timerValue[1]);
	printf("CPU->DEV : %10.2f\n",timerValue[2]);
	printf("DEV->CPU : %10.2f\n",timerValue[4]);
	printf("CPU->HD  : %10.2f\n",timerValue[3]);
	    printf("--------------------------------------------\n");
    printf("TOTAL    : %10.2f\n\n",timerValue[0]+timerValue[1]+timerValue[2]+timerValue[4]+timerValue[3]);
	
	
	cudaFree(ddata);
	cudaFree(ddata2);
		cudaFree(dpos2);
		cudaFree(dpos);
		cudaFree(dcum);
 
	return 0;
}

template <typename element>
bool GPUQSort<element>::errCheck(int e)
{
	if(e==cudaSuccess)
		return true;

	err = e;
	cudaFree(ddata);
	cudaFree(ddata2);
	return false;
}

template <typename element>
GPUQSort<element>::GPUQSort():init(false),workset(0),params(0),length(0),lqparams(0),dlqparams(0),
							  dhists(0),dlength(0),dparams(0)
{
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	if(!strcmp(deviceProp.name,"GeForce 8800 GTX"))
	{
		TK = 1.17125033316e-005f;
		TM = 52.855721393f;
		MK = 3.7480010661e-005f;
		MM = 476.338308458f;
		SK = 4.68500133262e-005f;
		SM = 211.422885572f;
	}
	else
	if(!strcmp(deviceProp.name,"GeForce 8600 GTS"))
	{
		TK = 0.0f;
		TM = 64.0f;
		MK = 0.0000951623403898f;
		MM = 476.338308458f;
		SK = 0.0000321583081317f;
		SM = 202.666666667f;
	}
	else
	{
		TK = 0;
		TM = 128;
		MK = 0;
		MM = 512;
		SK = 0;
		SM = 512;
	}

	if(cudaMallocHost((void**)&workset,MAXBLOCKS*2*sizeof(BlockSize<element>))!=cudaSuccess) return;
	if(cudaMallocHost((void**)&params,MAXBLOCKS*sizeof(Params<element>))!=cudaSuccess) return;
	if(cudaMallocHost((void**)&length,sizeof(Length<element>))!=cudaSuccess) return;
	if(cudaMallocHost((void**)&lqparams,MAXBLOCKS*sizeof(LQSortParams))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dlqparams,MAXBLOCKS*sizeof(LQSortParams))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dhists,sizeof(Hist))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dlength,sizeof(Length<element>))!=cudaSuccess) return;
	if(cudaMalloc((void**)&dparams,MAXBLOCKS*sizeof(Params<element>))!=cudaSuccess) return;

	init = true;
}

/**
* Returns the latest error message
* @returns the latest error message
*/
template <typename element>
const char* GPUQSort<element>::getErrorStr()
{
	return cudaGetErrorString((cudaError_t)err);
}

template <typename element>
GPUQSort<element>::~GPUQSort()
{
	cudaFreeHost(workset);
	cudaFreeHost(params);
	cudaFreeHost(length);
	cudaFreeHost(lqparams);
	cudaFree(dparams);
	cudaFree(dlqparams);
	cudaFree(dhists);
	cudaFree(dlength);
}

// Exported functions

char* expErrMsg = "No errors";

 GPUQSort<unsigned int>* s=0;

extern "C" 
int gpuqsort( struct Files *fList,
        struct params *p, unsigned int blockscount, unsigned int threads, unsigned int sbsize, unsigned int phase)
{
	if(s==0)
		s=new GPUQSort<unsigned int>();


	if(s->sort(fList,p, blockscount, threads, sbsize, phase)!=0)
	{
		expErrMsg = (char*)s->getErrorStr();
		return 1;
	}
	else
		return 0;
}


extern "C"
const char* getGPUSortErrorStr()
{
	return expErrMsg;
}
