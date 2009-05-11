
#ifndef PQSORTH
#define PQSORTH

#define MAXTHREADS 1024
#define MAXBLOCKS 2048
//#define SBSIZE 1024


extern "C" 
int gpuqsort( struct Files *fList,
        struct params *p, unsigned int blockscount = 0, unsigned int threads = 0, unsigned int sbsize = 0, unsigned int phase = 0);

extern "C"
const char* getGPUSortErrorStr();
// Float support removed due to some problems with CUDA 2.0 and templates
// Will be fixed
//extern "C" 
//DLLEXPORT int gpuqsortf(float* data, unsigned int size, double* timerValue=0);

/**
* Returns the latest error message
* @returns the latest error message
*/
//extern "C" DLLEXPORT const char* getGPUSortErrorStr();


template <typename element> struct BlockSize;
template <typename element> struct Params;
template <typename element> struct Length;
struct Hist;
struct LQSortParams;

template <typename element>
class GPUQSort
{
	element* ddata; 
	element* ddata2; 
        float* dcum;
	unsigned int * dpos; 
	unsigned int * dpos2; 
	struct Params<element>* params;
	struct Params<element>* dparams;

	LQSortParams* lqparams;
	LQSortParams* dlqparams;

	Hist* dhists;
	Length<element>* dlength;
	Length<element>* length;
	BlockSize<element>* workset;

	float TK,TM,MK,MM,SM,SK;

	int err;
	bool init;

	bool errCheck(int e);
public:
	GPUQSort();
	~GPUQSort();

	int sort(  struct Files *fList,
        struct params *p, unsigned int blockscount = 0, unsigned int threads = 0, unsigned int sbsize = 0, unsigned int phase = 0);
	const char* getErrorStr();
};


#endif
