LIB_SOURCES 	= gpuqsort.cu
QNORM_SOURCES 	= GPUQnorm5.cpp 
LIB_OBJ		= $(LIB_SOURCES:.cu=.o)
QNORM_OBJ	= $(QNORM_SOURCES:.cpp=.o)

all: gpuqsortlib QNORM

gpuqsortlib: $(LIB_OBJ)
	mkdir -p lib/linux/static
	nvcc -o lib/linux/static/libgpuqsort.a -lib $(LIB_OBJ)


QNORM: $(QNORM_OBJ) gpuqsortlib
	mkdir -p bin/linux
	nvcc -o bin/linux/GPUQnorm5 $(QNORM_OBJ) -lgpuqsort -lglut -L./lib/linux/static

%.o: QSORT/%.cu
	nvcc -O3 -I$(CUDA_SDK_DIR)/common/inc -Iinc -c $< -o $@

%.o: QNORM/%.cpp
	gcc -O3 -Iinc -c $< -o $@
