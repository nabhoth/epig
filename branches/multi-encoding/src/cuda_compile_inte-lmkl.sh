#!bin/sh

set LD_LIBRARY_PATH /opt/intel/cmkl/10.1.1.019/lib/em64t:/usr/local/cuda/lib:/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/lib:$LD_LIBRARY_PATH

/usr/local/cuda/bin/nvcc -I. -I.. -I/usr/local/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/opt/intel/cmkl/10.1.1.019/include -c -o tools.o tools.cc
/usr/local/cuda/bin/nvcc -I. -I.. -I/usr/local/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/opt/intel/cmkl/10.1.1.019/include -c -o main.o main.cc
/usr/local/cuda/bin/nvcc -I. -I.. -I/usr/local/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/opt/intel/cmkl/10.1.1.019/include -c -o cstGA.o cstGA.cc
/usr/local/cuda/bin/nvcc -I. -I.. -I/usr/local/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/opt/intel/cmkl/10.1.1.019/include -c -o Qe.o Qe.cc
/usr/local/cuda/bin/nvcc -I. -I.. -I/usr/local/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/opt/intel/cmkl/10.1.1.019/include -c -o ccuda_multi.o ccuda_multi.cu
/usr/local/cuda/bin/nvcc -g -O2 -o epig tools.o main.o cstGA.o Qe.o ccuda_multi.o -lm -lmkl -lmkl_blacs_ilp64 -lmkl_gf_ilp64 -lmkl_lapack -lguide -lpthread -lcublas -lcudart -L/home/kameyama/lukacm/Eclipse_Epig_CUBLAS/src -L/opt/intel/cmkl/10.1.1.019/lib/em64t -L/usr/local/cuda/lib -L/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/lib
