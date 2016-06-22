#!/bin/sh
ARCH=32
echo "Compiling for $ARCH bit architecture"


if test $ARCH -eq 32
then
INC='-I. -I.. -I/usr/local/cuda32/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/usr/local/include/gsl'
LIB='-lm -lgslcblas -lpthread -lcublas -lcudart -L/home/kameyama/lukacm/Eclipse_Epig_CUBLAS/src -L/usr/local/lib -L/usr/local/cuda32/cuda/lib -L/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/lib'
QMDD='-m32 -I/home/local_space/QMDD'
CC='/usr/local/cuda32/cuda/bin/nvcc'
echo 'Remeber to set LD_LIBRARY_PATH /usr/local/lib:/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/lib:/usr/local/cuda32/cuda/lib'
else
INC='-I. -I.. -I/usr/local/cuda/include -I/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/inc -I/usr/local/include/gsl'
LIB='-lm -lgslcblas -lpthread -lcublas -lcudart -L/home/kameyama/lukacm/Eclipse_Epig_CUBLAS/src -L/usr/local/lib -L/usr/local/cuda/lib -L/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/lib'
QMDD='-m64 -I/home/local_space/QMDD'
CC='/usr/local/cuda/bin/nvcc'
echo 'Remember to set :D_LIBRARY_PATH /usr/local/lib:/usr/local/cuda/lib:/home/kameyama/lukacm/NVIDIA_CUDA_SDK/common/lib:$LD_LIBRARY_PATH'
fi

$CC ${INC} ${QMDD} -c -o tools.o tools.cc
$CC ${INC} ${QMDD} -c -o main.o main.cc
$CC ${INC} ${QMDD} -c -o cstGA.o cstGA.cc
$CC ${INC} ${QMDD} -c -o Qe.o Qe.cc
$CC ${INC} ${QMDD} -c -o ccuda_multi.o ccuda_multi.cu
$CC -g -O2 -o epig_fsm tools.o main.o cstGA.o Qe.o ccuda_multi.o $LIB
