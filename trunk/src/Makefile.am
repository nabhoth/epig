SUFFIXES = .cu .cc .o .obj
bin_PROGRAMS = epig

if MSRMNT
epig_SPECS=cuda_msrmnt.cu
else
if SQNC
epig_SPECS=cuda_sqncs.cu
else
if ELEM
epig_SPECS=cuda_stdnrt.cu
endif
endif
endif

if QMDD
epig_SOURCES= $(epig_SPECS) standard_multi.cc tools.cc main.cc cstGA.cc Qe.cc qmddfunc.cc
else
if CUDA
epig_SOURCES= $(epig_SPECS) ccuda_multi.cu cuda_tools.cc main.cc cstGA.cc Qe.cc
else 
if GSLBLAS
epig_SOURCES= $(epig_SPECS) standard_multi.cc tools.cc main.cc cstGA.cc Qe.cc
endif
endif
endif

CFLAGS=
.cu.o:  
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $<

