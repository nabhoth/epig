SUFFIXES = .cu .cc .o .obj
bin_PROGRAMS = epig


if QMDD
epig_SOURCES=standard_multi.cc tools.cc main.cc cstGA.cc Qe.cc qmddfunc.cc #llist.cc
else 
if GSLBLAS
epig_SOURCES=standard_multi.cc tools.cc main.cc cstGA.cc Qe.cc #llist.cc
else
if CUDA
epig_SOURCES=ccuda_multi.cu tools.cc main.cc cstGA.cc Qe.cc# llist.cc
endif
endif
endif

CFLAGS=
.cu.o:  
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $<

