PACKAGE=epig
VERSION=0.1.0
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE



AC_PROG_CC
AC_PROG_CPP
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_MKDIR_P

AC_SUBST(CFLAGS, "")
AC_SUBST(CXXFLAGS, "")
AC_SUBST(LDFLAGS, "-L./")
AC_SUBST(BITAGE, "64")

AC_MSG_RESULT(${BITAGE});


dnl #CUDA_REQUIRED_VERSION=1.0.7
dnl #PKG_CHECK_MODULES(CUDA, cuda-1.0 >= $CUDA_REQUIRED_VERSION)
dnl #AC_SUBST(INTI_CFLAGS)
dnl #AC_SUBST(INTI_LIBS)
dnl #CXX='/usr/local/cuda32/cuda/bin/nvcc'

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lguide':
AC_CHECK_LIB([guide], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])



AC_MSG_CHECKING(for main in -lcudart);
AC_ARG_WITH(
	cudart,
	AC_HELP_STRING(
		[--with-cudart=ARG],
		[specify cuda location]
	),
	[
		wants_CUDA="yes"
		CUDA_PREFIX=$with_cudart
		if test ! -f "${CUDA_PREFIX}/lib/libcudart.so"; then
					AC_MSG_FAILURE([path to CUDA is not correct])
		fi
	],
	[wants_CUDA="no"]
)

if test "$wants_CUDA" = "yes" ; then
AC_SUBST(CUDA_PREFIX)
AC_SUBST(CXX,"${CUDA_PREFIX}/bin/nvcc")
AC_PROG_CXX
CUDA_LIBS=" -L${CUDA_PREFIX}/lib/cuda -lcudart"
CUDA_CFLAGS="-D__CUDA__ -D__CUDA_SDK__ -I${CUDA_PREFIX}/include"
LDFLAGS+=" $CUDA_LIBS"
CXXFLAGS+=" $CUDA_CFLAGS"
AC_SUBST(CUDA_LIBS)
AC_SUBST(CUDA_CFLAGS)
#else AC_MSG_FAILURE([CUDA is required when using EpiG])
fi
AC_MSG_RESULT($wants_CUDA);
AM_CONDITIONAL(CUDA, test "$wants_CUDA" = "yes")


#AC_MSG_CHECKING(for main in cudasdk);
#AC_ARG_WITH(
#	cudasdk,
#	AC_HELP_STRING(
#		[--with-cudasdk=ARG],
#		[specify cudasdk location]
#	),
#	[
#		wants_CUDASDK="yes"
#		CUDASDK_PREFIX=$with_cudasdk
#		if test ! -f "${CUDASDK_PREFIX}/common/inc/cutil.h"; then
#					AC_MSG_FAILURE([path to CUDASDK is not correct])
#					AC_MSG_FAILURE([CUDASDK is required when using CUDA])
#		fi
#	],
#	[wants_CUDASDK="no"]
#)
#
#
#if test "$wants_CUDASDK" = "yes" ; then
#AC_SUBST(CUDASDK_PREFIX)
#CUDASDK_LIBS="-L${CUDASDK_PREFIX}/common/lib"
#CUDASDK_CFLAGS="-D__CUDA__ -D__CUDA_SDK__ -I${CUDASDK_PREFIX}/common/inc"
#LDFLAGS+=" $CUDASDK_LIBS"
#CXXFLAGS+=" $CUDASDK_CFLAGS"
#AC_SUBST(CUDASDK_LIBS)
#AC_SUBST(CUDASDK_CFLAGS)
#fi
#AC_MSG_RESULT($wants_CUDASDK);
#AM_CONDITIONAL(CUDASDK, test "$wants_CUDASDK" = "yes")

AC_MSG_CHECKING(for main in -lgslcblas);
AC_ARG_WITH(
	gslcblas,
	AC_HELP_STRING(
		[--with-gslcblas=ARG],
		[specify gslcblas location]
	),
	[
		GSLBLAS_PREFIX=$with_gslcblas
		#GSLBLAS_PREFIX="$withval"
		wants_GSLBLAS="yes"
		if test ! -f "${GSLBLAS_PREFIX}/lib64/libgsl.so"; then
					AC_MSG_FAILURE([path to GSLBLAS is not correct])
		fi
	],
	[wants_GSLBLAS="no"]
)
if test "$wants_GSLBLAS" = "yes" ; then
AC_SUBST(GSLBLAS_PREFIX)
GSLBLAS_LIBS="-L${GSLBLAS_PREFIX}/lib64 -lgslcblas -lgsl"
AC_PROG_CXX
if test "$wants_CUDASDK" != "yes" ; then
GSLBLAS_CFLAGS="-D__GSLBLAS__ -I${GSLBLAS_PREFIX}/include/gsl"
else
GSLBLAS_CFLAGS="-I${GSLBLAS_PREFIX}/include/gsl"
fi
LDFLAGS+=" $GSLBLAS_LIBS"
CXXFLAGS+=" $GSLBLAS_CFLAGS"
AC_SUBST(GSLBLAS_LIBS)
AC_SUBST(GSLBLAS_CFLAGS)
fi
if test "$wants_CUDASDK" != "yes" ; then
AC_MSG_RESULT($wants_GSLBLAS);
AM_CONDITIONAL(GSLBLAS, test "$wants_GSLBLAS" = "yes")
else
AC_MSG_RESULT($wants_GSLBLAS);
AM_CONDITIONAL(GSLBLAS, test "$wants_GSLBLAS" != "yes")
fi
AC_MSG_CHECKING(for main in -lmpich);
AC_ARG_WITH(
	mpi,
	AC_HELP_STRING(
		[--with-mpi=ARG],
		[specify mpi location]
	),
	[
		wants_MPI="yes"
		MPI_PREFIX=$with_mpi
		if test ! -f "${MPI_PREFIX}/lib/libmpich.so"; then
					AC_MSG_FAILURE([path to MPI is not correct])
		fi
	],
	[wants_MPI="no"]
)
if test "$wants_MPI" = "yes" ; then
AC_SUBST(MPI_PREFIX)
AC_SUBST(CXX,"${MPI_PREFIX}bin/mpicc")
AC_PROG_CXX
MPI_LIBS="-L${MPI_PREFIX}/lib -limpich"
MPI_CFLAGS="-D__MPI__ -I${MPI_PREFIX}/include"
LDFLAGS+=" $MPI_LIBS"
CXXFLAGS+=" $MPI_CFLAGS"
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_CFLAGS)
fi
AC_MSG_RESULT($wants_MPI);
AM_CONDITIONAL(MPI, test "$wants_MPI" = "yes")

AC_MSG_CHECKING(for main in qmdd);
AC_ARG_WITH(
	qmdd,
	AC_HELP_STRING(
		[--with-qmdd=ARG],
		[specify qmdd location]
	),
	[
		wants_QMDD="yes"
		QMDD_PREFIX=$with_qmdd
		if test ! -f "${QMDD_PREFIX}/QMDDinclude.h"; then
					AC_MSG_FAILURE([path to QMDD is not correct])
		fi
	],
	[wants_QMDD="no"]
)
if test "$wants_QMDD" = "yes" ; then
AC_SUBST(QMDD_PREFIX)
QMDD_LIBS="-L${QMDD_PREFIX}"
QMDD_CFLAGS="-D__QMDD__ -m32 -I${QMDD_PREFIX}"
LDFLAGS+=" $QMDD_LIBS"
CXXFLAGS+=" $QMDD_CFLAGS"
AC_SUBST(QMDD_LIBS)
AC_SUBST(QMDD_CFLAGS)
else
CXXFLAGS+=" -m${BITAGE}"
fi
AC_MSG_RESULT($wants_QMDD);
AM_CONDITIONAL(QMDD, test "$wants_QMDD" = "yes")

AC_MSG_CHECKING(for main in measurement);
AC_ARG_WITH(measurement,
	AC_HELP_STRING(
		[--with-measurement],
	        [ specifies if measurement is to be used during evaluation of quantum circuit]
	),
    [wants_MSRMNT=yes],
    [wants_MSRMNT=no]
)

if test "$wants_MSRMNT" = "yes" ; then
ALGOPARAMS+=" -D__MSRMNT__"
fi
AC_MSG_RESULT($wants_MSRMNT);
AM_CONDITIONAL(MSRMNT, test "$wants_MSRMNT" = "yes")

AC_MSG_CHECKING(for main in sequence);
AC_ARG_WITH(sequence,
	AC_HELP_STRING(
		[--with-sequence],
	        [ specifies if sequence detector is to be evolved]
	),
    [wants_SQNC=yes],
    [wants_SQNC=no]
)

if test "$wants_SQNC" = "yes" ; then
ALGOPARAMS+=" -D__SQNC__"
fi
AC_MSG_RESULT($wants_SQNC);
AM_CONDITIONAL(SQNC, test "$wants_SQNC" = "yes")

AC_MSG_CHECKING(for main in polynomial);
AC_ARG_WITH(polynomial,
	AC_HELP_STRING(
		[--with-polynomial],
	        [ specifies if we are evolving quantum entangled states specified by a polynomial]
	),
    [wants_PLNM=yes],
    [wants_PLNM=no]
)

if test "$wants_PLNM" = "yes" ; then
ALGOPARAMS+=" -D__PLNM__"
fi
AC_MSG_RESULT($wants_PLNM);
AM_CONDITIONAL(PLNM, test "$wants_PLNM" = "yes")

AC_MSG_CHECKING(for main in elemnt);
AC_ARG_WITH(elemnt,
	AC_HELP_STRING(
		[--with-elem],
		[ specifies if we the evaluation is using element-to-element comparison for determining the fitness value (yes is default)]
	),
    [wants_ELEM=no],
    [wants_ELEM=yes]
)

if test "$wants_MSRMNT" = "no" ; then
wants_ELEM=yes
ALGOPARAMS+=" -D__ELEM__"
AC_MSG_RESULT($wants_ELEM);
AM_CONDITIONAL(ELEM, test "$wants_ELEM" = "yes")
else
wants_ELEM=no
AC_MSG_RESULT($wants_ELEM);
AM_CONDITIONAL(ELEM, test "$wants_ELEM" = "no")
fi



#####STDR - NOT Working now
#AC_MSG_CHECKING(for main in stdr);
#AC_ARG_WITH(stdr,
#	AC_HELP_STRING(
#		[--with-stdr],
#	        [ specifies if the standard compilation is required (enough of memory)]
#	),
#    [wants_STDR=yes],
#    [wants_STDR=no]
#)
#
#if test "$wants_STDR" = "yes" ; then
#CXXFLAGS+=" -D__STDR__"
#fi
#AC_MSG_RESULT($wants_STDR);
#AM_CONDITIONAL(STDR, test "$wants_STDR" = "yes")


CXXFLAGS+=" $ALGOPARAMS"

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([stdint.h stdlib.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_PID_T
AC_STRUCT_TM
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([getpagesize pow sqrt])

AC_CONFIG_FILES([Makefile
		src/Makefile])
AC_OUTPUT
