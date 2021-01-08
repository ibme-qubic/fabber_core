include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_core

USRINCFLAGS = -DFABBER_SRC_DIR="\"${PWD}\"" -DFABBER_BUILD_DIR="\"${PWD}\""
LIBS = -lfsl-newimage -lfsl-miscmaths -lfsl-utils \
       -lfsl-cprob -lfsl-NewNifti -lfsl-znz -ldl
TESTLIBS = -lgtest -lpthread

#
# Executables and libraries provided by this project
#

XFILES = fabber mvntool
SOFILES = libfsl-fabbercore.so libfsl-fabberexec.so
SCRIPTS = fabber_var

# Sets of objects separated into logical divisions

# Basic objects - things that have nothing directly to do with inference
BASICOBJS = tools.o rundata.o dist_mvn.o easylog.o fabber_capi.o version.o dist_gamma.o rundata_array.o

# Core objects - things that implement the framework for inference
COREOBJS =  noisemodel.o fwdmodel.o inference.o fwdmodel_linear.o fwdmodel_poly.o convergence.o motioncorr.o priors.o transforms.o

# Infernce methods
INFERENCEOBJS = inference_vb.o inference_nlls.o covariance_cache.o

# Noise models
NOISEOBJS = noisemodel_white.o noisemodel_ar.o

# Configuration
CONFIGOBJS = setup.o factories.o

# Library for executables
EXECOBJS = rundata_newimage.o fabber_core.o

# Executable main object
CLIENTOBJS =  fabber_main.o

# Unit tests
TESTOBJS = test/fabbertest.o test/test_inference.o test/test_priors.o test/test_vb.o test/test_convergence.o test/test_commandline.o test/test_rundata.o

# Everything together
OBJS = ${BASICOBJS} ${COREOBJS} ${INFERENCEOBJS} ${NOISEOBJS} ${CONFIGOBJS}

# For debugging:
#OPTFLAGS = -ggdb -Wall

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

# Targets

all: ${XFILES} ${SOFILES}

clean:
	${RM} -f /tmp/fslgrot *.o mvn_tool/*.o *.a *.so *.exe core depend.mk fabber_test

mvntool: ${OBJS} mvn_tool/mvntool.o rundata_newimage.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

# Build a fabber exectuable, this will have nothing but the generic models so it not practically useful for data analysis
fabber: ${OBJS} ${EXECOBJS} ${CLIENTOBJS}
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

# Library build
libfsl-fabbercore.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

libfsl-fabberexec.so : ${EXECOBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

# Unit tests
test: ${OBJS} ${EXECOBJS} ${CLIENTOBJS} ${TESTOBJS}
	${CXX} ${CXXFLAGS} -o fabber_test $^ ${LDFLAGS} ${TESTLIBS}

# DO NOT DELETE
