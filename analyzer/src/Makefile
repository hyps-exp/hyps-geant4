##	Makefile
### For Using CLHEP
G4INCLUDE = $(shell geant4-config --cflags)
G4LIBS    = $(shell geant4-config --libs)

### ROOT
RTINCLUDE       = $(shell root-config --cflags)
RTLIBS          = $(shell root-config --libs)

CFLAGS      = -g -O2 -std=c++17

# debug_flags += -DMemoryLeak

FLAGS = $(CFLAGS) $(DFLAGS) -I. $(G4INCLUDE) $(RTINCLUDE)

# $(debug_flags)

LIBOBJ = libSksAna.a
SKSLIB = -L. -lSksAna
LIBS = $(G4LIBS) $(RTLIBS)
LIBSONL = $(G4LIBS) $(RTLIBS)

C++ = g++

BINDIR = ../bin

SRCS    = $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

OBJ =   ConfMan.o RawData.o \
	HodoRawHit.o HodoAnalyzer.o HodoParamMan.o HodoPHCMan.o \
	Hodo2Hit.o Hodo1Hit.o\
	BH2Hit.o FiberHit.o FiberCluster.o\
	HodoCluster.o BH2Cluster.o \
	DCRawHit.o DCHit.o DCAnalyzer.o DCLocalTrack.o DCLTrackHit.o\
	TrackHit.o MWPCCluster.o TrackMaker.o BH1Filter.o BH2Filter.o\
	K18Track.o K18Track_BFT.o K18TrackFCN.o K18TransMatrix.o K18Particle.o Minuit.o\
	SksTrack.o SksParticle.o RungeKuttaUtilities.o MathTools.o\
	FieldMan.o SksFieldMap.o FieldElements.o \
	DCGeomMan.o DCGeomRecord.o DCTdcCalibMan.o DCDriftParamMan.o \
	DCTrackSearch.o DCPairHitCluster.o\
	Kinematics.o DCTrackUtilities.o \
        DebugTimer.o DebugCounter.o \
	FileManager.o CFTFiberHit.o \
	CFTFiberCluster.o CFTLocalTrack.o EnergyCorrection.o EvDispCFT.o \
	CFTParticle.o SksParticle.o MatrixTrigMan.o CFTTrackingEffMan.o
# 	ScalerAna.o\
# 	DSio.o CFile.o EvDisp.o   		\
# 	Sdc1HitCluster.o  D           \
# 	SimuData.o 	\
#

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: analysCFT analysKURAMA analysSigmaMinusProtonScat analysSigmaPlusProtonScat analysPPScat analysLambdaProtonScat analysXiMinusProtonScat

.cc.o:
	$(C++) $(FLAGS) -MMD -c $<


analysCFT: lib UserCFT.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserCFT.o 	\
	$(SKSLIB) $(LIBS)

analysKURAMA: lib UserKURAMA.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserKURAMA.o 	\
	$(SKSLIB) $(LIBS)

analysSigmaMinusProtonScat: lib UserSigmaMinusProtonScat.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSigmaMinusProtonScat.o 	\
	$(SKSLIB) $(LIBS)

analysSigmaPlusProtonScat: lib UserSigmaPlusProtonScat.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserSigmaPlusProtonScat.o 	\
	$(SKSLIB) $(LIBS)

analysPPScat: lib UserPPScat.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserPPScat.o 	\
	$(SKSLIB) $(LIBS)

analysLambdaProtonScat: lib UserLambdaProtonScat.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserLambdaProtonScat.o 	\
	$(SKSLIB) $(LIBS)

analysXiMinusProtonScat: lib UserXiMinusProtonScat.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserXiMinusProtonScat.o 	\
	$(SKSLIB) $(LIBS)


clean:
	rm -f *.o *.d *~ core* libSksAna.a
