CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA

all: selection.exe

selection.exe: selection.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTLIBS) selection.cc -o selection.exe
