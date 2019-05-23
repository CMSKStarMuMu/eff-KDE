#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore
ROOTCINT  := $(shell which rootcint)

#directories
SOURCEDIR  := ./src

#exe_files
EXECUTABLE0 := plotEff
EXECUTABLE1 := preComp_Integrals_MC
EXTRACLASS  := RooDataHist.cxx

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE1): $(EXECUTABLE1).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE0) $(EXECUTABLE1)
