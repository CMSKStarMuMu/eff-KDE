#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore
ROOTCINT  := $(shell which rootcint)

#directories
SOURCEDIR  := ./src
INCLUDEDIR := ./interface

#exe_files
EXECUTABLE0 := plotEff
EXECUTABLE3 := fit_genMC
EXTRACLASS := RooDataHist.cxx
CLASS2     := DecayRate
CLASSDICT  := AngDict

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(SOURCEDIR)/$(CLASS2).cc $(CLASSDICT).cc $(SOURCEDIR)/$(EXTRACLASS)

$(CLASSDICT): $(INCLUDEDIR)/$(CLASS2).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE3): $(EXECUTABLE3).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)


#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE0) $(EXECUTABLE3)
