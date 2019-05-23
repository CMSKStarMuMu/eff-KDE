#root_stuff (root libraries and needed root options)
ROOTLIBS  := $(shell root-config --glibs)
ROOTFLAGS := $(shell root-config --cflags --libs) -lRooFit -lRooFitCore -lMathMore
ROOTCINT  := $(shell which rootcint)

#directories
SOURCEDIR  := ./src
INCLUDEDIR := ./interface

#exe_files
EXECUTABLE0 := plotEff
EXECUTABLE1 := preComp_Integrals_MC
EXECUTABLE2 := fit_recoMC_singleComponent

EXTRACLASS := RooDataHist.cxx
CLASS0     := PdfRT
CLASS1     := PdfWT
CLASSDICT  := AngDict

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(SOURCEDIR)/$(CLASS0).cc $(SOURCEDIR)/$(CLASS1).cc $(CLASSDICT).cc $(SOURCEDIR)/$(EXTRACLASS)

$(CLASSDICT): $(INCLUDEDIR)/$(CLASS0).h $(INCLUDEDIR)/$(CLASS1).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE1): $(EXECUTABLE1).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE2): $(EXECUTABLE2).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)


#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE0) $(EXECUTABLE1)
