
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
EXECUTABLE2 := preComp_Integrals_MC_Swave
EXECUTABLE3 := fit_recoMC_singleComponent
EXECUTABLE4 := fit_genMC
EXECUTABLE5 := fit_recoMC_fullAngular
EXECUTABLE6 := composeEff_rooKeys_parSub
EXECUTABLE7 := composeToyEff_rooKeys_parSub
EXECUTABLE8 := plotToyEff
EXECUTABLE9 := preComp_toyIntegrals_MC
EXECUTABLE10 := plotEffComparison
EXECUTABLE11 := extractEff

EXTRACLASS  := RooDataHist.cxx
EXTRACLASS2 := RooNDKeysPdf.cxx
CLASS0     := PdfRT
CLASS1     := PdfWT
CLASS2     := DecayRate
CLASS3     := PdfSigAng
CLASS4     := ShapeSigAng
CLASSDICT  := AngDict

#compiling options
DEBUGFLAGS := -O3 -Wall -std=c++11
CXXFLAGS := $(DEBUGFLAGS) 

#compile class
LIBS := $(SOURCEDIR)/$(CLASS0).cc $(SOURCEDIR)/$(CLASS1).cc $(SOURCEDIR)/$(CLASS2).cc $(SOURCEDIR)/$(CLASS3).cc $(SOURCEDIR)/$(CLASS4).cc $(CLASSDICT).cc $(SOURCEDIR)/$(EXTRACLASS)

$(CLASSDICT): $(INCLUDEDIR)/$(CLASS0).h $(INCLUDEDIR)/$(CLASS1).h $(INCLUDEDIR)/$(CLASS2).h $(INCLUDEDIR)/$(CLASS3).h $(INCLUDEDIR)/$(CLASS4).h
	@echo "Generating dictionary $@ using rootcint ..."
	$(ROOTCINT) -f $@.cc -c $^

$(EXECUTABLE0): $(EXECUTABLE0).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE1): $(EXECUTABLE1).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE2): $(EXECUTABLE2).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE3): $(EXECUTABLE3).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE4): $(EXECUTABLE4).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE5): $(EXECUTABLE5).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)

$(EXECUTABLE6): $(EXECUTABLE6).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS2) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE7): $(EXECUTABLE7).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(SOURCEDIR)/$(EXTRACLASS2) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE8): $(EXECUTABLE8).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE9): $(EXECUTABLE9).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE10): $(EXECUTABLE10).cc 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SOURCEDIR)/$(EXTRACLASS) $(ROOTLIBS) $(ROOTFLAGS)

$(EXECUTABLE11): $(EXECUTABLE11).cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) $(ROOTLIBS) $(ROOTFLAGS) -I$(INCLUDEDIR)



#cleaning options
.PHONY: clean
clean:
	rm -f $(EXECUTABLE0) $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3) $(EXECUTABLE4) $(EXECUTABLE5)
