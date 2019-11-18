#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>

#include "PdfRT.h"
#include "PdfWT.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* c [4*nBins];


void simfit_recoMC_singleComponentBin(int q2Bin, int parity, int tagFlag, bool plot, bool save, std::vector<int> years)
{

  string shortString = Form("b%ip%it%i",q2Bin,parity,tagFlag);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency

  string datasetString = (tagFlag==1?"data_ctRECO":"data_wtRECO");
  datasetString = datasetString + Form((parity==1?"_ev_b%i":"_od_b%i"),q2Bin);
  string effString = Form((tagFlag==0?"effWHist_b%ip%i":"effCHist_b%ip%i"),q2Bin,parity);
  string intHistString = "MCint_"+shortString;
  string all_years = "";
  string year = ""; 

//   std::vector<string> filename_data;
  std::vector<TFile*> fin_data, fin_eff;
  std::vector<RooWorkspace*> wsp;
  std::vector<RooDataSet*> data;
  std::vector<RooAbsReal*> eff;
  std::vector<TH3D*> effHist;
  std::vector<TH1D*> intHist;
  std::vector< std::vector<double> > intVec(years.size(), std::vector<double>(0));
  std::vector<RooAbsPdf*> PDF_sig_ang_singleComponent (0);

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;


  RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl  = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1  = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2  = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3  = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    sample.defineType(("data"+year).c_str());
    all_years += year;
  }
  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/effDataset_b%i.root",years[iy], q2Bin); 

    // import data (or MC as data proxy)
    fin_data.push_back( TFile::Open( filename_data.c_str() ) );
    if ( !fin_data[iy] || !fin_data[iy]->IsOpen() ) {
      cout << "File not found: " << filename_data << endl;
      return;
    }
  
    wsp.push_back( (RooWorkspace*)fin_data[iy]->Get(Form("ws_b%ip%i", q2Bin, 1-parity ) ) );
    if ( !wsp[iy] || wsp[iy]->IsZombie() ) {
      cout<<"Workspace not found in file: "<<filename_data<<endl;
      return;
    }
  
    data.push_back( (RooDataSet*)wsp[iy]->data(datasetString.c_str())) ;
    if ( !data[iy] || data[iy]->IsZombie() ) {
      cout<<"Dataset "<<datasetString<<" not found in file: "<<filename_data<<endl;
      return;
    }

    // import KDE efficiency histograms and partial integral histograms
    string filename = "files_forSimFit/KDEeff_b";
    filename = filename + Form((parity==0 ? "%i_ev_%i.root" : "%i_od_%i.root"),q2Bin,years[iy]);
    fin_eff.push_back( new TFile( filename.c_str(), "READ" ));
    if ( !fin_eff[iy] || !fin_eff[iy]->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      return;
    }

    effHist.push_back( (TH3D*)fin_eff[iy]->Get(effString.c_str()));
    if ( !effHist[iy] || effHist[iy]->IsZombie() ) {
      cout<<"Efficiency histogram "<< effString <<" not found in file: "<< filename <<endl;
      return;
    }

    // create efficiency functions
    RooDataHist* effData = new RooDataHist(("effData_"+shortString+"_"+year).c_str(),"effData",vars,effHist[iy]);
    eff.push_back( new RooHistFunc(("eff_"+shortString+"_"+year).c_str(),
                                   ("eff"+year).c_str() ,
                                   vars,
                                   *effData,
                                   1));

    // import precomputed integrals and fill a std::vector
    intHist.push_back( (TH1D*)fin_eff[iy]->Get(intHistString.c_str()));
    intVec.push_back (vector<double> (0));
    if ( !intHist[iy] || intHist[iy]->IsZombie() ) {
      cout << "Integral histogram " << intHistString <<" not found in file: "<< filename << endl << "Abort" << endl;
      return;
    } else if ( strcmp( intHist[iy]->GetTitle(), effHist[iy]->GetTitle() ) ) {
    // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
      cout << "Integral histogram "<< intHistString << " is incoherent with efficiency " << effString << " in file: " << filename << endl;
      cout << "Efficiency conf: "  << effHist[iy]->GetTitle() << endl;
      cout << "Integral conf: "    << intHist[iy]->GetTitle() << endl << "Abort"<<endl;
      return;
    } 
    else {
      for (int i=1; i<=intHist[iy]->GetNbinsX(); ++i) {
        intVec[iy].push_back(intHist[iy]->GetBinContent(i));
      }
    }

    // define angular PDF for signal and only one tag component, using the custom class
    // efficiency function and integral values are passed as arguments
    if (tagFlag==1) PDF_sig_ang_singleComponent.push_back( new PdfRT(("PDF_sig_ang_singleComponent_"+shortString+"_"+year).c_str(),
                                                                     ("PDF_sig_ang_singleComponent_"+year).c_str(),
      							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p, *eff[iy], intVec[iy]));
    else            PDF_sig_ang_singleComponent.push_back( new PdfWT(("PDF_sig_ang_singleComponent_"+shortString+"_"+year).c_str(),
                                                                     ("PDF_sig_ang_singleComponent_"+year).c_str(),
  							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p, *eff[iy], intVec[iy]));
    
    // insert sample in the category map, to be imported in the combined dataset
    map.insert(map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year).c_str(), data[iy]));
    // associate model with the data
    simPdf->addPdf(*PDF_sig_ang_singleComponent[iy], ("data"+year).c_str());
  }

  for(auto it = map.cbegin(); it != map.cend(); ++it)
    std::cout << "dataset: " << it->first << ", with n entries: " << it->second->sumEntries() << "\n";

  // to start the fit, parameters are restored to the center of the parameter space
  Fl ->setVal(0.5);
  P1 ->setVal(0);
  P2 ->setVal(0);
  P3 ->setVal(0);
  P4p->setVal(0);
  P5p->setVal(0);
  P6p->setVal(0);
  P8p->setVal(0);

  // Construct combined dataset in (x,sample)
  RooDataSet combData ("combData", "combined data", 
                         vars,
                         Index(sample),
                         Import(map)); 


  // perform fit in two steps:
  // first with strategy=0 and no MINOS, to get close to the likilihood maximum
  simPdf->fitTo( combData,
                 Minimizer("Minuit2","migrad"), 
                 Extended(false), 
                 Timer(true),
                 NumCPU(1),
                 Hesse(true),
                 Strategy(0),
                 Offset(true));
  // second with full accuracy (strategy=2) and MINOS, to get the best result
  RooFitResult* fitResult = simPdf->fitTo( combData,
                                           Minimizer("Minuit2","migrad"),
                                           Extended(false), 
                                           Save(true),
//                                            Timer(true),
                                           NumCPU(1),
                                           Hesse(true),
                                           Strategy(2),
                                           Minos(true),
                                           Offset(true)
                                         );
  // The two step fit seemed necessary to get convergence in all q2 bins
  // if one observes that the convergence is obtained just running the second step, it is better to avoid the first step (to gain computing time)
  // Offset(true) parameter make the FCN value to be defined at 0 at the begin of the minimization
  // it is needed to get convergence in all q2 bins, and avoid the "machine precision reached" errors
  // which are due to too high FCN absolute values compared to the minimisation steps

  fitResult->Print("v");

  if (save) {
    // Save fit results in file
    TFile* fout = new TFile(("simFitResults/simFitResult_recoMC_singleComponent" + all_years + ".root").c_str(),"UPDATE");
    fitResult->Write(("simFitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  int confIndex = 2*nBins*parity + nBins*tagFlag + q2Bin;
  string longString  = (tagFlag==1?"Fit to correctly tagged events":"Fit to wrongly tagged events");
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,700);
  c[confIndex]->Divide(3, years.size());
  
   for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    RooPlot* xframe = ctK->frame(Title((longString+year).c_str()));
    RooPlot* yframe = ctL->frame(Title((longString+year).c_str()));
    RooPlot* zframe = phi->frame(Title((longString+year).c_str()));
    xframe->GetYaxis()->SetTitleOffset(1.8);
    yframe->GetYaxis()->SetTitleOffset(1.8);
    zframe->GetYaxis()->SetTitleOffset(1.8);
    xframe->SetMaximum(xframe->GetMaximum()*1.15);
    yframe->SetMaximum(yframe->GetMaximum()*1.15);
    zframe->SetMaximum(zframe->GetMaximum()*1.15);
    xframe->SetMinimum(0);
    yframe->SetMinimum(0);
    zframe->SetMinimum(0);
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);

    combData.plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year).c_str()), Name(("plData"+year).c_str()));
    combData.plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year).c_str()));
    combData.plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year).c_str()));

    simPdf->plotOn(xframe,Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1),Name(("plPDF"+year).c_str()));
    simPdf->plotOn(yframe,Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1));
    simPdf->plotOn(zframe,Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1));

    c[confIndex]->cd(iy*3+1);
    gPad->SetLeftMargin(0.19); 
    xframe->Draw();
    leg->Draw("same");
    c[confIndex]->cd(iy*3+2);
    gPad->SetLeftMargin(0.19); 
    yframe->Draw();
    leg->Draw("same");
    c[confIndex]->cd(iy*3+3);
    gPad->SetLeftMargin(0.19); 
    zframe->Draw();
    leg->SetTextSize(0.03);
    leg->AddEntry(xframe->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
    leg->AddEntry(xframe->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
    leg->Draw("same");


  }

  c[confIndex]->SaveAs( ("plotSimFit_d/simFitResult_recoMC_singleComponent_" + shortString + "_" + all_years + ".pdf").c_str() );
}


void simfit_recoMC_singleComponentBin2(int q2Bin, int parity, int tagFlag, bool plot, bool save, std::vector<int> years)
{
  if ( tagFlag==-1 )
    for (tagFlag=0; tagFlag<2; ++tagFlag)
      simfit_recoMC_singleComponentBin(q2Bin, parity, tagFlag, plot, save, years);
  else
    simfit_recoMC_singleComponentBin(q2Bin, parity, tagFlag, plot, save, years);
}

void simfit_recoMC_singleComponentBin1(int q2Bin, int parity, int tagFlag, bool plot, bool save, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_singleComponentBin2(q2Bin, parity, tagFlag, plot, save, years);
  else
    simfit_recoMC_singleComponentBin2(q2Bin, parity, tagFlag, plot, save, years);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively
  // tag format: [0] mistagged
  //             [1] correctly-tagged
  //             [-1] for each tag recursively

  int q2Bin   = -1;
  int parity  = -1; 
  int tagFlag = -1;

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);
  if ( argc >= 4 ) tagFlag = atoi(argv[3]);

  bool plot = true;
  bool save = true;

  if ( argc >= 5 && atoi(argv[4]) == 0 ) plot = false;
  if ( argc >= 6 && atoi(argv[5]) == 0 ) save = false;

  std::vector<int> years;
  if ( argc >= 7 && atoi(argv[6]) != 0 ) years.push_back(atoi(argv[6]));
  else {
    cout << " no specific years selected, using default: 2016";
    years.push_back(2016);
  }
  if ( argc >= 8 && atoi(argv[7]) != 0 ) years.push_back(atoi(argv[7]));
  if ( argc >= 9 && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;
  if ( tagFlag < -1 || tagFlag > 1      ) return 1;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;
  if ( tagFlag==-1 ) cout<<"Running both the tag conditions"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_recoMC_singleComponentBin1(q2Bin, parity, tagFlag, plot, save, years);
  else
    simfit_recoMC_singleComponentBin1(q2Bin, parity, tagFlag, plot, save, years);

  return 0;

}
