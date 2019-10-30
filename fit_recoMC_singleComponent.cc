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
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooNumIntConfig.h>

#include "PdfRT.h"
#include "PdfWT.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* c [4*nBins];

void fit_recoMC_singleComponentBin(int q2Bin, int parity, int tagFlag, bool plot, bool save)
{

  string shortString = Form("b%ip%it%i",q2Bin,parity,tagFlag);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency
  string filename_data = Form("effDataset_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,1-parity)); //complementary statistics
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename_data<<endl;
    return;
  }
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
    cout<<"Variables not found in file: "<<filename_data<<endl;
    return;
  }
  RooArgList vars (* ctK,* ctL,* phi);
  string datasetString = (tagFlag==1?"data_ctRECO":"data_wtRECO");
  datasetString = datasetString + Form((parity==1?"_ev_b%i":"_od_b%i"),q2Bin);
  RooDataSet* data = (RooDataSet*)wsp->data(datasetString.c_str());
  if ( !data || data->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename_data<<endl;
    return;
  }

  // import KDE efficiency histograms and partial integral histograms
  string filename = "files/KDEeff_b";
  filename = filename + Form((parity==0?"%i_ev.root":"%i_od.root"),q2Bin);
  TFile* fin = new TFile( filename.c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  string effString = Form((tagFlag==0?"effWHist_b%ip%i":"effCHist_b%ip%i"),q2Bin,parity);
  TH3D* effHist = (TH3D*)fin->Get(effString.c_str());
  if ( !effHist || effHist->IsZombie() ) {
    cout<<"Efficiency histogram "<<effString<<" not found in file: "<<filename<<endl;
    return;
  }
  // create efficiency functions
  RooDataHist* effData = new RooDataHist(("effData_"+shortString).c_str(),"effData",vars,effHist);
  RooAbsReal* eff = new RooHistFunc(("eff_"+shortString).c_str(),"eff",vars,*effData,1);
  // import precomputed integrals and fill a std::vector
  vector<double> intVec (0);
  string intHistString = "MCint_"+shortString;
  TH1D* intHist = (TH1D*)fin->Get(intHistString.c_str());
  if ( !intHist || intHist->IsZombie() ) {
    // cout<<"Integral histogram "<<intHistString<<" not found in file: "<<filename<<endl<<"Using rooFit integration"<<endl;
    // intVec.push_back(0);
    cout<<"Integral histogram "<<intHistString<<" not found in file: "<<filename<<endl<<"Abort"<<endl;
    return;
  } else if ( strcmp( intHist->GetTitle(), effHist->GetTitle() ) ) {
    // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
    // (usually when efficiencies have been reproduced and overwritten, but the integrals are still referring to the old version)
    cout<<"Integral histogram "<<intHistString<<" is incoherent with efficiency "<<effString<<" in file: "<<filename<<endl;
    cout<<"Efficiency conf: "<<effHist->GetTitle()<<endl;
    cout<<"Integral conf: "<<intHist->GetTitle()<<endl;
    // cout<<"Using rooFit integration"<<endl;
    // intVec.push_back(0);
    cout<<"Abort"<<endl;
    return;
  } else {
    for (int i=1; i<=intHist->GetNbinsX(); ++i) intVec.push_back(intHist->GetBinContent(i));
  }

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

  // define angular PDF for signal and only one tag component, using the custom class
  // efficiency function and integral values are passed as arguments
  RooAbsPdf* PDF_sig_ang_singleComponent = 0;
  if (tagFlag==1) PDF_sig_ang_singleComponent = new PdfRT(("PDF_sig_ang_singleComponent_"+shortString).c_str(),"PDF_sig_ang_singleComponent",
							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*eff,intVec);
  else            PDF_sig_ang_singleComponent = new PdfWT(("PDF_sig_ang_singleComponent_"+shortString).c_str(),"PDF_sig_ang_singleComponent",
							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*eff,intVec);

  // create a variable to compute negative log-likelihood of the PDF on the dataset
  // this will be computed with angular paramenters at three states:
  // results of the fit to GEN sample, center of the parameter space (initial fit condition), and best-fit values
  RooAbsReal* nll = PDF_sig_ang_singleComponent->createNLL(*data);
  double nllGen = 0;

  // import results of the fit to GEN sample
  // if it is not found, nllGen is left at default value
  TFile* finGen = TFile::Open("fitResults/fitResult_genMC.root");
  if ( !finGen || finGen->IsZombie() ) {
    cout<<"Missing gen file: fitResults/fitResult_genMC.root"<<endl;
  } else {
    RooFitResult* fitResultGen = (RooFitResult*)finGen->Get(Form("fitResult_b%ip%i",q2Bin,parity));
    if ( !fitResultGen || fitResultGen->IsZombie() ) {
      cout<<"No "<<Form("fitResult_b%ip%i",q2Bin,parity)<<" in file: fitResults/fitResult_genMC.root"<<endl;
    } else {
      // set parameters to the fitResult values
      Fl ->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("Fl"))->getValV());
      P1 ->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P1"))->getValV());
      P2 ->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P2"))->getValV());
      P3 ->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P3"))->getValV());
      P4p->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P4p"))->getValV());
      P5p->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P5p"))->getValV());
      P6p->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P6p"))->getValV());
      P8p->setVal(((RooRealVar*)fitResultGen->floatParsFinal().find("P8p"))->getValV());
      // compute NLL
      nllGen = nll->getVal();
    }
  }

  // to start the fit, parameters are restored to the center of the parameter space
  Fl ->setVal(0.5);
  P1 ->setVal(0);
  P2 ->setVal(0);
  P3 ->setVal(0);
  P4p->setVal(0);
  P5p->setVal(0);
  P6p->setVal(0);
  P8p->setVal(0);

  // fit initial NLL value
  double nllCenter = nll->getVal();

  // Customise rooFit integrator (here used only to plot fit projections)
  // PDF_sig_ang_singleComponent->getIntegratorConfig()->methodND().setLabel("RooMCIntegrator");

  // perform fit in two steps:
  // first with strategy=0 and no MINOS, to get close to the likilihood maximum
  PDF_sig_ang_singleComponent->fitTo(*data,Minimizer("Minuit2","migrad"),Timer(true),NumCPU(1),Hesse(true),Strategy(0),Offset(true));
  // second with full accuracy (strategy=2) and MINOS, to get the best result
  RooFitResult* fitResult = PDF_sig_ang_singleComponent->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true),NumCPU(1),Hesse(true),Strategy(2),Minos(true),Offset(true));
  // The two step fit seemed necessary to get convergence in all q2 bins
  // if one observes that the convergence is obtained just running the second step, it is better to avoid the first step (to gain computing time)
  // Offset(true) parameter make the FCN value to be defined at 0 at the begin of the minimization
  // it is needed to get convergence in all q2 bins, and avoid the "machine precision reached" errors
  // which are due to too high FCN absolute values compared to the minimisation steps

  fitResult->Print("v");

  // Get NLL at best fit conditions and (if GEN fitResult was found) save the GEN-RECO deltaNLL value in the fitResult title
  double nllReco = nll->getVal();
  fitResult->SetTitle(Form("deltaNLL=%.1f",(nllGen==0?0.0:nllGen-nllReco)));

  // Print NLL values and difference
  cout<<"NLL values: center="<<nllCenter<<" result="<<nllReco<<" gen="<<nllGen<<endl;
  if (nllGen!=0) cout<<"GEN-RECO deltaNLL="<<nllGen-nllReco<<endl;

  if (save) {
    // Save fit results in file
    TFile* fout = new TFile("fitResults/fitResult_recoMC_singleComponent.root","UPDATE");
    fitResult->Write(("fitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  int confIndex = 2*nBins*parity + nBins*tagFlag + q2Bin;
  string longString  = (tagFlag==1?"Fit to correctly tagged events":"Fit to wrongly tagged events");
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title(longString.c_str()));
  RooPlot* yframe = ctL->frame(Title(longString.c_str()));
  RooPlot* zframe = phi->frame(Title(longString.c_str()));
  data->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40),Name("plData"));
  data->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40));
  data->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40));
  PDF_sig_ang_singleComponent->plotOn(xframe,LineWidth(1),Name("plPDF"));
  PDF_sig_ang_singleComponent->plotOn(yframe,LineWidth(1));
  PDF_sig_ang_singleComponent->plotOn(zframe,LineWidth(1));
  xframe->GetYaxis()->SetTitleOffset(1.8);
  yframe->GetYaxis()->SetTitleOffset(1.8);
  zframe->GetYaxis()->SetTitleOffset(1.8);
  xframe->SetMaximum(xframe->GetMaximum()*1.15);
  yframe->SetMaximum(yframe->GetMaximum()*1.15);
  zframe->SetMaximum(zframe->GetMaximum()*1.15);
  xframe->SetMinimum(0);
  yframe->SetMinimum(0);
  zframe->SetMinimum(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(xframe->findObject("plData"),"Post-selection distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plPDF" ),"Decay rate x efficiency","l");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.19); 
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.19); 
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.19); 
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("plotFit_d/fitResult_recoMC_singleComponent_"+shortString+".pdf").c_str() );

}

void fit_recoMC_singleComponentBin2(int q2Bin, int parity, int tagFlag, bool plot, bool save)
{
  if ( tagFlag==-1 )
    for (tagFlag=0; tagFlag<2; ++tagFlag)
      fit_recoMC_singleComponentBin(q2Bin, parity, tagFlag, plot, save);
  else
    fit_recoMC_singleComponentBin(q2Bin, parity, tagFlag, plot, save);
}

void fit_recoMC_singleComponentBin1(int q2Bin, int parity, int tagFlag, bool plot, bool save)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      fit_recoMC_singleComponentBin2(q2Bin, parity, tagFlag, plot, save);
  else
    fit_recoMC_singleComponentBin2(q2Bin, parity, tagFlag, plot, save);
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

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;
  if ( tagFlag < -1 || tagFlag > 1      ) return 1;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;
  if ( tagFlag==-1 ) cout<<"Running both the tag conditions"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      fit_recoMC_singleComponentBin1(q2Bin, parity, tagFlag, plot, save);
  else
    fit_recoMC_singleComponentBin1(q2Bin, parity, tagFlag, plot, save);

  return 0;

}
