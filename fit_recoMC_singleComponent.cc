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

float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [4*nBins];

void fit_recoMC_singleComponentBin(int q2Bin, int parity, int tagFlag, bool plot, bool save)
{

  string shortString = Form("b%ip%it%i",q2Bin,parity,tagFlag);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  string filename_data = Form("effWeightedDataset_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%i",q2Bin));
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
  // import the complementary dataset, to fit statistically uncorrelated values
  datasetString = datasetString + Form((parity==1?"_ev_b%i":"_od_b%i"),q2Bin);
  RooDataSet* data = (RooDataSet*)wsp->data(datasetString.c_str());
  if ( !data || data->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename_data<<endl;
    return;
  }

  // import KDE efficiency histograms
  string filename = "/afs/cern.ch/user/a/aboletti/public/Run2-KstarMuMu/KDEeff/testVersion-integrals/";
  filename = filename + Form((parity==0?"KDEeff_b%i_ev.root":"KDEeff_b%i_od.root"),q2Bin);
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
  // import precomputed integrals
  vector<double> intVec (0);
  string intHistString = "MCint_"+shortString;
  TH1D* intHist = (TH1D*)fin->Get(intHistString.c_str());
  if ( !intHist || intHist->IsZombie() ) {
    cout<<"Integral histogram "<<intHistString<<" not found in file: "<<filename<<endl<<"Using rooFit integration"<<endl;
    intVec.push_back(0);
  } else if ( strcmp( intHist->GetTitle(), effHist->GetTitle() ) ) {
    cout<<"Integral histogram "<<intHistString<<" is incoherent with efficiency "<<effString<<" in file: "<<filename<<endl;
    cout<<"Efficiency conf: "<<effHist->GetTitle()<<endl;
    cout<<"Integral conf: "<<intHist->GetTitle()<<endl;
    cout<<"Using rooFit integration"<<endl;
    intVec.push_back(0);
  } else {
    for (int i=1; i<=intHist->GetNbinsX(); ++i) intVec.push_back(intHist->GetBinContent(i));
  }

  // define angular parameters
  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

  // Customise initial values of parameters
  // b3wt
  // Fl->setVal( 0.60);
  // P1->setVal(-0.16);
  // P2->setVal(-0.39);
  // P3->setVal(-0.01);
  // P4p->setVal(-0.88);
  // P5p->setVal( 0.75);
  // P6p->setVal(-0.03);
  // P8p->setVal(-0.06);
  // b0wt
  // Fl->setVal( 0.71);
  // P1->setVal(-0.02);
  // P2->setVal( 0.39);
  // P3->setVal(-0.01);
  // P4p->setVal( 0.10);
  // P5p->setVal(-0.35);
  // P6p->setVal( 0.04);
  // P8p->setVal( 0.01);

  RooAbsPdf* PDF_sig_ang_singleComponent = 0;
  if (tagFlag==1) PDF_sig_ang_singleComponent = new PdfRT(("PDF_sig_ang_singleComponent_"+shortString).c_str(),"PDF_sig_ang_singleComponent",
							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*eff,intVec);
  else            PDF_sig_ang_singleComponent = new PdfWT(("PDF_sig_ang_singleComponent_"+shortString).c_str(),"PDF_sig_ang_singleComponent",
							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*eff,intVec);

  // Customise rooFit integrator
  // PDF_sig_ang_singleComponent->getIntegratorConfig()->methodND().setLabel("RooMCIntegrator");

  RooFitResult* fitResult = PDF_sig_ang_singleComponent->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true),NumCPU(1),Hesse(true),Strategy(2),Minos(true)); 
  fitResult->Print("v");

  if (save) {
    TFile* fout = new TFile("fitResultWeighted_recoMC_singleComponent.root","UPDATE");
    fitResult->Write(("fitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  int confIndex = 2*nBins*parity + nBins*tagFlag + q2Bin;
  string longString  = (tagFlag==1?"Fit to correctly tagged events":"Fit to wrongly tagged events");
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit plojections
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

  c[confIndex]->SaveAs( ("fitResultWeighted_recoMC_singleComponent_"+shortString+".pdf").c_str() );

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
