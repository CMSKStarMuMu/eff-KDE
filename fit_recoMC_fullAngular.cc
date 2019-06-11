#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooNumIntConfig.h>

// #include "PdfRT.h"
// #include "PdfWT.h"
#include "PdfSigAng.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [2*nBins];

void fit_recoMC_fullAngularBin(int q2Bin, int parity, bool plot, bool save)
{

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  string filename_data = Form("effWeightedDataset_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,1-parity));
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
  // import the complementary datasets, to fit statistically uncorrelated values
  RooDataSet* dataCT = (RooDataSet*)wsp->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin));
  RooDataSet* dataWT = (RooDataSet*)wsp->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin));
  if ( !dataCT || dataCT->IsZombie() || !dataWT || dataWT->IsZombie() ) {
    cout<<"Datasets not found in file: "<<filename_data<<endl;
    return;
  }
  RooDataSet* data = new RooDataSet(*dataCT,("data_"+shortString).c_str());
  data->append(*dataWT);

  // import KDE efficiency histograms
  string filename = "/afs/cern.ch/user/a/aboletti/public/Run2-KstarMuMu/KDEeff/testVersion-integrals/";
  filename = filename + Form((parity==0?"KDEeff_b%i_ev.root":"KDEeff_b%i_od.root"),q2Bin);
  TFile* fin = new TFile( filename.c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  TH3D* effCHist = (TH3D*)fin->Get(Form("effCHist_b%ip%i",q2Bin,parity));
  TH3D* effWHist = (TH3D*)fin->Get(Form("effWHist_b%ip%i",q2Bin,parity));
  if ( !effCHist || effCHist->IsZombie() || !effWHist || effWHist->IsZombie() ) {
    cout<<"Efficiency histograms not found in file: "<<filename<<endl;
    return;
  }
  // create efficiency functions
  RooDataHist* effCData = new RooDataHist(("effCData_"+shortString).c_str(),"effCData",vars,effCHist);
  RooDataHist* effWData = new RooDataHist(("effWData_"+shortString).c_str(),"effWData",vars,effWHist);
  RooAbsReal* effC = new RooHistFunc(("effC_"+shortString).c_str(),"effC",vars,*effCData,1);
  RooAbsReal* effW = new RooHistFunc(("effW_"+shortString).c_str(),"effW",vars,*effWData,1);
  // import precomputed integrals
  vector<double> intCVec (0);
  vector<double> intWVec (0);
  TH1D* intCHist = (TH1D*)fin->Get(("MCint_"+shortString+"t1").c_str());
  TH1D* intWHist = (TH1D*)fin->Get(("MCint_"+shortString+"t0").c_str());
  if ( !intCHist || intCHist->IsZombie() || !intWHist || intWHist->IsZombie() ) {
    cout<<"Integral histograms not found in file: "<<filename<<endl<<"Using rooFit integration"<<endl;
    intCVec.push_back(0);
    intWVec.push_back(0);
  } else if ( strcmp( intCHist->GetTitle(), effCHist->GetTitle() ) || strcmp( intWHist->GetTitle(), effWHist->GetTitle() ) ) {
    cout<<"Integral histograms are incoherent with efficiency in file: "<<filename<<endl;
    cout<<"Efficiency (CT) conf: "<<effCHist->GetTitle()<<endl;
    cout<<"Integral (CT) conf: "<<intCHist->GetTitle()<<endl;
    cout<<"Efficiency (WT) conf: "<<effWHist->GetTitle()<<endl;
    cout<<"Integral (WT) conf: "<<intWHist->GetTitle()<<endl;
    cout<<"Using rooFit integration"<<endl;
    intCVec.push_back(0);
    intWVec.push_back(0);
  } else {
    for (int i=1; i<=intCHist->GetNbinsX(); ++i) intCVec.push_back(intCHist->GetBinContent(i));
    for (int i=1; i<=intWHist->GetNbinsX(); ++i) intWVec.push_back(intWHist->GetBinContent(i));
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
  RooRealVar* mFrac = new RooRealVar("mFrac","mistag fraction",1);
  mFrac->setConstant();

  RooAbsPdf* PDF_sig_ang_fullAngular = new PdfSigAng(("PDF_sig_ang_fullAngular_CT_"+shortString).c_str(),"PDF_sig_ang_fullAngular_CT",
						     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,*effC,*effW,intCVec,intWVec);

  // Customise rooFit integrator
  // PDF_sig_ang_fullAngular->getIntegratorConfig()->methodND().setLabel("RooMCIntegrator");

  RooFitResult* fitResult = PDF_sig_ang_fullAngular->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true),Hesse(true),Strategy(2),Minos(true),Offset(true),SumW2Error(true)); 
  fitResult->Print("v");

  if (save) {
    TFile* fout = new TFile("fitResultWeighted_recoMC_fullAngular.root","UPDATE");
    fitResult->Write(("fitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  int confIndex = nBins*parity + q2Bin;
  string longString  = "Fit to reconstructed events";
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
  PDF_sig_ang_fullAngular->plotOn(xframe,LineWidth(1),Name("plPDF"));
  PDF_sig_ang_fullAngular->plotOn(yframe,LineWidth(1));
  PDF_sig_ang_fullAngular->plotOn(zframe,LineWidth(1));
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

  c[confIndex]->SaveAs( ("fitResultWeighted_recoMC_fullAngular_"+shortString+".pdf").c_str() );

}

void fit_recoMC_fullAngularBin1(int q2Bin, int parity, bool plot, bool save)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      fit_recoMC_fullAngularBin(q2Bin, parity, plot, save);
  else
    fit_recoMC_fullAngularBin(q2Bin, parity, plot, save);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);

  bool plot = true;
  bool save = true;

  if ( argc >= 4 && atoi(argv[3]) == 0 ) plot = false;
  if ( argc >= 5 && atoi(argv[4]) == 0 ) save = false;

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      fit_recoMC_fullAngularBin1(q2Bin, parity, plot, save);
  else
    fit_recoMC_fullAngularBin1(q2Bin, parity, plot, save);

  return 0;

}
