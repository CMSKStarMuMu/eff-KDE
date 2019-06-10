#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>

#include "DecayRate.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [2*nBins];

void fit_genMCBin(int q2Bin, int parity, bool plot, bool save)
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
  // import the "other parity" dataset, to stay coherent with fit_recoMC notation
  string datasetString = Form((parity==1?"data_genDen_ev_b%i":"data_genDen_od_b%i"),q2Bin);
  RooDataSet* data = (RooDataSet*)wsp->data(datasetString.c_str());
  if ( !data || data->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename_data<<endl;
    return;
  }

  // define angular parameters
  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0,-1,1);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0,-1,1);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

  RooAbsPdf* PDF_sig_ang_decayRate = new DecayRate(("PDF_sig_ang_decayRate_"+shortString).c_str(),"PDF_sig_ang_decayRate",
						   *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  RooFitResult * fitResult = PDF_sig_ang_decayRate->fitTo(*data,Minimizer("Minuit2","migrad"),Save(true),Timer(true),NumCPU(2),Hesse(true),Strategy(2),Minos(false),SumW2Error(true)); 
  fitResult->Print("v");

  if (save) {
    TFile* fout = new TFile("fitResultWeighted_genMC.root","UPDATE");
    fitResult->Write(("fitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  int confIndex = nBins*parity + q2Bin;
  string longString  = "Fit to generation-level distributions";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit plojections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to GEN-level MC - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title(longString.c_str()));
  RooPlot* yframe = ctL->frame(Title(longString.c_str()));
  RooPlot* zframe = phi->frame(Title(longString.c_str()));
  data->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(100),Name("plData"));
  data->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(100));
  data->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(100));
  PDF_sig_ang_decayRate->plotOn(xframe,LineWidth(1),Name("plPDF"));
  PDF_sig_ang_decayRate->plotOn(yframe,LineWidth(1));
  PDF_sig_ang_decayRate->plotOn(zframe,LineWidth(1));
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
  leg->AddEntry(xframe->findObject("plData"),"Generation-level distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plPDF" ),"Angular decay rate" ,"l");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.15);
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.15);
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.15);
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("fitResultWeighted_genMC_"+shortString+".pdf").c_str() );

}

void fit_genMCBin1(int q2Bin, int parity, bool plot, bool save)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      fit_genMCBin(q2Bin, parity, plot, save);
  else
    fit_genMCBin(q2Bin, parity, plot, save);
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

  if ( argc >= 5 && atoi(argv[4]) == 0 ) plot = false;
  if ( argc >= 6 && atoi(argv[5]) == 0 ) save = false;

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      fit_genMCBin1(q2Bin, parity, plot, save);
  else
    fit_genMCBin1(q2Bin, parity, plot, save);

  return 0;

}
