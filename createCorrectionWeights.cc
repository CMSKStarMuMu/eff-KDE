#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TStyle.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>

#include "DecayRate.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* csx [2*nBins];
TCanvas* csy [2*nBins];
TCanvas* csz [2*nBins];
TCanvas* cp2 [2*nBins];
TCanvas* cp1 [2*nBins];

void createCorrectionWeightsBin(int q2Bin, int parity, bool plot, bool save)
{

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  int hBins = 20;

  // Load variables and dataset
  string filename_data = Form("effDataset_b%i.root",q2Bin);
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

  // Load angular parameters from fit results
  string filename = "fitResult_genMC.root";
  TFile* fin = TFile::Open( filename.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  string fitResultName = "fitResult_" + shortString;
  RooFitResult* fitResult = (RooFitResult*)fin->Get(fitResultName.c_str());
  if ( !fitResult || fitResult->IsZombie() ) {
    cout<<"Fit result "<<fitResultName<<" not found in file: "<<filename<<endl;
    return;
  }
  if ( fitResult->status()!=0 )
    cout<<"WARNING: using fit results from not converged fit"<<endl;
  RooRealVar* Fl  = (RooRealVar*)fitResult->floatParsFinal().find("Fl" );
  RooRealVar* P1  = (RooRealVar*)fitResult->floatParsFinal().find("P1" );
  RooRealVar* P2  = (RooRealVar*)fitResult->floatParsFinal().find("P2" );
  RooRealVar* P3  = (RooRealVar*)fitResult->floatParsFinal().find("P3" );
  RooRealVar* P4p = (RooRealVar*)fitResult->floatParsFinal().find("P4p");
  RooRealVar* P5p = (RooRealVar*)fitResult->floatParsFinal().find("P5p");
  RooRealVar* P6p = (RooRealVar*)fitResult->floatParsFinal().find("P6p");
  RooRealVar* P8p = (RooRealVar*)fitResult->floatParsFinal().find("P8p");

  // Define decay rate
  RooAbsPdf* PDF_sig_ang_decayRate = new DecayRate(("PDF_sig_ang_decayRate_"+shortString).c_str(),"PDF_sig_ang_decayRate",
						   *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // Create histogram with weights
  TH3D* hDecRate = (TH3D*)PDF_sig_ang_decayRate->createHistogram( ("hDecRate"+shortString).c_str(), *ctK, Binning(hBins,-1,1),
							   YVar( *ctL, Binning(hBins,-1,1) ),
							   ZVar( *phi, Binning(hBins,-TMath::Pi(),TMath::Pi()) ) );
  TH3D* hGenDist = (TH3D*)data                 ->createHistogram( ("hGenDist"+shortString).c_str(), *ctK, Binning(hBins,-1,1),
							   YVar( *ctL, Binning(hBins,-1,1) ),
							   ZVar( *phi, Binning(hBins,-TMath::Pi(),TMath::Pi()) ) );
  // TH3D* hWeight = new TH3D(*hDecRate);
  // hWeight->SetName(("hWeight_"+shortString).c_str());
  // hWeight->Divide(hGenDist);
  TH3D* hWeight = new TH3D( ("hWeight_"+shortString).c_str(), ("hWeight_"+shortString).c_str(), hBins,-1,1, hBins,-1,1, hBins,-TMath::Pi(),TMath::Pi() );
  for (int iBinX=1; iBinX<=hWeight->GetNbinsX(); ++iBinX)
    for (int iBinY=1; iBinY<=hWeight->GetNbinsY(); ++iBinY)
      for (int iBinZ=1; iBinZ<=hWeight->GetNbinsZ(); ++iBinZ) {
	hWeight->SetBinContent(iBinX,iBinY,iBinZ,hDecRate->GetBinContent(iBinX,iBinY,iBinZ)/hGenDist->GetBinContent(iBinX,iBinY,iBinZ));
	hWeight->SetBinError  (iBinX,iBinY,iBinZ,hWeight ->GetBinContent(iBinX,iBinY,iBinZ)/sqrt(hGenDist->GetBinContent(iBinX,iBinY,iBinZ)));
      }
  hWeight->Scale(hBins*hBins*hBins/hWeight->Integral());

  if (save) {
    TFile* fout = new TFile("correctionWeights.root","UPDATE");
    hWeight->Write(0,TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  gStyle->SetOptStat(0);

  int confIndex = nBins*parity + q2Bin;
  string longString  = "Fit to generation-level distributions";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);
  string confString = "MC-correction-weights_"+shortString;

  // Plot 1D slices of the KDE function and original dataset
  csx[confIndex] = new TCanvas(("csx"+shortString).c_str(),(longString+" - cos(theta_k) slices").c_str(),1500,1500) ;
  csy[confIndex] = new TCanvas(("csy"+shortString).c_str(),(longString+" - cos(theta_l) slices").c_str(),1500,1500) ;
  csz[confIndex] = new TCanvas(("csz"+shortString).c_str(),(longString+" - phi slices").c_str(),1500,1500) ;
  csx[confIndex]->Divide(5,5);
  csy[confIndex]->Divide(5,5);
  csz[confIndex]->Divide(5,5);
  TH1D* histSliceX [25];
  TH1D* histSliceY [25];
  TH1D* histSliceZ [25];

  // width of the slices in the hidden variables ("border" is half of it)
  double border = 0.025;
  // variables to be filled with global maxima
  double maxValX = 0;
  double maxValY = 0;
  double maxValZ = 0;

  // loop over slice grid
  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

      // central values and borders of the slices in the hidden variables
      double centA = -0.8 + 1.6*i/4;
      double centB = -0.8 + 1.6*j/4;

      int xBinHA = hWeight->GetXaxis()->FindBin(centA+border);
      int xBinLA = hWeight->GetXaxis()->FindBin(centA-border);
      int yBinHA = hWeight->GetYaxis()->FindBin(centA+border);
      int yBinLA = hWeight->GetYaxis()->FindBin(centA-border);
      int yBinHB = hWeight->GetYaxis()->FindBin(centB+border);
      int yBinLB = hWeight->GetYaxis()->FindBin(centB-border);
      int zBinHB = hWeight->GetZaxis()->FindBin((centB+border)*TMath::Pi());
      int zBinLB = hWeight->GetZaxis()->FindBin((centB-border)*TMath::Pi());
      
	
      // producing 1D slices of KDE descriptions
      histSliceX[5*i+j] = (TH1D*) hWeight->ProjectionX( Form("histSliceX_%i_%i_%s",i,j,shortString.c_str()),
							yBinLA, yBinHA, zBinLB, zBinHB );
      histSliceY[5*i+j] = (TH1D*) hWeight->ProjectionY( Form("histSliceY_%i_%i_%s",i,j,shortString.c_str()),
							xBinLA, xBinHA, zBinLB, zBinHB );
      histSliceZ[5*i+j] = (TH1D*) hWeight->ProjectionZ( Form("histSliceZ_%i_%i_%s",i,j,shortString.c_str()),
							xBinLA, xBinHA, yBinLB, yBinHB );
      histSliceX[5*i+j]->SetTitle( Form("%s - slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Weight",longString.c_str(),centA,centB*TMath::Pi()) );
      histSliceY[5*i+j]->SetTitle( Form("%s - slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Weight",longString.c_str(),centA,centB*TMath::Pi()) );
      histSliceZ[5*i+j]->SetTitle( Form("%s - slice ctK=%1.2f ctL=%1.2f;#phi;Weight"           ,longString.c_str(),centA,centB) );

      // scale slices to show averages instead of sums over hidden bins
      histSliceX[5*i+j]->Scale( 1.0 / (1+yBinHA-yBinLA) / (1+zBinHB-zBinLB) );
      histSliceY[5*i+j]->Scale( 1.0 / (1+xBinHA-xBinLA) / (1+zBinHB-zBinLB) );
      histSliceZ[5*i+j]->Scale( 1.0 / (1+xBinHA-xBinLA) / (1+yBinHB-yBinLB) );

      // plot in canvas
      csx[confIndex]->cd(5*j+i+1);
      // gPad->SetLeftMargin(0.18);
      histSliceX[5*i+j]->GetYaxis()->SetTitleOffset(1.7);
      histSliceX[5*i+j]->Draw();

      csy[confIndex]->cd(5*j+i+1);
      // gPad->SetLeftMargin(0.18);
      histSliceY[5*i+j]->GetYaxis()->SetTitleOffset(1.7);
      histSliceY[5*i+j]->Draw();

      csz[confIndex]->cd(5*j+i+1);
      // gPad->SetLeftMargin(0.18);
      histSliceZ[5*i+j]->GetYaxis()->SetTitleOffset(1.7);
      histSliceZ[5*i+j]->Draw();

      // checking maximum value
      if ( maxValX<histSliceX[5*i+j]->GetMaximum() ) maxValX = histSliceX[5*i+j]->GetMaximum();
      if ( maxValY<histSliceY[5*i+j]->GetMaximum() ) maxValY = histSliceY[5*i+j]->GetMaximum();
      if ( maxValZ<histSliceZ[5*i+j]->GetMaximum() ) maxValZ = histSliceZ[5*i+j]->GetMaximum();
      
    }

  // set uniform y-axis ranges
  for (int i=0; i<25; ++i) {
    histSliceX[i]->SetMaximum(maxValX);
    histSliceY[i]->SetMaximum(maxValY);
    histSliceZ[i]->SetMaximum(maxValZ);
  }
  csx[confIndex]->SaveAs( (confString+Form("_CTKslices_dp%i.pdf",(int)(border*200))).c_str() );
  csy[confIndex]->SaveAs( (confString+Form("_CTLslices_dp%i.pdf",(int)(border*200))).c_str() );
  csz[confIndex]->SaveAs( (confString+Form("_PHIslices_dp%i.pdf",(int)(border*200))).c_str() );

  // Plot 1D projections
  cp1[confIndex] = new TCanvas(("cp1"+shortString).c_str(),(longString+" - 1D Projections").c_str(),2000,700) ;
  cp1[confIndex]->Divide(3,1);
  TH1D* histProj1X = (TH1D*)hWeight->ProjectionX( Form("histProj1X_%s",shortString.c_str()) );
  TH1D* histProj1Y = (TH1D*)hWeight->ProjectionY( Form("histProj1Y_%s",shortString.c_str()) );
  TH1D* histProj1Z = (TH1D*)hWeight->ProjectionZ( Form("histProj1Z_%s",shortString.c_str()) );
  histProj1X->SetTitle( (longString+";cos(#theta_{K});Average weight").c_str() );
  histProj1Y->SetTitle( (longString+";cos(#theta_{L});Average weight").c_str() );
  histProj1Z->SetTitle( (longString+";#phi;Average weight").c_str() );
  histProj1X->Scale( 1.0 / hBins / hBins );
  histProj1Y->Scale( 1.0 / hBins / hBins );
  histProj1Z->Scale( 1.0 / hBins / hBins );
  cp1[confIndex]->cd(1);
  // gPad->SetLeftMargin(0.18);
  histProj1X->GetYaxis()->SetTitleOffset(1.7);
  histProj1X->Draw();
  cp1[confIndex]->cd(2);
  // gPad->SetLeftMargin(0.18);
  histProj1Y->GetYaxis()->SetTitleOffset(1.7);
  histProj1Y->Draw();
  cp1[confIndex]->cd(3);
  // gPad->SetLeftMargin(0.18);
  histProj1Z->GetYaxis()->SetTitleOffset(1.7);
  histProj1Z->Draw();
  cp1[confIndex]->SaveAs( (confString+"_Proj1D.pdf").c_str() );
    
  // Plot 2D projections
  cp2[confIndex] = new TCanvas(("cp2"+shortString).c_str(),(longString+" - 2D Projections").c_str(),2000,700) ;
  cp2[confIndex]->Divide(3,1);
  TH2D* histProj2XY = (TH2D*)hWeight->Project3D( "xy" );
  TH2D* histProj2XZ = (TH2D*)hWeight->Project3D( "xz" );
  TH2D* histProj2YZ = (TH2D*)hWeight->Project3D( "yz" );
  histProj2XY->SetName( Form("histProj2XY_%s",shortString.c_str()) );
  histProj2XZ->SetName( Form("histProj2XZ_%s",shortString.c_str()) );
  histProj2YZ->SetName( Form("histProj2YZ_%s",shortString.c_str()) );
  histProj2XY->SetTitle( Form("%s;cos(#theta_{L});cos(#theta_{K});Events",longString.c_str()) );
  histProj2XZ->SetTitle( Form("%s;#phi;cos(#theta_{K});;Events",longString.c_str()) );
  histProj2YZ->SetTitle( Form("%s;#phi;cos(#theta_{L});Events",longString.c_str()) );
  histProj2XY->GetXaxis()->SetTitleOffset(1.4);
  histProj2XY->GetYaxis()->SetTitleOffset(2);
  histProj2XZ->GetXaxis()->SetTitleOffset(1.4);
  histProj2XZ->GetYaxis()->SetTitleOffset(2);
  histProj2YZ->GetXaxis()->SetTitleOffset(1.4);
  histProj2YZ->GetYaxis()->SetTitleOffset(2);
  cp2[confIndex]->cd(1);
  histProj2XY->Draw("SURF3");
  cp2[confIndex]->cd(2);
  histProj2XZ->Draw("SURF3");
  cp2[confIndex]->cd(3);
  histProj2YZ->Draw("SURF3");
  cp2[confIndex]->SaveAs( (confString+"_Proj2D.pdf").c_str() );

}

void createCorrectionWeightsBin1(int q2Bin, int parity, bool plot, bool save)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      createCorrectionWeightsBin(q2Bin, parity, plot, save);
  else
    createCorrectionWeightsBin(q2Bin, parity, plot, save);
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
      createCorrectionWeightsBin1(q2Bin, parity, plot, save);
  else
    createCorrectionWeightsBin1(q2Bin, parity, plot, save);

  return 0;

}
