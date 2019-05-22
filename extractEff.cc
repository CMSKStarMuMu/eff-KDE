#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

TH3D* InvertHisto(TH3D* hin, string hname);

// q2-bin format: [0-8] for one bin
//                [-1] for each bin recursively
// parity format: [0] even
//                [1] odd
//                [-1] for each parity recursively

void extractEffBin(int q2Bin, int parity, float width0, float width1, float width2, float width3, float width4, float width5, int xbins, int ybins, int zbins)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  string confString = shortString+Form("_rooNDKeysPdf_w-%1.2f-%1.2f-%1.2f-%1.2f-%1.2f-%1.2f_bins-%i-%i-%i",width0,width1,width2,width3,width4,width5,xbins,ybins,zbins);

  TH3D* KDEhist [6];
  float width [6];
  width[0] = width0;
  width[1] = width1;
  width[2] = width2;
  width[3] = width3;
  width[4] = width4;
  width[5] = width5;

  // import histogram
  for (int effIndx=0; effIndx<6; ++effIndx) {
    string inFileName = Form((parity==0?"KDEhist_b%i_ev.root":"KDEhist_b%i_od.root"),q2Bin);
    TFile* fin = TFile::Open( inFileName.c_str() );
    if ( !fin || !fin->IsOpen() ) {
      cout<<"File not found: "<<inFileName<<endl;
      return;
    }
    string histoName = Form("hist_indx%i_%1.2f_%i_%i_%i",effIndx,width[effIndx],xbins,ybins,zbins);
    KDEhist[effIndx] = (TH3D*)fin->Get(histoName.c_str());
    if ( !KDEhist[effIndx] || KDEhist[effIndx]->IsZombie() ) {
      cout<<"Histogram "<<histoName<<" not found in file: "<<inFileName<<endl;
      return;
    }
  }

  // compose efficiencies histograms
  TH3D* factHist = (TH3D*)KDEhist[1]->Clone("factHist");
  factHist->Divide(KDEhist[0]);
  factHist->Divide(KDEhist[2]);
  TH3D* effCHist = (TH3D*)factHist->Clone(("effCHist_"+shortString).c_str());
  effCHist->Multiply(KDEhist[3]);
  TH3D* effWHist = InvertHisto(factHist,("effWHist_"+shortString).c_str());
  effWHist->Multiply(KDEhist[4]);
  TH3D* effHist = (TH3D*)effCHist->Clone(("effHist_"+shortString).c_str());
  effHist->Add(effWHist);
  TH3D* corrFracHist = (TH3D*)KDEhist[3]->Clone(("corrFracHist_"+shortString).c_str());
  corrFracHist->Divide(KDEhist[5]);
  TH3D* mistFracHist = (TH3D*)KDEhist[4]->Clone(("mistFracHist_"+shortString).c_str());
  mistFracHist->Divide(KDEhist[5]);

  // Check histograms to be in range
  // Usually the correct-tag and mistag fractions need some adjustment to be forced <1
  for (int iBinX=1; iBinX<=effHist->GetNbinsX(); ++iBinX)
    for (int iBinY=1; iBinY<=effHist->GetNbinsY(); ++iBinY)
      for (int iBinZ=1; iBinZ<=effHist->GetNbinsZ(); ++iBinZ) {
	double effVal  = effHist ->GetBinContent(iBinX,iBinY,iBinZ);
	double effCVal = effCHist->GetBinContent(iBinX,iBinY,iBinZ);
	double effWVal = effWHist->GetBinContent(iBinX,iBinY,iBinZ);
	double corrFracVal = corrFracHist->GetBinContent(iBinX,iBinY,iBinZ);
	double mistFracVal = mistFracHist->GetBinContent(iBinX,iBinY,iBinZ);
	if ( effVal <=0 ) { cout<<"ERROR: NEGATIVE EFF!"<<endl; return;}
	if ( effCVal<=0 ) { cout<<"ERROR: NEGATIVE EFF-CT!"<<endl; return;}
	if ( effWVal<=0 ) { cout<<"ERROR: NEGATIVE EFF-WT!"<<endl; return;}
	if ( corrFracVal<=0 ) { cout<<"ERROR: NEGATIVE CORR-FRAC!"<<endl; return;}
	if ( mistFracVal<=0 ) { cout<<"ERROR: NEGATIVE MIST-FRAC!"<<endl; return;}
	if ( effVal>1 ) {
	  cout<<"Warning: correcting eff value > 1"<<endl;
	  effHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  effHist->SetBinError(iBinX,iBinY,iBinZ,effHist->GetBinError(iBinX,iBinY,iBinZ)/effVal);
	}
	if ( effCVal>1 ) {
	  cout<<"Warning: correcting eff-ct value > 1"<<endl;
	  effCHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  effCHist->SetBinError(iBinX,iBinY,iBinZ,effCHist->GetBinError(iBinX,iBinY,iBinZ)/effCVal);
	}
	if ( effWVal>1 ) {
	  cout<<"Warning: correcting eff-wt value > 1"<<endl;
	  effWHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  effWHist->SetBinError(iBinX,iBinY,iBinZ,effWHist->GetBinError(iBinX,iBinY,iBinZ)/effWVal);
	}
	if ( corrFracVal>1 ) {
	  corrFracHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  corrFracHist->SetBinError(iBinX,iBinY,iBinZ,corrFracHist->GetBinError(iBinX,iBinY,iBinZ)/corrFracVal);
	}
	if ( mistFracVal>1 ) {
	  mistFracHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  mistFracHist->SetBinError(iBinX,iBinY,iBinZ,mistFracHist->GetBinError(iBinX,iBinY,iBinZ)/mistFracVal);
	}
      }

  // Save width and binning info in histograms' titles
  effHist->SetTitle(("eff-hist_"+confString).c_str());
  effCHist->SetTitle(("eff-ct-hist_"+confString).c_str());
  effWHist->SetTitle(("eff-wt-hist_"+confString).c_str());
  corrFracHist->SetTitle(("corrFrac-hist_"+confString).c_str());
  mistFracHist->SetTitle(("mistFrac-hist_"+confString).c_str());

  // save histograms in file
  TFile* fout = TFile::Open( Form((parity==0?"KDEeff_b%i_ev.root":"KDEeff_b%i_od.root"),q2Bin), "UPDATE" );
  effHist->Write(0,TObject::kWriteDelete);
  effCHist->Write(0,TObject::kWriteDelete);
  effWHist->Write(0,TObject::kWriteDelete);
  corrFracHist->Write(0,TObject::kWriteDelete);
  mistFracHist->Write(0,TObject::kWriteDelete);
  fout->Close();
  
}

TH3D* InvertHisto(TH3D* hin, string hname)
{
  // clone and reset the TH3
  TH3D* hout=(TH3D*)hin->Clone(hname.c_str());
  hout->Reset();

  // fill bin content and error one by one
  for (int ix=0; ix!=hin->GetXaxis()->GetNbins()+1;++ix) {
    int ix_inv = hin->GetXaxis()->GetNbins()+1-ix;
    for (int iy=0; iy!=hin->GetYaxis()->GetNbins()+1;++iy) {
      int iy_inv = hin->GetYaxis()->GetNbins()+1-iy;
      for (int iz=0; iz!=hin->GetZaxis()->GetNbins()+1;++iz) {
	int iz_inv = hin->GetZaxis()->GetNbins()+1-iz;
        hout->SetBinContent( ix_inv, iy_inv, iz_inv, hin->GetBinContent(ix,iy,iz) );
        hout->SetBinError  ( ix_inv, iy_inv, iz_inv, hin->GetBinError  (ix,iy,iz) );
      }
    }
  }

  return hout;

}

void extractEffBin1(int q2Bin, int parity, float width0, float width1, float width2, float width3, float width4, float width5, int xbins, int ybins, int zbins)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      extractEffBin(q2Bin, parity, width0, width1, width2, width3, width4, width5, xbins, ybins, zbins);
  else
    extractEffBin(q2Bin, parity, width0, width1, width2, width3, width4, width5, xbins, ybins, zbins);
}

void extractEff(int q2Bin = -1, int parity = -1, float width0 = 0.3, float width1 = 0.3, float width2 = 0.3, float width3 = 0.3, float width4 = 0.3, float width5 = 0.3, int xbins=50, int ybins = 0, int zbins = 0)
{

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( parity<-1 || parity>1 ) return;

  if ( width0<=0 ) return;
  if ( width1<=0 ) return;
  if ( width2<=0 ) return;
  if ( width3<=0 ) return;
  if ( width4<=0 ) return;
  if ( width5<=0 ) return;

  if ( xbins<1 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      extractEffBin1(q2Bin, parity, width0, width1, width2, width3, width4, width5, xbins, ybins, zbins);
  else
    extractEffBin1(q2Bin, parity, width0, width1, width2, width3, width4, width5, xbins, ybins, zbins);
  
}
