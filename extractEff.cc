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

void extractEffBin(int q2Bin, int parity, float width00, float width01, float width02, float width10, float width11, float width12, float width20, float width21, float width22, float width30, float width31, float width32, float width40, float width41, float width42, int xbins, int ybins, int zbins)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  bool doCT = true;
  bool doWT = true;

  // import the five efficiency terms
  TH3D* KDEhist [5];
  string histoName [5];
  histoName[0] = Form("hist_indx0_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",width00,width01,width02,xbins,ybins,zbins);
  histoName[1] = Form("hist_indx1_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",width10,width11,width12,xbins,ybins,zbins);
  histoName[2] = Form("hist_indx2_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",width20,width21,width22,xbins,ybins,zbins);
  histoName[3] = Form("hist_indx3_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",width30,width31,width32,xbins,ybins,zbins);
  histoName[4] = Form("hist_indx4_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",width40,width41,width42,xbins,ybins,zbins);

  string inFileName = Form((parity==0?"files/KDEhist_b%i_ev.root":"files/KDEhist_b%i_od.root"),q2Bin);
  TFile* fin = TFile::Open( inFileName.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<inFileName<<endl;
    return;
  }
  for (int effIndx=0; effIndx<5; ++effIndx) {
    KDEhist[effIndx] = (TH3D*)fin->Get(histoName[effIndx].c_str());
    if ( !KDEhist[effIndx] || KDEhist[effIndx]->IsZombie() ) {
      cout<<"Histogram "<<histoName[effIndx]<<" not found in file: "<<inFileName<<endl;
      // in case either correct-tag or wrong-tag numerator is not found, the efficiency is produced for the other
      if (effIndx==3) doCT = false;
      else if (effIndx==4) doWT = false;
      else return;
    }
  }
  if ( !doCT && !doWT ) return;

  // compose efficiencies histograms
  TH3D* factHist = (TH3D*)KDEhist[1]->Clone("factHist");
  factHist->Divide(KDEhist[0]);
  factHist->Divide(KDEhist[2]);
  TH3D* effCHist = (TH3D*)factHist->Clone(("effCHist_"+shortString).c_str());
  if (doCT) effCHist->Multiply(KDEhist[3]);
  // in wrong-tagged efficiency, all the terms but the RECO-numerator are flipped (each variable changes sign)
  TH3D* effWHist = InvertHisto(factHist,("effWHist_"+shortString).c_str());
  if (doWT) effWHist->Multiply(KDEhist[4]);

  // Check efficiency values (in each bin) to be in range [0;1]
  // in case of eff<=0, give an error (it should never happen, since the check is already run in the macro to merge histograms)
  // in case of eff>1, give a warning and set it to 1
  for (int iBinX=1; iBinX<=effCHist->GetNbinsX(); ++iBinX)
    for (int iBinY=1; iBinY<=effCHist->GetNbinsY(); ++iBinY)
      for (int iBinZ=1; iBinZ<=effCHist->GetNbinsZ(); ++iBinZ) {
	double effCVal = effCHist->GetBinContent(iBinX,iBinY,iBinZ);
	double effWVal = effWHist->GetBinContent(iBinX,iBinY,iBinZ);
	if ( doCT && effCVal<=0 ) { cout<<"ERROR: NEGATIVE EFF-CT!"<<endl; return;}
	if ( doWT && effWVal<=0 ) { cout<<"ERROR: NEGATIVE EFF-WT!"<<endl; return;}
	if ( doCT && effCVal>1 ) {
	  cout<<"Warning: correcting eff-ct value > 1"<<endl;
	  effCHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  effCHist->SetBinError(iBinX,iBinY,iBinZ,effCHist->GetBinError(iBinX,iBinY,iBinZ)/effCVal);
	}
	if ( doWT && effWVal>1 ) {
	  cout<<"Warning: correcting eff-wt value > 1"<<endl;
	  effWHist->SetBinContent(iBinX,iBinY,iBinZ,1);
	  effWHist->SetBinError(iBinX,iBinY,iBinZ,effWHist->GetBinError(iBinX,iBinY,iBinZ)/effWVal);
	}
      }

  // String tag to be saved as title of the efficiency histogram
  // this will transmit information about the KDE configuration used
  // (and make sure an efficiency is only used with its own precomputed integrals)
  string confString = Form("_KDEw_%1.1f-%1.1f-%1.1f_",width00,width01,width02);
  confString = confString+Form("%1.1f-%1.1f-%1.1f_",width10,width11,width12);
  confString = confString+Form("%1.1f-%1.1f-%1.1f_",width20,width21,width22);
  if (doCT) effCHist->SetTitle(Form("eff-ct-hist_%s%1.1f-%1.1f-%1.1f_bins-%i-%i-%i",confString,width30,width31,width32,xbins,ybins,zbins));
  if (doWT) effWHist->SetTitle(Form("eff-wt-hist_%s%1.1f-%1.1f-%1.1f_bins-%i-%i-%i",confString,width40,width41,width42,xbins,ybins,zbins));

  // save histograms in file
  // from this point on the names of files and objects will only contain information about bin number, parity, and tag condition
  // in this way, there in the rest of the code there is no need to specify the KDE configuration any time it is run
  TFile* fout = TFile::Open( Form((parity==0?"files/KDEeff_b%i_ev.root":"files/KDEeff_b%i_od.root"),q2Bin), "UPDATE" );
  if (doCT) effCHist->Write(0,TObject::kWriteDelete);
  if (doWT) effWHist->Write(0,TObject::kWriteDelete);
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

void extractEffBin1(int q2Bin, int parity, float width00, float width01, float width02, float width10, float width11, float width12, float width20, float width21, float width22, float width30, float width31, float width32, float width40, float width41, float width42, int xbins, int ybins, int zbins)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      extractEffBin(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins);
  else
    extractEffBin(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins);
}

void extractEff(int q2Bin, int parity,
		float width40, float width41, float width42,
		float width30, float width31, float width32,
		float width20, float width21, float width22,
		float width10, float width11, float width12,
		float width00, float width01, float width02,
		int xbins=50, int ybins = 0, int zbins = 0)
{

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( parity<-1 || parity>1 ) return;

  if ( width00<=0 ) return;
  if ( width01<=0 ) return;
  if ( width02<=0 ) return;
  if ( width10<=0 ) return;
  if ( width20<=0 ) return;
  if ( width30<=0 ) return;
  if ( width40<=0 ) return;
  if ( width11<=0 ) return;
  if ( width21<=0 ) return;
  if ( width31<=0 ) return;
  if ( width41<=0 ) return;
  if ( width12<=0 ) return;
  if ( width22<=0 ) return;
  if ( width32<=0 ) return;
  if ( width42<=0 ) return;

  if ( xbins<1 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      extractEffBin1(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins);
  else
    extractEffBin1(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins);
  
}
