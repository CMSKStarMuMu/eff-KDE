#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"

using namespace RooFit ;
using namespace std ;

bool safeFilling = false;

static const int nBins = 9;

TCanvas* csx [12*nBins];
TCanvas* csy [12*nBins];
TCanvas* csz [12*nBins];
TCanvas* cp1 [12*nBins];
TCanvas* cp2 [12*nBins];

// q2-bin format: [0-8] for one bin
//                [-1] for each bin recursively
// effIndx format: [0] GEN no-filter
//                 [1] GEN filtered
//                 [2] GEN filtered from full MC sample
//                 [3] correct-tag RECO candidates
//                 [4] wrong-tag RECO candidates
//                 [5] RECO candidates
//                 [-1] for each effIndx recursively
// parity format: [0] even
//                [1] odd
//                [-1] for each parity recursively

void mergeParSub_rooKeysBin(int q2Bin, int effIndx, int parity, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv)
{
  string shortString = Form("b%ie%ip%i",q2Bin,effIndx,parity);
  cout<<"Conf: "<<shortString<<endl;

  string datasetString = "data_";
  switch (effIndx) {
  case 0:  datasetString = datasetString+"genDen"; break;
  case 1:  datasetString = datasetString+"genNum"; break;
  case 2:  datasetString = datasetString+"den"   ; break;
  case 3:  datasetString = datasetString+"ctRECO"; break;
  default: datasetString = datasetString+"wtRECO";
  }
  datasetString = datasetString + Form((parity==0?"_ev_b%i":"_od_b%i"),q2Bin);

  string longString = Form(parity==0?" distribution - q2 bin %i (even)":" distribution - q2 bin %i (odd)",q2Bin);
  switch (effIndx) {
  case 0: longString = "GEN denominator"     +longString; break;
  case 1: longString = "GEN numerator"       +longString; break;
  case 2: longString = "Full denominator"    +longString; break;
  case 3: longString = "Selected correct-tag"+longString; break;
  case 4: longString = "Selected wrong-tag"  +longString; break;
  case 5: longString = "Selected"            +longString;
  }

  string confString = Form("tmp/KDEhist_%s_rooKeys_m_w0-%.2f_w1-%.2f_w2-%.2f_%i_%i_%i",shortString.c_str(),widthCTK,widthCTL,widthPHI,xbins,ybins,zbins);

  TH3D* KDEhist = new TH3D(Form("KDEhist_%s",shortString.c_str()),Form("KDEhist_%s",shortString.c_str()),xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  KDEhist->Sumw2();

  // import partial histograms
  string inFileName;
  int goodHistCnt = 0;
  for (int ndiv=0; ndiv<totdiv; ++ndiv) {
    inFileName = confString+Form("_%i-frac-%i.root",ndiv,totdiv);
    TFile* fin = TFile::Open( inFileName.c_str() );
    if ( !fin || !fin->IsOpen() ) {
      cout<<"File not found: "<<inFileName<<endl;
      continue;
    }
    TH3D* partHist = (TH3D*)fin->Get(("KDEhist_"+shortString+"__ctK_ctL_phi").c_str());
    if ( !partHist || partHist->IsZombie() ) {
      cout<<"Histogram not found in file: "<<inFileName<<endl;
      continue;
    }
    KDEhist->Add(partHist);
    ++goodHistCnt;
    delete partHist;
    fin->Close();
    delete fin;
  }
  if ( goodHistCnt!=totdiv ) {
    cout<<"Warning! Not all partial histograms found: "<<goodHistCnt<<" of "<<totdiv<<endl;
    if (safeFilling) return;
  }

  // check histogram
  double minVal = KDEhist->GetMinimum();
  if ( minVal<=0 ) {
    cout<<"Histogram has empty bins, this is bad. Abort!"<<endl;
    return;
  }

  // save histogram in file
  TFile* fout = TFile::Open( Form((parity==0?"files/KDEhist_b%i_ev.root":"files/KDEhist_b%i_od.root"),q2Bin), "UPDATE" );
  KDEhist->Write( Form("hist_indx%i_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",effIndx,widthCTK,widthCTL,widthPHI,xbins,ybins,zbins), TObject::kWriteDelete );
  fout->Close();
  
  // Remind user to delete partial files
  cout<<endl<<"Please, remove partial files running:\nrm "<<confString<<Form("_*-frac-%i.root",totdiv)<<endl<<endl;
  
}

void mergeParSub_rooKeysBin2(int q2Bin, int effIndx, int parity, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      mergeParSub_rooKeysBin(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv);
  else
    mergeParSub_rooKeysBin(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv);
}

void mergeParSub_rooKeysBin1(int q2Bin, int effIndx, int parity, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv)
{
  if ( effIndx==-1 )
    for (effIndx=0; effIndx<6; ++effIndx)
      mergeParSub_rooKeysBin2(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv);
  else
    mergeParSub_rooKeysBin2(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv);
}

void mergeParSub_rooKeys(int q2Bin = -1, int effIndx = -1, int parity = -1, float widthCTK = 0.3, float widthCTL = 0.3, float widthPHI = 0.3, int xbins=50, int ybins = 50, int zbins = 50, int totdiv = 50)
{

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( effIndx<-1 || effIndx>5 ) return;

  if ( parity<-1 || parity>1 ) return;

  if ( widthCTK<=0 ) return;
  if ( widthCTL<=0 ) return;
  if ( widthPHI<=0 ) return;

  if ( xbins<1 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( totdiv<1 ) return;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( effIndx==-1 ) cout<<"Running all the efficiency terms"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      mergeParSub_rooKeysBin1(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv);
  else
    mergeParSub_rooKeysBin1(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv);
  
}
