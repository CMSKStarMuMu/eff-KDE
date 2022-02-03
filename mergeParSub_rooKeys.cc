#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"

using namespace RooFit ;
using namespace std ;

// Use this flag to avoid writing a KDE histogram when there are failed parallel jobs
// (the efficiency would be correct in any case and could be used in a fit, but the result is reproducible only revoving exactly the same file(s))
bool safeFilling = true;

static const int nBins = 9;

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

void mergeParSub_rooKeysBin(int q2Bin, int effIndx, int parity, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv, int year, int vers)
{
  string shortString = Form("b%ie%ip%i",q2Bin,effIndx,parity);
  cout<<"Conf: "<<shortString<<endl;

  // string containing the path and begin of the filename of output from parallel jobs
  string confString = Form("/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta/tmp_v%i/KDEhistTheta_%s_rooKeys_m_w0-%.2f_w1-%.2f_w2-%.2f_%i_%i_%i",vers,shortString.c_str(),widthCTK,widthCTL,widthPHI,xbins,ybins,zbins);

  // full histogram to fill
  vector<Double_t> xboundaries (xbins+1);
  vector<Double_t> yboundaries (ybins+1);
  vector<Double_t> zboundaries (zbins+1);
  for (int i=0; i<=xbins; ++i)
    xboundaries[i] = TMath::ACos(1.0-2.0*i/xbins);
  for (int i=0; i<=ybins; ++i)
    yboundaries[i] = TMath::ACos(1.0-2.0*i/ybins);
  for (int i=0; i<=zbins; ++i)
    zboundaries[i] = 3.141593*(2.0*i/zbins-1.0);
  
  TH3D* KDEhist = new TH3D(Form("KDEhist_%s",shortString.c_str()),Form("KDEhist_%s",shortString.c_str()),xbins,&xboundaries[0],ybins,&yboundaries[0],zbins,&zboundaries[0]);
  KDEhist->Sumw2();
  // Double_t* xboundaries = new Double_t[xbins+1];
  // Double_t* yboundaries = new Double_t[ybins+1];
  // Double_t* zboundaries = new Double_t[zbins+1];
  // TH3D* KDEhist = 0;

  TH3D* KDEhistFlat = new TH3D(Form("KDEhistFlat_%s",shortString.c_str()),Form("KDEhist_%s",shortString.c_str()),xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  KDEhistFlat->Sumw2();

  // import partial histograms
  string inFileName;
  int goodHistCnt = 0;
  for (int ndiv=0; ndiv<totdiv; ++ndiv) {
    // add final part of filename and open file
    inFileName = confString+Form("_%i-frac-%i_%i.root",ndiv,totdiv,year);
    TFile* fin = TFile::Open( inFileName.c_str() );
    if ( !fin || !fin->IsOpen() ) {
      cout<<"File not found: "<<inFileName<<endl;
      continue;
    }
    // extract and check output histogram from parallel job
    TH3D* partHist = (TH3D*)fin->Get(("KDEhist_"+shortString+"__ctK_ctL_phi").c_str());
    if ( !partHist || partHist->IsZombie() ) {
      cout<<"Histogram not found in file: "<<inFileName<<endl;
      continue;
    }
    // cout<<"Part int "<<partHist->Integral()<<" "<<partHist->Integral("width")<<endl;
    // if (!KDEhist) {
    //   // KDEhist = (TH3D*)partHist->Clone(Form("KDEhist_%s",shortString.c_str()));
    //   // KDEhist->Reset();
    //   partHist->GetXaxis()->GetLowEdge(xboundaries);
    //   partHist->GetYaxis()->GetLowEdge(yboundaries);
    //   partHist->GetZaxis()->GetLowEdge(zboundaries);
    //   xboundaries[xbins]=-1*zboundaries[0];
    //   yboundaries[ybins]=-1*zboundaries[0];
    //   zboundaries[zbins]=-1*zboundaries[0];
    //   cout<<"x: "<<xboundaries[0]<<" "<<xboundaries[25]<<" "<<xboundaries[50]<<endl;
    //   cout<<"y: "<<yboundaries[0]<<" "<<yboundaries[25]<<" "<<yboundaries[50]<<endl;
    //   cout<<"z: "<<zboundaries[0]<<" "<<zboundaries[25]<<" "<<zboundaries[50]<<endl;
    //   KDEhist = new TH3D(Form("KDEhist_%s",shortString.c_str()),Form("KDEhist_%s",shortString.c_str()),xbins,xboundaries,ybins,yboundaries,zbins,zboundaries);
    //   KDEhist->Sumw2();
    // }
    // add it to the full histogram
    // KDEhist->Add(partHist);
    for (int ix=1; ix<=xbins; ++ix)
      for (int iy=1; iy<=ybins; ++iy)
	for (int iz=1; iz<=zbins; ++iz) {
	  int iBin =KDEhist    ->GetBin(ix,iy,iz);
	  int iBinF=KDEhistFlat->GetBin(1+xbins-ix,1+ybins-iy,iz);
	  double partVal = partHist->GetBinContent(ix,iy,iz);
	  double volRat = ( ( KDEhist->GetXaxis()->GetBinWidth(ix) *
			      KDEhist->GetYaxis()->GetBinWidth(iy) *
			      KDEhist->GetZaxis()->GetBinWidth(iz) ) /
			    ( KDEhistFlat->GetXaxis()->GetBinWidth(1+xbins-ix) *
			      KDEhistFlat->GetYaxis()->GetBinWidth(1+ybins-iy) *
			      KDEhistFlat->GetZaxis()->GetBinWidth(iz) ) );
	  KDEhist    ->AddBinContent(iBin ,partVal);
	  KDEhistFlat->AddBinContent(iBinF,partVal*volRat);
	}
    // cout<<"Hist int "<<KDEhist->Integral()<<" "<<KDEhist->Integral("width")<<endl;
    // cout<<"HiFl int "<<KDEhistFlat->Integral()<<" "<<KDEhistFlat->Integral("width")<<endl;
    // count partial histogram successfully merged
    ++goodHistCnt;
    delete partHist;
    fin->Close();
    delete fin;
  }
  if ( goodHistCnt!=totdiv ) {
    // if not all the partial histograms are found, give a warning and correct the normalisation
    // (or make the macro fail, in case the flag is activated)
    cout<<"Warning! Not all partial histograms found: "<<goodHistCnt<<" of "<<totdiv<<endl;
    if (safeFilling) return;
    KDEhist->Scale(1.0*totdiv/goodHistCnt);
  }

  KDEhistFlat->Scale(KDEhist->Integral()/KDEhistFlat->Integral());

  // check histogram against empty or negative bin contents, which would make the fit fail
  // (this is not expected from KDE description and should never happen,
  // since the tails of gaussian kernels are always at values greater than zero)
  double minVal = KDEhist->GetMinimum();
  if ( minVal<=0 ) {
    cout<<"Histogram has empty bins, this is bad. Abort!"<<endl;
    return;
  }

  // save histograms in files
  // to facilitate plotting same terms with different configurations (SF and sampling bins) for comparisons, thay are saved in the same file
  // to reduce the number of files produced to ~10/20, different terms of the same efficiency are saved in the same file
  // (this also allows a single extractEff call to access a single file)
  string dirName = "/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta/";
  TFile* fout = TFile::Open( Form((parity==0?"%sfiles/KDEhist_b%i_ev_%i_v%i.root":"%sfiles/KDEhist_b%i_od_%i_v%i.root"),dirName.c_str(),q2Bin,year,vers), "UPDATE" );
  KDEhist->Write( Form("histTheta_indx%i_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",effIndx,widthCTK,widthCTL,widthPHI,xbins,ybins,zbins), TObject::kWriteDelete );
  KDEhistFlat->Write( Form("hist_indx%i_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",effIndx,widthCTK,widthCTL,widthPHI,xbins,ybins,zbins), TObject::kWriteDelete );
  fout->Close();
  
  // Remind user to delete partial files
  cout<<endl<<"Please, remove partial files running:\nrm "<<confString<<Form("_*-frac-%i_%i.root",totdiv,year)<<endl<<endl;
  
}

void mergeParSub_rooKeysBin2(int q2Bin, int effIndx, int parity, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv, int year, int vers)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      mergeParSub_rooKeysBin(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv, year, vers);
  else
    mergeParSub_rooKeysBin(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv, year, vers);
}

void mergeParSub_rooKeysBin1(int q2Bin, int effIndx, int parity, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv, int year, int vers)
{
  if ( effIndx==-1 )
    for (effIndx=0; effIndx<6; ++effIndx)
      mergeParSub_rooKeysBin2(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv, year, vers);
  else
    mergeParSub_rooKeysBin2(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv, year, vers);
}

void mergeParSub_rooKeys(int q2Bin = -1, int effIndx = -1, int parity = -1, float widthCTK = 0.3, float widthCTL = 0.3, float widthPHI = 0.3, int xbins=50, int ybins = 50, int zbins = 50, int totdiv = 50, int year = 2016, int vers = -1)
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
      mergeParSub_rooKeysBin1(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv, year, vers);
  else
    mergeParSub_rooKeysBin1(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, totdiv, year, vers);
  
}
