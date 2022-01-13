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

void mergeToyEff_parSub(int q2Bin, int effIndx, int parity, int seed, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int totdiv, int year)
{

  if ( q2Bin<0 || q2Bin>=nBins ) return;

  if ( effIndx<0 || effIndx>5 ) return;

  if ( parity<0 || parity>1 ) return;

  if ( seed<0 ) return;

  if ( widthCTK<=0 ) return;
  if ( widthCTL<=0 ) return;
  if ( widthPHI<=0 ) return;

  if ( xbins<1 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( totdiv<1 ) return;

  string shortString = Form("b%ie%ip%i",q2Bin,effIndx,parity);
  cout<<"Conf: "<<shortString<<endl;

  // string containing the path and begin of the filename of output from parallel jobs
  string folder = "/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/";
  string confString = folder + Form("tmp_b%i_toy%i/KDEhist_%s_rooKeys_m_w0-%.2f_w1-%.2f_w2-%.2f_%i_%i_%i",q2Bin,seed,shortString.c_str(),widthCTK,widthCTL,widthPHI,xbins,ybins,zbins);

  // full histogram to fill
  TH3D* KDEhist = new TH3D(Form("KDEhist_%s",shortString.c_str()),Form("KDEhist_%s",shortString.c_str()),xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  KDEhist->Sumw2();

  // import partial histograms
  string inFileName;
  int goodHistCnt = 0;
  for (int ndiv=0; ndiv<totdiv; ++ndiv) {
    // add final part of filename and open file
    inFileName = confString+Form("_%i-frac-%i_%i_toy%i.root",ndiv,totdiv,year,seed);
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
    // add it to the full histogram
    KDEhist->Add(partHist);
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
    if (safeFilling) {
      fstream fout("resub-mergeToyEff.list");
      fout<<q2Bin<<","<<effIndx<<","<<parity<<","<<seed<<","<<widthCTK<<","<<widthCTL<<","<<widthPHI<<","<<xbins<<","<<ybins<<","<<zbins<<","<<totdiv<<","<<year<<endl;
      fout.close();
      return;
    }
    KDEhist->Scale(1.0*totdiv/goodHistCnt);
  }

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
  string foutname = folder + Form((parity==0?"files/KDEhist_b%i_ev_%i_toy%i.root":"files/KDEhist_b%i_od_%i_toy%i.root"),q2Bin,year,seed);
  TFile* fout = TFile::Open( foutname.c_str(), "UPDATE" );
  KDEhist->Write( Form("hist_indx%i_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",effIndx,widthCTK,widthCTL,widthPHI,xbins,ybins,zbins), TObject::kWriteDelete );
  fout->Close();
  
  // Remind user to delete partial files
  cout<<endl<<"Please, remove partial files running:\nrm "<<confString<<Form("_*-frac-%i_%i_toy%i.root",totdiv,year,seed)<<endl<<endl;
  
}
