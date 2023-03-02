#include <TFile.h>
#include <TChain.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TRandom3.h>
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
#include "ShapeSigAng.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

TH3D* InvertHisto(TH3D* hin, string hname);

// q2-bin format: [0-8] for one bin
//                [-1] for each bin recursively
// parity format: [0] even
//                [1] odd
//                [-1] for each parity recursively

void extractEffBin(int q2Bin, int parity, float width00, float width01, float width02, float width10, float width11, float width12, float width20, float width21, float width22, float width30, float width31, float width32, float width40, float width41, float width42, int xbins, int ybins, int zbins, int year, int vers)
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

  string inFileName = Form((parity==0?"files/KDEhist_b%i_ev_%i_v%i.root":"files/KDEhist_b%i_od_%i_v%i.root"),q2Bin,year,vers%10);
  TFile* fin = TFile::Open( inFileName.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<inFileName<<endl;
    return;
  }
  for (int effIndx=0; effIndx<5; ++effIndx) {
    if (effIndx==3 && vers>9) {
      inFileName = Form((parity==0?"files/KDEhist_b%i_ev_%i_v%i.root":"files/KDEhist_b%i_od_%i_v%i.root"),q2Bin,year,vers);
      fin = TFile::Open( inFileName.c_str() );
      if ( !fin || !fin->IsOpen() ) {
	cout<<"File not found: "<<inFileName<<endl;
	return;
      }
    }
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

  double MCmFrac = KDEhist[4]->Integral() / KDEhist[3]->Integral();
  cout<<"MC-based mFrac: "<<MCmFrac<<endl;

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
  if (doCT) effCHist->SetTitle(Form("eff-ct-hist_%s%1.1f-%1.1f-%1.1f_bins-%i-%i-%i",confString.c_str(),width30,width31,width32,xbins,ybins,zbins));
  if (doWT) effWHist->SetTitle(Form("eff-wt-hist_%s%1.1f-%1.1f-%1.1f_bins-%i-%i-%i",confString.c_str(),width40,width41,width42,xbins,ybins,zbins));

  vector<string> ParName = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};
  vector<RooRealVar*> Par(ParName.size(),0);

  string finGenName = "fitResult_genMC_penalty.root";
  cout<<"Opening file "<<finGenName<<endl;
  auto finGen = TFile::Open(finGenName.c_str());
  if ( !finGen || finGen->IsZombie() ) {
    cout<<"Missing gen file: "<<finGenName<<endl;
    return;
  }
  cout<<"Loading bin "<<q2Bin<<" results from file "<<finGenName<<endl;
  auto wspRes = (RooWorkspace*)finGen->Get(Form("ws_b%ip%i_s0_pow1.0",q2Bin,parity));
  if (!wspRes || wspRes->IsZombie()) {
    cout<<"Workspace not found: "<<Form("ws_b%ip%i_s0_pow1.0",q2Bin,parity)<<endl;
    return;
  }
  auto fitResultGen = (RooFitResult*)wspRes->obj("fitResult");
  if (!fitResultGen || fitResultGen->IsZombie() || fitResultGen->status()!=0 || fitResultGen->covQual()!=3) {
    cout<<"Non valid fit result in workspace: "<<Form("ws_b%ip%i_s0_pow1.0",q2Bin,parity)<<endl;
    return;
  }
  auto pars = fitResultGen->floatParsFinal();
  for (uint ParIndx=0; ParIndx<Par.size(); ++ParIndx)
    Par[ParIndx] = (RooRealVar*)pars.find(ParName[ParIndx].c_str());

  RooRealVar* ctK = new RooRealVar("ctK","cos(#theta_{K})",-1,1);
  RooRealVar* ctL = new RooRealVar("ctL","cos(#theta_{L})",-1,1);
  RooRealVar* phi = new RooRealVar("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgList vars (*ctK, *ctL, *phi);
  auto effCData = new RooDataHist(("effCData_"+shortString+Form("_%i",year)).c_str(),"effCData",vars,effCHist);
  auto effWData = new RooDataHist(("effWData_"+shortString+Form("_%i",year)).c_str(),"effWData",vars,effWHist);
  auto effC = new RooHistFunc(("effC_"+shortString+Form("_%i",year)).c_str(),
			      Form("effC%i",year),
			      vars,
			      *effCData,
			      1);
  auto effW = new RooHistFunc(("effW_"+shortString+Form("_%i",year)).c_str(),
			      Form("effW%i",year),
			      vars,
			      *effWData,
			      1);

  RooAbsReal* ang_ct = new ShapeSigAng( ("PDF_sig_ang_ct_"+shortString+Form("_%i",year)).c_str(),
					Form("PDF_sig_ang_ct_%i",year),
					*ctK,*ctL,*phi,
					*Par[0],*Par[1],*Par[2],*Par[3],*Par[4],*Par[5],*Par[6],*Par[7],
					*effC, vector<double>(0),
					true
					);
  RooAbsReal* ang_wt = new ShapeSigAng( ("PDF_sig_ang_wt_"+shortString+Form("_%i",year)).c_str(),
					Form("PDF_sig_ang_wt_%i",year),
					*ctK,*ctL,*phi,
					*Par[0],*Par[1],*Par[2],*Par[3],*Par[4],*Par[5],*Par[6],*Par[7],
					*effW, vector<double>(0),
					false
					);

  auto intC = ang_ct->createIntegral(vars);
  auto intW = ang_wt->createIntegral(vars);
  double mtf = intW->getVal()/intC->getVal();
  cout<<"Efficiency mFrac before correction: "<<mtf<<endl;
  effWHist->Scale(MCmFrac/mtf);

  // save histograms in file
  // from this point on the names of files and objects will only contain information about bin number, parity, and tag condition
  // in this way, there in the rest of the code there is no need to specify the KDE configuration any time it is run
  TFile* fout = TFile::Open( Form((parity==0?"files/KDEeff_b%i_ev_%i_v%i.root":"files/KDEeff_b%i_od_%i_v%i.root"),q2Bin,year,vers), "UPDATE" );
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

void extractEffBin1(int q2Bin, int parity, float width00, float width01, float width02, float width10, float width11, float width12, float width20, float width21, float width22, float width30, float width31, float width32, float width40, float width41, float width42, int xbins, int ybins, int zbins, int year, int vers)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      extractEffBin(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins, year, vers);
  else
    extractEffBin(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins, year, vers);
}

int main(int argc, char** argv)
{

  int q2Bin = -1;
  int parity = -1;
  float width40 = 0;
  float width41 = 0;
  float width42 = 0;
  float width30 = 0;
  float width31 = 0;
  float width32 = 0;
  float width20 = 0;
  float width21 = 0;
  float width22 = 0;
  float width10 = 0;
  float width11 = 0;
  float width12 = 0;
  float width00 = 0;
  float width01 = 0;
  float width02 = 0;
  int xbins = 50;
  int ybins = 0;
  int zbins = 0;
  int year = 2016;
  int vers = -1;

  if ( argc > 1 ) q2Bin = atoi(argv[1]);
  if ( argc > 2 ) parity = atoi(argv[2]);
  if ( argc > 3 ) width40 = atof(argv[3]);
  if ( argc > 4 ) width41 = atof(argv[4]);
  if ( argc > 5 ) width42 = atof(argv[5]);
  if ( argc > 6 ) width30 = atof(argv[6]);
  if ( argc > 7 ) width31 = atof(argv[7]);
  if ( argc > 8 ) width32 = atof(argv[8]);
  if ( argc > 9 ) width20 = atof(argv[9]);
  if ( argc > 10 ) width21 = atof(argv[10]);
  if ( argc > 11 ) width22 = atof(argv[11]);
  if ( argc > 12 ) width10 = atof(argv[12]);
  if ( argc > 13 ) width11 = atof(argv[13]);
  if ( argc > 14 ) width12 = atof(argv[14]);
  if ( argc > 15 ) width00 = atof(argv[15]);
  if ( argc > 16 ) width01 = atof(argv[16]);
  if ( argc > 17 ) width02 = atof(argv[17]);
  if ( argc > 18 ) xbins = atoi(argv[18]);
  if ( argc > 19 ) ybins = atoi(argv[19]);
  if ( argc > 20 ) zbins = atoi(argv[20]);
  if ( argc > 21 ) year = atoi(argv[21]);
  if ( argc > 22 ) vers = atoi(argv[22]);

  if ( q2Bin<-1 || q2Bin>=nBins ) return 1;

  if ( parity<-1 || parity>1 ) return 1;

  if ( width00<=0 ) return 1;
  if ( width01<=0 ) return 1;
  if ( width02<=0 ) return 1;
  if ( width10<=0 ) return 1;
  if ( width20<=0 ) return 1;
  if ( width30<=0 ) return 1;
  if ( width40<=0 ) return 1;
  if ( width11<=0 ) return 1;
  if ( width21<=0 ) return 1;
  if ( width31<=0 ) return 1;
  if ( width41<=0 ) return 1;
  if ( width12<=0 ) return 1;
  if ( width22<=0 ) return 1;
  if ( width32<=0 ) return 1;
  if ( width42<=0 ) return 1;

  if ( xbins<1 ) return 1;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      extractEffBin1(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins, year, vers);
  else
    extractEffBin1(q2Bin, parity, width00, width01, width02, width10, width11, width12, width20, width21, width22, width30, width31, width32, width40, width41, width42, xbins, ybins, zbins, year, vers);

  return 0;

}
