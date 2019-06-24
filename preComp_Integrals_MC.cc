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

using namespace RooFit;
using namespace std;

static const int nBins = 9;
static const int nFunc = 11;

void preComp_Integrals_MCBin(int q2Bin, int parity, int tagFlag, int cnt_hit, int seed, int altIndx)
{

  string shortString = Form("b%ip%it%i",q2Bin,parity,tagFlag);
  cout<<"Conf: "<<shortString<<endl;

  RooRealVar* ctK = new RooRealVar("ctK","cos(#theta_{K})",-1,1);
  RooRealVar* ctL = new RooRealVar("ctL","cos(#theta_{L})",-1,1);
  RooRealVar* phi = new RooRealVar("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgList vars (*ctK, *ctL, *phi);

  // import KDE efficiency histograms
  string filename = Form((parity==0?"KDEeff_b%i_ev_alt%i.root":"KDEeff_b%i_od_alt%i.root"),q2Bin,altIndx);
  TFile* fin = new TFile( filename.c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  TH3D* effHist = (TH3D*)fin->Get(Form((tagFlag==0?"effWHist_b%ip%i":"effCHist_b%ip%i"),q2Bin,parity));
  if ( !effHist || effHist->IsZombie() ) {
    cout<<"Efficiency histogram not found in file: "<<filename<<endl;
    return;
  }
  double maxEff = effHist->GetMaximum();
  // create efficiency functions
  RooDataHist* effData = new RooDataHist(("effData_"+shortString).c_str(),"effData",vars,effHist);
  RooAbsReal* eff = new RooHistFunc(("eff_"+shortString).c_str(),"eff",vars,*effData,1);

  // define partial PDF components
  RooFormulaVar* function [nFunc];
  function[0 ] = new RooFormulaVar(("function0" +shortString).c_str(),"function0" ,"1-ctK*ctK",vars);
  function[1 ] = new RooFormulaVar(("function1" +shortString).c_str(),"function1" ,"ctK*ctK",vars);
  function[2 ] = new RooFormulaVar(("function2" +shortString).c_str(),"function2" ,"(1-ctK*ctK)*(2*ctL*ctL-1)",vars);
  function[3 ] = new RooFormulaVar(("function3" +shortString).c_str(),"function3" ,"ctK*ctK*(2*ctL*ctL-1)",vars);
  function[4 ] = new RooFormulaVar(("function4" +shortString).c_str(),"function4" ,"(1-ctK*ctK)*(1-ctL*ctL)*cos(2*phi)",vars);
  function[5 ] = new RooFormulaVar(("function5" +shortString).c_str(),"function5" ,"4*ctK*ctL*sqrt((1-ctK*ctK)*(1-ctL*ctL))*cos(phi)",vars);
  function[6 ] = new RooFormulaVar(("function6" +shortString).c_str(),"function6" ,"2*ctK*sqrt((1-ctK*ctK)*(1-ctL*ctL))*cos(phi)",vars);
  function[7 ] = new RooFormulaVar(("function7" +shortString).c_str(),"function7" ,"2*ctK*sqrt((1-ctK*ctK)*(1-ctL*ctL))*sin(phi)",vars);
  function[8 ] = new RooFormulaVar(("function8" +shortString).c_str(),"function8" ,"4*ctK*ctL*sqrt((1-ctK*ctK)*(1-ctL*ctL))*sin(phi)",vars);
  function[9 ] = new RooFormulaVar(("function9" +shortString).c_str(),"function9" ,"(1-ctK*ctK)*ctL",vars);
  function[10] = new RooFormulaVar(("function10"+shortString).c_str(),"function10","(1-ctK*ctK)*(1-ctL*ctL)*sin(2*phi)",vars);

  // define counters
  int cnt_p [nFunc];
  int cnt_m [nFunc];
  for (int i=0; i<nFunc; ++i) cnt_p[i]=cnt_m[i]=0;

  // start the generation
  TRandom3 randGen (seed);
  double hVal, effVal, intVal;
  int iPoint, iFunc;
  double piVal = TMath::Pi();
  for (iPoint=0; cnt_p[1]<cnt_hit; ++iPoint) {
    ctK->setVal(2*randGen.Rndm()-1);
    ctL->setVal(2*randGen.Rndm()-1);
    phi->setVal((2*randGen.Rndm()-1)*piVal);

    hVal=randGen.Rndm()*maxEff;
    effVal=eff->getVal();
    if (effVal<hVal) continue;

    for (iFunc=0; iFunc<nFunc; ++iFunc) {
      intVal=effVal*function[iFunc]->getVal();
      if (intVal>hVal) ++cnt_p[iFunc];
      else if (-1*intVal>hVal) ++cnt_m[iFunc];
    }
  }

  // coumpute integral values and save all ingredients in histograms
  TH1D* hvolume  = new TH1D(("hvolume" +shortString).c_str(),effHist->GetTitle(),nFunc,-0.5,nFunc-0.5);
  TH1D* hcnt_p   = new TH1D(("hcnt_p"  +shortString).c_str(),"hcnt_p"  ,nFunc,-0.5,nFunc-0.5);
  TH1D* hcnt_m   = new TH1D(("hcnt_m"  +shortString).c_str(),"hcnt_m"  ,nFunc,-0.5,nFunc-0.5);
  TH1D* hcnt_tot = new TH1D(("hcnt_tot"+shortString).c_str(),"hcnt_tot",nFunc,-0.5,nFunc-0.5);
  double volume = maxEff*8*TMath::Pi();

  for (int i=0; i<nFunc; ++i) {

    hvolume ->SetBinContent(i+1,volume);
    hcnt_p  ->SetBinContent(i+1,cnt_p[i]);
    hcnt_m  ->SetBinContent(i+1,cnt_m[i]);
    hcnt_tot->SetBinContent(i+1,iPoint);

    // cout<<i<<"\t"<<cnt_p[i]<<"\t"<<cnt_m[i]<<"\t"
    // 	<<volume*(cnt_p[i]-cnt_m[i])/iPoint<<"\t"
    // 	<<volume*sqrt(1.0*(cnt_p[i]+cnt_m[i])*(iPoint-cnt_p[i]-cnt_m[i])/iPoint)/iPoint<<endl;

  }

  // save histograms to file
  TFile* fout = TFile::Open(Form("PreIntMC_%s_s%i_alt%i.root",shortString.c_str(),seed,altIndx),"RECREATE");
  fout->cd();
  hvolume ->Write();
  hcnt_p  ->Write();
  hcnt_m  ->Write();
  hcnt_tot->Write();
  fout->Close();

}

void preComp_Integrals_MCBin2(int q2Bin, int parity, int tagFlag, int cnt_hit, int seed, int altIndx)
{
  if ( tagFlag==-1 )
    for (tagFlag=0; tagFlag<2; ++tagFlag)
      preComp_Integrals_MCBin(q2Bin, parity, tagFlag, cnt_hit, seed, altIndx);
  else
    preComp_Integrals_MCBin(q2Bin, parity, tagFlag, cnt_hit, seed, altIndx);
}

void preComp_Integrals_MCBin1(int q2Bin, int parity, int tagFlag, int cnt_hit, int seed, int altIndx)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      preComp_Integrals_MCBin2(q2Bin, parity, tagFlag, cnt_hit, seed, altIndx);
  else
    preComp_Integrals_MCBin2(q2Bin, parity, tagFlag, cnt_hit, seed, altIndx);
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
  int cnt_hit = 1e6;
  int seed    = 1;
  int altIndx = 53324;

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);
  if ( argc >= 4 ) tagFlag = atoi(argv[3]);
  if ( argc >= 5 ) cnt_hit = atoi(argv[4]);
  if ( argc >= 6 ) seed    = atoi(argv[5]);
  if ( argc >= 7 ) altIndx = atoi(argv[6]);

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;
  if ( tagFlag < -1 || tagFlag > 1      ) return 1;

  if ( seed    < 0 ) return 1;
  if ( cnt_hit < 0 ) return 1;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      preComp_Integrals_MCBin1(q2Bin, parity, tagFlag, cnt_hit, seed, altIndx);
  else
    preComp_Integrals_MCBin1(q2Bin, parity, tagFlag, cnt_hit, seed, altIndx);

  return 0;

}
