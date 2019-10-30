#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;
static const int nFunc = 11;

// q2-bin format: [0-8] for one bin
//                [-1] for each bin recursively
// parity format: [0] even efficiency
//                [1] odd efficiency
//                [-1] for each parity recursively
// tag format: [0] mistagged
//             [1] correctly-tagged
//             [-1] for each tag recursively

void mergeParSub_preComp_Integrals_MCBin(int q2Bin, int parity, int tagFlag, int minSeed, int maxSeed)
{

  string shortString = Form("b%ip%it%i",q2Bin,parity,tagFlag);
  cout<<"Conf: "<<shortString<<endl;

  // Open efficiency file to append integral histogram
  string efffilename = Form((parity==0?"files/KDEeff_b%i_ev.root":"files/KDEeff_b%i_od.root"),q2Bin);
  TFile* fout = new TFile( efffilename.c_str(), "UPDATE" );
  if ( !fout || !fout->IsOpen() ) {
    cout<<"File not found: "<<efffilename<<endl;
    return;
  }

  // create empty histograms to fill with output of parallel jobs
  double volume = 0;
  TH1D* hcnt_p   = new TH1D("hcnt_p"  ,"hcnt_p"  ,nFunc,-0.5,nFunc-0.5);
  TH1D* hcnt_m   = new TH1D("hcnt_m"  ,"hcnt_m"  ,nFunc,-0.5,nFunc-0.5);
  TH1D* hcnt_tot = new TH1D("hcnt_tot","hcnt_tot",nFunc,-0.5,nFunc-0.5);

  string filename;
  string effTitleTag;

  for (int seed=minSeed; seed<=maxSeed; ++seed) {

    // open output file from parallel jobs
    filename = "tmp/PreIntMC_"+shortString+Form("_s%i.root",seed);
    TFile* fin = TFile::Open( filename.c_str() );
    if ( !fin || !fin->IsOpen() ) {
      continue;
    }

    // get from the first file the volume value and the eff_config tag
    // and check that other files have the same values
    TH1D* hvol = (TH1D*)fin->Get(("hvolume"+shortString).c_str());
    if (volume==0) {
      volume = hvol->GetBinContent(1);
      effTitleTag = hvol->GetTitle();
    } else {
      if (volume != hvol->GetBinContent(1)) {
	cout<<"Incoherent volume values in file "<<filename<<endl;
	continue;
      }
      if ( strcmp( effTitleTag.c_str(), hvol->GetTitle() ) != 0 ) {
	cout<<"Incoherent efficiency settings in file "<<filename<<endl;
	continue;
      }
    }

    // sum histogram with counter values
    hcnt_p  ->Add((TH1D*)fin->Get(("hcnt_p"  +shortString).c_str()));
    hcnt_m  ->Add((TH1D*)fin->Get(("hcnt_m"  +shortString).c_str()));
    hcnt_tot->Add((TH1D*)fin->Get(("hcnt_tot"+shortString).c_str()));

  }

  // Get and print total number of points generated
  double cnt_tot = hcnt_tot->GetBinContent(1);
  cout<<"Tot = "<<cnt_tot<<endl;

  // Compute integral values and binomial errors and save them in a histogram
  TH1D* integrals_MC = new TH1D(("MCint_"+shortString).c_str(),effTitleTag.c_str(),nFunc,-0.5,nFunc-0.5);
  for (int i=0; i<nFunc; ++i) {
    integrals_MC->SetBinContent(i+1,volume*(hcnt_p->GetBinContent(i+1)-hcnt_m->GetBinContent(i+1))/cnt_tot);
    integrals_MC->SetBinError  (i+1,volume*sqrt(1.0*(hcnt_p->GetBinContent(i+1)+hcnt_m->GetBinContent(i+1))*(cnt_tot-hcnt_p->GetBinContent(i+1)-hcnt_m->GetBinContent(i+1))/cnt_tot)/cnt_tot);

    // cout<<i<<"\t"<<hcnt_p->GetBinContent(i+1)<<"\t"<<hcnt_m->GetBinContent(i+1)<<"\t"
    // 	<<volume*(hcnt_p->GetBinContent(i+1)-hcnt_m->GetBinContent(i+1))/cnt_tot<<"\t"
    // 	<<volume*sqrt(1.0*(hcnt_p->GetBinContent(i+1)+hcnt_m->GetBinContent(i+1))*(cnt_tot-hcnt_p->GetBinContent(i+1)-hcnt_m->GetBinContent(i+1))/cnt_tot)/cnt_tot<<endl;

  }

  // append the histogram to the efficiency file
  fout->cd();
  integrals_MC->Write(0,TObject::kWriteDelete);
  fout->Close();

}

void mergeParSub_preComp_Integrals_MCBin2(int q2Bin, int parity, int tagFlag, int minSeed, int maxSeed)
{
  if ( tagFlag==-1 )
    for (tagFlag=0; tagFlag<2; ++tagFlag)
      mergeParSub_preComp_Integrals_MCBin(q2Bin, parity, tagFlag, minSeed, maxSeed);
  else
    mergeParSub_preComp_Integrals_MCBin(q2Bin, parity, tagFlag, minSeed, maxSeed);
}

void mergeParSub_preComp_Integrals_MCBin1(int q2Bin, int parity, int tagFlag, int minSeed, int maxSeed)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      mergeParSub_preComp_Integrals_MCBin2(q2Bin, parity, tagFlag, minSeed, maxSeed);
  else
    mergeParSub_preComp_Integrals_MCBin2(q2Bin, parity, tagFlag, minSeed, maxSeed);
}

void mergeParSub_preComp_Integrals_MC(int q2Bin, int parity, int tagFlag, int maxSeed, int minSeed=1)
{

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return;
  if ( parity  < -1 || parity  > 1      ) return;
  if ( tagFlag < -1 || tagFlag > 1      ) return;

  if ( minSeed < 1 ) return;
  if ( maxSeed < minSeed ) return;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;
  if ( tagFlag==-1 ) cout<<"Running both the tag conditions"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      mergeParSub_preComp_Integrals_MCBin1(q2Bin, parity, tagFlag, minSeed, maxSeed);
  else
    mergeParSub_preComp_Integrals_MCBin1(q2Bin, parity, tagFlag, minSeed, maxSeed);
  
}

