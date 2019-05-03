#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"
#include "RooNDKeysPdf.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

void composeEff_rooKeys_parSubBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins, float width, int ndiv, int totdiv);

void composeEff_rooKeys_parSub(int q2Bin, int tagFlag, int xbins=25, int ybins = 0, int zbins = 0, float width = 1, int ndiv = 0, int totdiv = 1)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( width <= 0 ) return;

  if ( ndiv<0 || ndiv>=totdiv ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Computing efficiency for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      composeEff_rooKeys_parSubBin(q2Bin, (tagFlag==1), xbins, ybins, zbins, width, ndiv, totdiv);
    }
    if (tagFlag == 2) {
      cout<<"Computing efficiency for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      composeEff_rooKeys_parSubBin(q2Bin, true,  xbins, ybins, zbins, width, ndiv, totdiv);
      composeEff_rooKeys_parSubBin(q2Bin, false, xbins, ybins, zbins, width, ndiv, totdiv);
    }
  }
  if (q2Bin == -1) {
    cout<<"Computing efficiency for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	composeEff_rooKeys_parSubBin(q2Bin, (tagFlag==1), xbins, ybins, zbins, width, ndiv, totdiv);
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	composeEff_rooKeys_parSubBin(q2Bin, true,  xbins, ybins, zbins, width, ndiv, totdiv);
	composeEff_rooKeys_parSubBin(q2Bin, false, xbins, ybins, zbins, width, ndiv, totdiv);
      }
    }
  }
  
}

void composeEff_rooKeys_parSubBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins, float width, int ndiv, int totdiv)
{

  string shortString = Form(tagFlag?"b%ict":"b%iwt",q2Bin);
  string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);

  // Load datasets
  TFile* fin = TFile::Open( ("effDataset_"+shortString+".root").c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: effDataset_"+shortString+".root"<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin->Get(("ws_"+shortString).c_str());
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: effDataset_"+shortString+".root"<<endl;
    return;
  }
  RooDataSet* totdata    = (RooDataSet*)wsp->data(("data_"   +shortString).c_str());
  RooDataSet* totnumData = (RooDataSet*)wsp->data(("numData_"+shortString).c_str());
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  RooArgSet vars (* ctK,* ctL,* phi);
  // split dataset according to job number
  double denN = totdata   ->numEntries();
  double numN = totnumData->numEntries();
  cout<<"Global efficiency: "<<numN<<" / "<<denN<<" = "<<numN/denN<<endl;
  cout<<"Partition "<<ndiv<<" of "<<totdiv
      <<": den "<<(int)((ndiv*denN)/totdiv)<<"->"<<((int)(((ndiv+1)*denN)/totdiv))-1
      << " num "<<(int)((ndiv*numN)/totdiv)<<"->"<<((int)(((ndiv+1)*numN)/totdiv))-1<<endl;
  RooAbsData* data    = totdata   ->reduce(EventRange((int)((ndiv*denN)/totdiv),((int)(((ndiv+1)*denN)/totdiv))-1));
  RooAbsData* numData = totnumData->reduce(EventRange((int)((ndiv*numN)/totdiv),((int)(((ndiv+1)*numN)/totdiv))-1));

  // create the KDE description of numerator and denominator datasets
  cout<<"Start creating denominator KDE"<<endl;
  TStopwatch t;
  t.Start();
  RooNDKeysPdf* denKDE = new RooNDKeysPdf("denKDE","denKDE",RooArgList(*ctK,*ctL,*phi),*data   ,"m",width);
  cout<<"Start creating numerator KDE"<<endl;
  RooNDKeysPdf* numKDE = new RooNDKeysPdf("numKDE","numKDE",RooArgList(*ctK,*ctL,*phi),*numData,"m",width);
  t.Stop();
  t.Print();

  // create numerator and denominator histograms
  cout<<"Start creating denominator histogram"<<endl;
  TStopwatch t1;
  t1.Start();
  TH3D* denHist = (TH3D*)denKDE->createHistogram( ("denHist"+shortString).c_str(),
  						   *ctK,     Binning(xbins,-1,1) ,
  						   YVar(*ctL,Binning(ybins,-1,1)),
  						   ZVar(*phi,Binning(zbins,-TMath::Pi(),TMath::Pi())),
						   Scaling(false) );
  cout<<"Start creating numerator histogram"<<endl;
  TH3D* numHist = (TH3D*)numKDE->createHistogram( ("numHist"+shortString).c_str(),
  						   *ctK,     Binning(xbins,-1,1) ,
  						   YVar(*ctL,Binning(ybins,-1,1)),
  						   ZVar(*phi,Binning(zbins,-TMath::Pi(),TMath::Pi())),
						   Scaling(false) );
  t1.Stop();
  t1.Print();
  denHist->Scale(   data->numEntries()/denHist->Integral());
  numHist->Scale(numData->numEntries()/numHist->Integral());

  // save histograms to file
  RooWorkspace *wsp_out = new RooWorkspace("ws","Workspace with efficiency parameterisation");
  wsp_out->import( *denHist );
  wsp_out->import( *numHist );
  wsp_out->import( vars, Silence() );
  wsp_out->writeToFile( ( "effKDE_"+shortString+Form("_rooKeys_mw%.2f_%i_%i_%i_%i-frac-%i.root",width,xbins,ybins,zbins,ndiv,totdiv)).c_str() );

}
