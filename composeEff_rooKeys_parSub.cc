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

void composeEff_rooKeys_parSub(int q2Bin, int effIndx, int parity, float widthCTK = 0.5, float widthCTL = 0.5, float widthPHI = 0.5, int xbins=50, int ybins = 0, int zbins = 0, int ndiv = 0, int totdiv = 1)
{
  // effIndx format: [0] GEN no-filter
  //                 [1] GEN filtered
  //                 [2] GEN filtered from full MC sample
  //                 [3] correct-tag RECO candidates
  //                 [4] wrong-tag RECO candidates
  //                 [5] RECO candidates
  // parity format: [0] even
  //                [1] odd

  if ( xbins<1 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( widthCTK<=0 ) return;
  if ( widthCTL<=0 ) return;
  if ( widthPHI<=0 ) return;

  if ( ndiv<0 || ndiv>=totdiv ) return;

  if ( q2Bin<0 || q2Bin>=nBins ) return;

  if ( effIndx<0 || effIndx>5 ) return;

  if ( parity<0 || parity>1 ) return;

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

  // Load variables and dataset
  string filename = Form("effDataset_b%i.root",q2Bin);
  TFile* fin = TFile::Open( filename.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin->Get(Form("ws_b%ip%i",q2Bin,parity));
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename<<endl;
    return;
  }
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
    cout<<"Variables not found in file: "<<filename<<endl;
    return;
  }
  RooDataSet* totdata = (RooDataSet*)wsp->data(datasetString.c_str());
  if ( !totdata || totdata->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename<<endl;
    return;
  }
  if ( effIndx==5 ) {
    RooDataSet* extradata = (RooDataSet*)wsp->data(Form((parity==0?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin));
    if ( !extradata || extradata->IsZombie() ) {
      cout<<"Dataset data_ctRECO_"<<(parity==0?"ev":"od")<<"_b"<<q2Bin<<" not found in file: "<<filename<<endl;
      return;
    }
    totdata->append(*extradata);
  }
  double normVal = 1;
  if ( effIndx==0 ) {
    TH1I* normHist = (TH1I*)fin->Get( Form("n_genDen_b%i",q2Bin) );
    if ( !normHist || normHist->IsZombie() ) cout<<"Histogram n_genDen_b"<<q2Bin<<" not found in file: "<<filename<<"\nUsing the statistics of the sample for function normalisation"<<endl;
    else normVal = normHist->GetBinContent(parity+1);
  }

  // split dataset according to job number
  // caveat: in the "EventRange" option, the end of the range is not included
  double totNev = totdata->numEntries();
  cout<<"Partition "<<ndiv<<" of "<<totdiv<<": "
      <<(int)((ndiv*totNev)/totdiv)<<"->"
      <<((int)(((ndiv+1)*totNev)/totdiv))-1<<" of total "<<totNev<<endl;
  RooAbsData* data = totdata->reduce(EventRange((int)((ndiv*totNev)/totdiv),(int)(((ndiv+1)*totNev)/totdiv)));

  // create vector containing width scale factors
  TVectorD rho (3);
  rho[0] = widthCTK * pow(data->sumEntries()/10000,1./7.);
  rho[1] = widthCTL * pow(data->sumEntries()/10000,1./7.);
  rho[2] = widthPHI * pow(data->sumEntries()/10000,1./7.);

  // create the KDE description of numerator and denominator datasets
  TStopwatch t;
  t.Start();
  RooNDKeysPdf* KDEfunc = new RooNDKeysPdf("KDEfunc","KDEfunc",RooArgList(*ctK,*ctL,*phi),*data,rho,"m",3,false);
  t.Stop();
  t.Print();

  // create numerator and denominator histograms
  TStopwatch t1;
  t1.Start();
  TH3D* KDEhist = (TH3D*)KDEfunc->createHistogram( ("KDEhist_"+shortString).c_str(),
  						   *ctK,     Binning(xbins,-1,1) ,
  						   YVar(*ctL,Binning(ybins,-1,1)),
  						   ZVar(*phi,Binning(zbins,-TMath::Pi(),TMath::Pi())),
						   Scaling(false) );
  t1.Stop();
  t1.Print();
  if ( effIndx==0 ) KDEhist->Scale(normVal/totdiv/KDEhist->Integral());
  else KDEhist->Scale(data->sumEntries()/KDEhist->Integral());

  // save histogram to file
  TFile* fout = TFile::Open(Form("KDEhist_%s_rooKeys_m_w0-%.2f_w1-%.2f_w2-%.2f_%i_%i_%i_%i-frac-%i.root",
				 shortString.c_str(),widthCTK,widthCTL,widthPHI,xbins,ybins,zbins,ndiv,totdiv),"RECREATE");
  fout->cd();
  KDEhist->Write();
  fout->Close();

}
