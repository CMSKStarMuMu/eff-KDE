#include "TFile.h"
#include "RooWorkspace.h"
#include "TStopwatch.h"
#include "TH3D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TVectorD.h"
#include "RooNDKeysPdf.h"
#include "TMath.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

void composeEff_rooKeys_parSub(int q2Bin, int effIndx, int parity, float widthCTK = 0.5, float widthCTL = 0.5, float widthPHI = 0.5, int xbins=50, int ybins = 0, int zbins = 0, int ndiv = 0, int totdiv = 1, int year=2016)
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
  string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE-Swave/effDataset_b%i_%i.root",q2Bin,year);
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
  // full reconstructed sample is obtained appending wrong-tag and correct-tag events
  if ( effIndx==5 ) {
    RooDataSet* extradata = (RooDataSet*)wsp->data(Form((parity==0?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin));
    if ( !extradata || extradata->IsZombie() ) {
      cout<<"Dataset data_ctRECO_"<<(parity==0?"ev":"od")<<"_b"<<q2Bin<<" not found in file: "<<filename<<endl;
      return;
    }
    totdata->append(*extradata);
  }

  // split dataset according to job number
  // caveat: in the "EventRange" option, the end of the range is not included
  double totNev = totdata->numEntries();
  cout<<"Partition "<<ndiv<<" of "<<totdiv<<": "
      <<(int)((ndiv*totNev)/totdiv)<<"->"
      <<((int)(((ndiv+1)*totNev)/totdiv))-1<<" of total "<<totNev<<endl;
  RooAbsData* data = totdata->reduce(EventRange((int)((ndiv*totNev)/totdiv),(int)(((ndiv+1)*totNev)/totdiv)));

  // create vector containing width scale factors, including a term to make the kernel width independent from the sample's statistics
  // (the inverse of this factor is used in the RooNDKeysPdf class definition)
  TVectorD rho (3);
  rho[0] = widthCTK * pow(data->sumEntries()/10000,1./7.);
  rho[1] = widthCTL * pow(data->sumEntries()/10000,1./7.);
  rho[2] = widthPHI * pow(data->sumEntries()/10000,1./7.);

  // create the KDE description
  new RooNDKeysPdf("KDEfunc","KDEfunc",RooArgList(*ctK,*ctL,*phi),*data,rho,"m",3,false);

}

int main(int argc, char** argv)
{

  int q2Bin = -1;
  int effIndx = 0;
  int parity = 0;
  float widthCTK = 0.5;
  float widthCTL = 0.5;
  float widthPHI = 0.5;
  int xbins = 50;
  int ybins = 50;
  int zbins = 50;
  int ndiv = 50;
  int totdiv = 1;
  int year = 2016;

  if ( argc >= 2 ) q2Bin = atoi(argv[1]);
  if ( argc >= 3 ) effIndx = atoi(argv[2]);
  if ( argc >= 4 ) parity = atoi(argv[3]);
  if ( argc >= 5 ) widthCTK = atof(argv[4]);
  if ( argc >= 6 ) widthCTL = atof(argv[5]);
  if ( argc >= 7 ) widthPHI = atof(argv[6]);
  if ( argc >= 8 ) xbins = atoi(argv[7]);
  if ( argc >= 9 ) ybins = atoi(argv[8]);
  if ( argc >= 10 ) zbins = atoi(argv[9]);
  if ( argc >= 11 ) ndiv = atoi(argv[10]);
  if ( argc >= 12 ) totdiv = atoi(argv[11]);
  if ( argc >= 13 ) year = atoi(argv[12]);

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  <  0 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      composeEff_rooKeys_parSub(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, ndiv, totdiv, year);
  else
    composeEff_rooKeys_parSub(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, ndiv, totdiv, year);
    
  return 0;

}
