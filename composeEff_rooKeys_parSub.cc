#include "TFile.h"
#include "RooWorkspace.h"
#include "TStopwatch.h"
#include "TH3D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TVectorD.h"
#include "RooNDKeysPdf.h"
#include "TMath.h"
#include "RooBinning.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

void composeEff_rooKeys_parSub(int q2Bin, int effIndx, int parity, float widthCTK = 0.5, float widthCTL = 0.5, float widthPHI = 0.5, int xbins=50, int ybins = 0, int zbins = 0, int ndiv = 0, int totdiv = 1, int year=2016, int vers=-1)
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

  string XGBstr = "";
  if (vers>9 && vers<100) XGBstr = Form("_XGBv%i",vers/10);
  else if (vers>99) XGBstr = Form("_TMVAv%i",(vers-100)/10);

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
  string fileDir = "/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta/";
  string filename = fileDir + Form("effDatasetTheta_b%i_%i%s.root",q2Bin,year,XGBstr.c_str());
  if (q2Bin==4 && effIndx==2)
    filename = fileDir + Form("effDatasetTheta_b%i_%i%s_den.root",q2Bin,year,XGBstr.c_str());
  TFile* fin = TFile::Open( filename.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin->Get(Form("ws_b%ip%i",q2Bin,parity));
  RooWorkspace* wsp2= (RooWorkspace*)fin->Get(Form("ws2_b%ip%i",q2Bin,parity));
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename<<endl;
    return;
  }
  RooRealVar* ctK = wsp->var("tK");
  RooRealVar* ctL = wsp->var("tL");
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
  if (effIndx==2 && wsp2) {
    totdata->append(*((RooDataSet*)wsp2->data(Form((parity==0?"data_den2_ev_b%i":"data_den2_od_b%i"),q2Bin))));
  }
  auto den2ds = (RooDataSet*)wsp->data(Form((parity==0?"data_den2_ev_b%i":"data_den2_od_b%i"),q2Bin));
  if (effIndx==2 && den2ds) totdata->append(*den2ds);
  // full reconstructed sample is obtained appending wrong-tag and correct-tag events
  if ( effIndx==5 ) {
    RooDataSet* extradata = (RooDataSet*)wsp->data(Form((parity==0?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin));
    if ( !extradata || extradata->IsZombie() ) {
      cout<<"Dataset data_ctRECO_"<<(parity==0?"ev":"od")<<"_b"<<q2Bin<<" not found in file: "<<filename<<endl;
      return;
    }
    totdata->append(*extradata);
  }
  // Get full number of events in GEN denominator sample (with no FSR veto)
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

  // create vector containing width scale factors, including a term to make the kernel width independent from the sample's statistics
  // (the inverse of this factor is used in the RooNDKeysPdf class definition)
  TVectorD rho (3);
  rho[0] = widthCTK * pow(data->sumEntries()/10000,1./7.);
  rho[1] = widthCTL * pow(data->sumEntries()/10000,1./7.);
  rho[2] = widthPHI * pow(data->sumEntries()/10000,1./7.);

  // create the KDE description
  TStopwatch t;
  t.Start();
  RooNDKeysPdf* KDEfunc = new RooNDKeysPdf("KDEfunc","KDEfunc",RooArgList(*ctK,*ctL,*phi),*data,rho,"m",3,false);
  t.Stop();
  t.Print();

  // sample the KDE function to save it in a file (this is the most time-consuming process)
  if (xbins!=ybins || xbins%2==1) {
    cout<<"xbins!=ybins or xbins%2==1 not implemented. Abort!"<<endl;
    return;
  }
  vector<Double_t> thetaCentre (xbins/2);
  vector<Double_t> xboundaries (xbins+3);
  // vector<Double_t> yboundaries (ybins+3);
  xboundaries[0] = 0;
  xboundaries[xbins+2] = TMath::Pi();
  xboundaries[xbins/2+1] = TMath::Pi()/2;
  for (int i=0; i<xbins/2; ++i) {
    thetaCentre[i] = TMath::ACos(-1.0*(2.0*i+1)/xbins);
    xboundaries[xbins/2+i+2] = 2*thetaCentre[i] - xboundaries[xbins/2+i+1];
    xboundaries[xbins/2-i] = TMath::Pi() - xboundaries[xbins/2+i+2];
  }
  // for (int i=0; i<=xbins; ++i)
  //   xboundaries[i] = TMath::ACos(1.0-2.0*i/xbins);
  // for (int i=0; i<=ybins; ++i)
  //   yboundaries[i] = TMath::ACos(1.0-2.0*i/ybins);
  RooBinning cosBinX ( xbins+2, &xboundaries[0], "cosBinX" );
  RooBinning cosBinY ( xbins+2, &xboundaries[0], "cosBinY" );

  TStopwatch t1;
  t1.Start();
  TH3D* KDEhist = (TH3D*)KDEfunc->createHistogram( ("KDEhist_"+shortString).c_str(),
  						   *ctK,     Binning(cosBinX) ,
  						   YVar(*ctL,Binning(cosBinY)),
  						   ZVar(*phi,Binning(zbins,-TMath::Pi(),TMath::Pi())),
						   Scaling(false) );
  t1.Stop();
  t1.Print();
  // scale the output to have the same normalisation as the input dataset (or as the full sample, with no FSR veto, for GEN-denominator)
  cout<<"Normalizations: numEntries = "<<data->numEntries()<<" sumEntries = "<<data->sumEntries()<<" raw hist integral = "<<KDEhist->Integral()<<endl;
  KDEhist->Scale(data->sumEntries()/KDEhist->Integral());
  if ( effIndx==0 ) KDEhist->Scale(1.0*normVal/totNev);

  // save histogram to file
  string foutName = fileDir + Form("tmp_v%i/KDEhistTheta_",vers) + shortString + Form("_rooKeys_m_w0-%.2f_w1-%.2f_w2-%.2f_%i_%i_%i_%i-frac-%i_%i.root",widthCTK,widthCTL,widthPHI,xbins,ybins,zbins,ndiv,totdiv,year);
  TFile* fout = TFile::Open(foutName.c_str(),"RECREATE");
  fout->cd();
  KDEhist->Write();
  fout->Close();

  delete KDEfunc;
  delete data;
  fin->Close();  

  return;

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
  int vers = -1;

  if ( argc > 1 ) q2Bin = atoi(argv[1]);
  if ( argc > 2 ) effIndx = atoi(argv[2]);
  if ( argc > 3 ) parity = atoi(argv[3]);
  if ( argc > 4 ) widthCTK = atof(argv[4]);
  if ( argc > 5 ) widthCTL = atof(argv[5]);
  if ( argc > 6 ) widthPHI = atof(argv[6]);
  if ( argc > 7 ) xbins = atoi(argv[7]);
  if ( argc > 8 ) ybins = atoi(argv[8]);
  if ( argc > 9 ) zbins = atoi(argv[9]);
  if ( argc > 10 ) ndiv = atoi(argv[10]);
  if ( argc > 11 ) totdiv = atoi(argv[11]);
  if ( argc > 12 ) year = atoi(argv[12]);
  if ( argc > 13 ) vers = atoi(argv[13]);

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  <  0 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      composeEff_rooKeys_parSub(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, ndiv, totdiv, year, vers);
  else
    composeEff_rooKeys_parSub(q2Bin, effIndx, parity, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, ndiv, totdiv, year, vers);
    
  return 0;

}
