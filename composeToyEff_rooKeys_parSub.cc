#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TStopwatch.h>

#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooRandom.h>
#include <RooNDKeysPdf.h>
#include <RooBinning.h>

using namespace RooFit;
using namespace std;

static const int nBins = 9;

void composeToyEff_rooKeys_parSub(int q2Bin, int effIndx, int parity, int seed, float widthCTK, float widthCTL, float widthPHI, int xbins, int ybins, int zbins, int ndiv, int totdiv, int year, int vers)
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

  // Load variables and dataset
  string fileDir = "/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta/";
  string filename_data = fileDir + Form("effDatasetTheta_b%i_%i.root",q2Bin,year);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,parity));
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename_data<<endl;
    return;
  }
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
    cout<<"Variables not found in file: "<<filename_data<<endl;
    return;
  }
  RooArgList vars (* ctK,* ctL,* phi);
  auto originData = (RooDataSet*)wsp->data(datasetString.c_str());
  if ( !originData || originData->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename_data<<endl;
    return;
  }
  // Get full number of events in GEN denominator sample (with no FSR veto)
  double normVal = 1;
  if ( effIndx==0 ) {
    TH1I* normHist = (TH1I*)fin_data->Get( Form("n_genDen_b%i",q2Bin) );
    if ( !normHist || normHist->IsZombie() ) cout<<"Histogram n_genDen_b"<<q2Bin<<" not found in file: "<<filename_data<<"\nUsing the statistics of the sample for function normalisation"<<endl;
    else normVal = normHist->GetBinContent(parity+1);
  }
  
  // import KDE histograms
  vector<TH3D*> KDEhists;
  vector<TString> KDEconfs;
  string inFileDir = "/lstore/cms/boletti/Run2-BdToKstarMuMu/eff-KDE-theta/files/";
  string inFileName = inFileDir + Form((parity==0?"KDEhist_b%i_ev_%i_v%i.root":"KDEhist_b%i_od_%i_v%i.root"),q2Bin,year,vers);
  TFile* fin = TFile::Open( inFileName.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<inFileName<<endl;
    return;
  }
  auto inFunc = (TH3D*)fin->Get( Form("histTheta_indx%i_w0-%1.2f_w1-%1.2f_w2-%1.2f_%i_%i_%i",effIndx,widthCTK,widthCTL,widthPHI,xbins,ybins,zbins) );

  auto inFunc_dh = new RooDataHist("inFunc_dh","inFunc_dh",vars,inFunc);
  auto inFunc_hp = new RooHistPdf("inFunc_hp","inFunc_hp",vars,*inFunc_dh,1);

  int nEv = originData->numEntries()/totdiv;
  RooRandom::randomGenerator()->SetSeed(1+ndiv+totdiv*seed);
  auto data = inFunc_hp->generate(vars,nEv,Name("data"));
  
  
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

  vector<Double_t> xboundaries (xbins+1);
  vector<Double_t> yboundaries (ybins+1);
  for (int i=0; i<=xbins; ++i)
    xboundaries[i] = TMath::ACos(1.0-2.0*i/xbins);
  for (int i=0; i<=ybins; ++i)
    yboundaries[i] = TMath::ACos(1.0-2.0*i/ybins);
  RooBinning cosBinX ( xbins, &xboundaries[0], "cosBinX" );
  RooBinning cosBinY ( ybins, &yboundaries[0], "cosBinY" );

  // sample the KDE function to save it in a file (this is the most time-consuming process)
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
  if ( effIndx==0 ) KDEhist->Scale(normVal/totdiv/KDEhist->Integral());
  else KDEhist->Scale(data->sumEntries()/KDEhist->Integral());

  // save histogram to file
  cout<<"out open"<<endl;
  string foutName = fileDir + Form("tmptoy_v%i/KDEhistTheta_",vers) + shortString + Form("_rooKeys_m_w0-%.2f_w1-%.2f_w2-%.2f_%i_%i_%i_%i-frac-%i_%i_toy%i.root",widthCTK,widthCTL,widthPHI,xbins,ybins,zbins,ndiv,totdiv,year,seed);
  TFile* fout = TFile::Open(foutName.c_str(),"RECREATE");
  fout->cd();
  cout<<"out write"<<endl;
  KDEhist->Write();
  cout<<"out close"<<endl;
  fout->Close();

  cout<<"del KDEfunc"<<endl;
  delete KDEfunc;
  cout<<"del data"<<endl;
  delete data;
  cout<<"in close"<<endl;
  fin->Close();  
  cout<<"data close"<<endl;
  fin_data->Close();  
  cout<<"done"<<endl;

  return;
  
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency vs. odd data
  //                [1] odd efficiency vs. even data
  //                [-1] for each parity recursively

  int q2Bin  = 0;
  int effIndx= 0;
  int parity = 0;
  int seed = 1;
  float widthCTK = 0.5;
  float widthCTL = 0.5;
  float widthPHI = 0.5;
  int xbins = 50;
  int ybins = 0;
  int zbins = 0;
  int ndiv = 0;
  int totdiv = 1;
  int year   = 2016;
  int vers   = -1;

  if ( argc > 1 ) q2Bin    = atoi(argv[1]);
  if ( argc > 2 ) effIndx  = atoi(argv[2]);
  if ( argc > 3 ) parity   = atoi(argv[3]);
  if ( argc > 4 ) widthCTK = atof(argv[4]);
  if ( argc > 5 ) widthCTL = atof(argv[5]);
  if ( argc > 6 ) widthPHI = atof(argv[6]);
  if ( argc > 7 ) xbins    = atoi(argv[7]);
  if ( argc > 8 ) ybins    = atoi(argv[8]);
  if ( argc > 9 ) zbins    = atoi(argv[9]);
  if ( argc > 11 ) totdiv  = atoi(argv[11]);
  if ( argc > 12 ) year    = atoi(argv[12]);
  if ( argc > 13 ) vers    = atoi(argv[13]);

  if ( argc > 10 ) {
    ndiv = atoi(argv[10]) % totdiv;
    seed = atoi(argv[10]) / totdiv;
  }

  if ( seed<0 ) return 1;

  if ( xbins<1 ) return 1;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( widthCTK<=0 ) return 1;
  if ( widthCTL<=0 ) return 1;
  if ( widthPHI<=0 ) return 1;

  if ( ndiv<0 || ndiv>=totdiv ) return 1;

  if ( q2Bin<0 || q2Bin>=nBins ) return 1;

  if ( effIndx<0 || effIndx>5 ) return 1;

  if ( parity<0 || parity>1 ) return 1;

  composeToyEff_rooKeys_parSub(q2Bin, effIndx, parity, seed, widthCTK, widthCTL, widthPHI, xbins, ybins, zbins, ndiv, totdiv, year, vers);

  return 0;

}
