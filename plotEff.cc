#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TStyle.h>

#include <RooRealVar.h>
#include <RooAbsReal.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* csxC [2*nBins];
TCanvas* csyC [2*nBins];
TCanvas* cszC [2*nBins];
TCanvas* cp1C [2*nBins];
TCanvas* cp2C [2*nBins];
TCanvas* cctC [2*nBins];
TCanvas* csxW [2*nBins];
TCanvas* csyW [2*nBins];
TCanvas* cszW [2*nBins];
TCanvas* cp1W [2*nBins];
TCanvas* cp2W [2*nBins];
TCanvas* cctW [2*nBins];

TH1D* Invert1Dhist(TH1D* hin, string hname);

void plotEffBin(int q2Bin, int parity, bool doClosure)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  bool doCT = true;
  bool doWT = true;

  int confIndex = nBins*parity + q2Bin;

  int distBins = 10;

  // Load variables and dataset
  string filename_data = Form("effDataset_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  // import the both datasets
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,1-parity));
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename_data<<endl;
    return;
  }
  RooWorkspace* wsp_corr = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,parity));
  if ( !wsp_corr || wsp_corr->IsZombie() ) {
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
  string datasetString = Form((parity==1?"_ev_b%i":"_od_b%i"),q2Bin);
  RooDataSet* data_genDen = (RooDataSet*)wsp->data(("data_genDen"+datasetString).c_str());
  RooDataSet* data_genNum = (RooDataSet*)wsp->data(("data_genNum"+datasetString).c_str());
  RooDataSet* data_den    = (RooDataSet*)wsp->data(("data_den"   +datasetString).c_str());
  RooDataSet* data_ctRECO = (RooDataSet*)wsp->data(("data_ctRECO"+datasetString).c_str());
  RooDataSet* data_wtRECO = (RooDataSet*)wsp->data(("data_wtRECO"+datasetString).c_str());
  datasetString = Form((parity==0?"_ev_b%i":"_od_b%i"),q2Bin);
  RooDataSet* data_genDen_corr = (RooDataSet*)wsp_corr->data(("data_genDen"+datasetString).c_str());
  RooDataSet* data_genNum_corr = (RooDataSet*)wsp_corr->data(("data_genNum"+datasetString).c_str());
  RooDataSet* data_den_corr    = (RooDataSet*)wsp_corr->data(("data_den"   +datasetString).c_str());
  RooDataSet* data_ctRECO_corr = (RooDataSet*)wsp_corr->data(("data_ctRECO"+datasetString).c_str());
  RooDataSet* data_wtRECO_corr = (RooDataSet*)wsp_corr->data(("data_wtRECO"+datasetString).c_str());
  
  // import KDE efficiency histograms
  string filename = Form((parity==0?"files/KDEeff_b%i_ev.root":"files/KDEeff_b%i_od.root"),q2Bin);
  TFile* fin = new TFile( filename.c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  TH3D* effCHist = (TH3D*)fin->Get(("effCHist_"+shortString).c_str());
  TH3D* effWHist = (TH3D*)fin->Get(("effWHist_"+shortString).c_str());
  if ( !effCHist || effCHist->IsZombie() ) {
    cout<<"Correct-tag efficiency histograms not found in file: "<<filename<<endl;
    doCT = false;
  }
  if ( !effWHist || effWHist->IsZombie() ) {
    cout<<"Wrong-tag efficiency histograms not found in file: "<<filename<<endl;
    doWT = false;
  }
  if ( !doCT && !doWT ) return;
  // create efficiency functions
  RooDataHist* effCData = 0;
  RooDataHist* effWData = 0;
  RooAbsReal* effC = 0;
  RooAbsReal* effW = 0;
  if (doCT) {
    effCData = new RooDataHist("effCData","effCData",vars,effCHist);
    effC = new RooHistFunc("effC","effC",vars,*effCData,1);
  }
  if (doWT) {
    effWData = new RooDataHist("effWData","effWData",vars,effWHist);
    effW = new RooHistFunc("effW","effW",vars,*effWData,1);
  }

  string confString = "plotEff_d/effKDE_test_"+shortString;

  gStyle->SetOptStat(0);

  // Plot 1D slices of the efficiency function and binned efficiency
  if (doCT) {
    csxC[confIndex] = new TCanvas(("csxC"+shortString).c_str(),(shortString+"_effC_ctK").c_str(),1500,1500) ;
    csyC[confIndex] = new TCanvas(("csyC"+shortString).c_str(),(shortString+"_effC_ctL").c_str(),1500,1500) ;
    cszC[confIndex] = new TCanvas(("cszC"+shortString).c_str(),(shortString+"_effC_phi").c_str(),1500,1500) ;
    csxC[confIndex]->Divide(5,5);
    csyC[confIndex]->Divide(5,5);
    cszC[confIndex]->Divide(5,5);
  }
  if (doWT) {
    csxW[confIndex] = new TCanvas(("csxW"+shortString).c_str(),(shortString+"_effW_ctK").c_str(),1500,1500) ;
    csyW[confIndex] = new TCanvas(("csyW"+shortString).c_str(),(shortString+"_effW_ctL").c_str(),1500,1500) ;
    cszW[confIndex] = new TCanvas(("cszW"+shortString).c_str(),(shortString+"_effW_phi").c_str(),1500,1500) ;
    csxW[confIndex]->Divide(5,5);
    csyW[confIndex]->Divide(5,5);
    cszW[confIndex]->Divide(5,5);
  }
  vector <TH1D*> effCSliceX;
  vector <TH1D*> effCSliceY;
  vector <TH1D*> effCSliceZ;
  vector <TH1D*> effWSliceX;
  vector <TH1D*> effWSliceY;
  vector <TH1D*> effWSliceZ;
  vector <TH1D*> effCSliceX_corr;
  vector <TH1D*> effCSliceY_corr;
  vector <TH1D*> effCSliceZ_corr;
  vector <TH1D*> effWSliceX_corr;
  vector <TH1D*> effWSliceY_corr;
  vector <TH1D*> effWSliceZ_corr;
  vector <RooPlot*> fsxC;
  vector <RooPlot*> fsxW;
  vector <RooPlot*> fsyC;
  vector <RooPlot*> fsyW;
  vector <RooPlot*> fszC;
  vector <RooPlot*> fszW;

  // TLegend* leg = new TLegend (0.35,0.8,0.9,0.9);

  // width of the slices in the hidden variables ("border" is half of it)
  double border = 0.04;

  // variables to be filled with global efficiency maximum
  double maxEffCX = 0;
  double maxEffWX = 0;
  double maxEffCY = 0;
  double maxEffWY = 0;
  double maxEffCZ = 0;
  double maxEffWZ = 0;

  // loop over slice grid
  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

      // central values and borders of the slices in the hidden variables
      double centA = -0.8 + 1.6*i/4;
      double centB = -0.8 + 1.6*j/4;
      string cutX = Form("fabs(ctL-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centA,border,TMath::Pi(),centB,border);
      string cutY = Form("fabs(ctK-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centA,border,TMath::Pi(),centB,border);
      string cutZ = Form("fabs(ctK-(%1.1f))<%1.3f && fabs(ctL-(%1.1f))<%1.3f",centA,border,centB,border);
      // cout<<cutX<<endl<<cutY<<endl<<cutZ<<endl;

      // slicing distributions
      TH1D* genDenDistX = (TH1D*)data_genDen->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_genDen_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* genNumDistX = (TH1D*)data_genNum->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_genNum_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* denDistX    = (TH1D*)data_den   ->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_den_%i_%i_%s"   ,i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* genDenDistY = (TH1D*)data_genDen->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_genDen_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* genNumDistY = (TH1D*)data_genNum->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_genNum_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* denDistY    = (TH1D*)data_den   ->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_den_%i_%i_%s"   ,i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* genDenDistZ = (TH1D*)data_genDen->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_genDen_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* genNumDistZ = (TH1D*)data_genNum->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_genNum_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* denDistZ    = (TH1D*)data_den   ->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_den_%i_%i_%s"   ,i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* ctRECODistX = 0;
      TH1D* wtRECODistX = 0;
      TH1D* ctRECODistY = 0;
      TH1D* wtRECODistY = 0;
      TH1D* ctRECODistZ = 0;
      TH1D* wtRECODistZ = 0;
      if (doCT) {
	ctRECODistX = (TH1D*)data_ctRECO->reduce(RooArgSet(*ctK),cutX.c_str())
	  ->createHistogram( Form("distX_ctRECO_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
	ctRECODistY = (TH1D*)data_ctRECO->reduce(RooArgSet(*ctL),cutY.c_str())
	  ->createHistogram( Form("distY_ctRECO_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
	ctRECODistZ = (TH1D*)data_ctRECO->reduce(RooArgSet(*phi),cutZ.c_str())
	  ->createHistogram( Form("distZ_ctRECO_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      }
      if (doWT) {
	wtRECODistX = (TH1D*)data_wtRECO->reduce(RooArgSet(*ctK),cutX.c_str())
	  ->createHistogram( Form("distX_wtRECO_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
	wtRECODistY = (TH1D*)data_wtRECO->reduce(RooArgSet(*ctL),cutY.c_str())
	  ->createHistogram( Form("distY_wtRECO_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
	wtRECODistZ = (TH1D*)data_wtRECO->reduce(RooArgSet(*phi),cutZ.c_str())
	  ->createHistogram( Form("distZ_wtRECO_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      }
      // slicing correlated distributions
      TH1D* genDenDistX_corr = (TH1D*)data_genDen_corr->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_genDen_corr_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* genNumDistX_corr = (TH1D*)data_genNum_corr->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_genNum_corr_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* denDistX_corr    = (TH1D*)data_den_corr   ->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_den_corr_%i_%i_%s"   ,i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* genDenDistY_corr = (TH1D*)data_genDen_corr->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_genDen_corr_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* genNumDistY_corr = (TH1D*)data_genNum_corr->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_genNum_corr_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* denDistY_corr    = (TH1D*)data_den_corr   ->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_den_corr_%i_%i_%s"   ,i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* genDenDistZ_corr = (TH1D*)data_genDen_corr->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_genDen_corr_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* genNumDistZ_corr = (TH1D*)data_genNum_corr->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_genNum_corr_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* denDistZ_corr    = (TH1D*)data_den_corr   ->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_den_corr_%i_%i_%s"   ,i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* ctRECODistX_corr = 0;
      TH1D* wtRECODistX_corr = 0;
      TH1D* ctRECODistY_corr = 0;
      TH1D* wtRECODistY_corr = 0;
      TH1D* ctRECODistZ_corr = 0;
      TH1D* wtRECODistZ_corr = 0;
      if (doCT) {
	ctRECODistX_corr = (TH1D*)data_ctRECO_corr->reduce(RooArgSet(*ctK),cutX.c_str())
	  ->createHistogram( Form("distX_ctRECO_corr_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
	ctRECODistY_corr = (TH1D*)data_ctRECO_corr->reduce(RooArgSet(*ctL),cutY.c_str())
	  ->createHistogram( Form("distY_ctRECO_corr_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
	ctRECODistZ_corr = (TH1D*)data_ctRECO_corr->reduce(RooArgSet(*phi),cutZ.c_str())
	  ->createHistogram( Form("distZ_ctRECO_corr_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      }
      if (doWT) {
	wtRECODistX_corr = (TH1D*)data_wtRECO_corr->reduce(RooArgSet(*ctK),cutX.c_str())
	  ->createHistogram( Form("distX_wtRECO_corr_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
	wtRECODistY_corr = (TH1D*)data_wtRECO_corr->reduce(RooArgSet(*ctL),cutY.c_str())
	  ->createHistogram( Form("distY_wtRECO_corr_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
	wtRECODistZ_corr = (TH1D*)data_wtRECO_corr->reduce(RooArgSet(*phi),cutZ.c_str())
	  ->createHistogram( Form("distZ_wtRECO_corr_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      }
      // composing binned efficiencies from sliced distributions
      TH1D* factDistX = (TH1D*)genNumDistX->Clone(Form("factDistX_%i_%i_%s",i,j,shortString.c_str()));
      factDistX->Divide(genDenDistX);
      factDistX->Divide(denDistX);
      TH1D* effCDistX = (TH1D*)factDistX->Clone(Form("effCDistX_%i_%i_%s",i,j,shortString.c_str()));
      if (doCT) effCDistX->Multiply(ctRECODistX);
      TH1D* effWDistX = Invert1Dhist(factDistX,Form("effWDistX_%i_%i_%s",i,j,shortString.c_str()));
      if (doWT) effWDistX->Multiply(wtRECODistX);
      TH1D* factDistY = (TH1D*)genNumDistY->Clone(Form("factDistY_%i_%i_%s",i,j,shortString.c_str()));
      factDistY->Divide(genDenDistY);
      factDistY->Divide(denDistY);
      TH1D* effCDistY = (TH1D*)factDistY->Clone(Form("effCDistY_%i_%i_%s",i,j,shortString.c_str()));
      if (doCT) effCDistY->Multiply(ctRECODistY);
      TH1D* effWDistY = Invert1Dhist(factDistY,Form("effWDistY_%i_%i_%s",i,j,shortString.c_str()));
      if (doWT) effWDistY->Multiply(wtRECODistY);
      TH1D* factDistZ = (TH1D*)genNumDistZ->Clone(Form("factDistZ_%i_%i_%s",i,j,shortString.c_str()));
      factDistZ->Divide(genDenDistZ);
      factDistZ->Divide(denDistZ);
      TH1D* effCDistZ = (TH1D*)factDistZ->Clone(Form("effCDistZ_%i_%i_%s",i,j,shortString.c_str()));
      if (doCT) effCDistZ->Multiply(ctRECODistZ);
      TH1D* effWDistZ = Invert1Dhist(factDistZ,Form("effWDistZ_%i_%i_%s",i,j,shortString.c_str()));
      if (doWT) effWDistZ->Multiply(wtRECODistZ);
      // composing binned efficiencies from sliced distributions of correlated dataset
      TH1D* factDistX_corr = (TH1D*)genNumDistX_corr->Clone(Form("factDistX_corr_%i_%i_%s",i,j,shortString.c_str()));
      factDistX_corr->Divide(genDenDistX_corr);
      factDistX_corr->Divide(denDistX_corr);
      TH1D* effCDistX_corr = (TH1D*)factDistX_corr->Clone(Form("effCDistX_corr_%i_%i_%s",i,j,shortString.c_str()));
      if (doCT) effCDistX_corr->Multiply(ctRECODistX_corr);
      TH1D* effWDistX_corr = Invert1Dhist(factDistX_corr,Form("effWDistX_corr_%i_%i_%s",i,j,shortString.c_str()));
      if (doWT) effWDistX_corr->Multiply(wtRECODistX_corr);
      TH1D* factDistY_corr = (TH1D*)genNumDistY_corr->Clone(Form("factDistY_corr_%i_%i_%s",i,j,shortString.c_str()));
      factDistY_corr->Divide(genDenDistY_corr);
      factDistY_corr->Divide(denDistY_corr);
      TH1D* effCDistY_corr = (TH1D*)factDistY_corr->Clone(Form("effCDistY_corr_%i_%i_%s",i,j,shortString.c_str()));
      if (doCT) effCDistY_corr->Multiply(ctRECODistY_corr);
      TH1D* effWDistY_corr = Invert1Dhist(factDistY_corr,Form("effWDistY_corr_%i_%i_%s",i,j,shortString.c_str()));
      if (doWT) effWDistY_corr->Multiply(wtRECODistY_corr);
      TH1D* factDistZ_corr = (TH1D*)genNumDistZ_corr->Clone(Form("factDistZ_corr_%i_%i_%s",i,j,shortString.c_str()));
      factDistZ_corr->Divide(genDenDistZ_corr);
      factDistZ_corr->Divide(denDistZ_corr);
      TH1D* effCDistZ_corr = (TH1D*)factDistZ_corr->Clone(Form("effCDistZ_corr_%i_%i_%s",i,j,shortString.c_str()));
      if (doCT) effCDistZ_corr->Multiply(ctRECODistZ_corr);
      TH1D* effWDistZ_corr = Invert1Dhist(factDistZ_corr,Form("effWDistZ_corr_%i_%i_%s",i,j,shortString.c_str()));
      if (doWT) effWDistZ_corr->Multiply(wtRECODistZ_corr);
      if (doCT) {
	// set titles
	effCDistX->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effCDistY->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effCDistZ->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
	// set minimum value
	effCDistX->SetMinimum(0);
	effCDistY->SetMinimum(0);
	effCDistZ->SetMinimum(0);
	// set color
	effCDistX_corr->SetLineColor(418);
	effCDistY_corr->SetLineColor(418);
	effCDistZ_corr->SetLineColor(418);
	// fill vector of slices
	effCSliceX.push_back( effCDistX );
	effCSliceY.push_back( effCDistY );
	effCSliceZ.push_back( effCDistZ );
	effCSliceX_corr.push_back( effCDistX_corr );
	effCSliceY_corr.push_back( effCDistY_corr );
	effCSliceZ_corr.push_back( effCDistZ_corr );
      }
      if (doWT) {
	effWDistX->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effWDistY->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effWDistZ->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
	// set minimum value
	effWDistX->SetMinimum(0);
	effWDistY->SetMinimum(0);
	effWDistZ->SetMinimum(0);
	// set color
	effWDistX_corr->SetLineColor(418);
	effWDistY_corr->SetLineColor(418);
	effWDistZ_corr->SetLineColor(418);
	// fill vector of slices
	effWSliceX.push_back( effWDistX );
	effWSliceY.push_back( effWDistY );
	effWSliceZ.push_back( effWDistZ );
	effWSliceX_corr.push_back( effWDistX_corr );
	effWSliceY_corr.push_back( effWDistY_corr );
	effWSliceZ_corr.push_back( effWDistZ_corr );
      }

      // producing 1D slices of efficiency description 
      ctL->setVal(centA);
      phi->setVal(centB*TMath::Pi());
      if (doCT) {
	fsxC.push_back( ctK->frame( Name( Form("fsxC_%i_%i_%s",i,j,shortString.c_str()) ) ) );
	effC->plotOn( fsxC.back(), LineColor(kRed), Name(Form("effCX_%i_%i_%s",i,j,shortString.c_str()))) ;
      }
      if (doWT) {
	fsxW.push_back( ctK->frame( Name( Form("fsxW_%i_%i_%s",i,j,shortString.c_str()) ) ) );
	effW->plotOn( fsxW.back(), LineColor(kRed), Name(Form("effWX_%i_%i_%s",i,j,shortString.c_str()))) ;
      }
      ctK->setVal(centA);
      phi->setVal(centB*TMath::Pi());
      if (doCT) {
	fsyC.push_back( ctL->frame( Name( Form("fsyC_%i_%i_%s",i,j,shortString.c_str()) ) ) );
	effC->plotOn( fsyC.back(), LineColor(kRed), Name(Form("effCY_%i_%i_%s",i,j,shortString.c_str()))) ;
      }
      if (doWT) {
	fsyW.push_back( ctL->frame( Name( Form("fsyW_%i_%i_%s",i,j,shortString.c_str()) ) ) );
	effW->plotOn( fsyW.back(), LineColor(kRed), Name(Form("effWY_%i_%i_%s",i,j,shortString.c_str()))) ;
      }
      ctK->setVal(centA);
      ctL->setVal(centB);
      if (doCT) {
	fszC.push_back( phi->frame( Name( Form("fszC_%i_%i_%s",i,j,shortString.c_str()) ) ) );
	effC->plotOn( fszC.back(), LineColor(kRed), Name(Form("effCZ_%i_%i_%s",i,j,shortString.c_str()))) ;
      }
      if (doWT) {
	fszW.push_back( phi->frame( Name( Form("fszW_%i_%i_%s",i,j,shortString.c_str()) ) ) );
	effW->plotOn( fszW.back(), LineColor(kRed), Name(Form("effWZ_%i_%i_%s",i,j,shortString.c_str()))) ;
      }

      // if (i+j==0) {
      // 	leg->AddEntry(effHistsX.back(),"Binned efficiency" ,"lep");
      // 	leg->AddEntry(xframes.back()->findObject(Form("effx_%i_%i",i,j)),"KDE efficiency","l");
      // }
      // leg->Draw("same");

      if (doCT) {
	// plot ctK slices
	csxC[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effCSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
	effCSliceX.back()->Draw();
	fsxC.back()->Draw("same");
	effCSliceX.back()->Draw("same");
	effCSliceX_corr.back()->Draw("same");
	// plot ctL slices
	csyC[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effCSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
	effCSliceY.back()->Draw();
	fsyC.back()->Draw("same");
	effCSliceY.back()->Draw("same");
	effCSliceY_corr.back()->Draw("same");
	// plot phi slices
	cszC[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effCSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
	effCSliceZ.back()->Draw();
	fszC.back()->Draw("same");
	effCSliceZ.back()->Draw("same");
	effCSliceZ_corr.back()->Draw("same");
      }

      if (doWT) {
	// plot ctK slices
	csxW[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effWSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
	effWSliceX.back()->Draw();
	fsxW.back()->Draw("same");
	effWSliceX.back()->Draw("same");
	effWSliceX_corr.back()->Draw("same");
	// plot ctL slices
	csyW[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effWSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
	effWSliceY.back()->Draw();
	fsyW.back()->Draw("same");
	effWSliceY.back()->Draw("same");
	effWSliceY_corr.back()->Draw("same");
	// plot phi slices
	cszW[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effWSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
	effWSliceZ.back()->Draw();
	fszW.back()->Draw("same");
	effWSliceZ.back()->Draw("same");
	effWSliceZ_corr.back()->Draw("same");
      }

      // checking maximum values
      if (doCT) {
	if ( maxEffCX<effCSliceX.back()->GetMaximum() ) maxEffCX = effCSliceX.back()->GetMaximum();
	if ( maxEffCY<effCSliceY.back()->GetMaximum() ) maxEffCY = effCSliceY.back()->GetMaximum();
	if ( maxEffCZ<effCSliceZ.back()->GetMaximum() ) maxEffCZ = effCSliceZ.back()->GetMaximum();
	if ( maxEffCX<effCSliceX_corr.back()->GetMaximum() ) maxEffCX = effCSliceX_corr.back()->GetMaximum();
	if ( maxEffCY<effCSliceY_corr.back()->GetMaximum() ) maxEffCY = effCSliceY_corr.back()->GetMaximum();
	if ( maxEffCZ<effCSliceZ_corr.back()->GetMaximum() ) maxEffCZ = effCSliceZ_corr.back()->GetMaximum();
      }
      if (doWT) {
	if ( maxEffWX<effWSliceX.back()->GetMaximum() ) maxEffWX = effWSliceX.back()->GetMaximum();
	if ( maxEffWY<effWSliceY.back()->GetMaximum() ) maxEffWY = effWSliceY.back()->GetMaximum();
	if ( maxEffWZ<effWSliceZ.back()->GetMaximum() ) maxEffWZ = effWSliceZ.back()->GetMaximum();
	if ( maxEffWX<effWSliceX_corr.back()->GetMaximum() ) maxEffWX = effWSliceX_corr.back()->GetMaximum();
	if ( maxEffWY<effWSliceY_corr.back()->GetMaximum() ) maxEffWY = effWSliceY_corr.back()->GetMaximum();
	if ( maxEffWZ<effWSliceZ_corr.back()->GetMaximum() ) maxEffWZ = effWSliceZ_corr.back()->GetMaximum();
      }

    }    

  // set uniform y-axis ranges
  if (doCT) {
    for (vector<TH1D*>::iterator hist = effCSliceX.begin(); hist != effCSliceX.end(); ++hist) (*hist)->SetMaximum(maxEffCX*1.1);
    for (vector<TH1D*>::iterator hist = effCSliceY.begin(); hist != effCSliceY.end(); ++hist) (*hist)->SetMaximum(maxEffCY*1.1);
    for (vector<TH1D*>::iterator hist = effCSliceZ.begin(); hist != effCSliceZ.end(); ++hist) (*hist)->SetMaximum(maxEffCZ*1.1);
  }
  if (doWT) {
    for (vector<TH1D*>::iterator hist = effWSliceX.begin(); hist != effWSliceX.end(); ++hist) (*hist)->SetMaximum(maxEffWX*1.1);
    for (vector<TH1D*>::iterator hist = effWSliceY.begin(); hist != effWSliceY.end(); ++hist) (*hist)->SetMaximum(maxEffWY*1.1);
    for (vector<TH1D*>::iterator hist = effWSliceZ.begin(); hist != effWSliceZ.end(); ++hist) (*hist)->SetMaximum(maxEffWZ*1.1);
  }

  if (doCT) {
    csxC[confIndex]->SaveAs( (confString+Form("_eff-ct_CTKslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
    csyC[confIndex]->SaveAs( (confString+Form("_eff-ct_CTLslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
    cszC[confIndex]->SaveAs( (confString+Form("_eff-ct_PHIslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  }
  if (doWT) {
    csxW[confIndex]->SaveAs( (confString+Form("_eff-wt_CTKslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
    csyW[confIndex]->SaveAs( (confString+Form("_eff-wt_CTLslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
    cszW[confIndex]->SaveAs( (confString+Form("_eff-wt_PHIslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  }

  // Plot projections of efficiency functions
  TH1D* effCProj_x = 0;
  TH1D* effCProj_y = 0;
  TH1D* effCProj_z = 0;
  TH2D* effCProj_xy = 0;
  TH2D* effCProj_xz = 0;
  TH2D* effCProj_yz = 0;
  TH1D* effWProj_x = 0;
  TH1D* effWProj_y = 0;
  TH1D* effWProj_z = 0;
  TH2D* effWProj_xy = 0;
  TH2D* effWProj_xz = 0;
  TH2D* effWProj_yz = 0;
  if (doCT) {
    effCProj_x  = effCHist->ProjectionX();
    effCProj_y  = effCHist->ProjectionY();
    effCProj_z  = effCHist->ProjectionZ();
    effCProj_xy = (TH2D*)effCHist->Project3D("xy");
    effCProj_xz = (TH2D*)effCHist->Project3D("xz");
    effCProj_yz = (TH2D*)effCHist->Project3D("yz");
  }
  if (doWT) {
    effWProj_x  = effWHist->ProjectionX();
    effWProj_y  = effWHist->ProjectionY();
    effWProj_z  = effWHist->ProjectionZ();
    effWProj_xy = (TH2D*)effWHist->Project3D("xy");
    effWProj_xz = (TH2D*)effWHist->Project3D("xz");
    effWProj_yz = (TH2D*)effWHist->Project3D("yz");
  }
  // Plot 1D projections
  if (doCT) {
    cp1C[confIndex] = new TCanvas(Form("cp1C_%s",shortString.c_str()),Form("cp1C_%s",shortString.c_str()),1800,800);
    effCProj_x->SetTitle( Form("Correct-tag efficiency projection (q2-bin %i, %s);cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
    effCProj_y->SetTitle( Form("Correct-tag efficiency projection (q2-bin %i, %s);cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
    effCProj_z->SetTitle( Form("Correct-tag efficiency projection (q2-bin %i, %s);#phi",q2Bin,(parity==0?"even":"odd")) );
    effCProj_x->Scale(1.0/2500);
    effCProj_y->Scale(1.0/2500);
    effCProj_z->Scale(1.0/2500);
    effCProj_x->SetMinimum(0.0);
    effCProj_y->SetMinimum(0.0);
    effCProj_z->SetMinimum(0.0);
    cp1C[confIndex]->Divide(3,1);
    cp1C[confIndex]->cd(1);
    effCProj_x->Draw("LHIST");
    cp1C[confIndex]->cd(2);
    effCProj_y->Draw("LHIST");
    cp1C[confIndex]->cd(3);
    effCProj_z->Draw("LHIST");
    cp1C[confIndex]->SaveAs( (confString+"_eff-ct_1DProj.pdf").c_str() );
  }
  if (doWT) {
    cp1W[confIndex] = new TCanvas(Form("cp1W_%s",shortString.c_str()),Form("cp1W_%s",shortString.c_str()),1800,800);
    effWProj_x->SetTitle( Form("Wrong-tag efficiency projection (q2-bin %i, %s);cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
    effWProj_y->SetTitle( Form("Wrong-tag efficiency projection (q2-bin %i, %s);cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
    effWProj_z->SetTitle( Form("Wrong-tag efficiency projection (q2-bin %i, %s);#phi",q2Bin,(parity==0?"even":"odd")) );
    effWProj_x->Scale(1.0/2500);
    effWProj_y->Scale(1.0/2500);
    effWProj_z->Scale(1.0/2500);
    effWProj_x->SetMinimum(0.0);
    effWProj_y->SetMinimum(0.0);
    effWProj_z->SetMinimum(0.0);
    cp1W[confIndex]->Divide(3,1);
    cp1W[confIndex]->cd(1);
    effWProj_x->Draw("LHIST");
    cp1W[confIndex]->cd(2);
    effWProj_y->Draw("LHIST");
    cp1W[confIndex]->cd(3);
    effWProj_z->Draw("LHIST");
    cp1W[confIndex]->SaveAs( (confString+"_eff-wt_1DProj.pdf").c_str() );
  }
  // 2D projections
  if (doCT) {
    cp2C[confIndex] = new TCanvas(Form("cp2C_%s",shortString.c_str()),Form("cp2C_%s",shortString.c_str()),1800,800);
    effCProj_xy->SetTitle( Form("Correct-tag efficiency projection (q2-bin %i, %s);cos(#theta_{L});cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
    effCProj_xz->SetTitle( Form("Correct-tag efficiency projection (q2-bin %i, %s);#phi;cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
    effCProj_yz->SetTitle( Form("Correct-tag efficiency projection (q2-bin %i, %s);#phi;cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
    effCProj_xy->Scale(1.0/50);
    effCProj_xz->Scale(1.0/50);
    effCProj_yz->Scale(1.0/50);
    effCProj_xy->GetXaxis()->SetTitleOffset(1.4);
    effCProj_xz->GetXaxis()->SetTitleOffset(1.4);
    effCProj_yz->GetXaxis()->SetTitleOffset(1.4);
    effCProj_xy->GetYaxis()->SetTitleOffset(2);
    effCProj_xz->GetYaxis()->SetTitleOffset(2);
    effCProj_yz->GetYaxis()->SetTitleOffset(2);
    effCProj_xy->SetMinimum(0.0);
    effCProj_xz->SetMinimum(0.0);
    effCProj_yz->SetMinimum(0.0);
    cp2C[confIndex]->Divide(3,1);
    cp2C[confIndex]->cd(1);
    effCProj_xy->Draw("SURF3");
    cp2C[confIndex]->cd(2);
    effCProj_xz->Draw("SURF3");
    cp2C[confIndex]->cd(3);
    effCProj_yz->Draw("SURF3");
    cp2C[confIndex]->SaveAs( (confString+"_eff-ct_2DProj.pdf").c_str() );
  }
  if (doWT) {
    cp2W[confIndex] = new TCanvas(Form("cp2W_%s",shortString.c_str()),Form("cp2W_%s",shortString.c_str()),1800,800);
    effWProj_xy->SetTitle( Form("Wrong-tag efficiency projection (q2-bin %i, %s);cos(#theta_{L});cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
    effWProj_xz->SetTitle( Form("Wrong-tag efficiency projection (q2-bin %i, %s);#phi;cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
    effWProj_yz->SetTitle( Form("Wrong-tag efficiency projection (q2-bin %i, %s);#phi;cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
    effWProj_xy->Scale(1.0/50);
    effWProj_xz->Scale(1.0/50);
    effWProj_yz->Scale(1.0/50);
    effWProj_xy->GetXaxis()->SetTitleOffset(1.4);
    effWProj_xz->GetXaxis()->SetTitleOffset(1.4);
    effWProj_yz->GetXaxis()->SetTitleOffset(1.4);
    effWProj_xy->GetYaxis()->SetTitleOffset(2);
    effWProj_xz->GetYaxis()->SetTitleOffset(2);
    effWProj_yz->GetYaxis()->SetTitleOffset(2);
    effWProj_xy->SetMinimum(0.0);
    effWProj_xz->SetMinimum(0.0);
    effWProj_yz->SetMinimum(0.0);
    cp2W[confIndex]->Divide(3,1);
    cp2W[confIndex]->cd(1);
    effWProj_xy->Draw("SURF3");
    cp2W[confIndex]->cd(2);
    effWProj_xz->Draw("SURF3");
    cp2W[confIndex]->cd(3);
    effWProj_yz->Draw("SURF3");
    cp2W[confIndex]->SaveAs( (confString+"_eff-wt_2DProj.pdf").c_str() );
  }

  if (!doClosure) return;

  int nbins_ct = 30;

  // create and fill GEN histograms
  TH1D* hXctGEN_effC = 0;
  TH1D* hXctGEN_effW = 0;
  TH1D* hYctGEN_effC = 0;
  TH1D* hYctGEN_effW = 0;
  TH1D* hZctGEN_effC = 0;
  TH1D* hZctGEN_effW = 0;
  if (doCT) {
    hXctGEN_effC = new TH1D( Form("hXctGEN_effC_%s",shortString.c_str()), Form("hXctGEN_effC_%s",shortString.c_str()), nbins_ct, -1, 1 );
    hYctGEN_effC = new TH1D( Form("hYctGEN_effC_%s",shortString.c_str()), Form("hYctGEN_effC_%s",shortString.c_str()), nbins_ct, -1, 1 );
    hZctGEN_effC = new TH1D( Form("hZctGEN_effC_%s",shortString.c_str()), Form("hZctGEN_effC_%s",shortString.c_str()), nbins_ct, -TMath::Pi(), TMath::Pi() );
    hXctGEN_effC->Sumw2();
    hYctGEN_effC->Sumw2();
    hZctGEN_effC->Sumw2();
  }
  if (doWT) {
    hXctGEN_effW = new TH1D( Form("hXctGEN_effW_%s",shortString.c_str()), Form("hXctGEN_effW_%s",shortString.c_str()), nbins_ct, -1, 1 );
    hYctGEN_effW = new TH1D( Form("hYctGEN_effW_%s",shortString.c_str()), Form("hYctGEN_effW_%s",shortString.c_str()), nbins_ct, -1, 1 );
    hZctGEN_effW = new TH1D( Form("hZctGEN_effW_%s",shortString.c_str()), Form("hZctGEN_effW_%s",shortString.c_str()), nbins_ct, -TMath::Pi(), TMath::Pi() );
    hXctGEN_effW->Sumw2();
    hYctGEN_effW->Sumw2();
    hZctGEN_effW->Sumw2();
  }
  double gen_ctK, gen_ctL, gen_phi, wei;
  for (int iEv=0; iEv<data_genDen->numEntries(); ++iEv) {
    const RooArgSet* set = data_genDen->get(iEv);
    gen_ctK = ((RooRealVar*)set->find("ctK"))->getVal();
    gen_ctL = ((RooRealVar*)set->find("ctL"))->getVal();
    gen_phi = ((RooRealVar*)set->find("phi"))->getVal();
    if (doCT) {
      ctK->setVal(gen_ctK);
      ctL->setVal(gen_ctL);
      phi->setVal(gen_phi);
      wei = effC->getVal();
      hXctGEN_effC->Fill( gen_ctK, wei );
      hYctGEN_effC->Fill( gen_ctL, wei );
      hZctGEN_effC->Fill( gen_phi, wei );
    }
    if (doWT) {
      ctK->setVal(-1*gen_ctK);
      ctL->setVal(-1*gen_ctL);
      phi->setVal(-1*gen_phi);
      wei = effW->getVal();
      hXctGEN_effW->Fill( -1*gen_ctK, wei );
      hYctGEN_effW->Fill( -1*gen_ctL, wei );
      hZctGEN_effW->Fill( -1*gen_phi, wei );
    }
  }

  // create RECO histograms
  TH1D* hXctREC_effC = 0;
  TH1D* hXctREC_effW = 0;
  TH1D* hYctREC_effC = 0;
  TH1D* hYctREC_effW = 0;
  TH1D* hZctREC_effC = 0;
  TH1D* hZctREC_effW = 0;
  if (doCT) {
    hXctREC_effC = (TH1D*)data_ctRECO->createHistogram( Form("hXctREC_effC_%s",shortString.c_str()), *ctK, Binning(nbins_ct,-1,1) );
    hYctREC_effC = (TH1D*)data_ctRECO->createHistogram( Form("hYctREC_effC_%s",shortString.c_str()), *ctL, Binning(nbins_ct,-1,1) );
    hZctREC_effC = (TH1D*)data_ctRECO->createHistogram( Form("hZctREC_effC_%s",shortString.c_str()), *phi, Binning(nbins_ct,-TMath::Pi(),TMath::Pi()) );
  }
  if (doWT) {
    hXctREC_effW = (TH1D*)data_wtRECO->createHistogram( Form("hXctREC_effW_%s",shortString.c_str()), *ctK, Binning(nbins_ct,-1,1) );
    hYctREC_effW = (TH1D*)data_wtRECO->createHistogram( Form("hYctREC_effW_%s",shortString.c_str()), *ctL, Binning(nbins_ct,-1,1) );
    hZctREC_effW = (TH1D*)data_wtRECO->createHistogram( Form("hZctREC_effW_%s",shortString.c_str()), *phi, Binning(nbins_ct,-TMath::Pi(),TMath::Pi()) );
  }

  // prepare for plot
  if (doCT) {
    hXctGEN_effC->SetTitle( Form("Closure test of correct-tag efficiency (q2-bin %i, %s);cos(#theta_{K});Events",q2Bin,(parity==0?"even":"odd")) );
    hYctGEN_effC->SetTitle( Form("Closure test of correct-tag efficiency (q2-bin %i, %s);cos(#theta_{L});Events",q2Bin,(parity==0?"even":"odd")) );
    hZctGEN_effC->SetTitle( Form("Closure test of correct-tag efficiency (q2-bin %i, %s);#phi;Events",q2Bin,(parity==0?"even":"odd")) );
    hXctGEN_effC->Scale( hXctREC_effC->Integral() / hXctGEN_effC->Integral() );
    hYctGEN_effC->Scale( hYctREC_effC->Integral() / hYctGEN_effC->Integral() );
    hZctGEN_effC->Scale( hZctREC_effC->Integral() / hZctGEN_effC->Integral() );
    hXctGEN_effC->SetMinimum(0);
    hYctGEN_effC->SetMinimum(0);
    hZctGEN_effC->SetMinimum(0);
    hXctGEN_effC->SetLineColor(kRed+1);
    hYctGEN_effC->SetLineColor(kRed+1);
    hZctGEN_effC->SetLineColor(kRed+1);
    hXctGEN_effC->SetLineWidth(2);
    hYctGEN_effC->SetLineWidth(2);
    hZctGEN_effC->SetLineWidth(2);
    hXctREC_effC->SetLineWidth(2);
    hYctREC_effC->SetLineWidth(2);
    hZctREC_effC->SetLineWidth(2);
  }
  if (doWT) {
    hXctGEN_effW->SetTitle( Form(  "Closure test of wrong-tag efficiency (q2-bin %i, %s);cos(#theta_{K});Events",q2Bin,(parity==0?"even":"odd")) );
    hYctGEN_effW->SetTitle( Form(  "Closure test of wrong-tag efficiency (q2-bin %i, %s);cos(#theta_{L});Events",q2Bin,(parity==0?"even":"odd")) );
    hZctGEN_effW->SetTitle( Form(  "Closure test of wrong-tag efficiency (q2-bin %i, %s);#phi;Events",q2Bin,(parity==0?"even":"odd")) );
    hXctGEN_effW->Scale( hXctREC_effW->Integral() / hXctGEN_effW->Integral() );
    hYctGEN_effW->Scale( hYctREC_effW->Integral() / hYctGEN_effW->Integral() );
    hZctGEN_effW->Scale( hZctREC_effW->Integral() / hZctGEN_effW->Integral() );
    hXctGEN_effW->SetMinimum(0);
    hYctGEN_effW->SetMinimum(0);
    hZctGEN_effW->SetMinimum(0);
    hXctGEN_effW->SetLineColor(kRed+1);
    hYctGEN_effW->SetLineColor(kRed+1);
    hZctGEN_effW->SetLineColor(kRed+1);
    hXctGEN_effW->SetLineWidth(2);
    hYctGEN_effW->SetLineWidth(2);
    hZctGEN_effW->SetLineWidth(2);
    hXctREC_effW->SetLineWidth(2);
    hYctREC_effW->SetLineWidth(2);
    hZctREC_effW->SetLineWidth(2);
  }
  // Plot closure test
  gPad->SetLeftMargin(0.17); 
  if (doCT) {
    cctC[confIndex] = new TCanvas(("ccteffC_"+shortString).c_str(),("ccteffC_"+shortString).c_str(),2000,700);
    cctC[confIndex]->Divide(3,1);
    cctC[confIndex]->cd(1);
    hXctGEN_effC->Draw();
    hXctREC_effC->Draw("same");
    cctC[confIndex]->cd(2);
    hYctGEN_effC->Draw();
    hYctREC_effC->Draw("same");
    cctC[confIndex]->cd(3);
    hZctGEN_effC->Draw();
    hZctREC_effC->Draw("same");
    cctC[confIndex]->SaveAs( (confString+"_eff-ct_ClosureTest.pdf").c_str() );
  }
  if (doWT) {
    cctW[confIndex] = new TCanvas(("ccteffW_"+shortString).c_str(),("ccteffW_"+shortString).c_str(),2000,700);
    // TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
    // xframe->GetYaxis()->SetTitleOffset(1.6);
    // yframe->GetYaxis()->SetTitleOffset(1.6);
    // zframe->GetYaxis()->SetTitleOffset(1.6);
    // xframe->SetMaximum(xframe->GetMaximum()*1.15);
    // yframe->SetMaximum(yframe->GetMaximum()*1.15);
    // zframe->SetMaximum(zframe->GetMaximum()*1.15);
    // leg->SetTextSize(0.03);
    // leg->AddEntry(xframe->findObject("plNumDist"),"Post-selection RECO distribution" ,"lep");
    // leg->AddEntry(xframe->findObject("plDenDist"),"Efficiency-corrected GEN distribution" ,"lep");
    // leg->Draw("same");
    cctW[confIndex]->Divide(3,1);
    cctW[confIndex]->cd(1);
    hXctGEN_effW->Draw();
    hXctREC_effW->Draw("same");
    cctW[confIndex]->cd(2);
    hYctGEN_effW->Draw();
    hYctREC_effW->Draw("same");
    cctW[confIndex]->cd(3);
    hZctGEN_effW->Draw();
    hZctREC_effW->Draw("same");
    cctW[confIndex]->SaveAs( (confString+"_eff-wt_ClosureTest.pdf").c_str() );
  }

}

TH1D* Invert1Dhist(TH1D* hin, string hname)
{
  // clone and reset the TH1
  TH1D* hout=(TH1D*)hin->Clone(hname.c_str());
  hout->Reset();

  // fill bin content and error one by one
  for (int ix=0; ix!=hin->GetXaxis()->GetNbins()+1;++ix) {
    int ix_inv = hin->GetXaxis()->GetNbins()+1-ix;
    hout->SetBinContent( ix_inv, hin->GetBinContent(ix) );
    hout->SetBinError  ( ix_inv, hin->GetBinError  (ix) );
  }

  return hout;

}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency vs. odd data
  //                [1] odd efficiency vs. even data
  //                [-1] for each parity recursively

  int q2Bin  = -1;
  int parity = -1; 

  if ( argc >= 2 ) q2Bin  = atoi(argv[1]);
  if ( argc >= 3 ) parity = atoi(argv[2]);

  if ( q2Bin  < -1 || q2Bin  >= nBins ) return 1;
  if ( parity < -1 || parity > 1      ) return 1;

  bool doClosure = true;
  if ( argc >= 4 && atoi(argv[3]) == 0 ) doClosure = false;

  if ( q2Bin > -1 ) {
    if ( parity > -1 ) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<(parity==1?" - odd":" - even")<<" events"<<endl;
      plotEffBin( q2Bin, parity, doClosure );
    } else {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<" - both event parities"<<endl;
      plotEffBin( q2Bin, 0, doClosure );
      plotEffBin( q2Bin, 1, doClosure );
    }
  } else {
    cout<<"Plotting efficiency for all q2 bins - "<<(parity==1?"odd events":(parity==0?"even events":"both event parities"))<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if ( parity > -1 ) plotEffBin( q2Bin, parity, doClosure );
      else {
	plotEffBin( q2Bin, 0, doClosure );
	plotEffBin( q2Bin, 1, doClosure );
      }
    }
  }

  return 0;

}
