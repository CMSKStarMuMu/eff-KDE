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

TCanvas* csx [2*nBins];
TCanvas* csy [2*nBins];
TCanvas* csz [2*nBins];
TCanvas* cp1 [2*nBins];
TCanvas* cp2 [2*nBins];
TCanvas* cct [2*nBins];
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
TCanvas* csxcf [2*nBins];
TCanvas* csycf [2*nBins];
TCanvas* cszcf [2*nBins];
TCanvas* cp1cf [2*nBins];
TCanvas* cp2cf [2*nBins];
TCanvas* cctcf [2*nBins];
TCanvas* csxmf [2*nBins];
TCanvas* csymf [2*nBins];
TCanvas* cszmf [2*nBins];
TCanvas* cp1mf [2*nBins];
TCanvas* cp2mf [2*nBins];
TCanvas* cctmf [2*nBins];

TH1D* Invert1Dhist(TH1D* hin, string hname);

void plotEffBin(int q2Bin, int parity, bool doClosure)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  int confIndex = nBins*parity + q2Bin;

  int distBins = 10;

  // Load variables and dataset
  string filename_data = Form("/eos/user/a/aboletti/BdToKstarMuMu/datasets/PUweight/effDataset_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  // import the complementary dataset, to compare statistically uncorrelated values
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,1-parity));
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
  string datasetString = Form((parity==1?"_ev_b%i":"_od_b%i"),q2Bin);
  RooDataSet* data_genDen = (RooDataSet*)wsp->data(("data_genDen"+datasetString).c_str());
  RooDataSet* data_genNum = (RooDataSet*)wsp->data(("data_genNum"+datasetString).c_str());
  RooDataSet* data_den    = (RooDataSet*)wsp->data(("data_den"   +datasetString).c_str());
  RooDataSet* data_ctRECO = (RooDataSet*)wsp->data(("data_ctRECO"+datasetString).c_str());
  RooDataSet* data_wtRECO = (RooDataSet*)wsp->data(("data_wtRECO"+datasetString).c_str());
  
  // import KDE efficiency histograms
  string filename = Form((parity==0?"files/KDEeff_b%i_ev.root":"files/KDEeff_b%i_od.root"),q2Bin);
  TFile* fin = new TFile( filename.c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  TH3D* effHist = (TH3D*)fin->Get(("effHist_"+shortString).c_str());
  TH3D* effCHist = (TH3D*)fin->Get(("effCHist_"+shortString).c_str());
  TH3D* effWHist = (TH3D*)fin->Get(("effWHist_"+shortString).c_str());
  TH3D* corrFracHist = (TH3D*)fin->Get(("corrFracHist_"+shortString).c_str());
  TH3D* mistFracHist = (TH3D*)fin->Get(("mistFracHist_"+shortString).c_str());
  if ( !effHist || !effCHist || !effWHist || !corrFracHist || !mistFracHist ||
       effHist->IsZombie() || effCHist->IsZombie() || effWHist->IsZombie() || corrFracHist->IsZombie() || mistFracHist->IsZombie() ) {
    cout<<"Efficiency histograms not found in file: "<<filename<<endl;
    return;
  }
  // create efficiency functions
  RooDataHist* effData = new RooDataHist("effData","effData",vars,effHist);
  RooDataHist* effCData = new RooDataHist("effCData","effCData",vars,effCHist);
  RooDataHist* effWData = new RooDataHist("effWData","effWData",vars,effWHist);
  RooDataHist* corrFracData = new RooDataHist("corrFracData","corrFracData",vars,corrFracHist);
  RooDataHist* mistFracData = new RooDataHist("mistFracData","mistFracData",vars,mistFracHist);
  RooAbsReal* eff = new RooHistFunc("eff","eff",vars,*effData,1);
  RooAbsReal* effC = new RooHistFunc("effC","effC",vars,*effCData,1);
  RooAbsReal* effW = new RooHistFunc("effW","effW",vars,*effWData,1);
  RooAbsReal* corrFrac = new RooHistFunc("corrFrac","corrFrac",vars,*corrFracData,1);
  RooAbsReal* mistFrac = new RooHistFunc("mistFrac","mistFrac",vars,*mistFracData,1);

  string confString = "plotEff_d/effKDE_test_"+shortString;

  gStyle->SetOptStat(0);

  // Plot 1D slices of the efficiency function and binned efficiency
  csx[confIndex] = new TCanvas(("csx"+shortString).c_str(),(shortString+"_eff_ctK").c_str(),1500,1500) ;
  csy[confIndex] = new TCanvas(("csy"+shortString).c_str(),(shortString+"_eff_ctL").c_str(),1500,1500) ;
  csz[confIndex] = new TCanvas(("csz"+shortString).c_str(),(shortString+"_eff_phi").c_str(),1500,1500) ;
  csxC[confIndex] = new TCanvas(("csxC"+shortString).c_str(),(shortString+"_effC_ctK").c_str(),1500,1500) ;
  csyC[confIndex] = new TCanvas(("csyC"+shortString).c_str(),(shortString+"_effC_ctL").c_str(),1500,1500) ;
  cszC[confIndex] = new TCanvas(("cszC"+shortString).c_str(),(shortString+"_effC_phi").c_str(),1500,1500) ;
  csxW[confIndex] = new TCanvas(("csxW"+shortString).c_str(),(shortString+"_effW_ctK").c_str(),1500,1500) ;
  csyW[confIndex] = new TCanvas(("csyW"+shortString).c_str(),(shortString+"_effW_ctL").c_str(),1500,1500) ;
  cszW[confIndex] = new TCanvas(("cszW"+shortString).c_str(),(shortString+"_effW_phi").c_str(),1500,1500) ;
  csxcf[confIndex] = new TCanvas(("csxcf"+shortString).c_str(),(shortString+"_cff_ctK").c_str(),1500,1500) ;
  csycf[confIndex] = new TCanvas(("csycf"+shortString).c_str(),(shortString+"_cf_ctL").c_str(),1500,1500) ;
  cszcf[confIndex] = new TCanvas(("cszcf"+shortString).c_str(),(shortString+"_cf_phi").c_str(),1500,1500) ;
  csxmf[confIndex] = new TCanvas(("csxmf"+shortString).c_str(),(shortString+"_mf_ctK").c_str(),1500,1500) ;
  csymf[confIndex] = new TCanvas(("csymf"+shortString).c_str(),(shortString+"_mf_ctL").c_str(),1500,1500) ;
  cszmf[confIndex] = new TCanvas(("cszmf"+shortString).c_str(),(shortString+"_mf_phi").c_str(),1500,1500) ;
  csx[confIndex]->Divide(5,5);
  csy[confIndex]->Divide(5,5);
  csz[confIndex]->Divide(5,5);
  csxC[confIndex]->Divide(5,5);
  csyC[confIndex]->Divide(5,5);
  cszC[confIndex]->Divide(5,5);
  csxW[confIndex]->Divide(5,5);
  csyW[confIndex]->Divide(5,5);
  cszW[confIndex]->Divide(5,5);
  csxcf[confIndex]->Divide(5,5);
  csycf[confIndex]->Divide(5,5);
  cszcf[confIndex]->Divide(5,5);
  csxmf[confIndex]->Divide(5,5);
  csymf[confIndex]->Divide(5,5);
  cszmf[confIndex]->Divide(5,5);
  vector <TH1D*> effSliceX;
  vector <TH1D*> effSliceY;
  vector <TH1D*> effSliceZ;
  vector <TH1D*> effCSliceX;
  vector <TH1D*> effCSliceY;
  vector <TH1D*> effCSliceZ;
  vector <TH1D*> effWSliceX;
  vector <TH1D*> effWSliceY;
  vector <TH1D*> effWSliceZ;
  vector <TH1D*> mfSliceX;
  vector <TH1D*> mfSliceY;
  vector <TH1D*> mfSliceZ;
  vector <TH1D*> cfSliceX;
  vector <TH1D*> cfSliceY;
  vector <TH1D*> cfSliceZ;
  vector <RooPlot*> fsx;
  vector <RooPlot*> fsxC;
  vector <RooPlot*> fsxW;
  vector <RooPlot*> fsxcf;
  vector <RooPlot*> fsxmf;
  vector <RooPlot*> fsy;
  vector <RooPlot*> fsyC;
  vector <RooPlot*> fsyW;
  vector <RooPlot*> fsycf;
  vector <RooPlot*> fsymf;
  vector <RooPlot*> fsz;
  vector <RooPlot*> fszC;
  vector <RooPlot*> fszW;
  vector <RooPlot*> fszcf;
  vector <RooPlot*> fszmf;

  // TLegend* leg = new TLegend (0.35,0.8,0.9,0.9);

  // width of the slices in the hidden variables ("border" is half of it)
  double border = 0.04;

  // variables to be filled with global efficiency maximum
  double maxEffX = 0;
  double maxEffCX = 0;
  double maxEffWX = 0;
  double maxCFX = 0;
  double maxMFX = 0;
  double maxEffY = 0;
  double maxEffCY = 0;
  double maxEffWY = 0;
  double maxCFY = 0;
  double maxMFY = 0;
  double maxEffZ = 0;
  double maxEffCZ = 0;
  double maxEffWZ = 0;
  double maxCFZ = 0;
  double maxMFZ = 0;

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
      TH1D* ctRECODistX = (TH1D*)data_ctRECO->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_ctRECO_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* wtRECODistX = (TH1D*)data_wtRECO->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distX_wtRECO_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(distBins,-1,1) );
      TH1D* genDenDistY = (TH1D*)data_genDen->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_genDen_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* genNumDistY = (TH1D*)data_genNum->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_genNum_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* denDistY    = (TH1D*)data_den   ->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_den_%i_%i_%s"   ,i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* ctRECODistY = (TH1D*)data_ctRECO->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_ctRECO_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* wtRECODistY = (TH1D*)data_wtRECO->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distY_wtRECO_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(distBins,-1,1) );
      TH1D* genDenDistZ = (TH1D*)data_genDen->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_genDen_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* genNumDistZ = (TH1D*)data_genNum->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_genNum_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* denDistZ    = (TH1D*)data_den   ->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_den_%i_%i_%s"   ,i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* ctRECODistZ = (TH1D*)data_ctRECO->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_ctRECO_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* wtRECODistZ = (TH1D*)data_wtRECO->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distZ_wtRECO_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(distBins,-TMath::Pi(),TMath::Pi()) );
      TH1D* RECODistX = (TH1D*)ctRECODistX->Clone(Form("distX_RECO_%i_%i_%s",i,j,shortString.c_str()));
      RECODistX->Add(wtRECODistX);
      TH1D* RECODistY = (TH1D*)ctRECODistY->Clone(Form("distY_RECO_%i_%i_%s",i,j,shortString.c_str()));
      RECODistY->Add(wtRECODistY);
      TH1D* RECODistZ = (TH1D*)ctRECODistZ->Clone(Form("distZ_RECO_%i_%i_%s",i,j,shortString.c_str()));
      RECODistZ->Add(wtRECODistZ);
      // composing binned efficiencies from sliced distributions
      TH1D* factDistX = (TH1D*)genNumDistX->Clone(Form("factDistX_%i_%i_%s",i,j,shortString.c_str()));
      factDistX->Divide(genDenDistX);
      factDistX->Divide(denDistX);
      TH1D* effCDistX = (TH1D*)factDistX->Clone(Form("effCDistX_%i_%i_%s",i,j,shortString.c_str()));
      effCDistX->Multiply(ctRECODistX);
      TH1D* effWDistX = Invert1Dhist(factDistX,Form("effWDistX_%i_%i_%s",i,j,shortString.c_str()));
      effWDistX->Multiply(wtRECODistX);
      TH1D* effDistX = (TH1D*)effCDistX->Clone(Form("effDistX_%i_%i_%s",i,j,shortString.c_str()));
      effDistX->Add(effWDistX);
      TH1D* corrFracDistX = (TH1D*)ctRECODistX->Clone(Form("corrFracDistX_%i_%i_%s",i,j,shortString.c_str()));
      corrFracDistX->Divide(RECODistX);
      TH1D* mistFracDistX = (TH1D*)wtRECODistX->Clone(Form("mistFracDistX_%i_%i_%s",i,j,shortString.c_str()));
      mistFracDistX->Divide(RECODistX);
      TH1D* factDistY = (TH1D*)genNumDistY->Clone(Form("factDistY_%i_%i_%s",i,j,shortString.c_str()));
      factDistY->Divide(genDenDistY);
      factDistY->Divide(denDistY);
      TH1D* effCDistY = (TH1D*)factDistY->Clone(Form("effCDistY_%i_%i_%s",i,j,shortString.c_str()));
      effCDistY->Multiply(ctRECODistY);
      TH1D* effWDistY = Invert1Dhist(factDistY,Form("effWDistY_%i_%i_%s",i,j,shortString.c_str()));
      effWDistY->Multiply(wtRECODistY);
      TH1D* effDistY = (TH1D*)effCDistY->Clone(Form("effDistY_%i_%i_%s",i,j,shortString.c_str()));
      effDistY->Add(effWDistY);
      TH1D* corrFracDistY = (TH1D*)ctRECODistY->Clone(Form("corrFracDistY_%i_%i_%s",i,j,shortString.c_str()));
      corrFracDistY->Divide(RECODistY);
      TH1D* mistFracDistY = (TH1D*)wtRECODistY->Clone(Form("mistFracDistY_%i_%i_%s",i,j,shortString.c_str()));
      mistFracDistY->Divide(RECODistY);
      TH1D* factDistZ = (TH1D*)genNumDistZ->Clone(Form("factDistZ_%i_%i_%s",i,j,shortString.c_str()));
      factDistZ->Divide(genDenDistZ);
      factDistZ->Divide(denDistZ);
      TH1D* effCDistZ = (TH1D*)factDistZ->Clone(Form("effCDistZ_%i_%i_%s",i,j,shortString.c_str()));
      effCDistZ->Multiply(ctRECODistZ);
      TH1D* effWDistZ = Invert1Dhist(factDistZ,Form("effWDistZ_%i_%i_%s",i,j,shortString.c_str()));
      effWDistZ->Multiply(wtRECODistZ);
      TH1D* effDistZ = (TH1D*)effCDistZ->Clone(Form("effDistZ_%i_%i_%s",i,j,shortString.c_str()));
      effDistZ->Add(effWDistZ);
      TH1D* corrFracDistZ = (TH1D*)ctRECODistZ->Clone(Form("corrFracDistZ_%i_%i_%s",i,j,shortString.c_str()));
      corrFracDistZ->Divide(RECODistZ);
      TH1D* mistFracDistZ = (TH1D*)wtRECODistZ->Clone(Form("mistFracDistZ_%i_%i_%s",i,j,shortString.c_str()));
      mistFracDistZ->Divide(RECODistZ);
      // set titles
      effDistX->SetTitle( Form("Full efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      effDistY->SetTitle( Form("Full efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      effDistZ->SetTitle( Form("Full efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      effCDistX->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      effCDistY->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      effCDistZ->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      effWDistX->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      effWDistY->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      effWDistZ->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      corrFracDistX->SetTitle( Form("Correct-tag fraction (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Correct-tag fraction",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      corrFracDistY->SetTitle( Form("Correct-tag fraction (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Correct-tag fraction",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      corrFracDistZ->SetTitle( Form("Correct-tag fraction (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Correct-tag fraction",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      mistFracDistX->SetTitle( Form("Mistag fraction (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Mistag fraction",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      mistFracDistY->SetTitle( Form("Mistag fraction (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Mistag fraction",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
      mistFracDistZ->SetTitle( Form("Mistag fraction (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Mistag fraction",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      // set minimum value
      effDistX->SetMinimum(0);
      effCDistX->SetMinimum(0);
      effWDistX->SetMinimum(0);
      corrFracDistX->SetMinimum(0);
      mistFracDistX->SetMinimum(0);
      effDistY->SetMinimum(0);
      effCDistY->SetMinimum(0);
      effWDistY->SetMinimum(0);
      corrFracDistY->SetMinimum(0);
      mistFracDistY->SetMinimum(0);
      effDistZ->SetMinimum(0);
      effCDistZ->SetMinimum(0);
      effWDistZ->SetMinimum(0);
      corrFracDistZ->SetMinimum(0);
      mistFracDistZ->SetMinimum(0);
      // fill vector of slices
      effSliceX.push_back( effDistX );
      effCSliceX.push_back( effCDistX );
      effWSliceX.push_back( effWDistX );
      cfSliceX.push_back( corrFracDistX );
      mfSliceX.push_back( mistFracDistX );
      effSliceY.push_back( effDistY );
      effCSliceY.push_back( effCDistY );
      effWSliceY.push_back( effWDistY );
      cfSliceY.push_back( corrFracDistY );
      mfSliceY.push_back( mistFracDistY );
      effSliceZ.push_back( effDistZ );
      effCSliceZ.push_back( effCDistZ );
      effWSliceZ.push_back( effWDistZ );
      cfSliceZ.push_back( corrFracDistZ );
      mfSliceZ.push_back( mistFracDistZ );

      // producing 1D slices of efficiency description 
      ctL->setVal(centA);
      phi->setVal(centB*TMath::Pi());
      fsx  .push_back( ctK->frame( Name( Form("fsx_%i_%i_%s"  ,i,j,shortString.c_str()) ) ) );
      fsxC .push_back( ctK->frame( Name( Form("fsxC_%i_%i_%s" ,i,j,shortString.c_str()) ) ) );
      fsxW .push_back( ctK->frame( Name( Form("fsxW_%i_%i_%s" ,i,j,shortString.c_str()) ) ) );
      fsxcf.push_back( ctK->frame( Name( Form("fsxcf_%i_%i_%s",i,j,shortString.c_str()) ) ) );
      fsxmf.push_back( ctK->frame( Name( Form("fsxmf_%i_%i_%s",i,j,shortString.c_str()) ) ) );
      eff     ->plotOn( fsx  .back(), LineColor(kRed), Name(Form("effX_%i_%i_%s"     ,i,j,shortString.c_str()))) ;
      effC    ->plotOn( fsxC .back(), LineColor(kRed), Name(Form("effCX_%i_%i_%s"    ,i,j,shortString.c_str()))) ;
      effW    ->plotOn( fsxW .back(), LineColor(kRed), Name(Form("effWX_%i_%i_%s"    ,i,j,shortString.c_str()))) ;
      corrFrac->plotOn( fsxcf.back(), LineColor(kRed), Name(Form("corrFracX_%i_%i_%s",i,j,shortString.c_str()))) ;
      mistFrac->plotOn( fsxmf.back(), LineColor(kRed), Name(Form("mistFracX_%i_%i_%s",i,j,shortString.c_str()))) ;
      ctK->setVal(centA);
      phi->setVal(centB*TMath::Pi());
      fsy  .push_back( ctL->frame( Name( Form("fsy_%i_%i_%s"  ,i,j,shortString.c_str()) ) ) );
      fsyC .push_back( ctL->frame( Name( Form("fsyC_%i_%i_%s" ,i,j,shortString.c_str()) ) ) );
      fsyW .push_back( ctL->frame( Name( Form("fsyW_%i_%i_%s" ,i,j,shortString.c_str()) ) ) );
      fsycf.push_back( ctL->frame( Name( Form("fsycf_%i_%i_%s",i,j,shortString.c_str()) ) ) );
      fsymf.push_back( ctL->frame( Name( Form("fsymf_%i_%i_%s",i,j,shortString.c_str()) ) ) );
      eff     ->plotOn( fsy  .back(), LineColor(kRed), Name(Form("effY_%i_%i_%s"     ,i,j,shortString.c_str()))) ;
      effC    ->plotOn( fsyC .back(), LineColor(kRed), Name(Form("effCY_%i_%i_%s"    ,i,j,shortString.c_str()))) ;
      effW    ->plotOn( fsyW .back(), LineColor(kRed), Name(Form("effWY_%i_%i_%s"    ,i,j,shortString.c_str()))) ;
      corrFrac->plotOn( fsycf.back(), LineColor(kRed), Name(Form("corrFracY_%i_%i_%s",i,j,shortString.c_str()))) ;
      mistFrac->plotOn( fsymf.back(), LineColor(kRed), Name(Form("mistFracY_%i_%i_%s",i,j,shortString.c_str()))) ;
      ctK->setVal(centA);
      ctL->setVal(centB);
      fsz  .push_back( phi->frame( Name( Form("fsz_%i_%i_%s"  ,i,j,shortString.c_str()) ) ) );
      fszC .push_back( phi->frame( Name( Form("fszC_%i_%i_%s" ,i,j,shortString.c_str()) ) ) );
      fszW .push_back( phi->frame( Name( Form("fszW_%i_%i_%s" ,i,j,shortString.c_str()) ) ) );
      fszcf.push_back( phi->frame( Name( Form("fszcf_%i_%i_%s",i,j,shortString.c_str()) ) ) );
      fszmf.push_back( phi->frame( Name( Form("fszmf_%i_%i_%s",i,j,shortString.c_str()) ) ) );
      eff     ->plotOn( fsz  .back(), LineColor(kRed), Name(Form("effZ_%i_%i_%s"     ,i,j,shortString.c_str()))) ;
      effC    ->plotOn( fszC .back(), LineColor(kRed), Name(Form("effCZ_%i_%i_%s"    ,i,j,shortString.c_str()))) ;
      effW    ->plotOn( fszW .back(), LineColor(kRed), Name(Form("effWZ_%i_%i_%s"    ,i,j,shortString.c_str()))) ;
      corrFrac->plotOn( fszcf.back(), LineColor(kRed), Name(Form("corrFracZ_%i_%i_%s",i,j,shortString.c_str()))) ;
      mistFrac->plotOn( fszmf.back(), LineColor(kRed), Name(Form("mistFracZ_%i_%i_%s",i,j,shortString.c_str()))) ;

      // plot ctK slices
      csx[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
      effSliceX.back()->Draw();
      fsx.back()->Draw("same");
      effSliceX.back()->Draw("same");

      // if (i+j==0) {
      // 	leg->AddEntry(effHistsX.back(),"Binned efficiency" ,"lep");
      // 	leg->AddEntry(xframes.back()->findObject(Form("effx_%i_%i",i,j)),"KDE efficiency","l");
      // }
      // leg->Draw("same");

      csxC[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effCSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
      effCSliceX.back()->Draw();
      fsxC.back()->Draw("same");
      effCSliceX.back()->Draw("same");

      csxW[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effWSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
      effWSliceX.back()->Draw();
      fsxW.back()->Draw("same");
      effWSliceX.back()->Draw("same");

      csxcf[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      cfSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
      cfSliceX.back()->Draw();
      fsxcf.back()->Draw("same");
      cfSliceX.back()->Draw("same");

      csxmf[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      mfSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
      mfSliceX.back()->Draw();
      fsxmf.back()->Draw("same");
      mfSliceX.back()->Draw("same");

      // plot ctL slices
      csy[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
      effSliceY.back()->Draw();
      fsy.back()->Draw("same");
      effSliceY.back()->Draw("same");

      csyC[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effCSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
      effCSliceY.back()->Draw();
      fsyC.back()->Draw("same");
      effCSliceY.back()->Draw("same");

      csyW[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effWSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
      effWSliceY.back()->Draw();
      fsyW.back()->Draw("same");
      effWSliceY.back()->Draw("same");

      csycf[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      cfSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
      cfSliceY.back()->Draw();
      fsycf.back()->Draw("same");
      cfSliceY.back()->Draw("same");

      csymf[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      mfSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
      mfSliceY.back()->Draw();
      fsymf.back()->Draw("same");
      mfSliceY.back()->Draw("same");

      // plot phi slices
      csz[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
      effSliceZ.back()->Draw();
      fsz.back()->Draw("same");
      effSliceZ.back()->Draw("same");

      cszC[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effCSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
      effCSliceZ.back()->Draw();
      fszC.back()->Draw("same");
      effCSliceZ.back()->Draw("same");

      cszW[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      effWSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
      effWSliceZ.back()->Draw();
      fszW.back()->Draw("same");
      effWSliceZ.back()->Draw("same");

      cszcf[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      cfSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
      cfSliceZ.back()->Draw();
      fszcf.back()->Draw("same");
      cfSliceZ.back()->Draw("same");

      cszmf[confIndex]->cd(5*j+i+1);
      gPad->SetLeftMargin(0.18);
      mfSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
      mfSliceZ.back()->Draw();
      fszmf.back()->Draw("same");
      mfSliceZ.back()->Draw("same");

      // checking maximum values
      if ( maxEffX <effSliceX .back()->GetMaximum() ) maxEffX  = effSliceX .back()->GetMaximum();
      if ( maxEffCX<effCSliceX.back()->GetMaximum() ) maxEffCX = effCSliceX.back()->GetMaximum();
      if ( maxEffWX<effWSliceX.back()->GetMaximum() ) maxEffWX = effWSliceX.back()->GetMaximum();
      if ( maxCFX  <cfSliceX  .back()->GetMaximum() ) maxCFX   = cfSliceX  .back()->GetMaximum();
      if ( maxMFX  <mfSliceX  .back()->GetMaximum() ) maxMFX   = mfSliceX  .back()->GetMaximum();
      if ( maxEffY <effSliceY .back()->GetMaximum() ) maxEffY  = effSliceY .back()->GetMaximum();
      if ( maxEffCY<effCSliceY.back()->GetMaximum() ) maxEffCY = effCSliceY.back()->GetMaximum();
      if ( maxEffWY<effWSliceY.back()->GetMaximum() ) maxEffWY = effWSliceY.back()->GetMaximum();
      if ( maxCFY  <cfSliceY  .back()->GetMaximum() ) maxCFY   = cfSliceY  .back()->GetMaximum();
      if ( maxMFY  <mfSliceY  .back()->GetMaximum() ) maxMFY   = mfSliceY  .back()->GetMaximum();
      if ( maxEffZ <effSliceZ .back()->GetMaximum() ) maxEffZ  = effSliceZ .back()->GetMaximum();
      if ( maxEffCZ<effCSliceZ.back()->GetMaximum() ) maxEffCZ = effCSliceZ.back()->GetMaximum();
      if ( maxEffWZ<effWSliceZ.back()->GetMaximum() ) maxEffWZ = effWSliceZ.back()->GetMaximum();
      if ( maxCFZ  <cfSliceZ  .back()->GetMaximum() ) maxCFZ   = cfSliceZ  .back()->GetMaximum();
      if ( maxMFZ  <mfSliceZ  .back()->GetMaximum() ) maxMFZ   = mfSliceZ  .back()->GetMaximum();

    }    

  // set uniform y-axis ranges
  for (vector<TH1D*>::iterator hist = effSliceX .begin(); hist != effSliceX .end(); ++hist) (*hist)->SetMaximum(maxEffX *1.1);
  for (vector<TH1D*>::iterator hist = effCSliceX.begin(); hist != effCSliceX.end(); ++hist) (*hist)->SetMaximum(maxEffCX*1.1);
  for (vector<TH1D*>::iterator hist = effWSliceX.begin(); hist != effWSliceX.end(); ++hist) (*hist)->SetMaximum(maxEffWX*1.1);
  for (vector<TH1D*>::iterator hist = cfSliceX  .begin(); hist != cfSliceX  .end(); ++hist) (*hist)->SetMaximum(maxCFX  *1.1);
  for (vector<TH1D*>::iterator hist = mfSliceX  .begin(); hist != mfSliceX  .end(); ++hist) (*hist)->SetMaximum(maxMFX  *1.1);
  for (vector<TH1D*>::iterator hist = effSliceY .begin(); hist != effSliceY .end(); ++hist) (*hist)->SetMaximum(maxEffY *1.1);
  for (vector<TH1D*>::iterator hist = effCSliceY.begin(); hist != effCSliceY.end(); ++hist) (*hist)->SetMaximum(maxEffCY*1.1);
  for (vector<TH1D*>::iterator hist = effWSliceY.begin(); hist != effWSliceY.end(); ++hist) (*hist)->SetMaximum(maxEffWY*1.1);
  for (vector<TH1D*>::iterator hist = cfSliceY  .begin(); hist != cfSliceY  .end(); ++hist) (*hist)->SetMaximum(maxCFY  *1.1);
  for (vector<TH1D*>::iterator hist = mfSliceY  .begin(); hist != mfSliceY  .end(); ++hist) (*hist)->SetMaximum(maxMFY  *1.1);
  for (vector<TH1D*>::iterator hist = effSliceZ .begin(); hist != effSliceZ .end(); ++hist) (*hist)->SetMaximum(maxEffZ *1.1);
  for (vector<TH1D*>::iterator hist = effCSliceZ.begin(); hist != effCSliceZ.end(); ++hist) (*hist)->SetMaximum(maxEffCZ*1.1);
  for (vector<TH1D*>::iterator hist = effWSliceZ.begin(); hist != effWSliceZ.end(); ++hist) (*hist)->SetMaximum(maxEffWZ*1.1);
  for (vector<TH1D*>::iterator hist = cfSliceZ  .begin(); hist != cfSliceZ  .end(); ++hist) (*hist)->SetMaximum(maxCFZ  *1.1);
  for (vector<TH1D*>::iterator hist = mfSliceZ  .begin(); hist != mfSliceZ  .end(); ++hist) (*hist)->SetMaximum(maxMFZ  *1.1);

  csx  [confIndex]->SaveAs( (confString+Form("_eff_CTKslices_comp_dp%i.pdf"    ,(int)(border*200))).c_str() );
  csxC [confIndex]->SaveAs( (confString+Form("_eff-ct_CTKslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  csxW [confIndex]->SaveAs( (confString+Form("_eff-wt_CTKslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  csxmf[confIndex]->SaveAs( (confString+Form("_mt-frac_CTKslices_comp_dp%i.pdf",(int)(border*200))).c_str() );
  csxcf[confIndex]->SaveAs( (confString+Form("_ct-frac_CTKslices_comp_dp%i.pdf",(int)(border*200))).c_str() );
  csy  [confIndex]->SaveAs( (confString+Form("_eff_CTLslices_comp_dp%i.pdf"    ,(int)(border*200))).c_str() );
  csyC [confIndex]->SaveAs( (confString+Form("_eff-ct_CTLslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  csyW [confIndex]->SaveAs( (confString+Form("_eff-wt_CTLslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  csymf[confIndex]->SaveAs( (confString+Form("_mt-frac_CTLslices_comp_dp%i.pdf",(int)(border*200))).c_str() );
  csycf[confIndex]->SaveAs( (confString+Form("_ct-frac_CTLslices_comp_dp%i.pdf",(int)(border*200))).c_str() );
  csz  [confIndex]->SaveAs( (confString+Form("_eff_PHIslices_comp_dp%i.pdf"    ,(int)(border*200))).c_str() );
  cszC [confIndex]->SaveAs( (confString+Form("_eff-ct_PHIslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  cszW [confIndex]->SaveAs( (confString+Form("_eff-wt_PHIslices_comp_dp%i.pdf" ,(int)(border*200))).c_str() );
  cszmf[confIndex]->SaveAs( (confString+Form("_mt-frac_PHIslices_comp_dp%i.pdf",(int)(border*200))).c_str() );
  cszcf[confIndex]->SaveAs( (confString+Form("_ct-frac_PHIslices_comp_dp%i.pdf",(int)(border*200))).c_str() );

  // Plot projections of efficiency functions
  // test: method 1
  // auto effProj_xyz = (TH3D*)effC->createHistogram( Form("effProj_xyz_%s"  ,shortString.c_str()), *ctK, Binning(50,-1,1), YVar(*ctL, Binning(50,-1,1)), ZVar(*phi, Binning(50,-TMath::Pi(),TMath::Pi())), Extended(kFALSE) );
  // auto effProj_x  = effProj_xyz->ProjectionX();
  // auto effProj_y  = effProj_xyz->ProjectionY();
  // auto effProj_z  = effProj_xyz->ProjectionZ();
  // auto effProj_xy = effProj_xyz->Project3D("xy");
  // auto effProj_xz = effProj_xyz->Project3D("xz");
  // auto effProj_yz = effProj_xyz->Project3D("yz");
  // // test: method 2
  // auto effCProj_x  = effCHist->ProjectionX();
  // auto effCProj_y  = effCHist->ProjectionY();
  // auto effCProj_z  = effCHist->ProjectionZ();
  // auto effCProj_xy = effCHist->Project3D("xy");
  // auto effCProj_xz = effCHist->Project3D("xz");
  // auto effCProj_yz = effCHist->Project3D("yz");
  // // test: method 4
  // RooPlot* cfProj_x = ctK->frame( Name( Form("cfProjx_%s"  ,shortString.c_str()) ) );
  // RooPlot* cfProj_y = ctL->frame( Name( Form("cfProjy_%s"  ,shortString.c_str()) ) );
  // RooPlot* cfProj_z = phi->frame( Name( Form("cfProjz_%s"  ,shortString.c_str()) ) );
  // effC->plotOn( cfProj_x, Project( RooArgSet(*ctL,*phi) ) );
  // effC->plotOn( cfProj_y, Project( RooArgSet(*ctK,*phi) ) );
  // effC->plotOn( cfProj_z, Project( RooArgSet(*ctK,*ctL) ) );
  // auto cfProj_xy = (TH2D*)effC->createHistogram( Form("cfProjxy_%s"  ,shortString.c_str()), *ctK, Binning(50,-1,1), YVar(*ctL, Binning(50,-1,1)), Scaling(kFALSE), Extended(kFALSE) );
  // auto cfProj_xz = (TH2D*)effC->createHistogram( Form("cfProjxy_%s"  ,shortString.c_str()), *ctK, Binning(50,-1,1), YVar(*phi, Binning(50,-TMath::Pi(),TMath::Pi())), ConditionalObservables(RooArgSet(*ctL)), Scaling(kFALSE), Extended(kFALSE) );
  // auto cfProj_yz = (TH2D*)effC->createHistogram( Form("cfProjxy_%s"  ,shortString.c_str()), *ctL, Binning(50,-1,1), YVar(*phi, Binning(50,-TMath::Pi(),TMath::Pi())), ConditionalObservables(RooArgSet(*ctK)), Scaling(kFALSE), Extended(kFALSE) );
  // // test: method 5
  // auto mfProj_x  = (TH1D*)effC->createHistogram( Form("mfProjx_%s"  ,shortString.c_str()), *ctK, Binning(50,-1,1), Extended(kFALSE) );
  // auto mfProj_y  = (TH1D*)effC->createHistogram( Form("mfProjy_%s"  ,shortString.c_str()), *ctL, Binning(50,-1,1), ConditionalObservables(RooArgSet(*ctK,*phi)), Extended(kFALSE) );
  // auto mfProj_z  = (TH1D*)effC->createHistogram( Form("mfProjz_%s"  ,shortString.c_str()), *phi, Binning(50,-TMath::Pi(),TMath::Pi()), Scaling(kFALSE), Extended(kFALSE) );
  // auto mfProj_xy = (TH2D*)effC->createHistogram( Form("mfProjxy_%s"  ,shortString.c_str()), *ctK, Binning(50,-1,1), YVar(*ctL, Binning(50,-1,1)), Extended(kFALSE) );
  // auto mfProj_xz = (TH2D*)effC->createHistogram( Form("mfProjxy_%s"  ,shortString.c_str()), *ctK, Binning(50,-1,1), YVar(*phi, Binning(50,-TMath::Pi(),TMath::Pi())), ConditionalObservables(RooArgSet(*ctL)), Extended(kFALSE) );
  // auto mfProj_yz = (TH2D*)effC->createHistogram( Form("mfProjxy_%s"  ,shortString.c_str()), *ctL, Binning(50,-1,1), YVar(*phi, Binning(50,-TMath::Pi(),TMath::Pi())), ConditionalObservables(RooArgSet(*ctK)), Extended(kFALSE) );
  // USING METHOD 2
  auto effProj_x  = effHist->ProjectionX();
  auto effProj_y  = effHist->ProjectionY();
  auto effProj_z  = effHist->ProjectionZ();
  auto effProj_xy = effHist->Project3D("xy");
  auto effProj_xz = effHist->Project3D("xz");
  auto effProj_yz = effHist->Project3D("yz");
  auto effCProj_x  = effCHist->ProjectionX();
  auto effCProj_y  = effCHist->ProjectionY();
  auto effCProj_z  = effCHist->ProjectionZ();
  auto effCProj_xy = effCHist->Project3D("xy");
  auto effCProj_xz = effCHist->Project3D("xz");
  auto effCProj_yz = effCHist->Project3D("yz");
  auto effWProj_x  = effWHist->ProjectionX();
  auto effWProj_y  = effWHist->ProjectionY();
  auto effWProj_z  = effWHist->ProjectionZ();
  auto effWProj_xy = effWHist->Project3D("xy");
  auto effWProj_xz = effWHist->Project3D("xz");
  auto effWProj_yz = effWHist->Project3D("yz");
  auto cfProj_x  = corrFracHist->ProjectionX();
  auto cfProj_y  = corrFracHist->ProjectionY();
  auto cfProj_z  = corrFracHist->ProjectionZ();
  auto cfProj_xy = corrFracHist->Project3D("xy");
  auto cfProj_xz = corrFracHist->Project3D("xz");
  auto cfProj_yz = corrFracHist->Project3D("yz");
  auto mfProj_x  = mistFracHist->ProjectionX();
  auto mfProj_y  = mistFracHist->ProjectionY();
  auto mfProj_z  = mistFracHist->ProjectionZ();
  auto mfProj_xy = mistFracHist->Project3D("xy");
  auto mfProj_xz = mistFracHist->Project3D("xz");
  auto mfProj_yz = mistFracHist->Project3D("yz");
  // Plot 1D projections
  cp1[confIndex] = new TCanvas(Form("cp1_%s",shortString.c_str()),Form("cp1_%s",shortString.c_str()),1800,800);
  effProj_x->SetTitle( Form("Full efficiency projection (q2-bin %i, %s);cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  effProj_y->SetTitle( Form("Full efficiency projection (q2-bin %i, %s);cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
  effProj_z->SetTitle( Form("Full efficiency projection (q2-bin %i, %s);#phi",q2Bin,(parity==0?"even":"odd")) );
  effProj_x->Scale(1.0/2500);
  effProj_y->Scale(1.0/2500);
  effProj_z->Scale(1.0/2500);
  effProj_x->SetMinimum(0.0);
  effProj_y->SetMinimum(0.0);
  effProj_z->SetMinimum(0.0);
  cp1[confIndex]->Divide(3,1);
  cp1[confIndex]->cd(1);
  effProj_x->Draw("LHIST");
  cp1[confIndex]->cd(2);
  effProj_y->Draw("LHIST");
  cp1[confIndex]->cd(3);
  effProj_z->Draw("LHIST");
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
  cp1cf[confIndex] = new TCanvas(Form("cp1cf_%s",shortString.c_str()),Form("cp1cf_%s",shortString.c_str()),1800,800);
  cfProj_x->SetTitle( Form("Correct-tag fraction projection (q2-bin %i, %s);cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  cfProj_y->SetTitle( Form("Correct-tag fraction projection (q2-bin %i, %s);cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
  cfProj_z->SetTitle( Form("Correct-tag fraction projection (q2-bin %i, %s);#phi",q2Bin,(parity==0?"even":"odd")) );
  cfProj_x->Scale(1.0/2500);
  cfProj_y->Scale(1.0/2500);
  cfProj_z->Scale(1.0/2500);
  cfProj_x->SetMinimum(0.0);
  cfProj_y->SetMinimum(0.0);
  cfProj_z->SetMinimum(0.0);
  cp1cf[confIndex]->Divide(3,1);
  cp1cf[confIndex]->cd(1);
  cfProj_x->Draw("LHIST");
  cp1cf[confIndex]->cd(2);
  cfProj_y->Draw("LHIST");
  cp1cf[confIndex]->cd(3);
  cfProj_z->Draw("LHIST");
  cp1mf[confIndex] = new TCanvas(Form("cp1mf_%s",shortString.c_str()),Form("cp1mf_%s",shortString.c_str()),1800,800);
  mfProj_x->SetTitle( Form("Mistag fraction projection (q2-bin %i, %s);cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  mfProj_y->SetTitle( Form("Mistag fraction projection (q2-bin %i, %s);cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
  mfProj_z->SetTitle( Form("Mistag fraction projection (q2-bin %i, %s);#phi",q2Bin,(parity==0?"even":"odd")) );
  mfProj_x->Scale(1.0/2500);
  mfProj_y->Scale(1.0/2500);
  mfProj_z->Scale(1.0/2500);
  mfProj_x->SetMinimum(0.0);
  mfProj_y->SetMinimum(0.0);
  mfProj_z->SetMinimum(0.0);
  cp1mf[confIndex]->Divide(3,1);
  cp1mf[confIndex]->cd(1);
  mfProj_x->Draw("LHIST");
  cp1mf[confIndex]->cd(2);
  mfProj_y->Draw("LHIST");
  cp1mf[confIndex]->cd(3);
  mfProj_z->Draw("LHIST");
  // h3_x->SetLineColor(2);
  // h3_y->SetLineColor(2);
  // h3_z->SetLineColor(2);
  // 2D projections
  cp2[confIndex] = new TCanvas(Form("cp2_%s",shortString.c_str()),Form("cp2_%s",shortString.c_str()),1800,800);
  effProj_xy->SetTitle( Form("Full efficiency projection (q2-bin %i, %s);cos(#theta_{L});cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  effProj_xz->SetTitle( Form("Full efficiency projection (q2-bin %i, %s);#phi;cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  effProj_yz->SetTitle( Form("Full efficiency projection (q2-bin %i, %s);#phi;cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
  effProj_xy->Scale(1.0/50);
  effProj_xz->Scale(1.0/50);
  effProj_yz->Scale(1.0/50);
  effProj_xy->GetXaxis()->SetTitleOffset(1.4);
  effProj_xz->GetXaxis()->SetTitleOffset(1.4);
  effProj_yz->GetXaxis()->SetTitleOffset(1.4);
  effProj_xy->GetYaxis()->SetTitleOffset(2);
  effProj_xz->GetYaxis()->SetTitleOffset(2);
  effProj_yz->GetYaxis()->SetTitleOffset(2);
  effProj_xy->SetMinimum(0.0);
  effProj_xz->SetMinimum(0.0);
  effProj_yz->SetMinimum(0.0);
  cp2[confIndex]->Divide(3,1);
  cp2[confIndex]->cd(1);
  effProj_xy->Draw("SURF3");
  cp2[confIndex]->cd(2);
  effProj_xz->Draw("SURF3");
  cp2[confIndex]->cd(3);
  effProj_yz->Draw("SURF3");
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
  cp2cf[confIndex] = new TCanvas(Form("cp2cf_%s",shortString.c_str()),Form("cp2cf_%s",shortString.c_str()),1800,800);
  cfProj_xy->SetTitle( Form("Correct-tag fraction projection (q2-bin %i, %s);cos(#theta_{L});cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  cfProj_xz->SetTitle( Form("Correct-tag fraction projection (q2-bin %i, %s);#phi;cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  cfProj_yz->SetTitle( Form("Correct-tag fraction projection (q2-bin %i, %s);#phi;cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
  cfProj_xy->Scale(1.0/50);
  cfProj_xz->Scale(1.0/50);
  cfProj_yz->Scale(1.0/50);
  cfProj_xy->GetXaxis()->SetTitleOffset(1.4);
  cfProj_xz->GetXaxis()->SetTitleOffset(1.4);
  cfProj_yz->GetXaxis()->SetTitleOffset(1.4);
  cfProj_xy->GetYaxis()->SetTitleOffset(2);
  cfProj_xz->GetYaxis()->SetTitleOffset(2);
  cfProj_yz->GetYaxis()->SetTitleOffset(2);
  cfProj_xy->SetMinimum(0.0);
  cfProj_xz->SetMinimum(0.0);
  cfProj_yz->SetMinimum(0.0);
  cp2cf[confIndex]->Divide(3,1);
  cp2cf[confIndex]->cd(1);
  cfProj_xy->Draw("SURF3");
  cp2cf[confIndex]->cd(2);
  cfProj_xz->Draw("SURF3");
  cp2cf[confIndex]->cd(3);
  cfProj_yz->Draw("SURF3");
  cp2mf[confIndex] = new TCanvas(Form("cp2mf_%s",shortString.c_str()),Form("cp2mf_%s",shortString.c_str()),1800,800);
  mfProj_xy->SetTitle( Form("Mistag fraction projection (q2-bin %i, %s);cos(#theta_{L});cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  mfProj_xz->SetTitle( Form("Mistag fraction projection (q2-bin %i, %s);#phi;cos(#theta_{K})",q2Bin,(parity==0?"even":"odd")) );
  mfProj_yz->SetTitle( Form("Mistag fraction projection (q2-bin %i, %s);#phi;cos(#theta_{L})",q2Bin,(parity==0?"even":"odd")) );
  mfProj_xy->Scale(1.0/50);
  mfProj_xz->Scale(1.0/50);
  mfProj_yz->Scale(1.0/50);
  mfProj_xy->GetXaxis()->SetTitleOffset(1.4);
  mfProj_xz->GetXaxis()->SetTitleOffset(1.4);
  mfProj_yz->GetXaxis()->SetTitleOffset(1.4);
  mfProj_xy->GetYaxis()->SetTitleOffset(2);
  mfProj_xz->GetYaxis()->SetTitleOffset(2);
  mfProj_yz->GetYaxis()->SetTitleOffset(2);
  mfProj_xy->SetMinimum(0.0);
  mfProj_xz->SetMinimum(0.0);
  mfProj_yz->SetMinimum(0.0);
  cp2mf[confIndex]->Divide(3,1);
  cp2mf[confIndex]->cd(1);
  mfProj_xy->Draw("SURF3");
  cp2mf[confIndex]->cd(2);
  mfProj_xz->Draw("SURF3");
  cp2mf[confIndex]->cd(3);
  mfProj_yz->Draw("SURF3");

  cp1  [confIndex]->SaveAs( (confString+"_eff_1DProj.pdf"    ).c_str() );    
  cp1C [confIndex]->SaveAs( (confString+"_eff-ct_1DProj.pdf" ).c_str() );    
  cp1W [confIndex]->SaveAs( (confString+"_eff-wt_1DProj.pdf" ).c_str() );    
  cp1cf[confIndex]->SaveAs( (confString+"_ct-frac_1DProj.pdf").c_str() );    
  cp1mf[confIndex]->SaveAs( (confString+"_mt-frac_1DProj.pdf").c_str() );    
  cp2  [confIndex]->SaveAs( (confString+"_eff_2DProj.pdf"    ).c_str() );    
  cp2C [confIndex]->SaveAs( (confString+"_eff-ct_2DProj.pdf" ).c_str() );    
  cp2W [confIndex]->SaveAs( (confString+"_eff-wt_2DProj.pdf" ).c_str() );    
  cp2cf[confIndex]->SaveAs( (confString+"_ct-frac_2DProj.pdf").c_str() );    
  cp2mf[confIndex]->SaveAs( (confString+"_mt-frac_2DProj.pdf").c_str() );    

  if (!doClosure) return;

  int nbins_ct = 30;

  // create and fill GEN histograms
  TH1D* hXctGEN_effC = new TH1D( Form("hXctGEN_effC_%s",shortString.c_str()), Form("hXctGEN_effC_%s",shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hXctGEN_cf   = new TH1D( Form("hXctGEN_effW_%s",shortString.c_str()), Form("hXctGEN_effW_%s",shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hXctGEN_effW = new TH1D( Form("hXctGEN_cf_%s"  ,shortString.c_str()), Form("hXctGEN_cf_%s"  ,shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hXctGEN_mf   = new TH1D( Form("hXctGEN_mf_%s"  ,shortString.c_str()), Form("hXctGEN_mf_%s"  ,shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hYctGEN_effC = new TH1D( Form("hYctGEN_effC_%s",shortString.c_str()), Form("hYctGEN_effC_%s",shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hYctGEN_cf   = new TH1D( Form("hYctGEN_effW_%s",shortString.c_str()), Form("hYctGEN_effW_%s",shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hYctGEN_effW = new TH1D( Form("hYctGEN_cf_%s"  ,shortString.c_str()), Form("hYctGEN_cf_%s"  ,shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hYctGEN_mf   = new TH1D( Form("hYctGEN_mf_%s"  ,shortString.c_str()), Form("hYctGEN_mf_%s"  ,shortString.c_str()), nbins_ct, -1, 1 );
  TH1D* hZctGEN_effC = new TH1D( Form("hZctGEN_effC_%s",shortString.c_str()), Form("hZctGEN_effC_%s",shortString.c_str()), nbins_ct, -TMath::Pi(), TMath::Pi() );
  TH1D* hZctGEN_cf   = new TH1D( Form("hZctGEN_effW_%s",shortString.c_str()), Form("hZctGEN_effW_%s",shortString.c_str()), nbins_ct, -TMath::Pi(), TMath::Pi() );
  TH1D* hZctGEN_effW = new TH1D( Form("hZctGEN_cf_%s"  ,shortString.c_str()), Form("hZctGEN_cf_%s"  ,shortString.c_str()), nbins_ct, -TMath::Pi(), TMath::Pi() );
  TH1D* hZctGEN_mf   = new TH1D( Form("hZctGEN_mf_%s"  ,shortString.c_str()), Form("hZctGEN_mf_%s"  ,shortString.c_str()), nbins_ct, -TMath::Pi(), TMath::Pi() );
  hXctGEN_effC->Sumw2();
  hXctGEN_effW->Sumw2();
  hXctGEN_cf  ->Sumw2();
  hXctGEN_mf  ->Sumw2();
  hYctGEN_effC->Sumw2();
  hYctGEN_effW->Sumw2();
  hYctGEN_cf  ->Sumw2();
  hYctGEN_mf  ->Sumw2();
  hZctGEN_effC->Sumw2();
  hZctGEN_effW->Sumw2();
  hZctGEN_cf  ->Sumw2();
  hZctGEN_mf  ->Sumw2();
  for (int iEv=0; iEv<data_genDen->numEntries(); ++iEv) {
    const RooArgSet* set = data_genDen->get(iEv);
    double gen_ctK = ((RooRealVar*)set->find("ctK"))->getVal();
    double gen_ctL = ((RooRealVar*)set->find("ctL"))->getVal();
    double gen_phi = ((RooRealVar*)set->find("phi"))->getVal();
    ctK->setVal(gen_ctK);
    ctL->setVal(gen_ctL);
    phi->setVal(gen_phi);
    double wei = effC->getVal();
    hXctGEN_effC->Fill( gen_ctK, wei );
    hYctGEN_effC->Fill( gen_ctL, wei );
    hZctGEN_effC->Fill( gen_phi, wei );
    ctK->setVal(-1*gen_ctK);
    ctL->setVal(-1*gen_ctL);
    phi->setVal(-1*gen_phi);
    wei = effW->getVal();
    hXctGEN_effW->Fill( -1*gen_ctK, wei );
    hYctGEN_effW->Fill( -1*gen_ctL, wei );
    hZctGEN_effW->Fill( -1*gen_phi, wei );
  }
  for (int iEv=0; iEv<data_ctRECO->numEntries(); ++iEv) {
    const RooArgSet* set = data_ctRECO->get(iEv);
    double reco_ctK = ((RooRealVar*)set->find("ctK"))->getVal();
    double reco_ctL = ((RooRealVar*)set->find("ctL"))->getVal();
    double reco_phi = ((RooRealVar*)set->find("phi"))->getVal();
    ctK->setVal(reco_ctK);
    ctL->setVal(reco_ctL);
    phi->setVal(reco_phi);
    double wei = corrFrac->getVal() * data_ctRECO->weight();
    hXctGEN_cf->Fill( reco_ctK, wei );
    hYctGEN_cf->Fill( reco_ctL, wei );
    hZctGEN_cf->Fill( reco_phi, wei );
    wei = mistFrac->getVal() * data_ctRECO->weight();
    hXctGEN_mf->Fill( reco_ctK, wei );
    hYctGEN_mf->Fill( reco_ctL, wei );
    hZctGEN_mf->Fill( reco_phi, wei );
  }
  for (int iEv=0; iEv<data_wtRECO->numEntries(); ++iEv) {
    const RooArgSet* set = data_wtRECO->get(iEv);
    double reco_ctK = ((RooRealVar*)set->find("ctK"))->getVal();
    double reco_ctL = ((RooRealVar*)set->find("ctL"))->getVal();
    double reco_phi = ((RooRealVar*)set->find("phi"))->getVal();
    ctK->setVal(reco_ctK);
    ctL->setVal(reco_ctL);
    phi->setVal(reco_phi);
    double wei = corrFrac->getVal() * data_wtRECO->weight();
    hXctGEN_cf->Fill( reco_ctK, wei );
    hYctGEN_cf->Fill( reco_ctL, wei );
    hZctGEN_cf->Fill( reco_phi, wei );
    wei = mistFrac->getVal() * data_wtRECO->weight();
    hXctGEN_mf->Fill( reco_ctK, wei );
    hYctGEN_mf->Fill( reco_ctL, wei );
    hZctGEN_mf->Fill( reco_phi, wei );
  }

  // create RECO histograms
  TH1D* hXctREC_effC = (TH1D*)data_ctRECO->createHistogram( Form("hXctREC_effC_%s",shortString.c_str()), *ctK, Binning(nbins_ct,-1,1) );
  TH1D* hXctREC_effW = (TH1D*)data_wtRECO->createHistogram( Form("hXctREC_effW_%s",shortString.c_str()), *ctK, Binning(nbins_ct,-1,1) );
  TH1D* hYctREC_effC = (TH1D*)data_ctRECO->createHistogram( Form("hYctREC_effC_%s",shortString.c_str()), *ctL, Binning(nbins_ct,-1,1) );
  TH1D* hYctREC_effW = (TH1D*)data_wtRECO->createHistogram( Form("hYctREC_effW_%s",shortString.c_str()), *ctL, Binning(nbins_ct,-1,1) );
  TH1D* hZctREC_effC = (TH1D*)data_ctRECO->createHistogram( Form("hZctREC_effC_%s",shortString.c_str()), *phi, Binning(nbins_ct,-TMath::Pi(),TMath::Pi()) );
  TH1D* hZctREC_effW = (TH1D*)data_wtRECO->createHistogram( Form("hZctREC_effW_%s",shortString.c_str()), *phi, Binning(nbins_ct,-TMath::Pi(),TMath::Pi()) );

  // prepare for plot
  hXctGEN_effC->SetTitle( Form("Closure test of correct-tag efficiency (q2-bin %i, %s);cos(#theta_{K});Events",q2Bin,(parity==0?"even":"odd")) );
  hXctGEN_effW->SetTitle( Form(  "Closure test of wrong-tag efficiency (q2-bin %i, %s);cos(#theta_{K});Events",q2Bin,(parity==0?"even":"odd")) );
  hXctGEN_cf  ->SetTitle( Form(  "Closure test of correct-tag fraction (q2-bin %i, %s);cos(#theta_{K});Events",q2Bin,(parity==0?"even":"odd")) );
  hXctGEN_mf  ->SetTitle( Form(       "Closure test of mistag fraction (q2-bin %i, %s);cos(#theta_{K});Events",q2Bin,(parity==0?"even":"odd")) );
  hYctGEN_effC->SetTitle( Form("Closure test of correct-tag efficiency (q2-bin %i, %s);cos(#theta_{L});Events",q2Bin,(parity==0?"even":"odd")) );
  hYctGEN_effW->SetTitle( Form(  "Closure test of wrong-tag efficiency (q2-bin %i, %s);cos(#theta_{L});Events",q2Bin,(parity==0?"even":"odd")) );
  hYctGEN_cf  ->SetTitle( Form(  "Closure test of correct-tag fraction (q2-bin %i, %s);cos(#theta_{L});Events",q2Bin,(parity==0?"even":"odd")) );
  hYctGEN_mf  ->SetTitle( Form(       "Closure test of mistag fraction (q2-bin %i, %s);cos(#theta_{L});Events",q2Bin,(parity==0?"even":"odd")) );
  hZctGEN_effC->SetTitle( Form("Closure test of correct-tag efficiency (q2-bin %i, %s);#phi;Events",q2Bin,(parity==0?"even":"odd")) );
  hZctGEN_effW->SetTitle( Form(  "Closure test of wrong-tag efficiency (q2-bin %i, %s);#phi;Events",q2Bin,(parity==0?"even":"odd")) );
  hZctGEN_cf  ->SetTitle( Form(  "Closure test of correct-tag fraction (q2-bin %i, %s);#phi;Events",q2Bin,(parity==0?"even":"odd")) );
  hZctGEN_mf  ->SetTitle( Form(       "Closure test of mistag fraction (q2-bin %i, %s);#phi;Events",q2Bin,(parity==0?"even":"odd")) );
  hXctGEN_effC->Scale( hXctREC_effC->Integral() / hXctGEN_effC->Integral() );
  hXctGEN_effW->Scale( hXctREC_effW->Integral() / hXctGEN_effW->Integral() );
  hYctGEN_effC->Scale( hYctREC_effC->Integral() / hYctGEN_effC->Integral() );
  hYctGEN_effW->Scale( hYctREC_effW->Integral() / hYctGEN_effW->Integral() );
  hZctGEN_effC->Scale( hZctREC_effC->Integral() / hZctGEN_effC->Integral() );
  hZctGEN_effW->Scale( hZctREC_effW->Integral() / hZctGEN_effW->Integral() );
  hXctGEN_effC->SetMinimum(0);
  hXctGEN_effW->SetMinimum(0);
  hXctGEN_cf  ->SetMinimum(0);
  hXctGEN_mf  ->SetMinimum(0);
  hYctGEN_effC->SetMinimum(0);
  hYctGEN_effW->SetMinimum(0);
  hYctGEN_cf  ->SetMinimum(0);
  hYctGEN_mf  ->SetMinimum(0);
  hZctGEN_effC->SetMinimum(0);
  hZctGEN_effW->SetMinimum(0);
  hZctGEN_cf  ->SetMinimum(0);
  hZctGEN_mf  ->SetMinimum(0);
  hXctGEN_effC->SetLineColor(kRed+1);
  hXctGEN_effW->SetLineColor(kRed+1);
  hXctGEN_cf->SetLineColor(kRed+1);
  hXctGEN_mf->SetLineColor(kRed+1);
  hYctGEN_effC->SetLineColor(kRed+1);
  hYctGEN_effW->SetLineColor(kRed+1);
  hYctGEN_cf->SetLineColor(kRed+1);
  hYctGEN_mf->SetLineColor(kRed+1);
  hZctGEN_effC->SetLineColor(kRed+1);
  hZctGEN_effW->SetLineColor(kRed+1);
  hZctGEN_cf->SetLineColor(kRed+1);
  hZctGEN_mf->SetLineColor(kRed+1);
  hXctGEN_effC->SetLineWidth(2);
  hXctGEN_effW->SetLineWidth(2);
  hXctGEN_cf->SetLineWidth(2);
  hXctGEN_mf->SetLineWidth(2);
  hYctGEN_effC->SetLineWidth(2);
  hYctGEN_effW->SetLineWidth(2);
  hYctGEN_cf->SetLineWidth(2);
  hYctGEN_mf->SetLineWidth(2);
  hZctGEN_effC->SetLineWidth(2);
  hZctGEN_effW->SetLineWidth(2);
  hZctGEN_cf->SetLineWidth(2);
  hZctGEN_mf->SetLineWidth(2);
  hXctREC_effC->SetLineWidth(2);
  hXctREC_effW->SetLineWidth(2);
  hYctREC_effC->SetLineWidth(2);
  hYctREC_effW->SetLineWidth(2);
  hZctREC_effC->SetLineWidth(2);
  hZctREC_effW->SetLineWidth(2);


  // Plot closure test
  cctC [confIndex] = new TCanvas(("ccteffC_"+shortString).c_str(),("ccteffC_"+shortString).c_str(),2000,700);
  cctW [confIndex] = new TCanvas(("ccteffW_"+shortString).c_str(),("ccteffW_"+shortString).c_str(),2000,700);
  cctcf[confIndex] = new TCanvas(("cctcf_"  +shortString).c_str(),("cctcf_"  +shortString).c_str(),2000,700);
  cctmf[confIndex] = new TCanvas(("cctmf_"  +shortString).c_str(),("cctmf_"  +shortString).c_str(),2000,700);
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

  gPad->SetLeftMargin(0.17); 

  cctC [confIndex]->Divide(3,1);
  cctW [confIndex]->Divide(3,1);
  cctcf[confIndex]->Divide(3,1);
  cctmf[confIndex]->Divide(3,1);
  cctC[confIndex]->cd(1);
  hXctGEN_effC->Draw();
  hXctREC_effC->Draw("same");
  cctW[confIndex]->cd(1); 
  hXctGEN_effW->Draw();
  hXctREC_effW->Draw("same");
  cctcf[confIndex]->cd(1);
  hXctGEN_cf->Draw();
  hXctREC_effC->Draw("same");
  cctmf[confIndex]->cd(1);
  hXctGEN_mf->Draw();
  hXctREC_effW->Draw("same");
  cctC[confIndex]->cd(2);
  hYctGEN_effC->Draw();
  hYctREC_effC->Draw("same");
  cctW[confIndex]->cd(2); 
  hYctGEN_effW->Draw();
  hYctREC_effW->Draw("same");
  cctcf[confIndex]->cd(2);
  hYctGEN_cf->Draw();
  hYctREC_effC->Draw("same");
  cctmf[confIndex]->cd(2);
  hYctGEN_mf->Draw();
  hYctREC_effW->Draw("same");
  cctC[confIndex]->cd(3);
  hZctGEN_effC->Draw();
  hZctREC_effC->Draw("same");
  cctW[confIndex]->cd(3); 
  hZctGEN_effW->Draw();
  hZctREC_effW->Draw("same");
  cctcf[confIndex]->cd(3);
  hZctGEN_cf->Draw();
  hZctREC_effC->Draw("same");
  cctmf[confIndex]->cd(3);
  hZctGEN_mf->Draw();
  hZctREC_effW->Draw("same");

  cctC [confIndex]->SaveAs( (confString+"_eff-ct_ClosureTest.pdf" ).c_str() );    
  cctW [confIndex]->SaveAs( (confString+"_eff-wt_ClosureTest.pdf" ).c_str() );    
  cctcf[confIndex]->SaveAs( (confString+"_ct-frac_ClosureTest.pdf").c_str() );    
  cctmf[confIndex]->SaveAs( (confString+"_mt-frac_ClosureTest.pdf").c_str() );    

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
