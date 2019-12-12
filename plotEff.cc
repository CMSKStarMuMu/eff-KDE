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

TCanvas* csC [3][2*nBins];
TCanvas* csW [3][2*nBins];

TCanvas* cp1C [2*nBins];
TCanvas* cp2C [2*nBins];
TCanvas* cctC [2*nBins];
TCanvas* cp1W [2*nBins];
TCanvas* cp2W [2*nBins];
TCanvas* cctW [2*nBins];


TH1D* Invert1Dhist(TH1D* hin, string hname);

void plotEffBin(int q2Bin, int parity, bool doClosure, int year)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  string confString = "plotEff_d/effKDE_test_"+shortString;
  gStyle->SetOptStat(0);

  bool doCT = true;
  bool doWT = true;

  int confIndex = nBins*parity + q2Bin;

  int distBins = 10;

  const char *varCoord[3];
  varCoord[0] = "X";
  varCoord[1] = "Y";
  varCoord[2] = "Z";

  const char *varNames[3];
  varNames[0] = "ctK";
  varNames[1] = "ctL";
  varNames[2] = "phi";

  // Load variables and dataset
  string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/effDataset_b%i_%i.root",year,q2Bin,year);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  // import both datasets (correlated and uncorrelated statistics)
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
  std::vector<RooRealVar*> vars_vec;
  vars_vec.push_back(ctK);
  vars_vec.push_back(ctL);
  vars_vec.push_back(phi);
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
 
  // import number of total genDen statistics, before FSR-veto
  // normVal = num GEN events tot / num GEN events w/o FSR ( > 1)
  double normVal = 1;
  double normVal_corr = 1;
  TH1I* normHist = (TH1I*)fin_data->Get( Form("n_genDen_b%i",q2Bin) );
  if ( !normHist || normHist->IsZombie() ) cout<<"Histogram n_genDen_b"<<q2Bin<<" not found in file: "<<filename_data<<"\nThe normalisation of the genDen histograms will not be corrected to pre-FSR-veto values\nand the efficiency histograms and functions will have different scaling"<<endl;
  else {
    normVal      = 1.0 * normHist->GetBinContent(2-parity) / data_genDen     ->sumEntries(); // opposite parity wrt efficiency
    normVal_corr = 1.0 * normHist->GetBinContent(parity+1) / data_genDen_corr->sumEntries(); // coherent parity
  }
  
  // import KDE efficiency histograms
  string filename = Form((parity==0?"files/KDEeff_b%i_ev_%i.root":"files/KDEeff_b%i_od_%i.root"),q2Bin, year);
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
    // sara: interpolation order set to 1 -> is it ok? enough?
    effC = new RooHistFunc("effC","effC",vars,*effCData,1);
  }
  if (doWT) {
    effWData = new RooDataHist("effWData","effWData",vars,effWHist);
    effW = new RooHistFunc("effW","effW",vars,*effWData,1);
  }


  // Plot 1D slices of the efficiency function and binned efficiency
  for (int ivar=0; ivar<3; ivar++){
  if (doCT) {
      csC[ivar][confIndex] = new TCanvas((Form("cs%sC",varCoord[ivar])+shortString).c_str(),(shortString+Form("_effC_%s",varNames[ivar])).c_str(),1500,1500) ;
      csC[ivar][confIndex]->Divide(5,5);
    }
  if (doWT) {
      csW[ivar][confIndex] = new TCanvas((Form("cs%sW",varCoord[ivar])+shortString).c_str(),(shortString+Form("_effW_%s",varNames[ivar])).c_str(),1500,1500) ;
      csW[ivar][confIndex]->Divide(5,5);
    }
  }

  vector <vector <TH1D*>>    effCSlice(3, std::vector<TH1D*>(0));
  vector <vector <TH1D*>>    effWSlice(3, std::vector<TH1D*>(0));
  vector <vector <TH1D*>>    effCSlice_corr(3, std::vector<TH1D*>(0));
  vector <vector <TH1D*>>    effWSlice_corr(3, std::vector<TH1D*>(0));
  vector <vector <RooPlot*>> fsC(3,std::vector<RooPlot*>(0) );
  vector <vector <RooPlot*>> fsW(3,std::vector<RooPlot*>(0) );

  // TLegend* leg = new TLegend (0.35,0.8,0.9,0.9);

  // width of the slices in the hidden variables ("border" is half of it)
  double border = 0.04;

  // variables to be filled with global efficiency maximum
  double maxEffC[3] = {0,0,0};
  double maxEffW[3] = {0,0,0};

  // loop over slice grid
  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

      std::vector<string> ctCuts, wtCuts;
      // central values and borders of the slices in the hidden variables
      double centA = -0.8 + 1.6*i/4;
      double centB = -0.8 + 1.6*j/4;
      ctCuts.push_back( Form("fabs(ctL-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centA,border,TMath::Pi(),centB,border)); // cutX
      ctCuts.push_back( Form("fabs(ctK-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centA,border,TMath::Pi(),centB,border));  //cutY
      ctCuts.push_back( Form("fabs(ctK-(%1.1f))<%1.3f && fabs(ctL-(%1.1f))<%1.3f",centA,border,centB,border));  // cytZ

      double centAwt = 0.8 - 1.6*i/4;
      double centBwt = 0.8 - 1.6*j/4;
      wtCuts.push_back( Form("fabs(ctL-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centAwt,border,TMath::Pi(),centBwt,border));
      wtCuts.push_back( Form("fabs(ctK-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centAwt,border,TMath::Pi(),centBwt,border));
      wtCuts.push_back( Form("fabs(ctK-(%1.1f))<%1.3f && fabs(ctL-(%1.1f))<%1.3f"        ,centAwt,border,centBwt,border));

      // slicing distributions
      TH1D* ctGENDenDist[3] = {0,0,0};
      TH1D* ctGENNumDist[3] = {0,0,0};
      TH1D* ctDenDist[3]    = {0,0,0};
      TH1D* ctRECODist[3]   = {0,0,0};

      TH1D* wtGENDenDist[3] = {0,0,0};
      TH1D* wtGENNumDist[3] = {0,0,0};
      TH1D* wtDenDist[3]    = {0,0,0};
      TH1D* wtRECODist[3]   = {0,0,0};

      // slicing distributions for the correlated dataset
      TH1D* ctGENDenDist_corr[3] = {0,0,0};
      TH1D* ctGENNumDist_corr[3] = {0,0,0};
      TH1D* ctDenDist_corr[3]    = {0,0,0};
      TH1D* ctRECODist_corr[3]   = {0,0,0};

      TH1D* wtGENDenDist_corr[3] = {0,0,0};
      TH1D* wtGENNumDist_corr[3] = {0,0,0};
      TH1D* wtDenDist_corr[3]    = {0,0,0};
      TH1D* wtRECODist_corr[3]   = {0,0,0};
      
      RooCmdArg bins = RooFit::Binning(distBins,-1,1);
      
      for (int ivar=0; ivar < 3; ivar++){
        if (ivar == 2) bins = RooFit::Binning(distBins, -TMath::Pi(), TMath::Pi());
        if(doCT){
          ctGENDenDist[ivar] = (TH1D*) data_genDen->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str()) 
            -> createHistogram( Form("dist%s_ctGenDen_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          ctGENNumDist[ivar] = (TH1D*) data_genNum->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    -> createHistogram( Form("dist%s_ctGenNum_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          ctDenDist[ivar]    = (TH1D*) data_den   ->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_ctDen_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()),  *vars_vec[ivar], bins);
   	  ctRECODist[ivar]   = (TH1D*) data_ctRECO->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_ctRECO_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);

          // correlated dataset
          ctGENDenDist_corr[ivar] = (TH1D*) data_genDen_corr->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str()) 
            -> createHistogram( Form("dist%s_ctGenDen_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          ctGENNumDist_corr[ivar] = (TH1D*) data_genNum_corr->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    -> createHistogram( Form("dist%s_ctGenNum_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          ctDenDist_corr[ivar]    = (TH1D*) data_den_corr   ->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_ctDen_corr_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()),  *vars_vec[ivar], bins);
   	  ctRECODist_corr[ivar]   = (TH1D*) data_ctRECO_corr->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_ctRECO_corr_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()),  *vars_vec[ivar], bins);

          // Scale genDen slices (currently with FSR-veto applied)
          // to match the total statistics of genDen (needed to obtain a correct absolute value in efficiency histograms and comparable with the efficiency function)
          ctGENDenDist[ivar]     ->Scale(normVal);
          ctGENDenDist_corr[ivar]->Scale(normVal_corr);
        }

        if (doWT) {
          wtGENDenDist[ivar] = (TH1D*) data_genDen->reduce(RooArgSet(vars[ivar]), wtCuts[ivar].c_str()) 
            -> createHistogram( Form("dist%s_wtGenDen_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          wtGENNumDist[ivar] = (TH1D*) data_genNum->reduce(RooArgSet(vars[ivar]), wtCuts[ivar].c_str())
	    -> createHistogram( Form("dist%s_wtGenNum_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          wtDenDist[ivar]    = (TH1D*) data_den   ->reduce(RooArgSet(vars[ivar]), wtCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_wtDen_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()),  *vars_vec[ivar], bins);
          // wt reco distribution uses cut with non-inverted sign
   	  wtRECODist[ivar]   = (TH1D*) data_wtRECO->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_wtRECO_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);

          // correlated dataset
          wtGENDenDist_corr[ivar] = (TH1D*) data_genDen_corr->reduce(RooArgSet(vars[ivar]), wtCuts[ivar].c_str()) 
            -> createHistogram( Form("dist%s_wtGenDen_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          wtGENNumDist_corr[ivar] = (TH1D*) data_genNum_corr->reduce(RooArgSet(vars[ivar]), wtCuts[ivar].c_str())
	    -> createHistogram( Form("dist%s_wtGenNum_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);
          wtDenDist_corr[ivar]    = (TH1D*) data_den_corr   ->reduce(RooArgSet(vars[ivar]), wtCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_wtDen_corr_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()),  *vars_vec[ivar], bins);
          // wt reco distribution uses cut with non-inverted sign
   	  wtRECODist_corr[ivar]   = (TH1D*) data_wtRECO_corr->reduce(RooArgSet(vars[ivar]), ctCuts[ivar].c_str())
	    ->createHistogram( Form("dist%s_wtRECO_corr_%i_%i_%s"   ,varCoord[ivar],i,j,shortString.c_str()), *vars_vec[ivar], bins);

          // Scale genDen slices (currently with FSR-veto applied)
          wtGENDenDist[ivar]     ->Scale(normVal);
          wtGENDenDist_corr[ivar]->Scale(normVal_corr);
        }
      }
      
      TH1D* effCDist[3] = {0,0,0};
      TH1D* effWDist[3] = {0,0,0};
      TH1D* effCDist_corr[3] = {0,0,0};
      TH1D* effWDist_corr[3] = {0,0,0};
      
      // composing binned efficiencies from sliced distributions
      for (int ivar=0; ivar< 3; ivar++){
        TH1D* ctDist = (TH1D*)ctGENNumDist[ivar]->Clone(Form("ctDist%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        ctDist->Divide(ctGENDenDist[ivar]);
        ctDist->Divide(ctDenDist[ivar]);
        effCDist[ivar] = (TH1D*)ctDist->Clone(Form("effCDist%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        if (doCT) effCDist[ivar]->Multiply(ctRECODist[ivar]);
      
        TH1D* wtDist = (TH1D*)wtGENNumDist[ivar]->Clone(Form("wtDist%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        wtDist->Divide(wtGENDenDist[ivar]);
        wtDist->Divide(wtDenDist[ivar]);
        effWDist[ivar] = Invert1Dhist(wtDist,Form("effWDist%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        if (doWT) effWDist[ivar]->Multiply(wtRECODist[ivar]);


        TH1D* ctDist_corr = (TH1D*)ctGENNumDist_corr[ivar]->Clone(Form("ctDist%s_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        ctDist_corr->Divide(ctGENDenDist_corr[ivar]);
        ctDist_corr->Divide(ctDenDist_corr[ivar]);
        effCDist_corr[ivar] = (TH1D*)ctDist_corr->Clone(Form("effCDist%s_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        if (doCT) effCDist_corr[ivar]->Multiply(ctRECODist_corr[ivar]);

        TH1D* wtDist_corr = (TH1D*)wtGENNumDist_corr[ivar]->Clone(Form("wtDist%s_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        wtDist_corr->Divide(wtGENDenDist_corr[ivar]);
        wtDist_corr->Divide(wtDenDist_corr[ivar]);
        effWDist_corr[ivar] = Invert1Dhist(wtDist_corr,Form("effWDist%s_corr_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()));
        if (doWT) effWDist_corr[ivar]->Multiply(wtRECODist_corr[ivar]);
      }

      // set histo properties
      if (doCT) {
	// set titles (to be improved)
	effCDist[0]->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effCDist[1]->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effCDist[2]->SetTitle( Form("Correct-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      }
      if (doWT) {
	effWDist[0]->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effWDist[1]->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB*TMath::Pi()) );
	effWDist[2]->SetTitle( Form("Wrong-tag efficiency (q2-bin %i, %s) slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",q2Bin,(parity==0?"even":"odd"),centA,centB) );
      }

      for (int ivar = 0; ivar < 3; ivar++){
        if (doCT) {
  	  effCDist[ivar]->SetMinimum(0);
          effCDist_corr[ivar]->SetLineColor(418);
  	  effCSlice[ivar].push_back( effCDist[ivar] );
  	  effCSlice_corr[ivar].push_back( effCDist_corr[ivar] );
        }
        if (doWT) {
  	  effWDist[ivar]->SetMinimum(0);
  	  effWDist_corr[ivar]->SetLineColor(418);
  	  effWSlice[ivar].push_back( effWDist[ivar] );
  	  effWSlice_corr[ivar].push_back( effWDist_corr[ivar] );
  	}
      }


      // producing 1D slices of efficiency description 
      for (int ivar = 0; ivar < 3; ivar++){
        if (ivar==0){ 
          ctL->setVal(centA);
          phi->setVal(centB*TMath::Pi());
        } else if (ivar==1){ 
          ctK->setVal(centA);
          phi->setVal(centB*TMath::Pi());
        } else if (ivar==2){ 
          ctK->setVal(centA); 
          ctL->setVal(centB);
        }
        if (doCT) {
  	  fsC[ivar].push_back( vars_vec[ivar]->frame( Name( Form("fs%sC_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()) ) ) );
  	  effC->plotOn( fsC[ivar].back(), LineColor(kRed), Name(Form("effC%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()))) ;
        }
        if (doWT) {
  	  fsW[ivar].push_back( vars_vec[ivar]->frame( Name( Form("fs%sW_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()) ) ) );
  	  effW->plotOn( fsW[ivar].back(), LineColor(kRed), Name(Form("effW%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()))) ;
        }
      }

      // if (i+j==0) {
      // 	leg->AddEntry(effHistsX.back(),"Binned efficiency" ,"lep");
      // 	leg->AddEntry(xframes.back()->findObject(Form("effx_%i_%i",i,j)),"KDE efficiency","l");
      // }
      // leg->Draw("same");

      for (int ivar=0; ivar< 3; ivar++){
        if (doCT) {
	  csC[ivar][confIndex]->cd(5*j+i+1);
	  gPad->SetLeftMargin(0.18);
	  effCSlice[ivar].back()->GetYaxis()->SetTitleOffset(1.7);
          effCSlice[ivar].back()->Draw();
          fsC[ivar].back()->Draw("same");
          effCSlice[ivar].back()->Draw("same");
          effCSlice_corr[ivar].back()->Draw("same");

          // checking maximum values
          if ( maxEffC[ivar]<effCSlice[ivar].back()     ->GetMaximum() ) maxEffC[ivar] = effCSlice[ivar].back()->GetMaximum();
          if ( maxEffC[ivar]<effCSlice_corr[ivar].back()->GetMaximum() ) maxEffC[ivar] = effCSlice_corr[ivar].back()->GetMaximum();
        }
        if (doWT) {
	  csW[ivar][confIndex]->cd(5*j+i+1);
	  gPad->SetLeftMargin(0.18);
	  effWSlice[ivar].back()->GetYaxis()->SetTitleOffset(1.7);
          effWSlice[ivar].back()->Draw();
          fsW[ivar].back()->Draw("same");
          effWSlice[ivar].back()->Draw("same");
          effWSlice_corr[ivar].back()->Draw("same");

          // checking maximum values
          if ( maxEffW[ivar]<effWSlice[ivar].back()     ->GetMaximum() ) maxEffW[ivar] = effWSlice[ivar].back()->GetMaximum();
          if ( maxEffW[ivar]<effWSlice_corr[ivar].back()->GetMaximum() ) maxEffW[ivar] = effWSlice_corr[ivar].back()->GetMaximum();
        }
      }
  } // end of slicing

  //set  uniform y-axis ranges and save
  for (int ivar=0; ivar< 3; ivar++){
    if (doCT) {
      for (vector<TH1D*>::iterator hist = effCSlice[ivar].begin(); hist != effCSlice[ivar].end(); ++hist) 
        (*hist)->SetMaximum(maxEffC[ivar]*1.1);
      csC[ivar][confIndex]->SaveAs( (confString+Form("_eff-ct_%sslices_comp_dp%i_%i.pdf",varNames[ivar] ,(int)(border*200), year)).c_str() );
    }
    if (doWT) {
      for (vector<TH1D*>::iterator hist = effWSlice[ivar].begin(); hist != effWSlice[ivar].end(); ++hist) 
        (*hist)->SetMaximum(maxEffW[ivar]*1.1);
      csW[ivar][confIndex]->SaveAs( (confString+Form("_eff-wt_%sslices_comp_dp%i_%i.pdf",varNames[ivar] ,(int)(border*200), year)).c_str() );
    }
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
    cp1C[confIndex]->SaveAs( (confString+Form("_eff-ct_1DProj_%i.pdf",year)).c_str() );
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
    cp1W[confIndex]->SaveAs( (confString+Form("_eff-wt_1DProj_%i.pdf",year)).c_str() );
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
    cp2C[confIndex]->SaveAs( (confString+Form("_eff-ct_2DProj_%i.pdf",year)).c_str() );
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
    cp2W[confIndex]->SaveAs( (confString+Form("_eff-wt_2DProj_%i.pdf", year)).c_str() );
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
    cctC[confIndex]->SaveAs( (confString+Form("_eff-ct_ClosureTest_%i.pdf",year)).c_str() );
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
    cctW[confIndex]->SaveAs( (confString+Form("_eff-wt_ClosureTest_%i.pdf",year)).c_str() );
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
  int year   = 2016;

  if ( argc >= 2 ) q2Bin  = atoi(argv[1]);
  if ( argc >= 3 ) parity = atoi(argv[2]);

  if ( q2Bin  < -1 || q2Bin  >= nBins ) return 1;
  if ( parity < -1 || parity > 1      ) return 1;

  bool doClosure = true;
  if ( argc >= 4 && atoi(argv[3]) == 0 ) doClosure = false;
  if ( argc >= 5 ) year = atoi(argv[4]);

  if ( q2Bin > -1 ) {
    if ( parity > -1 ) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<(parity==1?" - odd":" - even")<<" events"<<endl;
      plotEffBin( q2Bin, parity, doClosure, year );
    } else {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<" - both event parities"<<endl;
      plotEffBin( q2Bin, 0, doClosure, year );
      plotEffBin( q2Bin, 1, doClosure, year );
    }
  } else {
    cout<<"Plotting efficiency for all q2 bins - "<<(parity==1?"odd events":(parity==0?"even events":"both event parities"))<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if ( parity > -1 ) plotEffBin( q2Bin, parity, doClosure, year );
      else {
	plotEffBin( q2Bin, 0, doClosure, year );
	plotEffBin( q2Bin, 1, doClosure, year );
      }
    }
  }

  return 0;

}
