using namespace RooFit;
using namespace std;

/*
code to plot results from the Simultaneous fit of multiple years and compare them to
- GEN level results
- results of the fit to the single year datasets
  (should be provided)
*/


static const int nBins = 8; // was 8
// float binBorders [nBins+1] = { 1, 2};
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
// static const int nBins = 9;
// float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

static const int nPars = 8;
string ParName  [nPars] = { "Fl", "P1", "P2", "P3", "P4p", "P5p", "P6p", "P8p" };
string ParTitle [nPars] = { "F_{L}", "P_{1}", "P_{2}", "P_{3}", "P'_{4}", "P'_{5}", "P'_{6}", "P'_{8}" };
// double ParMin [nPars] = {0.4, -0.1, -0.2,-0.2, -0.5, -0.1, -0.1, -0.5};
// double ParMax [nPars] = {0.6,  0.1,  0.2, 0.2,  0.5,  0.1,  0.1,  0.5};
double ParMin [nPars] = {0,-1,-0.5,-0.5,-sqrt(2),-sqrt(2),-sqrt(2),-sqrt(2)};
double ParMax [nPars] = {1,1,0.5,0.5,sqrt(2),sqrt(2),sqrt(2),sqrt(2)};
float diffMax = 0.0999;

EColor colorlist[] = {kBlue, kRed, kOrange, kGreen, kMagenta};

TCanvas* c[nPars];


void fillVectors(RooFitResult* fitResult, int ParIndx, double Res[][nBins],  double ErrH[][nBins],  double ErrL[][nBins], 
                                                       double Diff[][nBins], double DErrH[][nBins], double DErrL[][nBins], 
                                                       double genRes[nBins], double genErrH[nBins], double genErrL[nBins], 
                                                       int year, int ibin){
  RooRealVar* Par = (RooRealVar*)fitResult->floatParsFinal().find(ParName[ParIndx].c_str());
  Res[year][ibin]  = Par->getValV();
  ErrH[year][ibin] = Par->getErrorHi();
  ErrL[year][ibin] = -1*Par->getErrorLo();
  Diff[year][ibin]  = Res[year][ibin]-genRes[ibin];
  DErrH[year][ibin] = sqrt(ErrH[year][ibin]*ErrH[year][ibin]+genErrL[ibin]*genErrL[ibin]);
  DErrL[year][ibin] = sqrt(ErrL[year][ibin]*ErrL[year][ibin]+genErrH[ibin]*genErrH[ibin]);
}

void fillVectors(RooFitResult* fitResult, int ParIndx, double Res[nBins], double ErrH[nBins], double ErrL[nBins], int ibin){
  RooRealVar* Par = (RooRealVar*)fitResult->floatParsFinal().find(ParName[ParIndx].c_str());
  Res[ibin]  = Par->getValV();
  ErrH[ibin] = Par->getErrorHi();
  ErrL[ibin] = -1*Par->getErrorLo();
}

void setGraphProperties(TGraphAsymmErrors* Gr, int icolor, string year, int ParIndx, string type){

  int addc = 0;
  if (type.find("WT") != std::string::npos) addc = 1;
  if (type.find("CT") != std::string::npos) addc = 2;

  Gr -> SetLineColor(colorlist[icolor] + addc);
  Gr -> SetMarkerColor(colorlist[icolor] + addc);
  Gr -> SetName(Form("Gr%s%i_%s" , type.c_str(), ParIndx, year.c_str()) );
  Gr -> SetLineWidth(2);
}

void setAxisProperties(TH1F* auxE2){
  auxE2->SetStats(kFALSE);
  auxE2->SetLineColor(1);
  auxE2->GetXaxis()->SetTitle("q^{2} (GeV^{2})");
  auxE2->GetXaxis()->SetTitleSize(0.12);
  auxE2->GetXaxis()->SetTitleOffset(0.95);
  auxE2->GetXaxis()->SetLabelSize( 0.10);
  auxE2->GetXaxis()->SetTickLength(0.1);
  auxE2->GetYaxis()->SetTitle("(RECO - GEN)");
  auxE2->GetYaxis()->SetTitleSize(0.10);
  auxE2->GetYaxis()->SetTitleOffset(0.45);
  auxE2->GetYaxis()->SetRangeUser(-1*diffMax,diffMax);
  auxE2->GetYaxis()->SetLabelSize(0.09);
  auxE2->GetYaxis()->SetNdivisions(505);
}

void plotFitResultsBin(int parity, int ParIndx, bool plotCT, bool plotWT, bool plotRECO, std::vector<string> years)
{


  double q2Val [nBins];
  double q2ErrH [nBins];
  double q2ErrL [nBins];

  for (int i=0; i<nBins; ++i) {
    q2Val [i] = 0.5 * (binBorders[i+1]+binBorders[i]);
    q2ErrH[i] = 0.5 * (binBorders[i+1]-binBorders[i]);
    q2ErrL[i] = 0.5 * (binBorders[i+1]-binBorders[i]);
  }

  double genRes [nBins];
  double genErrH [nBins];
  double genErrL [nBins];
  
  // n_years = number of datasets considered + 1*simultaneous dataset
  string combDataString = "";
  if (years.size() > 1){
    for (unsigned int iy = 0; iy < years.size() ; iy++) {
      combDataString += years[iy];    
    }
    years.push_back(combDataString.c_str());
  }
  unsigned int n_years = years.size();

  double ctRes  [n_years][nBins];
  double wtRes  [n_years][nBins];
  double Res    [n_years][nBins];

  double ctDiff  [n_years][nBins];
  double wtDiff  [n_years][nBins];
  double Diff    [n_years][nBins];

  double ctErrH  [n_years][nBins];
  double wtErrH  [n_years][nBins];
  double ErrH    [n_years][nBins];
  double ctErrL  [n_years][nBins];
  double wtErrL  [n_years][nBins];
  double ErrL    [n_years][nBins];

  double ctDiffErrH  [n_years][nBins];
  double wtDiffErrH  [n_years][nBins];
  double DiffErrH    [n_years][nBins];
  double ctDiffErrL  [n_years][nBins];
  double wtDiffErrL  [n_years][nBins];
  double DiffErrL    [n_years][nBins];
  
  std::vector<TGraphAsymmErrors*> GrCT, GrWT, Gr, GrDiffCT, GrDiffWT, GrDiff;

  TFile* finGen = TFile::Open("../eff-KDE/fitResults/fitResult_genMC.root");
  if ( !finGen || finGen->IsZombie() ) {
    cout<<"Missing gen file: fitResults/fitResult_genMC.root"<<endl;
    return;
  }

  // first fill for gen results
  for (int i=0; i<nBins; ++i) {
    RooFitResult* fitResultGen = (RooFitResult*)finGen->Get(Form("fitResult_b%ip%i",i,parity));
    if (fitResultGen && !fitResultGen->IsZombie() && fitResultGen->status()==0) {
      fillVectors(fitResultGen, ParIndx, genRes, genErrH, genErrL, i);
    }
  }
  finGen->Close();
  cout << genRes[0] << endl;
  
  string year = ""; 
  for (unsigned int iy = 0; iy < n_years; iy++) {
    year.clear(); year.assign(years[iy]);

    std::fill(std::begin(ctRes[iy]), std::end(ctRes[iy]), 0);
    std::fill(std::begin(wtRes[iy]), std::end(wtRes[iy]), 0);
    std::fill(std::begin(Res[iy]),   std::end(Res[iy]),   0);
//     genRes[i] = ctRes1[i] = ctRes2[i] = ctRes[i] = wtRes[i]  = ctDiff1[i] = ctDiff2[i] = ctDiff[i] =wtDiff[i] = Res[i] = Diff[i] = -2;
//     genErrH[i] = genErrL[i] = ctErrH1[i] = ctErrH2[i] = ctErrH[i] = ctErrL1[i] = ctErrL2[i] = ctErrL[i] = wtErrH[i] = wtErrL[i] = 0;
//     ctDiffErrH1[i] =  ctDiffErrH2[i] = ctDiffErrH[i] = ctDiffErrL1[i] = ctDiffErrL2[i] = ctDiffErrL[i] = wtDiffErrH[i] = wtDiffErrL[i] = 0;
//     ErrH[i] = ErrL[i] = DiffErrH[i] = DiffErrL[i] = 0;

    TFile* finReco;
    if (plotCT || plotWT) finReco = TFile::Open(("simFitResults/simFitResult_recoMC_singleComponent" + year + ".root").c_str());
    TFile* finFullReco;
    if (plotRECO) finFullReco = TFile::Open(("fitResults/fitResult_recoMC_fullAngular" + year + ".root").c_str());

    for (int i=0; i<nBins; ++i) {
      if (i==4 || i==6 || i==8) continue;

      // fill for CT events
      if ( plotCT && finReco && !finReco->IsZombie() ) {
         
	RooFitResult* fitResultCT = (RooFitResult*)finReco->Get(Form("simFitResult_b%ip%it1",i,parity));
	if (fitResultCT && !fitResultCT->IsZombie() && fitResultCT->status()==0) {
          fillVectors(fitResultCT, ParIndx, ctRes, ctErrH, ctErrL, ctDiff, ctDiffErrH, ctDiffErrL, genRes, genErrH, genErrL,  iy, i);
	}
      }    

      // fill for WT events
      if ( plotWT && finReco && !finReco->IsZombie() ) {
	RooFitResult* fitResultWT = (RooFitResult*)finReco->Get(Form("simFitResult_b%ip%it0",i,parity));
	if (fitResultWT && !fitResultWT->IsZombie() && fitResultWT->status()==0) {
          fillVectors(fitResultWT, ParIndx, wtRes, wtErrH, wtErrL, wtDiff, wtDiffErrH, wtDiffErrL, genRes, genErrH, genErrL,  iy, i);
	}
      }    

      // fill for full angular fit
      if ( plotRECO && finFullReco && !finFullReco->IsZombie() ) {
        RooFitResult* fitResult = (RooFitResult*)finFullReco->Get(Form("FitResult_b%ip%i",i,parity));
        if (fitResult && !fitResult->IsZombie() && fitResult->status()==0) {
          fillVectors(fitResult, ParIndx, Res, ErrH, ErrL, Diff, DiffErrH, DiffErrL, genRes, genErrH, genErrL,  iy, i);
        }
      }

    }
    
    // close files 
    if ( (plotCT || plotWT) && finReco ) finReco->Close();
    if ( plotRECO && finFullReco ) finFullReco->Close();
    
    // create TGraphs with fit results
    GrCT.push_back( new TGraphAsymmErrors(nBins, q2Val, ctRes[iy], q2ErrH, q2ErrL, ctErrL[iy], ctErrH[iy]) );
    setGraphProperties(GrCT[iy], iy, year, ParIndx, "CT");

    GrWT.push_back( new TGraphAsymmErrors(nBins, q2Val, wtRes[iy], q2ErrH, q2ErrL, wtErrL[iy], wtErrH[iy]) );
    setGraphProperties(GrWT[iy], iy, year, ParIndx, "WT");

    Gr  .push_back( new TGraphAsymmErrors(nBins, q2Val, Res[iy], q2ErrH, q2ErrL, ErrL[iy], ErrH[iy]) );
    setGraphProperties(Gr[iy], iy, year, ParIndx, "Full");
    
    // create TGraphs with difference to GEN Results
    GrDiffCT.push_back( new TGraphAsymmErrors(nBins, q2Val, ctDiff[iy], q2ErrH, q2ErrL, ctDiffErrL[iy], ctDiffErrH[iy]) );
    setGraphProperties(GrDiffCT[iy], iy, year, ParIndx, "DiffCT");

    GrDiffWT.push_back( new TGraphAsymmErrors(nBins, q2Val, wtDiff[iy], q2ErrH, q2ErrL, wtDiffErrL[iy], wtDiffErrH[iy]) );
    setGraphProperties(GrDiffWT[iy], iy, year, ParIndx, "DiffWT");

    GrDiff.push_back( new TGraphAsymmErrors(nBins, q2Val, Diff[iy], q2ErrH, q2ErrL, DiffErrL[iy], DiffErrH[iy]) );
    setGraphProperties(GrDiff[iy], iy, year, ParIndx, "Diff");
    
  }


  TGraphAsymmErrors* GrGen  = new TGraphAsymmErrors(nBins,q2Val,genRes ,q2ErrH,q2ErrL,genErrL,genErrH);
  GrGen->SetName(Form("GrGen%i",ParIndx));

  GrGen->SetTitle( ("Fit result comparison of "+ParTitle[ParIndx]+" parameter").c_str() );
  GrGen->GetYaxis()->SetTitle( ParTitle[ParIndx].c_str() );
  GrGen->GetYaxis()->SetTitleSize(20);
  GrGen->GetYaxis()->SetTitleFont(43);
  GrGen->GetYaxis()->SetTitleOffset(1.55);

  GrGen->GetYaxis()->SetRangeUser(ParMin[ParIndx],ParMax[ParIndx]);
  GrGen->SetLineWidth(2);

  // Grey bands for resonant regions
  double ResX [2] = {0.5*(binBorders[5]+binBorders[4]),0.5*(binBorders[7]+binBorders[6])};
  double ResXe[2] = {0.5*(binBorders[5]-binBorders[4]),0.5*(binBorders[7]-binBorders[6])};
  double ResY [2] = {0.5*(ParMax[ParIndx]+ParMin[ParIndx])-0.002*(ParMax[ParIndx]-ParMin[ParIndx]),
		     0.5*(ParMax[ParIndx]+ParMin[ParIndx])-0.002*(ParMax[ParIndx]-ParMin[ParIndx])};
  double ResYe[2] = {0.498*(ParMax[ParIndx]-ParMin[ParIndx]),0.498*(ParMax[ParIndx]-ParMin[ParIndx])};
  double ResD [2] = {-0.01*diffMax,-0.01*diffMax};
  double ResDe[2] = { 0.99*diffMax, 0.99*diffMax};
  TGraphErrors *resCover = new TGraphErrors(2,ResX,ResY,ResXe,ResYe);
  resCover->SetName(Form("resCover%i",ParIndx));
  resCover->SetFillColor(18);
  resCover->SetFillStyle(1001);
  TGraphErrors *resDiffCover = new TGraphErrors(2,ResX,ResD,ResXe,ResDe);
  resDiffCover->SetName(Form("resDiffCover%i",ParIndx));
  resDiffCover->SetFillColor(18);
  resDiffCover->SetFillStyle(1001);

  // Legend
  TLegend *leg;
  if (ParIndx==4 || ParIndx==5) leg = new TLegend(0.15,0.65,0.4,0.85);
  else if (ParIndx==2) leg = new TLegend(0.48,0.1,0.9,0.3);
  else leg = new TLegend(0.15,0.1,0.4,0.3);
  leg->SetName(Form("leg%i",ParIndx));
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.032);
  leg->AddEntry(GrGen,"Fit to generation-level events","lep");

  // Zero line
  TLine *line = new TLine(GrGen->GetXaxis()->GetXmin(),0,GrGen->GetXaxis()->GetXmax(),0);
  line->SetLineColor(14);
  line->SetLineStyle(7);

  c[ParIndx] = new TCanvas(Form("c%i",ParIndx),Form("c%i",ParIndx),800,800);
  c[ParIndx]->cd();
  TPad *pad1 = new TPad(Form("pad1_%i",ParIndx), "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  
     
  GrGen->Draw("AP");
  for (unsigned int iy = 0; iy < years.size() ; iy++) {
    if (plotCT)  {
      GrCT[iy] -> Draw("P");    
      leg->AddEntry(GrCT[iy] , Form("Sim Fit to correctly tagged events, %s", years[iy].c_str()), "lep");
    }  
    if (plotWT) {
      GrWT[iy] -> Draw("P");    
      leg->AddEntry(GrWT[iy] , Form("Sim Fit to wrongly tagged events, %s", years[iy].c_str()), "lep");
    }  
    if (plotRECO){
      Gr[iy]   -> Draw("P");    
      leg->AddEntry(Gr[iy],    Form("Fit to reconstructed events, %s", years[iy].c_str()), "lep");
    }  
  }
  
  resCover->Draw("e2");
  leg->Draw();

  GrGen->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( 0, 
                             ParMin[ParIndx],
                             0, 
                             ParMax[ParIndx], 
                             ParMin[ParIndx], 
                             ParMax[ParIndx], 
                             510, 
                             "");
  axis->SetName(Form("leg%i",ParIndx));
  axis->SetLabelFont(43);
  axis->SetLabelSize(15);
  axis->Draw();


  // plot difference wrt GEN results
  c[ParIndx]->cd();
  TPad *pad2 = new TPad(Form("pad2_%i",ParIndx), "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->Draw();
  pad2->cd();

  // first create axis
  TH1F* auxE2 = new TH1F("auxE2", "", nBins, GrGen->GetXaxis()->GetXmin(), GrGen->GetXaxis()->GetXmax());
  setAxisProperties(auxE2);
  auxE2->Draw();

  for (unsigned int iy = 0; iy < years.size() ; iy++) {
    if (plotCT)  GrDiffCT[iy] -> Draw("P");    
    if (plotWT)  GrDiffWT[iy] -> Draw("P");    
    if (plotRECO)  GrDiff[iy]   -> Draw("P");    
  }
  line->Draw();
  resDiffCover->Draw("e2");

  string confString = "plotSimFit_d/fitResult_";
  if (plotCT) confString = confString + "ctRes_";
  if (plotWT) confString = confString + "wtRes_";
  if (plotRECO) confString = confString + "recoRes_";
  c[ParIndx]->SaveAs( (confString+ParName[ParIndx]+".pdf").c_str() );

}

void plotSimFitResults(int parity, int ParIndx = -1, bool plotCT = true, bool plotWT = false, bool plotRECO = false)
{

  if ( parity<0 || parity>1 ) return;

  if ( ParIndx<-1 || ParIndx>=nPars ) return;

  if ( ParIndx==-1 ) cout<<"Running all the parameters"<<endl;

  std::vector<string> years;
  years.push_back("2016");
  years.push_back("2017");

  if ( ParIndx==-1 )
    for (ParIndx=0; ParIndx<nPars; ++ParIndx)
      plotFitResultsBin(parity, ParIndx, plotCT, plotWT, plotRECO, years);
  else
    plotFitResultsBin(parity, ParIndx, plotCT, plotWT, plotRECO, years);

}
