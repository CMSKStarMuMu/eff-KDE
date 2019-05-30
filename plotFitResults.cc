using namespace RooFit;
using namespace std;

static const int nBins = 8;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
// static const int nBins = 9;
// float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

static const int nPars = 8;
string ParName  [nPars] = { "Fl", "P1", "P2", "P3", "P4p", "P5p", "P6p", "P8p" };
string ParTitle [nPars] = { "F_{L}", "P_{1}", "P_{2}", "P_{3}", "P'_{4}", "P'_{5}", "P'_{6}", "P'_{8}" };
double ParMin [nPars] = {0,-1,-0.5,-0.5,-sqrt(2),-sqrt(2),-sqrt(2),-sqrt(2)};
double ParMax [nPars] = {1,1,0.5,0.5,sqrt(2),sqrt(2),sqrt(2),sqrt(2)};

TCanvas* c[nPars];

void plotFitResultsBin(int parity, int ParIndx, bool plotCT, bool plotWT, bool plotRECO)
{

  float diffMax = 0.0999;

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
  double ctRes [nBins];
  double wtRes [nBins];
  double Res [nBins];
  double ctDiff [nBins];
  double wtDiff [nBins];
  double Diff [nBins];
  double ctErrH [nBins];
  double wtErrH [nBins];
  double ErrH [nBins];
  double ctErrL [nBins];
  double wtErrL [nBins];
  double ErrL [nBins];
  double ctDiffErrH [nBins];
  double wtDiffErrH [nBins];
  double DiffErrH [nBins];
  double ctDiffErrL [nBins];
  double wtDiffErrL [nBins];
  double DiffErrL [nBins];

  TFile* finGen = TFile::Open("fitResult_genMC.root");
  TFile* finReco;
  TFile* finFullReco;
  if (plotCT || plotWT) finReco = TFile::Open("fitResult_recoMC_singleComponent.root");
  if (plotRECO) finFullReco = TFile::Open("fitResult_recoMC_fullAngular.root");

  for (int i=0; i<nBins; ++i) {

    if (i==4 || i==6 || i==8) {
      genRes[i]=ctRes[i]=wtRes[i]=ctDiff[i]=wtDiff[i]=Res[i]=Diff[i]=-2;
      genErrH[i]=genErrL[i]=ctErrH[i]=ctErrL[i]=wtErrH[i]=wtErrL[i]=0;
      ctDiffErrH[i]=ctDiffErrL[i]=wtDiffErrH[i]=wtDiffErrL[i]=0;
      ErrH[i]=ErrL[i]=DiffErrH[i]=DiffErrL[i]=0;
      continue;
    }

    RooFitResult* fitResultGen = (RooFitResult*)finGen->Get(Form("fitResult_b%ip%i",i,parity));
    RooRealVar* ParGen = (RooRealVar*)fitResultGen->floatParsFinal().find(ParName[ParIndx].c_str());
    genRes[i] = ParGen->getValV();
    genErrH[i] = ParGen->getErrorHi();
    genErrL[i] = -1*ParGen->getErrorLo();

    if (plotCT) {
      RooFitResult* fitResultCT = (RooFitResult*)finReco->Get(Form("fitResult_b%ip%it1",i,parity));
      RooRealVar* ParCT  = (RooRealVar*)fitResultCT ->floatParsFinal().find(ParName[ParIndx].c_str());
      ctRes[i] = ParCT->getValV();
      ctErrH[i] = ParCT->getErrorHi();
      ctErrL[i] = -1*ParCT->getErrorLo();
      ctDiff[i] = ctRes[i]-genRes[i];
      ctDiffErrH[i] = sqrt(ctErrH[i]*ctErrH[i]+genErrL[i]*genErrL[i]);
      ctDiffErrL[i] = sqrt(ctErrL[i]*ctErrL[i]+genErrH[i]*genErrH[i]);
    } else {
      ctRes[i] = ctDiff[i] = -2;
      ctErrH[i] = ctErrL[i] = ctDiffErrH[i] = ctDiffErrL[i] = 0;
    }

    if (plotWT) {
      RooFitResult* fitResultWT = (RooFitResult*)finReco->Get(Form("fitResult_b%ip%it0",i,parity));
      RooRealVar* ParWT  = (RooRealVar*)fitResultWT ->floatParsFinal().find(ParName[ParIndx].c_str());
      wtRes[i] = ParWT->getValV();
      wtErrH[i] = ParWT->getErrorHi();
      wtErrL[i] = -1*ParWT->getErrorLo();
      wtDiff[i] = wtRes[i]-genRes[i];
      wtDiffErrH[i] = sqrt(wtErrH[i]*wtErrH[i]+genErrL[i]*genErrL[i]);
      wtDiffErrL[i] = sqrt(wtErrL[i]*wtErrL[i]+genErrH[i]*genErrH[i]);
    } else {
      wtRes[i] = wtDiff[i] = -2;
      wtErrH[i] = wtErrL[i] = wtDiffErrH[i] = wtDiffErrL[i] = 0;
    }

    if (plotRECO) {
      RooFitResult* fitResult = (RooFitResult*)finFullReco->Get(Form("fitResult_b%ip%i",i,parity));
      RooRealVar* Par    = (RooRealVar*)fitResult   ->floatParsFinal().find(ParName[ParIndx].c_str());
      Res[i] = Par->getValV();
      ErrH[i] = Par->getErrorHi();
      ErrL[i] = -1*Par->getErrorLo();
      Diff[i] = Res[i]-genRes[i];
      DiffErrH[i] = sqrt(ErrH[i]*ErrH[i]+genErrL[i]*genErrL[i]);
      DiffErrL[i] = sqrt(ErrL[i]*ErrL[i]+genErrH[i]*genErrH[i]);
    } else {
      Res[i] = Diff[i] = -2;
      ErrH[i] = ErrL[i] = DiffErrH[i] = DiffErrL[i] = 0;
    }

  }

  TGraphAsymmErrors* GrGen = new TGraphAsymmErrors(nBins,q2Val,genRes,q2ErrH,q2ErrL,genErrL,genErrH);
  TGraphAsymmErrors* GrCT  = new TGraphAsymmErrors(nBins,q2Val,ctRes ,q2ErrH,q2ErrL,ctErrL,ctErrH);
  TGraphAsymmErrors* GrWT  = new TGraphAsymmErrors(nBins,q2Val,wtRes ,q2ErrH,q2ErrL,wtErrL,wtErrH);
  TGraphAsymmErrors* Gr    = new TGraphAsymmErrors(nBins,q2Val,Res   ,q2ErrH,q2ErrL,ErrL  ,ErrH  );
  TGraphAsymmErrors* GrDiffCT = new TGraphAsymmErrors(nBins,q2Val,ctDiff,q2ErrH,q2ErrL,ctDiffErrL,ctDiffErrH);
  TGraphAsymmErrors* GrDiffWT = new TGraphAsymmErrors(nBins,q2Val,wtDiff,q2ErrH,q2ErrL,wtDiffErrL,wtDiffErrH);
  TGraphAsymmErrors* GrDiff   = new TGraphAsymmErrors(nBins,q2Val,Diff  ,q2ErrH,q2ErrL,DiffErrL  ,DiffErrH  );
  GrGen->SetName(Form("GrGen%i",ParIndx));
  GrCT ->SetName(Form("GrCT%i" ,ParIndx));
  GrWT ->SetName(Form("GrWT%i" ,ParIndx));
  Gr   ->SetName(Form("Gr%i"   ,ParIndx));
  GrDiffCT->SetName(Form("GrDiffCT%i",ParIndx));
  GrDiffWT->SetName(Form("GrDiffWT%i",ParIndx));
  GrDiff  ->SetName(Form("GrDiff%i"  ,ParIndx));

  GrGen->SetTitle( ("Fit result comparison of "+ParTitle[ParIndx]+" parameter").c_str() );
  GrGen->GetYaxis()->SetTitle( ParTitle[ParIndx].c_str() );
  GrGen->GetYaxis()->SetTitleSize(20);
  GrGen->GetYaxis()->SetTitleFont(43);
  GrGen->GetYaxis()->SetTitleOffset(1.55);

  GrDiffWT->SetTitle( "" );
  GrDiffWT->GetYaxis()->SetTitle( "(RECO - GEN)" );
  GrDiffWT->GetYaxis()->SetNdivisions(505);
  GrDiffWT->GetYaxis()->SetTitleSize(15);
  GrDiffWT->GetYaxis()->SetTitleFont(43);
  GrDiffWT->GetYaxis()->SetTitleOffset(1.85);
  GrDiffWT->GetYaxis()->SetLabelFont(43);
  GrDiffWT->GetYaxis()->SetLabelSize(15);
  GrDiffWT->GetXaxis()->SetTitle( "q^{2} (GeV^{2})" );
  GrDiffWT->GetXaxis()->SetTitleSize(20);
  GrDiffWT->GetXaxis()->SetTitleFont(43);
  GrDiffWT->GetXaxis()->SetTitleOffset(3.);
  GrDiffWT->GetXaxis()->SetLabelFont(43);
  GrDiffWT->GetXaxis()->SetLabelSize(15);

  GrGen->GetYaxis()->SetRangeUser(ParMin[ParIndx],ParMax[ParIndx]);
  GrDiffWT->GetYaxis()->SetRangeUser(-1*diffMax,diffMax);

  Gr->SetLineColor(kGreen+2);
  GrDiff->SetLineColor(kGreen+2);
  GrCT->SetLineColor(kBlue);
  GrDiffCT->SetLineColor(kBlue);
  GrWT->SetLineColor(kRed+1);
  GrDiffWT->SetLineColor(kRed+1);
  GrGen->SetLineWidth(2);
  GrCT->SetLineWidth(2);
  GrWT->SetLineWidth(2);
  Gr->SetLineWidth(2);
  GrDiffCT->SetLineWidth(2);
  GrDiffWT->SetLineWidth(2);
  GrDiff->SetLineWidth(2);

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
  if (ParIndx==4 || ParIndx==5) leg = new TLegend(0.15,0.7,0.4,0.85);
  else if (ParIndx==2) leg = new TLegend(0.48,0.1,0.9,0.2);
  else leg = new TLegend(0.15,0.1,0.4,0.2);
  leg->SetName(Form("leg%i",ParIndx));
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.032);
  leg->AddEntry(GrGen,"Fit to generation-level events","lep");
  if (plotCT) leg->AddEntry(GrCT,"Fit to correctly tagged events","lep");
  if (plotWT) leg->AddEntry(GrWT,"Fit to wrongly tagged events","lep");
  if (plotRECO) leg->AddEntry(Gr,"Fit to reconstructed events","lep");

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
  GrWT->Draw("P");
  GrCT->Draw("P");
  Gr->Draw("P");
  resCover->Draw("e2");
  leg->Draw();

  GrGen->GetYaxis()->SetLabelSize(0.);
  TGaxis *axis = new TGaxis( 0, ((int)(ParMin[ParIndx]*10))/10.0+0.1, 0, ParMax[ParIndx], ((int)(ParMin[ParIndx]*10))/10.0+0.1, ParMax[ParIndx], 510, "");
  axis->SetName(Form("leg%i",ParIndx));
  axis->SetLabelFont(43);
  axis->SetLabelSize(15);
  axis->Draw();

  c[ParIndx]->cd();
  TPad *pad2 = new TPad(Form("pad2_%i",ParIndx), "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();

  GrDiffWT->Draw("AP");
  line->Draw();
  GrDiffWT->Draw("P");
  GrDiffCT->Draw("P");
  GrDiff->Draw("P");
  resDiffCover->Draw("e2");

  string confString = "fitResult_";
  if (plotCT) confString = confString + "ctRes_";
  if (plotWT) confString = confString + "wtRes_";
  if (plotRECO) confString = confString + "recoRes_";
  c[ParIndx]->SaveAs( (confString+ParName[ParIndx]+".pdf").c_str() );

}

void plotFitResults(int parity, int ParIndx = -1, bool plotCT = true, bool plotWT = true, bool plotRECO = true)
{

  if ( parity<0 || parity>1 ) return;

  if ( ParIndx<-1 || ParIndx>1 ) return;

  if ( ParIndx==-1 ) cout<<"Running all the parameters"<<endl;

  if ( ParIndx==-1 )
    for (ParIndx=0; ParIndx<nPars; ++ParIndx)
      plotFitResultsBin(parity, ParIndx, plotCT, plotWT, plotRECO);
  else
    plotFitResultsBin(parity, ParIndx, plotCT, plotWT, plotRECO);

}
