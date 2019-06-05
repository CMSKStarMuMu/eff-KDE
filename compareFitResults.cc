using namespace RooFit;
using namespace std;

static const int nBins = 8;

static const int nPars = 8;
string ParName  [nPars] = { "Fl", "P1", "P2", "P3", "P4p", "P5p", "P6p", "P8p" };
const char* ParTitle [nPars] = { "F_{L}", "P_{1}", "P_{2}", "P_{3}", "P'_{4}", "P'_{5}", "P'_{6}", "P'_{8}" };

static const int nAlt = 9;
int AltIndx [nAlt] = { 23324, 53324, 32324, 35324, 33524, 33334, 33344, 33323, 33325 };
int AltColor [nAlt] = { 633, 417, 879, 857, 839, 402, 887, 801, 921 };
const char* AltLabel [nAlt] = { "denGen SF: 0.3->0.2",
				"denGen SF: 0.3->0.5",
				"numGen SF: 0.3->0.2",
				"numGen SF: 0.3->0.5",
				"Den SF: 0.3->0.5",
				"ctRECO SF: 0.2->0.3",
				"ctRECO SF: 0.2->0.4",
				"wtRECO SF: 0.4->0.3",
				"wtRECO SF: 0.4->0.5" };

TCanvas* c[3];

void compareFitResults(int q2Bin, int parity, bool plotCT = true, bool plotWT = true, bool plotRECO = true)
{

  if ( q2Bin<0 || q2Bin>=nBins ) return;
  if ( parity<0 || parity>1 ) return;

  string shortString = Form("b%ip%i",q2Bin,parity);

  // maximum bias value plotted
  float diffMax = 0.1499;
  // float diffMax = 0.0799;

  // Histograms to plot axes
  TH1D * HistLab = new TH1D("HistLab","Fit bias;;(RECO - GEN)",nPars,0,nPars);

  // set of arreys for x values
  vector<double*> binC;

  // ser of arrays for y values and errors
  vector<double*> ctDiff (0);
  vector<double*> wtDiff (0);
  vector<double*> Diff (0);
  vector<double*> ctDiffErrH (0);
  vector<double*> wtDiffErrH (0);
  vector<double*> DiffErrH (0);
  vector<double*> ctDiffErrL (0);
  vector<double*> wtDiffErrL (0);
  vector<double*> DiffErrL (0);

  TFile* finGen = TFile::Open("fitResult_genMC.root");
  if ( !finGen || finGen->IsZombie() ) {
    cout<<"Missing gen file: fitResult_genMC.root"<<endl;
    return;
  }
  TFile* finReco[nAlt+1];
  TFile* finFullReco[nAlt+1];

  for (int iAlt=0; iAlt<=nAlt; ++iAlt) {

    // reference values
    double genRes[nPars];
    double genErrH[nPars];
    double genErrL[nPars];

    binC.push_back( new double[nPars] );
    double offset = ((iAlt+1)/2) * 0.5 / ((nAlt+3)/2);
    if (iAlt%2==0) offset = -1*offset;
    for (int ParIndx=0; ParIndx<nPars; ++ParIndx)
      binC.back()[ParIndx] = HistLab->GetBinCenter(ParIndx+1) + offset;

    ctDiff.push_back( new double[nPars] );
    wtDiff.push_back( new double[nPars] );
    Diff.push_back( new double[nPars] );
    ctDiffErrH.push_back( new double[nPars] );
    wtDiffErrH.push_back( new double[nPars] );
    DiffErrH.push_back( new double[nPars] );
    ctDiffErrL.push_back( new double[nPars] );
    wtDiffErrL.push_back( new double[nPars] );
    DiffErrL.push_back( new double[nPars] );

    // open files with RECO fit results
    if (plotCT || plotWT) {
      if (iAlt==0) finReco[iAlt] = TFile::Open("fitResult_recoMC_singleComponent.root");
      else finReco[iAlt] = TFile::Open(Form("fitResult_recoMC_singleComponent_alt%i.root",AltIndx[iAlt-1]));
    }
    if (plotRECO) {
      if (iAlt==0) finFullReco[iAlt] = TFile::Open("fitResult_recoMC_fullAngular.root");
      else finFullReco[iAlt] = TFile::Open(Form("fitResult_recoMC_fullAngular_alt%i.root",AltIndx[iAlt-1]));
    }

    // Fill reference values
    RooFitResult* fitResultGen = (RooFitResult*)finGen->Get(("fitResult_"+shortString).c_str());
    for (int ParIndx=0; ParIndx<nPars; ++ParIndx) {
      RooRealVar* ParGen = (RooRealVar*)fitResultGen->floatParsFinal().find(ParName[ParIndx].c_str());
      genRes [ParIndx] = ParGen->getValV();
      genErrH[ParIndx] = ParGen->getErrorHi();
      genErrL[ParIndx] = -1*ParGen->getErrorLo();
    }

    // Fill default values
    for (int ParIndx=0; ParIndx<nPars; ++ParIndx) {
      ctDiff.back()[ParIndx] = -2;
      ctDiffErrH.back()[ParIndx] = ctDiffErrL.back()[ParIndx] = 0;
      wtDiff.back()[ParIndx] = -2;
      wtDiffErrH.back()[ParIndx] = wtDiffErrL.back()[ParIndx] = 0;
      Diff.back()[ParIndx] = -2;
      DiffErrH.back()[ParIndx] = DiffErrL.back()[ParIndx] = 0;
    }

    if ( plotCT && finReco[iAlt] && !finReco[iAlt]->IsZombie() ) {
      // Fill values of correct-tag fit 
      RooFitResult* fitResultCT = (RooFitResult*)finReco[iAlt]->Get(("fitResult_"+shortString+"t1").c_str());
      if (fitResultCT && !fitResultCT->IsZombie()) for (int ParIndx=0; ParIndx<nPars; ++ParIndx) {
	RooRealVar* ParCT = (RooRealVar*)fitResultCT->floatParsFinal().find(ParName[ParIndx].c_str());
	ctDiff.back()[ParIndx] = ParCT->getValV()-genRes[ParIndx];
	ctDiffErrH.back()[ParIndx] = sqrt(ParCT->getErrorHi()*ParCT->getErrorHi()+genErrL[ParIndx]*genErrL[ParIndx]);
	ctDiffErrL.back()[ParIndx] = sqrt(ParCT->getErrorLo()*ParCT->getErrorLo()+genErrH[ParIndx]*genErrH[ParIndx]);
      }
    }

    if ( plotWT && finReco[iAlt] && !finReco[iAlt]->IsZombie() ) {
      // Fill values of wrong-tag fit 
      RooFitResult* fitResultWT = (RooFitResult*)finReco[iAlt]->Get(("fitResult_"+shortString+"t0").c_str());
      if (fitResultWT && !fitResultWT->IsZombie()) for (int ParIndx=0; ParIndx<nPars; ++ParIndx) {
	RooRealVar* ParWT = (RooRealVar*)fitResultWT->floatParsFinal().find(ParName[ParIndx].c_str());
	wtDiff.back()[ParIndx] = ParWT->getValV()-genRes[ParIndx];
	wtDiffErrH.back()[ParIndx] = sqrt(ParWT->getErrorHi()*ParWT->getErrorHi()+genErrL[ParIndx]*genErrL[ParIndx]);
	wtDiffErrL.back()[ParIndx] = sqrt(ParWT->getErrorLo()*ParWT->getErrorLo()+genErrH[ParIndx]*genErrH[ParIndx]);
      }
    }

    if ( plotRECO && finFullReco[iAlt] && !finFullReco[iAlt]->IsZombie() ) {
      // Fill values of RECO fit 
      RooFitResult* fitResult = (RooFitResult*)finFullReco[iAlt]->Get(("fitResult_"+shortString).c_str());
      if (fitResult && !fitResult->IsZombie()) for (int ParIndx=0; ParIndx<nPars; ++ParIndx) {
	RooRealVar* Par = (RooRealVar*)fitResult->floatParsFinal().find(ParName[ParIndx].c_str());
	Diff.back()[ParIndx] = Par->getValV()-genRes[ParIndx];
	DiffErrH.back()[ParIndx] = sqrt(Par->getErrorHi()*Par->getErrorHi()+genErrL[ParIndx]*genErrL[ParIndx]);
	DiffErrL.back()[ParIndx] = sqrt(Par->getErrorLo()*Par->getErrorLo()+genErrH[ParIndx]*genErrH[ParIndx]);
      }
    }

    // close files 
    if ( (plotCT || plotWT) && finReco[iAlt] ) finReco[iAlt]->Close();
    if ( plotRECO && finFullReco[iAlt] ) finFullReco[iAlt]->Close();

  }

  double binEL [nPars];
  double binEH [nPars];
  for (int i=0; i<nPars; ++i) {
    HistLab->GetXaxis()->SetBinLabel(i+1,ParTitle[i]);
    binEL[i] = binEH[i] = 0;
  }

  TGraphAsymmErrors* GrDiffCT[nAlt+1];
  TGraphAsymmErrors* GrDiffWT[nAlt+1];
  TGraphAsymmErrors* GrDiff[nAlt+1];

  for (int iAlt=0; iAlt<=nAlt; ++iAlt) {
    GrDiffCT[iAlt] = new TGraphAsymmErrors(nPars,binC[iAlt],ctDiff[iAlt],binEL,binEH,ctDiffErrL[iAlt],ctDiffErrH[iAlt]);
    GrDiffWT[iAlt] = new TGraphAsymmErrors(nPars,binC[iAlt],wtDiff[iAlt],binEL,binEH,wtDiffErrL[iAlt],wtDiffErrH[iAlt]);
    GrDiff  [iAlt] = new TGraphAsymmErrors(nPars,binC[iAlt],Diff  [iAlt],binEL,binEH,DiffErrL  [iAlt],DiffErrH  [iAlt]);
    GrDiffCT[iAlt]->SetName(Form("GrDiffCT%i",iAlt));
    GrDiffWT[iAlt]->SetName(Form("GrDiffWT%i",iAlt));
    GrDiff  [iAlt]->SetName(Form("GrDiff%i"  ,iAlt));
    GrDiffCT[iAlt]->SetMarkerStyle(20);
    GrDiffWT[iAlt]->SetMarkerStyle(20);
    GrDiff  [iAlt]->SetMarkerStyle(20);
    GrDiffCT[iAlt]->SetLineWidth(2);
    GrDiffWT[iAlt]->SetLineWidth(2);
    GrDiff  [iAlt]->SetLineWidth(2);
    if (iAlt>0) {
      GrDiffCT[iAlt]->SetLineColor(AltColor[iAlt-1]);
      GrDiffWT[iAlt]->SetLineColor(AltColor[iAlt-1]);
      GrDiff  [iAlt]->SetLineColor(AltColor[iAlt-1]);
      GrDiffCT[iAlt]->SetMarkerColor(AltColor[iAlt-1]);
      GrDiffWT[iAlt]->SetMarkerColor(AltColor[iAlt-1]);
      GrDiff  [iAlt]->SetMarkerColor(AltColor[iAlt-1]);
    }
  }

  // Setting of axes
  HistLab->SetMinimum(-1*diffMax);
  HistLab->SetMaximum(   diffMax);
  HistLab->SetStats(false);
  HistLab->GetYaxis()->SetTitleSize(20);
  HistLab->GetYaxis()->SetTitleFont(43);
  HistLab->GetYaxis()->SetTitleOffset(1.75);
  HistLab->GetYaxis()->SetLabelFont(43);
  HistLab->GetYaxis()->SetLabelSize(20);
  HistLab->GetXaxis()->SetLabelFont(43);
  HistLab->GetXaxis()->SetLabelSize(35);

  // Legend
  TLegend *leg = new TLegend(0.15,0.89-0.025*(nAlt+1),0.4,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.025);
  leg->AddEntry(GrDiff[0],"legacy","lep");
  for (int iAlt=1; iAlt<=nAlt; ++iAlt)
    leg->AddEntry(GrDiff[iAlt],AltLabel[iAlt-1],"lep");

  // Zero line
  TLine *line = new TLine(0,0,nPars,0);
  line->SetLineColor(14);
  line->SetLineStyle(7);

  // Title
  TLatex* TitleTex = new TLatex(0.5,0.91,"");
  TitleTex->SetTextFont(43);
  TitleTex->SetTextSize(32);
  TitleTex->SetTextAlign(20);
  TitleTex->SetNDC(true);
  string longString = Form(parity==0?"q2 bin %i (even eff)":"q2 bin %i (odd eff)",q2Bin);

  if (plotCT) {
    c[0] = new TCanvas("c0","c0",800,800);
    c[0]->cd();
    HistLab->Draw("AXIS");
    line->Draw();
    for (int iAlt=0; iAlt<=nAlt; ++iAlt) GrDiffCT[iAlt]->Draw("p");
    leg->Draw();
    TitleTex->DrawLatex(0.5,0.91,("Bias in fit to correct-tag events - "+longString).c_str());
    c[0]->SaveAs( ("fitResult_alternativeComp_"+shortString+"_ctRECO.pdf").c_str() );
  }

  if (plotWT) {
    c[1] = new TCanvas("c1","c1",800,800);
    HistLab->Draw("AXIS");
    line->Draw();
    for (int iAlt=0; iAlt<=nAlt; ++iAlt) GrDiffWT[iAlt]->Draw("p");
    leg->Draw();
    TitleTex->DrawLatex(0.5,0.91,("Bias in fit to wrong-tag events - "+longString).c_str());
    c[1]->SaveAs( ("fitResult_alternativeComp_"+shortString+"_wtRECO.pdf").c_str() );
  }
  
  if (plotRECO) {
    c[2] = new TCanvas("c2","c2",800,800);
    HistLab->Draw("AXIS");
    line->Draw();
    for (int iAlt=0; iAlt<=nAlt; ++iAlt) GrDiff[iAlt]->Draw("p");
    leg->Draw();
    TitleTex->DrawLatex(0.5,0.91,("Bias in fit to RECO events - "+longString).c_str());
    c[2]->SaveAs( ("fitResult_alternativeComp_"+shortString+"_RECO.pdf").c_str() );
  }

  // free memory
  for (int iAlt=0; iAlt<=nAlt; ++iAlt) {
    delete ctDiff[iAlt];
    delete wtDiff[iAlt];
    delete Diff[iAlt];
    delete ctDiffErrH[iAlt];
    delete wtDiffErrH[iAlt];
    delete DiffErrH[iAlt];
    delete ctDiffErrL[iAlt];
    delete wtDiffErrL[iAlt];
    delete DiffErrL[iAlt];
  }

}
