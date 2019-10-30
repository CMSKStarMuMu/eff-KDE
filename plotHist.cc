#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

TCanvas* csx [12*nBins];
TCanvas* csy [12*nBins];
TCanvas* csz [12*nBins];
TCanvas* cp1 [12*nBins];
TCanvas* cp2 [12*nBins];

// q2-bin format: [0-8] for one bin
//                [-1] for each bin recursively
// effIndx format: [0] GEN no-filter
//                 [1] GEN filtered
//                 [2] GEN filtered from full MC sample
//                 [3] correct-tag RECO candidates
//                 [4] wrong-tag RECO candidates
//                 [5] RECO candidates
//                 [-1] for each effIndx recursively
// parity format: [0] even
//                [1] odd
//                [-1] for each parity recursively

void plotHistBin(int q2Bin, int effIndx, int parity)
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

  string longString = Form(parity==0?" distribution - q2 bin %i (even)":" distribution - q2 bin %i (odd)",q2Bin);
  switch (effIndx) {
  case 0: longString = "GEN denominator"     +longString; break;
  case 1: longString = "GEN numerator"       +longString; break;
  case 2: longString = "Full denominator"    +longString; break;
  case 3: longString = "Selected correct-tag"+longString; break;
  case 4: longString = "Selected wrong-tag"  +longString; break;
  case 5: longString = "Selected"            +longString;
  }

  string confString = "plotHist_d/KDEhist_"+shortString;

  // Load variables and dataset
  string filename = Form("/eos/user/a/aboletti/BdToKstarMuMu/datasets/PUweight/effDataset_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  RooWorkspace* wsp_data = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,parity));
  if ( !wsp_data || wsp_data->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename<<endl;
    return;
  }
  RooRealVar* ctK = wsp_data->var("ctK");
  RooRealVar* ctL = wsp_data->var("ctL");
  RooRealVar* phi = wsp_data->var("phi");
  if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
    cout<<"Variables not found in file: "<<filename<<endl;
    return;
  }
  RooDataSet* data = (RooDataSet*)wsp_data->data(datasetString.c_str());
  if ( !data || data->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename<<endl;
    return;
  }
  if ( effIndx==5 ) {
    RooDataSet* extradata = (RooDataSet*)wsp_data->data(Form((parity==0?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin));
    if ( !extradata || extradata->IsZombie() ) {
      cout<<"Dataset data_ctRECO_"<<(parity==0?"ev":"od")<<"_b"<<q2Bin<<" not found in file: "<<filename<<endl;
      return;
    }
    data->append(*extradata);
  }
 
  // import KDE histograms
  vector<TH3D*> KDEhists;
  vector<TString> KDEconfs;
  string inFileName = Form((parity==0?"files/KDEhist_b%i_ev.root":"files/KDEhist_b%i_od.root"),q2Bin);
  TFile* fin = TFile::Open( inFileName.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<inFileName<<endl;
    return;
  }
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *) next())) {
    TString s = key->GetName();
    if (!s.Contains( Form("indx%i",effIndx) )) continue;
    cout<<s<<endl;
    KDEhists.push_back( (TH3D*)key->ReadObj() );
    KDEconfs.push_back( s(11,s.Length()-11) );
    KDEhists.back()->SetName(s.Data());
    // cout<<KDEhists.back()->GetName()<<"\t"<<KDEconfs.back()<<endl;
  }

  gStyle->SetOptStat(0);

  int confIndex = 6*nBins*parity + nBins*effIndx + q2Bin;

  // Plot 1D slices of the KDE function and original dataset
  csx[confIndex] = new TCanvas(("csx"+shortString).c_str(),(longString+" - cos(theta_k) slices").c_str(),1500,1500) ;
  csy[confIndex] = new TCanvas(("csy"+shortString).c_str(),(longString+" - cos(theta_l) slices").c_str(),1500,1500) ;
  csz[confIndex] = new TCanvas(("csz"+shortString).c_str(),(longString+" - phi slices").c_str(),1500,1500) ;
  csx[confIndex]->Divide(5,5);
  csy[confIndex]->Divide(5,5);
  csz[confIndex]->Divide(5,5);
  vector <TH1D*> histSliceX [25];
  vector <TH1D*> histSliceY [25];
  vector <TH1D*> histSliceZ [25];
  TH1D* distSliceX [25];
  TH1D* distSliceY [25];
  TH1D* distSliceZ [25];

  TLegend legSl (0.42,0.75,0.9,0.9);

  // width of the slices in the hidden variables ("border" is half of it)
  double border = 0.025;
  // variables to be filled with global maxima
  double maxValX = 0;
  double maxValY = 0;
  double maxValZ = 0;

  // loop over slice grid
  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

      // central values and borders of the slices in the hidden variables
      double centA = -0.8 + 1.6*i/4;
      double centB = -0.8 + 1.6*j/4;
      string cutX = Form("fabs(ctL-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centA,border,TMath::Pi(),centB,border);
      string cutY = Form("fabs(ctK-(%1.1f))<%1.3f && fabs((phi/%1.6f)-(%1.1f))<%1.3f",centA,border,TMath::Pi(),centB,border);
      string cutZ = Form("fabs(ctK-(%1.1f))<%1.3f && fabs(ctL-(%1.1f))<%1.3f",centA,border,centB,border);

      // slicing distribution
      distSliceX[5*i+j] = (TH1D*)data->reduce(RooArgSet(*ctK),cutX.c_str())
	->createHistogram( Form("distSliceX_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(20,-1,1) );
      distSliceY[5*i+j] = (TH1D*)data->reduce(RooArgSet(*ctL),cutY.c_str())
	->createHistogram( Form("distSliceY_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(20,-1,1) );
      distSliceZ[5*i+j] = (TH1D*)data->reduce(RooArgSet(*phi),cutZ.c_str())
	->createHistogram( Form("distSliceZ_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(20,-TMath::Pi(),TMath::Pi()) );
      distSliceX[5*i+j]->SetTitle( Form("%s - slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Events",longString.c_str(),centA,centB*TMath::Pi()) );
      distSliceY[5*i+j]->SetTitle( Form("%s - slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Events",longString.c_str(),centA,centB*TMath::Pi()) );
      distSliceZ[5*i+j]->SetTitle( Form("%s - slice ctK=%1.2f ctL=%1.2f;#phi;Events"           ,longString.c_str(),centA,centB) );
      if (i+j==0) legSl.AddEntry(distSliceX[0],"Event distribution","lep");
	
      // producing 1D slices of KDE descriptions
      for (int iHist=0; iHist<KDEhists.size(); ++iHist) {
	if ( iHist >= KDEconfs.size() ) {
	  cout<<"ERROR: histo vector and name vector have different sizes"<<endl;
	  return;
	}
	TH3D* KDEhist = KDEhists[iHist];
	histSliceX[5*i+j].push_back( (TH1D*) KDEhist->ProjectionX( Form("histSliceX_%i_%i_%s_%s",i,j,shortString.c_str(),KDEconfs[iHist].Data()),
								      KDEhist->GetYaxis()->FindBin(centA            ),
								      KDEhist->GetYaxis()->FindBin(centA            ),
								      KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()),
								      KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()) ) ); 
	histSliceY[5*i+j].push_back( (TH1D*) KDEhist->ProjectionY( Form("histSliceY_%i_%i_%s_%s",i,j,shortString.c_str(),KDEconfs[iHist].Data()),
								      KDEhist->GetXaxis()->FindBin(centA            ),
								      KDEhist->GetXaxis()->FindBin(centA            ),
								      KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()),
								      KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()) ) );
	histSliceZ[5*i+j].push_back( (TH1D*) KDEhist->ProjectionZ( Form("histSliceZ_%i_%i_%s_%s",i,j,shortString.c_str(),KDEconfs[iHist].Data()),
								      KDEhist->GetXaxis()->FindBin(centA),
								      KDEhist->GetXaxis()->FindBin(centA),
								      KDEhist->GetYaxis()->FindBin(centB),
								      KDEhist->GetYaxis()->FindBin(centB) ) );

	// scale KDE slices to match distribution density
	histSliceX[5*i+j].back()->Scale( 2*border * 2*border*TMath::Pi() / KDEhist->GetYaxis()->GetBinWidth(1) / KDEhist->GetZaxis()->GetBinWidth(1) );
	histSliceY[5*i+j].back()->Scale( 2*border * 2*border*TMath::Pi() / KDEhist->GetXaxis()->GetBinWidth(1) / KDEhist->GetZaxis()->GetBinWidth(1) );
	histSliceZ[5*i+j].back()->Scale( 2*border * 2*border             / KDEhist->GetXaxis()->GetBinWidth(1) / KDEhist->GetYaxis()->GetBinWidth(1) );
	histSliceX[5*i+j].back()->Scale( 1.0 * histSliceX[5*i+j].back()->GetNbinsX() / distSliceX[5*i+j]->GetNbinsX() );
	histSliceY[5*i+j].back()->Scale( 1.0 * histSliceY[5*i+j].back()->GetNbinsX() / distSliceY[5*i+j]->GetNbinsX() );
	histSliceZ[5*i+j].back()->Scale( 1.0 * histSliceZ[5*i+j].back()->GetNbinsX() / distSliceZ[5*i+j]->GetNbinsX() );

	if (i+j==0) legSl.AddEntry(histSliceX[0].back(),Form("KDE %s",KDEconfs[iHist].Data()),"l");

      }

      // plot in canvas
      csx[confIndex]->cd(5*j+i+1);
      // gPad->SetLeftMargin(0.18);
      distSliceX[5*i+j]->GetYaxis()->SetTitleOffset(1.4);
      distSliceX[5*i+j]->Draw();
      for (int iHist=0; iHist<KDEhists.size(); ++iHist) {
	histSliceX[5*i+j].at(iHist)->SetLineWidth(2);
	histSliceX[5*i+j].at(iHist)->SetLineColor(1+iHist);
	histSliceX[5*i+j].at(iHist)->Draw("sameLHIST");
      }
      distSliceX[5*i+j]->Draw("same");
      legSl.Draw("same");

      csy[confIndex]->cd(5*j+i+1);
      // gPad->SetLeftMargin(0.18);
      distSliceY[5*i+j]->GetYaxis()->SetTitleOffset(1.4);
      distSliceY[5*i+j]->Draw();
      for (int iHist=0; iHist<KDEhists.size(); ++iHist) {
	(histSliceY[5*i+j])[iHist]->SetLineWidth(2);
	(histSliceY[5*i+j])[iHist]->SetLineColor(1+iHist);
	(histSliceY[5*i+j])[iHist]->Draw("sameLHIST");
      }
      distSliceY[5*i+j]->Draw("same");
      legSl.Draw("same");

      csz[confIndex]->cd(5*j+i+1);
      // gPad->SetLeftMargin(0.18);
      distSliceZ[5*i+j]->GetYaxis()->SetTitleOffset(1.4);
      distSliceZ[5*i+j]->Draw();
      for (int iHist=0; iHist<KDEhists.size(); ++iHist) {
	(histSliceZ[5*i+j])[iHist]->SetLineWidth(2);
	(histSliceZ[5*i+j])[iHist]->SetLineColor(1+iHist);
	(histSliceZ[5*i+j])[iHist]->Draw("sameLHIST");
      }
      distSliceZ[5*i+j]->Draw("same");
      legSl.Draw("same");

      // checking maximum value
      if ( maxValX<distSliceX[5*i+j]->GetMaximum() ) maxValX = distSliceX[5*i+j]->GetMaximum();
      if ( maxValY<distSliceY[5*i+j]->GetMaximum() ) maxValY = distSliceY[5*i+j]->GetMaximum();
      if ( maxValZ<distSliceZ[5*i+j]->GetMaximum() ) maxValZ = distSliceZ[5*i+j]->GetMaximum();
      // if ( maxValX<histSliceX[5*i+j]->GetMaximum() ) maxValX = histSliceX[5*i+j]->GetMaximum();
      // if ( maxValY<histSliceY[5*i+j]->GetMaximum() ) maxValY = histSliceY[5*i+j]->GetMaximum();
      // if ( maxValZ<histSliceZ[5*i+j]->GetMaximum() ) maxValZ = histSliceZ[5*i+j]->GetMaximum();
      
    }
  // set uniform y-axis ranges
  for (int i=0; i<25; ++i) {
    distSliceX[i]->GetYaxis()->SetRangeUser(0,maxValX*1.25);
    distSliceY[i]->GetYaxis()->SetRangeUser(0,maxValY*1.25);
    distSliceZ[i]->GetYaxis()->SetRangeUser(0,maxValZ*1.25);
  }
  csx[confIndex]->SaveAs( (confString+Form("_CTKslices_dp%i.pdf",(int)(border*200))).c_str() );
  csy[confIndex]->SaveAs( (confString+Form("_CTLslices_dp%i.pdf",(int)(border*200))).c_str() );
  csz[confIndex]->SaveAs( (confString+Form("_PHIslices_dp%i.pdf",(int)(border*200))).c_str() );

  // Plot 1D projections
  cp1[confIndex] = new TCanvas(("cp1"+shortString).c_str(),(longString+" - 1D Projections").c_str(),2000,700) ;
  cp1[confIndex]->Divide(3,1);
  double x1leg = 0.11;
  double x2leg = 0.11;
  double x3leg = 0.11;
  double y1leg = 0.1;
  double y2leg = 0.1;
  double y3leg = 0.1;
  if (effIndx==1 || effIndx==2 || effIndx==3) {
    x3leg = 0.45;
    if (q2Bin<4) {
      x1leg = 0.45;
      y1leg = 0.78;
      x2leg = 0.28;
    } else
      x2leg = 0.45;
  }
  if (effIndx==4) {
    x3leg = 0.45;
    if (q2Bin<4) {
      x1leg = 0.28;
      x2leg = 0.28;
    } else {
      x2leg = 0.45;
      y2leg = 0.78;
    }
  }
  TLegend leg1Pr (x1leg,y1leg,x1leg+0.45,y1leg+0.12);
  TLegend leg2Pr (x2leg,y2leg,x2leg+0.45,y2leg+0.12);
  TLegend leg3Pr (x3leg,y3leg,x3leg+0.45,y3leg+0.12);

  TH1D* distProj1X = (TH1D*)data->createHistogram( Form("distProj1X_%s",shortString.c_str()), *ctK, Binning(50,-1,1) );
  TH1D* distProj1Y = (TH1D*)data->createHistogram( Form("distProj1Y_%s",shortString.c_str()), *ctL, Binning(50,-1,1) );
  TH1D* distProj1Z = (TH1D*)data->createHistogram( Form("distProj1Z_%s",shortString.c_str()), *phi, Binning(50,-TMath::Pi(),TMath::Pi()) );
  distProj1X->SetTitle( Form("%s;cos(#theta_{K});Events",longString.c_str()) );
  distProj1Y->SetTitle( Form("%s;cos(#theta_{L});Events",longString.c_str()) );
  distProj1Z->SetTitle( Form("%s;#phi;Events",longString.c_str()) );
  leg1Pr.AddEntry(distProj1X,"Event distribution","lep");
  leg2Pr.AddEntry(distProj1X,"Event distribution","lep");
  leg3Pr.AddEntry(distProj1X,"Event distribution","lep");
  vector <TH1D*> histProj1X;
  vector <TH1D*> histProj1Y;
  vector <TH1D*> histProj1Z;

  // Compute 1D cumulative distributions of dataset, used for Kolmogorov-Smirnov test
  auto lowXbinEdges = new double[KDEhists[0]->GetNbinsX()];
  auto lowYbinEdges = new double[KDEhists[0]->GetNbinsY()];
  auto lowZbinEdges = new double[KDEhists[0]->GetNbinsZ()];
  KDEhists[0]->GetXaxis()->GetLowEdge(lowXbinEdges);
  KDEhists[0]->GetYaxis()->GetLowEdge(lowYbinEdges);
  KDEhists[0]->GetZaxis()->GetLowEdge(lowZbinEdges);
  auto dataCumulX = new double[KDEhists[0]->GetNbinsX()];
  auto dataCumulY = new double[KDEhists[0]->GetNbinsY()];
  auto dataCumulZ = new double[KDEhists[0]->GetNbinsZ()];
  double totalData = data->sumEntries();
  for (int iBin=0; iBin<KDEhists[0]->GetNbinsX(); ++iBin) dataCumulX[iBin] = data->sumEntries(Form("ctK<%f",lowXbinEdges[iBin])) / totalData;
  for (int iBin=0; iBin<KDEhists[0]->GetNbinsY(); ++iBin) dataCumulY[iBin] = data->sumEntries(Form("ctL<%f",lowYbinEdges[iBin])) / totalData;
  for (int iBin=0; iBin<KDEhists[0]->GetNbinsZ(); ++iBin) dataCumulZ[iBin] = data->sumEntries(Form("phi<%f",lowZbinEdges[iBin])) / totalData;

  for (int iHist=0; iHist<KDEhists.size(); ++iHist) {
    if ( iHist >= KDEconfs.size() ) {
      cout<<"ERROR: histo vector and name vector have different sizes"<<endl;
      return;
    }
    TH3D* KDEhist = KDEhists[iHist];
    histProj1X.push_back( (TH1D*)KDEhist->ProjectionX( Form("histProj1X_%s_%s",shortString.c_str(),KDEconfs[iHist].Data()) ) );
    histProj1Y.push_back( (TH1D*)KDEhist->ProjectionY( Form("histProj1Y_%s_%s",shortString.c_str(),KDEconfs[iHist].Data()) ) );
    histProj1Z.push_back( (TH1D*)KDEhist->ProjectionZ( Form("histProj1Z_%s_%s",shortString.c_str(),KDEconfs[iHist].Data()) ) );
    histProj1X.back()->Scale( 1.0 * histProj1X.back()->GetNbinsX() / distProj1X->GetNbinsX() );
    histProj1Y.back()->Scale( 1.0 * histProj1Y.back()->GetNbinsX() / distProj1Y->GetNbinsX() );
    histProj1Z.back()->Scale( 1.0 * histProj1Z.back()->GetNbinsX() / distProj1Z->GetNbinsX() );

    // Compute Kolmogorov-Smirnov test on 1D projections to quantitatively compare performances
    auto partIntX = histProj1X.back()->GetIntegral();
    auto partIntY = histProj1Y.back()->GetIntegral();
    auto partIntZ = histProj1Z.back()->GetIntegral();
    double xKSmax = 0;
    double yKSmax = 0;
    double zKSmax = 0;
    double diff;
    for (int iBin=1; iBin<histProj1X.back()->GetNbinsX(); ++iBin) {
      diff = fabs( partIntX[iBin] - dataCumulX[iBin] );
      if ( xKSmax<diff ) xKSmax = diff;
    }
    for (int iBin=1; iBin<histProj1Y.back()->GetNbinsX(); ++iBin) {
      diff = fabs( partIntY[iBin] - dataCumulY[iBin] );
      if ( yKSmax<diff ) yKSmax = diff;
    }
    for (int iBin=1; iBin<histProj1Z.back()->GetNbinsX(); ++iBin) {
      diff = fabs( partIntZ[iBin] - dataCumulZ[iBin] );
      if ( zKSmax<diff ) zKSmax = diff;
    }
    leg1Pr.AddEntry(histProj1X.back(),Form("KDE %s - KS=%1.2f",KDEconfs[iHist].Data(),xKSmax*sqrt(data->numEntries())),"l");
    leg2Pr.AddEntry(histProj1X.back(),Form("KDE %s - KS=%1.2f",KDEconfs[iHist].Data(),yKSmax*sqrt(data->numEntries())),"l");
    leg3Pr.AddEntry(histProj1X.back(),Form("KDE %s - KS=%1.2f",KDEconfs[iHist].Data(),zKSmax*sqrt(data->numEntries())),"l");
  }

  cp1[confIndex]->cd(1);
  gPad->SetLeftMargin(0.11);
  distProj1X->GetYaxis()->SetTitleOffset(1.7);
  distProj1X->SetMinimum(0);
  distProj1X->Draw();
  for (int iHist=0; iHist<histProj1X.size(); ++iHist) {
    histProj1X[iHist]->SetLineWidth(2);
    histProj1X[iHist]->SetLineColor(1+iHist);
    histProj1X[iHist]->Draw("sameLHIST");
  }
  distProj1X->Draw("same");
  leg1Pr.Draw("same");
  cp1[confIndex]->cd(2);
  gPad->SetLeftMargin(0.11);
  distProj1Y->GetYaxis()->SetTitleOffset(1.7);
  distProj1Y->SetMinimum(0);
  distProj1Y->Draw();
  for (int iHist=0; iHist<histProj1Y.size(); ++iHist) {
    histProj1Y[iHist]->SetLineWidth(2);
    histProj1Y[iHist]->SetLineColor(1+iHist);
    histProj1Y[iHist]->Draw("sameLHIST");
  }
  distProj1Y->Draw("same");
  leg2Pr.Draw("same");
  cp1[confIndex]->cd(3);
  gPad->SetLeftMargin(0.11);
  distProj1Z->GetYaxis()->SetTitleOffset(1.7);
  distProj1Z->SetMinimum(0);
  distProj1Z->Draw();
  for (int iHist=0; iHist<histProj1Z.size(); ++iHist) {
    histProj1Z[iHist]->SetLineWidth(2);
    histProj1Z[iHist]->SetLineColor(1+iHist);
    histProj1Z[iHist]->Draw("sameLHIST");
  }
  distProj1Z->Draw("same");
  leg3Pr.Draw("same");
  cp1[confIndex]->SaveAs( (confString+"_Proj1D.pdf").c_str() );
    
  // Plot 2D projections
  cp2[confIndex] = new TCanvas(("cp2"+shortString).c_str(),(longString+" - 2D Projections").c_str(),1500,500*KDEhists.size()) ;
  cp2[confIndex]->Divide(3,KDEhists.size());
  for (int iHist=0; iHist<KDEhists.size(); ++iHist) {
    if ( iHist >= KDEconfs.size() ) {
      cout<<"ERROR: histo vector and name vector have different sizes"<<endl;
      return;
    }
    TH3D* KDEhist = KDEhists[iHist];
    TH2D* histProj2XY = (TH2D*)KDEhist->Project3D( "xy" );
    TH2D* histProj2XZ = (TH2D*)KDEhist->Project3D( "xz" );
    TH2D* histProj2YZ = (TH2D*)KDEhist->Project3D( "yz" );
    histProj2XY->SetName( Form("histProj2XY_%s_%s",shortString.c_str(),KDEconfs[iHist].Data()) );
    histProj2XZ->SetName( Form("histProj2XZ_%s_%s",shortString.c_str(),KDEconfs[iHist].Data()) );
    histProj2YZ->SetName( Form("histProj2YZ_%s_%s",shortString.c_str(),KDEconfs[iHist].Data()) );
    histProj2XY->SetTitle( Form("%s %s;cos(#theta_{L});cos(#theta_{K});",longString.c_str(),KDEconfs[iHist].Data()) );
    histProj2XZ->SetTitle( Form("%s %s;#phi;cos(#theta_{K});",longString.c_str(),KDEconfs[iHist].Data()) );
    histProj2YZ->SetTitle( Form("%s %s;#phi;cos(#theta_{L});",longString.c_str(),KDEconfs[iHist].Data()) );
    histProj2XY->GetXaxis()->SetTitleOffset(1.4);
    histProj2XY->GetYaxis()->SetTitleOffset(2);
    histProj2XZ->GetXaxis()->SetTitleOffset(1.4);
    histProj2XZ->GetYaxis()->SetTitleOffset(2);
    histProj2YZ->GetXaxis()->SetTitleOffset(1.4);
    histProj2YZ->GetYaxis()->SetTitleOffset(2);
    cp2[confIndex]->cd(1+3*iHist);
    histProj2XY->SetMinimum(0);
    histProj2XY->Draw("SURF3");
    cp2[confIndex]->cd(2+3*iHist);
    histProj2XZ->SetMinimum(0);
    histProj2XZ->Draw("SURF3");
    cp2[confIndex]->cd(3+3*iHist);
    histProj2YZ->SetMinimum(0);
    histProj2YZ->Draw("SURF3");
  }
  cp2[confIndex]->SaveAs( (confString+"_Proj2D.pdf").c_str() );

}

void plotHistBin2(int q2Bin, int effIndx, int parity)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      plotHistBin(q2Bin, effIndx, parity);
  else
    plotHistBin(q2Bin, effIndx, parity);
}

void plotHistBin1(int q2Bin, int effIndx, int parity)
{
  if ( effIndx==-1 )
    for (effIndx=0; effIndx<6; ++effIndx)
      plotHistBin2(q2Bin, effIndx, parity);
  else
    plotHistBin2(q2Bin, effIndx, parity);
}

void plotHist(int q2Bin = -1, int effIndx = -1, int parity = -1)
{

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( effIndx<-1 || effIndx>5 ) return;

  if ( parity<-1 || parity>1 ) return;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( effIndx==-1 ) cout<<"Running all the efficiency terms"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      plotHistBin1(q2Bin, effIndx, parity);
  else
    plotHistBin1(q2Bin, effIndx, parity);
  
}
