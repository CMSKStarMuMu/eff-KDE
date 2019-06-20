#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

bool safeFilling = false;

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

void mergeParSub_rooKeysBin(int q2Bin, int effIndx, int parity, float width, int xbins, int ybins, int zbins, int totdiv, bool plot)
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

  // string confString = Form("KDEhist_%s_rooKeys_mw%.2f_%i_%i_%i",shortString.c_str(),width,xbins,ybins,zbins);
  string confString = Form("/eos/user/a/aboletti/BdToKstarMuMu/eff-KDE/KDEhist_%s_rooKeys_mw%.2f_%i_%i_%i",shortString.c_str(),width,xbins,ybins,zbins);

  TH3D* KDEhist = new TH3D(Form("KDEhist_%s",shortString.c_str()),Form("KDEhist_%s",shortString.c_str()),xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  KDEhist->Sumw2();

  // import partial histograms
  string inFileName;
  int goodHistCnt = 0;
  for (int ndiv=0; ndiv<totdiv; ++ndiv) {
    inFileName = confString+Form("_%i-frac-%i.root",ndiv,totdiv);
    TFile* fin = TFile::Open( inFileName.c_str() );
    if ( !fin || !fin->IsOpen() ) {
      cout<<"File not found: "<<inFileName<<endl;
      continue;
    }
    TH3D* partHist = (TH3D*)fin->Get(("KDEhist_"+shortString+"__ctK_ctL_phi").c_str());
    if ( !partHist || partHist->IsZombie() ) {
      cout<<"Histogram not found in file: "<<inFileName<<endl;
      continue;
    }
    KDEhist->Add(partHist);
    ++goodHistCnt;
    delete partHist;
    fin->Close();
    delete fin;
  }
  if ( goodHistCnt!=totdiv ) {
    cout<<"Warning! Not all partial histograms found: "<<goodHistCnt<<" of "<<totdiv<<endl;
    if (safeFilling) return;
  }

  // check histogram
  double minVal = KDEhist->GetMinimum();
  if ( minVal<=0 ) {
    cout<<"Histogram has empty bins, this is bad. Abort!"<<endl;
    return;
  }

  // save histogram in map to file
  TFile* fout = TFile::Open( Form((parity==0?"files/KDEhist_b%i_ev.root":"files/KDEhist_b%i_od.root"),q2Bin), "UPDATE" );
  KDEhist->Write( Form("hist_indx%i_%1.2f_%i_%i_%i",effIndx,width,xbins,ybins,zbins), TObject::kWriteDelete );
  fout->Close();

  if (plot) {
    gStyle->SetOptStat(0);

    int confIndex = 6*nBins*parity + nBins*effIndx + q2Bin;

    // Load variables and dataset
    string filename = Form("effDataset_b%i.root",q2Bin);
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

    // Plot 1D slices of the KDE function and original dataset
    vector <TH1D*> histSliceX;
    vector <TH1D*> histSliceY;
    vector <TH1D*> histSliceZ;
    vector <TH1D*> distSliceX;
    vector <TH1D*> distSliceY;
    vector <TH1D*> distSliceZ;
    csx[confIndex] = new TCanvas(("csx"+shortString).c_str(),(longString+" - cos(theta_k) slices").c_str(),1500,1500) ;
    csy[confIndex] = new TCanvas(("csy"+shortString).c_str(),(longString+" - cos(theta_l) slices").c_str(),1500,1500) ;
    csz[confIndex] = new TCanvas(("csz"+shortString).c_str(),(longString+" - phi slices").c_str(),1500,1500) ;
    csx[confIndex]->Divide(5,5);
    csy[confIndex]->Divide(5,5);
    csz[confIndex]->Divide(5,5);

    // width of the slices in the hidden variables ("border" is half of it)
    double border = 0.05;
    // variables to be filled with global maxima
    double maxValX = 0;
    double maxValY = 0;
    double maxValZ = 0;

    // loop over slice grid
    for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

	// central values and borders of the slices in the hidden variables
	double centA = -0.8 + 1.6*i/4;
	double centB = -0.8 + 1.6*j/4;
	string cutX = Form("fabs(ctL-%1.1f)<%1.3f && fabs((phi/%1.6f)-%1.1f)<%1.3f",centA,border,TMath::Pi(),centB,border);
	string cutY = Form("fabs(ctK-%1.1f)<%1.3f && fabs((phi/%1.6f)-%1.1f)<%1.3f",centA,border,TMath::Pi(),centB,border);
	string cutZ = Form("fabs(ctK-%1.1f)<%1.3f && fabs(ctL-%1.1f)<%1.3f",centA,border,centB,border);

	// slicing distribution
	distSliceX.push_back( (TH1D*)data->reduce(RooArgSet(*ctK),cutX.c_str())
			      ->createHistogram( Form("distSliceX_%i_%i_%s",i,j,shortString.c_str()), *ctK, Binning(20,-1,1) ) );
	distSliceY.push_back( (TH1D*)data->reduce(RooArgSet(*ctL),cutY.c_str())
			      ->createHistogram( Form("distSliceY_%i_%i_%s",i,j,shortString.c_str()), *ctL, Binning(20,-1,1) ) );
	distSliceZ.push_back( (TH1D*)data->reduce(RooArgSet(*phi),cutZ.c_str())
			      ->createHistogram( Form("distSliceZ_%i_%i_%s",i,j,shortString.c_str()), *phi, Binning(20,-TMath::Pi(),TMath::Pi()) ) );
	distSliceX.back()->SetTitle( Form("%s - slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Events",longString.c_str(),centA,centB*TMath::Pi()) );
	distSliceY.back()->SetTitle( Form("%s - slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Events",longString.c_str(),centA,centB*TMath::Pi()) );
	distSliceZ.back()->SetTitle( Form("%s - slice ctK=%1.2f ctL=%1.2f;#phi;Events"           ,longString.c_str(),centA,centB) );
	
	// producing 1D slices of KDE description
	histSliceX.push_back( (TH1D*)KDEhist->ProjectionX( Form("histSliceX_%i_%i_%s",i,j,shortString.c_str()),
						    KDEhist->GetYaxis()->FindBin(centA            ), KDEhist->GetYaxis()->FindBin(centA            ),
						    KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()), KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()) ) ); 
	histSliceY.push_back( (TH1D*)KDEhist->ProjectionY( Form("histSliceY_%i_%i_%s",i,j,shortString.c_str()),
						    KDEhist->GetXaxis()->FindBin(centA            ), KDEhist->GetXaxis()->FindBin(centA            ),
						    KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()), KDEhist->GetZaxis()->FindBin(centB*TMath::Pi()) ) );
	histSliceZ.push_back( (TH1D*)KDEhist->ProjectionZ( Form("histSliceZ_%i_%i_%s",i,j,shortString.c_str()),
						    KDEhist->GetXaxis()->FindBin(centA), KDEhist->GetXaxis()->FindBin(centA),
						    KDEhist->GetYaxis()->FindBin(centB), KDEhist->GetYaxis()->FindBin(centB) ) );

	// scale KDE slices to match distribution density
	histSliceX.back()->Scale( 2*border * 2*border*TMath::Pi() / KDEhist->GetYaxis()->GetBinWidth(1) / KDEhist->GetZaxis()->GetBinWidth(1) );
	histSliceY.back()->Scale( 2*border * 2*border*TMath::Pi() / KDEhist->GetXaxis()->GetBinWidth(1) / KDEhist->GetZaxis()->GetBinWidth(1) );
	histSliceZ.back()->Scale( 2*border * 2*border             / KDEhist->GetXaxis()->GetBinWidth(1) / KDEhist->GetYaxis()->GetBinWidth(1) );
	histSliceX.back()->Scale( 1.0 * histSliceX.back()->GetNbinsX() / distSliceX.back()->GetNbinsX() );
	histSliceY.back()->Scale( 1.0 * histSliceY.back()->GetNbinsX() / distSliceY.back()->GetNbinsX() );
	histSliceZ.back()->Scale( 1.0 * histSliceZ.back()->GetNbinsX() / distSliceZ.back()->GetNbinsX() );
	histSliceX.back()->SetLineWidth(2);
	histSliceY.back()->SetLineWidth(2);
	histSliceZ.back()->SetLineWidth(2);
	histSliceX.back()->SetLineColor(kRed+1);
	histSliceY.back()->SetLineColor(kRed+1);
	histSliceZ.back()->SetLineColor(kRed+1);

	// plot in canvas
	csx[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	distSliceX.back()->GetYaxis()->SetTitleOffset(1.7);
	distSliceX.back()->Draw();
	histSliceX.back()->Draw("sameLHIST");

	csy[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	distSliceY.back()->GetYaxis()->SetTitleOffset(1.7);
	distSliceY.back()->Draw();
	histSliceY.back()->Draw("sameLHIST");

	csz[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	distSliceZ.back()->GetYaxis()->SetTitleOffset(1.7);
	distSliceZ.back()->Draw();
	histSliceZ.back()->Draw("sameLHIST");

	// checking maximum value
	if ( maxValX<distSliceX.back()->GetMaximum() ) maxValX = distSliceX.back()->GetMaximum();
	if ( maxValY<distSliceY.back()->GetMaximum() ) maxValY = distSliceY.back()->GetMaximum();
	if ( maxValZ<distSliceZ.back()->GetMaximum() ) maxValZ = distSliceZ.back()->GetMaximum();
	if ( maxValX<histSliceX.back()->GetMaximum() ) maxValX = histSliceX.back()->GetMaximum();
	if ( maxValY<histSliceY.back()->GetMaximum() ) maxValY = histSliceY.back()->GetMaximum();
	if ( maxValZ<histSliceZ.back()->GetMaximum() ) maxValZ = histSliceZ.back()->GetMaximum();

      }

    // set uniform y-axis ranges
    for (int i=0; i<distSliceX.size(); ++i) distSliceX[i]->GetYaxis()->SetRangeUser(0,maxValX*1.1);
    for (int i=0; i<distSliceY.size(); ++i) distSliceY[i]->GetYaxis()->SetRangeUser(0,maxValY*1.1);
    for (int i=0; i<distSliceZ.size(); ++i) distSliceZ[i]->GetYaxis()->SetRangeUser(0,maxValZ*1.1);

    csx[confIndex]->SaveAs( (confString+Form("_CTKslices_dp%i.pdf",(int)(border*200))).c_str() );
    csy[confIndex]->SaveAs( (confString+Form("_CTLslices_dp%i.pdf",(int)(border*200))).c_str() );
    csz[confIndex]->SaveAs( (confString+Form("_PHIslices_dp%i.pdf",(int)(border*200))).c_str() );

    // Plot 1D projections
    cp1[confIndex] = new TCanvas(("cp1"+shortString).c_str(),(longString+" - 1D Projections").c_str(),2000,700) ;
    cp1[confIndex]->Divide(3,1);
    TH1D* distProj1X = (TH1D*)data->createHistogram( Form("distProj1X_%s",shortString.c_str()), *ctK, Binning(20,-1,1) );
    TH1D* distProj1Y = (TH1D*)data->createHistogram( Form("distProj1Y_%s",shortString.c_str()), *ctL, Binning(20,-1,1) );
    TH1D* distProj1Z = (TH1D*)data->createHistogram( Form("distProj1Z_%s",shortString.c_str()), *phi, Binning(20,-TMath::Pi(),TMath::Pi()) );
    distProj1X->SetTitle( Form("%s;cos(#theta_{K});Events",longString.c_str()) );
    distProj1Y->SetTitle( Form("%s;cos(#theta_{L});Events",longString.c_str()) );
    distProj1Z->SetTitle( Form("%s;#phi;Events",longString.c_str()) );
    TH1D* histProj1X = (TH1D*)KDEhist->ProjectionX( Form("histProj1X_%s",shortString.c_str()) );
    TH1D* histProj1Y = (TH1D*)KDEhist->ProjectionY( Form("histProj1Y_%s",shortString.c_str()) );
    TH1D* histProj1Z = (TH1D*)KDEhist->ProjectionZ( Form("histProj1Z_%s",shortString.c_str()) );
    histProj1X->Scale( 1.0 * histProj1X->GetNbinsX() / distProj1X->GetNbinsX() );
    histProj1Y->Scale( 1.0 * histProj1Y->GetNbinsX() / distProj1Y->GetNbinsX() );
    histProj1Z->Scale( 1.0 * histProj1Z->GetNbinsX() / distProj1Z->GetNbinsX() );
    histProj1X->SetLineWidth(2);
    histProj1Y->SetLineWidth(2);
    histProj1Z->SetLineWidth(2);
    histProj1X->SetLineColor(kRed+1);
    histProj1Y->SetLineColor(kRed+1);
    histProj1Z->SetLineColor(kRed+1);
    cp1[confIndex]->cd(1);
    gPad->SetLeftMargin(0.18);
    distProj1X->GetYaxis()->SetTitleOffset(1.7);
    distProj1X->SetMinimum(0);
    distProj1X->Draw();
    histProj1X->Draw("sameLHIST");
    cp1[confIndex]->cd(2);
    gPad->SetLeftMargin(0.18);
    distProj1Y->GetYaxis()->SetTitleOffset(1.7);
    distProj1Y->SetMinimum(0);
    distProj1Y->Draw();
    histProj1Y->Draw("sameLHIST");
    cp1[confIndex]->cd(3);
    gPad->SetLeftMargin(0.18);
    distProj1Z->GetYaxis()->SetTitleOffset(1.7);
    distProj1Z->SetMinimum(0);
    distProj1Z->Draw();
    histProj1Z->Draw("sameLHIST");
    cp1[confIndex]->SaveAs( (confString+"_Proj1D.pdf").c_str() );
    
    // Plot 2D projections
    cp2[confIndex] = new TCanvas(("cp2"+shortString).c_str(),(longString+" - 2D Projections").c_str(),2000,700) ;
    cp2[confIndex]->Divide(3,1);
    TH2D* histProj2XY = (TH2D*)KDEhist->Project3D( "xy" );
    TH2D* histProj2XZ = (TH2D*)KDEhist->Project3D( "xz" );
    TH2D* histProj2YZ = (TH2D*)KDEhist->Project3D( "yz" );
    histProj2XY->SetName( Form("histProj2XY_%s",shortString.c_str()) );
    histProj2XZ->SetName( Form("histProj2XZ_%s",shortString.c_str()) );
    histProj2YZ->SetName( Form("histProj2YZ_%s",shortString.c_str()) );
    histProj2XY->SetTitle( Form("%s;cos(#theta_{K});cos(#theta_{K});Events",longString.c_str()) );
    histProj2XZ->SetTitle( Form("%s;cos(#theta_{K});#phi;Events",longString.c_str()) );
    histProj2YZ->SetTitle( Form("%s;cos(#theta_{L});#phi;Events",longString.c_str()) );
    cp2[confIndex]->cd(1);
    histProj2XY->SetMinimum(0);
    histProj2XY->Draw("SURF3");
    cp2[confIndex]->cd(2);
    histProj2XZ->SetMinimum(0);
    histProj2XZ->Draw("SURF3");
    cp2[confIndex]->cd(3);
    histProj2YZ->SetMinimum(0);
    histProj2YZ->Draw("SURF3");
    cp2[confIndex]->SaveAs( (confString+"_Proj2D.pdf").c_str() );

  }
  
  // Remind user to delete partial files
  cout<<endl<<"Please, remove partial files running:\nrm "<<confString<<Form("_*-frac-%i.root",totdiv)<<endl<<endl;
  
}

void mergeParSub_rooKeysBin2(int q2Bin, int effIndx, int parity, float width, int xbins, int ybins, int zbins, int totdiv, bool plot)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      mergeParSub_rooKeysBin(q2Bin, effIndx, parity, width, xbins, ybins, zbins, totdiv, plot);
  else
    mergeParSub_rooKeysBin(q2Bin, effIndx, parity, width, xbins, ybins, zbins, totdiv, plot);
}

void mergeParSub_rooKeysBin1(int q2Bin, int effIndx, int parity, float width, int xbins, int ybins, int zbins, int totdiv, bool plot)
{
  if ( effIndx==-1 )
    for (effIndx=0; effIndx<6; ++effIndx)
      mergeParSub_rooKeysBin2(q2Bin, effIndx, parity, width, xbins, ybins, zbins, totdiv, plot);
  else
    mergeParSub_rooKeysBin2(q2Bin, effIndx, parity, width, xbins, ybins, zbins, totdiv, plot);
}

void mergeParSub_rooKeys(int q2Bin = -1, int effIndx = -1, int parity = -1, float width = 0.3, int xbins=50, int ybins = 50, int zbins = 50, int totdiv = 50, bool plot = false)
{

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( effIndx<-1 || effIndx>5 ) return;

  if ( parity<-1 || parity>1 ) return;

  if ( width<=0 ) return;

  if ( xbins<1 ) return;
  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;

  if ( totdiv<1 ) return;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( effIndx==-1 ) cout<<"Running all the efficiency terms"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      mergeParSub_rooKeysBin1(q2Bin, effIndx, parity, width, xbins, ybins, zbins, totdiv, plot);
  else
    mergeParSub_rooKeysBin1(q2Bin, effIndx, parity, width, xbins, ybins, zbins, totdiv, plot);
  
}
