#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

void plotEff_customKDE_uniformBin(int q2Bin, bool tagFlag, int kernel, int xbins, int ybins, int zbins, bool doMirror);

TCanvas* c  [2*nBins];
TCanvas* c1 [2*nBins];
TCanvas* c3 [2*nBins];

void plotEff_customKDE_uniform(int q2Bin, int tagFlag, int kernel = 1000, int doMirror = 1, int xbins=25, int ybins = 0, int zbins = 0)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( doMirror<0 || doMirror>1 ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      plotEff_customKDE_uniformBin(q2Bin, (tagFlag==1), kernel, xbins, ybins, zbins, (doMirror==1));
    }
    if (tagFlag == 2) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      plotEff_customKDE_uniformBin(q2Bin, true,  kernel, xbins, ybins, zbins, (doMirror==1));
      plotEff_customKDE_uniformBin(q2Bin, false, kernel, xbins, ybins, zbins, (doMirror==1));
    }
  }
  if (q2Bin == -1) {
    cout<<"Plotting efficiency for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	plotEff_customKDE_uniformBin(q2Bin, (tagFlag==1), kernel, xbins, ybins, zbins, (doMirror==1));
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	plotEff_customKDE_uniformBin(q2Bin, true,  kernel, xbins, ybins, zbins, (doMirror==1));
	plotEff_customKDE_uniformBin(q2Bin, false, kernel, xbins, ybins, zbins, (doMirror==1));
      }
    }
  }
  
}

void plotEff_customKDE_uniformBin(int q2Bin, bool tagFlag, int kernel, int xbins, int ybins, int zbins, bool doMirror)
{

  string shortString = Form(tagFlag?"b%ict":"b%iwt",q2Bin);
  string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);
  string confString  = "effKDEun_"+shortString+Form((doMirror?"_kn%i_%i_%i_%i":"_NoMirr_kn%i_%i_%i_%i"),kernel,xbins,ybins,zbins);
  int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // Load datasets
  TFile* fin_data = TFile::Open( ("effDataset_"+shortString+".root").c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: effDataset_"+shortString+".root"<<endl;
    return;
  }
  RooWorkspace* ws_data = (RooWorkspace*)fin_data->Get(("ws_"+shortString).c_str());
  if ( !ws_data || ws_data->IsZombie() ) {
    cout<<"Workspace not found in file: effDataset_"+shortString+".root"<<endl;
    return;
  }
  RooDataSet* data    = (RooDataSet*)ws_data->data(("data_"   +shortString).c_str());
  RooDataSet* numData = (RooDataSet*)ws_data->data(("numData_"+shortString).c_str());
  RooRealVar* ctK = ws_data->var("ctK");
  RooRealVar* ctL = ws_data->var("ctL");
  RooRealVar* phi = ws_data->var("phi");
  RooArgList vars (* ctK,* ctL,* phi);

  // open file with efficiency and import efficiency function
  TFile* fin = new TFile( (confString+".root").c_str(), "READ" );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<confString<<".root"<<endl;
    return;
  }
  RooWorkspace* ws = (RooWorkspace*)fin->Get("ws");
  if ( !ws || ws->IsZombie() ) {
    cout<<"Workspace not found in file: "<<confString<<".root"<<endl;
    return;
  } 
  TH3D* denHist = (TH3D*)ws->obj("denHist");
  TH3D* numHist = (TH3D*)ws->obj("numHist");
  if ( !denHist || !numHist || denHist->IsZombie() || numHist->IsZombie() ) {
    cout<<"Efficiency histograms not found in file: "<<confString<<".root"<<endl;
    return;
  } 
  TH3D* effHist = new TH3D(*numHist);
  effHist->Divide(denHist);
  RooDataHist* effData = new RooDataHist("effData","effData",vars,effHist);
  RooAbsReal* eff = new RooHistFunc("eff","eff",vars,*effData);

  // Plot efficiency projections
  gStyle->SetOptStat(0);
  auto h3_xyz = (TH3D*)eff->createHistogram("ctK,ctL,phi",50,50,50);
  auto h3_x  = h3_xyz->ProjectionX();
  auto h3_y  = h3_xyz->ProjectionY();
  auto h3_z  = h3_xyz->ProjectionZ();
  auto h3_xy = h3_xyz->Project3D("xy");
  auto h3_xz = h3_xyz->Project3D("xz");
  auto h3_yz = h3_xyz->Project3D("yz");

  // 1D projections
  c1[confIndex] = new TCanvas(Form("c1-%i",confIndex),("Efficiency 1D projections - "+longString).c_str(),1800,800);
  h3_x->SetTitle(("Efficiency cos(#theta_{K}) projection - "+longString).c_str());
  h3_y->SetTitle(("Efficiency cos(#theta_{L}) projection - "+longString).c_str());
  h3_z->SetTitle(("Efficiency #phi projection - "           +longString).c_str());
  h3_x->SetLineColor(2);
  h3_y->SetLineColor(2);
  h3_z->SetLineColor(2);
  h3_x->SetMinimum(0.0);
  h3_y->SetMinimum(0.0);
  h3_z->SetMinimum(0.0);
  c1[confIndex]->Divide(3,1);
  c1[confIndex]->cd(1);
  h3_x->Draw();
  c1[confIndex]->cd(2);
  h3_y->Draw();
  c1[confIndex]->cd(3);
  h3_z->Draw();
  c1[confIndex]->SaveAs( (confString+"_1DProj.pdf").c_str() );

  // 2D projections
  c3[confIndex] = new TCanvas(Form("c3-%i",confIndex),("Efficiency 2D projections - "+longString).c_str(),1800,800);
  h3_xy->SetTitle(("Efficiency cos(#theta_{K})/cos(#theta_{L}) projection - "+longString).c_str());
  h3_xz->SetTitle(("Efficiency cos(#theta_{K})/#phi projection - "           +longString).c_str());
  h3_yz->SetTitle(("Efficiency cos(#theta_{L})/#phi projection - "           +longString).c_str());
  h3_xy->GetXaxis()->SetTitleOffset(1.4);
  h3_xy->GetYaxis()->SetTitleOffset(2);
  h3_xz->GetXaxis()->SetTitleOffset(1.4);
  h3_xz->GetYaxis()->SetTitleOffset(2);
  h3_yz->GetXaxis()->SetTitleOffset(1.4);
  h3_yz->GetYaxis()->SetTitleOffset(2);
  h3_xy->SetMinimum(0.0);
  h3_xz->SetMinimum(0.0);
  h3_yz->SetMinimum(0.0);
  c3[confIndex]->Divide(3,1);
  c3[confIndex]->cd(1);
  h3_xy->Draw("SURF3");
  c3[confIndex]->cd(2);
  h3_xz->Draw("SURF3");
  c3[confIndex]->cd(3);
  h3_yz->Draw("SURF3");
  c3[confIndex]->SaveAs( (confString+"_2DProj.pdf").c_str() );

  // create the weighted dataset for GEN events
  RooAbsReal* effVal = (RooAbsReal*)data->addColumn(*eff);
  RooDataSet* wdata = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,eff->GetName());

  // Plot projections for closure test (RECO vs. eff*GEN)
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Closure test - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title((longString+" cos(#theta_{K}) distributions").c_str()));
  RooPlot* yframe = ctL->frame(Title((longString+" cos(#theta_{L}) distributions").c_str()));
  RooPlot* zframe = phi->frame(Title((longString+" #phi distributions").c_str()));
  wdata->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2),Name("plDenDist"));
  wdata->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
  wdata->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),DataError(RooAbsData::SumW2));
  numData->plotOn(xframe,Binning(30),Name("plNumDist"));
  numData->plotOn(yframe,Binning(30));
  numData->plotOn(zframe,Binning(30));
  xframe->GetYaxis()->SetTitleOffset(1.6);
  yframe->GetYaxis()->SetTitleOffset(1.6);
  zframe->GetYaxis()->SetTitleOffset(1.6);
  xframe->SetMaximum(xframe->GetMaximum()*1.15);
  yframe->SetMaximum(yframe->GetMaximum()*1.15);
  zframe->SetMaximum(zframe->GetMaximum()*1.15);
  leg->SetTextSize(0.03);
  leg->AddEntry(xframe->findObject("plNumDist"),"Post-selection RECO distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plDenDist"),"Efficiency-corrected GEN distribution" ,"lep");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.17); 
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.17); 
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.17); 
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("closure_"+confString+".pdf").c_str() );

}

