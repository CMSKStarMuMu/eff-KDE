#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

TCanvas* cx1 [2*nBins];
TCanvas* cy1 [2*nBins];
TCanvas* cz1 [2*nBins];

void composeEff_customKDE_uniformBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins, int kernel, bool plot);

void fillHists(TH3D* denHist, TH3D* numHist, RooDataSet* data, RooDataSet* numData, int nev);

void composeEff_customKDE_uniform(int q2Bin, int tagFlag, int kernel = 100, int xbins=50, int ybins = 0, int zbins = 0, bool plot = true)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Computing efficiency for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      composeEff_customKDE_uniformBin(q2Bin, (tagFlag==1), xbins, ybins, zbins, kernel, plot);
    }
    if (tagFlag == 2) {
      cout<<"Computing efficiency for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      composeEff_customKDE_uniformBin(q2Bin, true,  xbins, ybins, zbins, kernel, plot);
      composeEff_customKDE_uniformBin(q2Bin, false, xbins, ybins, zbins, kernel, plot);
    }
  }
  if (q2Bin == -1) {
    cout<<"Computing efficiency for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	composeEff_customKDE_uniformBin(q2Bin, (tagFlag==1), xbins, ybins, zbins, kernel, plot);
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	composeEff_customKDE_uniformBin(q2Bin, true,  xbins, ybins, zbins, kernel, plot);
	composeEff_customKDE_uniformBin(q2Bin, false, xbins, ybins, zbins, kernel, plot);
      }
    }
  }
  
}

void composeEff_customKDE_uniformBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins, int kernel, bool plot)
{

  string shortString = Form(tagFlag?"b%ict":"b%iwt",q2Bin);
  string longString  = Form(tagFlag?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",q2Bin);
  int confIndex = (tagFlag?q2Bin:q2Bin+nBins);

  // Load datasets
  TFile* fin = TFile::Open( ("effDataset_"+shortString+".root").c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: effDataset_"+shortString+".root"<<endl;
    return;
  }
  RooWorkspace* wsp = (RooWorkspace*)fin->Get(("ws_"+shortString).c_str());
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: effDataset_"+shortString+".root"<<endl;
    return;
  }
  RooDataSet* data    = (RooDataSet*)wsp->data(("data_"   +shortString).c_str());
  RooDataSet* numData = (RooDataSet*)wsp->data(("numData_"+shortString).c_str());
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  RooArgSet vars (* ctK,* ctL,* phi);

  // create numerator and denominator histograms
  // TH3D* denHist = (TH3D*)data   ->createHistogram( ("denHist"+shortString).c_str(),
  // 						   *ctK,     Binning(xbins,-1,1) ,
  // 						   YVar(*ctL,Binning(ybins,-1,1)),
  // 						   ZVar(*phi,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
  // TH3D* numHist = (TH3D*)numData->createHistogram( ("numHist"+shortString).c_str(),
  // 						   *ctK,     Binning(xbins,-1,1) ,
  // 						   YVar(*ctL,Binning(ybins,-1,1)),
  // 						   ZVar(*phi,Binning(zbins,-TMath::Pi(),TMath::Pi())) );
  TH3D* denHist = new TH3D("denHist","denHist",xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  TH3D* numHist = new TH3D("numHist","numHist",xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  fillHists(denHist,numHist,data,numData,kernel);
  denHist->Sumw2();
  numHist->Sumw2();

  // save histograms to file
  RooWorkspace *wsp_out = new RooWorkspace("ws","Workspace with efficiency parameterisation");
  wsp_out->import( *denHist );
  wsp_out->import( *numHist );
  wsp_out->import( vars, Silence() );
  wsp_out->writeToFile( ( "effKDE_"+shortString+Form("_k%i_%i_%i_%i.root",kernel,xbins,ybins,zbins)).c_str() );

  if (plot) {
    // Plot 1D slices of the efficiency function and binned efficiency
    vector <TEfficiency*> effHistsX; 
    vector <TEfficiency*> effHistsY;
    vector <TEfficiency*> effHistsZ;
    cx1[confIndex] = new TCanvas(("cx1"+shortString).c_str(),("KDE efficiency - "+longString+" - cos(theta_k) slices").c_str(),1500,1500) ;
    cy1[confIndex] = new TCanvas(("cy1"+shortString).c_str(),("KDE efficiency - "+longString+" - cos(theta_l) slices").c_str(),1500,1500) ;
    cz1[confIndex] = new TCanvas(("cz1"+shortString).c_str(),("KDE efficiency - "+longString+" - phi slices"         ).c_str(),1500,1500) ;
    cx1[confIndex]->Divide(5,5);
    cy1[confIndex]->Divide(5,5);
    cz1[confIndex]->Divide(5,5);

    // width of the slices in the hidden variables ("border" is half of it)
    double border = 0.005;
    // variables to be filled with global efficiency maximum
    double maxEffX = 0;
    double maxEffY = 0;
    double maxEffZ = 0;

    // loop over slice grid
    for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {

	// central values and borders of the slices in the hidden variables
	double centA = -0.8 + 1.6*i/4;
	double centB = -0.8 + 1.6*j/4;
	double lowA  = TMath::Max( centA - border,  1e-4-1 );
	double lowB  = TMath::Max( centB - border,  1e-4-1 );
	double highA = TMath::Min( centA + border, -1e-4+1 );
	double highB = TMath::Min( centB + border, -1e-4+1 );

	// slicing num and den distributions    
	auto numProjX = numHist->ProjectionX("numProjX", 
					     numHist->GetYaxis()->FindBin(lowA            ), numHist->GetYaxis()->FindBin(highA            ),
					     numHist->GetZaxis()->FindBin(lowB*TMath::Pi()), numHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto numProjY = numHist->ProjectionY("numProjY", 
					     numHist->GetXaxis()->FindBin(lowA            ), numHist->GetXaxis()->FindBin(highA            ),
					     numHist->GetZaxis()->FindBin(lowB*TMath::Pi()), numHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto numProjZ = numHist->ProjectionZ("numProjZ", 
					     numHist->GetXaxis()->FindBin(lowA            ), numHist->GetXaxis()->FindBin(highA            ),
					     numHist->GetYaxis()->FindBin(lowB            ), numHist->GetYaxis()->FindBin(highB            ),"e");
	auto denProjX = denHist->ProjectionX("denProjX", 
					     denHist->GetYaxis()->FindBin(lowA            ), denHist->GetYaxis()->FindBin(highA            ),
					     denHist->GetZaxis()->FindBin(lowB*TMath::Pi()), denHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto denProjY = denHist->ProjectionY("denProjY", 
					     denHist->GetXaxis()->FindBin(lowA            ), denHist->GetXaxis()->FindBin(highA            ),
					     denHist->GetZaxis()->FindBin(lowB*TMath::Pi()), denHist->GetZaxis()->FindBin(highB*TMath::Pi()),"e");
	auto denProjZ = denHist->ProjectionZ("denProjZ", 
					     denHist->GetXaxis()->FindBin(lowA            ), denHist->GetXaxis()->FindBin(highA            ),
					     denHist->GetYaxis()->FindBin(lowB            ), denHist->GetYaxis()->FindBin(highB            ),"e");

	// producing 1D efficiencies from the slices
	effHistsX.push_back( new TEfficiency(*numProjX,*denProjX) );
	effHistsX.back()->SetName( Form("effHistX_%i_%i",i,j) );
	effHistsX.back()->SetTitle( ("Efficiency - "+longString+Form(" - slice ctL=%1.2f phi=%1.2f;cos(#theta_{K});Efficiency",centA,centB*TMath::Pi())).c_str() );
    
	effHistsY.push_back( new TEfficiency(*numProjY,*denProjY) );
	effHistsY.back()->SetName( Form("effHistY_%i_%i",i,j) );
	effHistsY.back()->SetTitle( ("Efficiency - "+longString+Form(" - slice ctK=%1.2f phi=%1.2f;cos(#theta_{L});Efficiency",centA,centB*TMath::Pi())).c_str() );

	effHistsZ.push_back( new TEfficiency(*numProjZ,*denProjZ) );
	effHistsZ.back()->SetName( Form("effHistZ_%i_%i",i,j) );
	effHistsZ.back()->SetTitle( ("Efficiency - "+longString+Form(" - slice ctK=%1.2f ctL=%1.2f;#phi;Efficiency",centA,centB)).c_str() );

	// plot in canvas
	cx1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effHistsX.back()->Draw();
	cx1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphx = effHistsX.back()->GetPaintedGraph(); 
	graphx->SetMinimum(0);
	auto effValsX = graphx->GetY();
	for (int iBin=0; iBin<graphx->GetN(); ++iBin) if (maxEffX<effValsX[iBin]) maxEffX = effValsX[iBin];
	graphx->GetYaxis()->SetTitleOffset(1.7);
	cx1[confIndex]->cd(5*j+i+1)->Update();

	cy1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effHistsY.back()->Draw();
	cy1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphy = effHistsY.back()->GetPaintedGraph(); 
	graphy->SetMinimum(0);
	auto effValsY = graphy->GetY();
	for (int iBin=0; iBin<graphy->GetN(); ++iBin) if (maxEffY<effValsY[iBin]) maxEffY = effValsY[iBin];
	graphy->GetYaxis()->SetTitleOffset(1.7);
	cy1[confIndex]->cd(5*j+i+1)->Update();

	cz1[confIndex]->cd(5*j+i+1);
	gPad->SetLeftMargin(0.18);
	effHistsZ.back()->Draw();
	cz1[confIndex]->cd(5*j+i+1)->Update(); 
	auto graphz = effHistsZ.back()->GetPaintedGraph(); 
	graphz->SetMinimum(0);
	auto effValsZ = graphz->GetY();
	for (int iBin=0; iBin<graphz->GetN(); ++iBin) if (maxEffZ<effValsZ[iBin]) maxEffZ = effValsZ[iBin];
	graphz->GetYaxis()->SetTitleOffset(1.7);
	cz1[confIndex]->cd(5*j+i+1)->Update();

      }    

    // set uniform y-axis ranges
    for (int i=0; i<effHistsX.size(); ++i) (effHistsX[i]->GetPaintedGraph())->SetMaximum(maxEffX*1.1);
    for (int i=0; i<effHistsY.size(); ++i) (effHistsY[i]->GetPaintedGraph())->SetMaximum(maxEffY*1.1);
    for (int i=0; i<effHistsZ.size(); ++i) (effHistsZ[i]->GetPaintedGraph())->SetMaximum(maxEffZ*1.1);
    for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {
	cx1[confIndex]->cd(5*j+i+1)->Update();
	cy1[confIndex]->cd(5*j+i+1)->Update();
	cz1[confIndex]->cd(5*j+i+1)->Update();
      }	

    cx1[confIndex]->SaveAs( ("effKDE_"+shortString+Form("_k%i_%i_%i_%i_CTKslices_dp%i.pdf",kernel,xbins,ybins,zbins,(int)(border*200))).c_str() );
    cy1[confIndex]->SaveAs( ("effKDE_"+shortString+Form("_k%i_%i_%i_%i_CTLslices_dp%i.pdf",kernel,xbins,ybins,zbins,(int)(border*200))).c_str() );
    cz1[confIndex]->SaveAs( ("effKDE_"+shortString+Form("_k%i_%i_%i_%i_PHIslices_dp%i.pdf",kernel,xbins,ybins,zbins,(int)(border*200))).c_str() );
  }
  
}

void fillHists(TH3D* denHist, TH3D* numHist, RooDataSet* data, RooDataSet* numData, int nev)
{
  int nBins = denHist->GetNbinsX() * denHist->GetNbinsY() * denHist->GetNbinsZ();
  double AvgRadSq = TMath::Power( 6.0*nev/TMath::Pi()/TMath::Pi()/data->sumEntries() , 2.0/3 );

  vector < TH1F* > distNum;
  vector < TH1F* > distDen;

  int iBin;

  for (iBin=0; iBin<nBins; ++iBin) {
    distNum.push_back( new TH1F(Form("distNum%i",iBin),Form("distNum%i",iBin),200,0,AvgRadSq*20) );
    distDen.push_back( new TH1F(Form("distDen%i",iBin),Form("distDen%i",iBin),200,0,AvgRadSq*20) );
  }

  double* xBinCenter = new double[denHist->GetNbinsX()];
  double* yBinCenter = new double[denHist->GetNbinsY()];
  double* zBinCenter = new double[denHist->GetNbinsZ()];
  denHist->GetXaxis()->GetCenter(xBinCenter);
  denHist->GetYaxis()->GetCenter(yBinCenter);
  denHist->GetZaxis()->GetCenter(zBinCenter);

  int iBinX, iBinY, iBinZ;
  double xVal, yVal, zVal;

  cout<<"Denominator preparation"<<endl;
  int counter=0;
  for (int iEv=0; iEv<data->sumEntries(); ++iEv) {
    if ( iEv > counter ) {
      cout<<counter*100/data->sumEntries()<<"%"<<endl;
      counter+=data->sumEntries()/10;
    }
    
    const RooArgSet *iPoint = data->get(iEv);
    xVal = iPoint->getRealValue("ctK");
    yVal = iPoint->getRealValue("ctL");
    zVal = iPoint->getRealValue("phi");

    for (iBinX=1; iBinX<=denHist->GetNbinsX(); ++iBinX)
      for(iBinY=1; iBinY<=denHist->GetNbinsY(); ++iBinY)
	for(iBinZ=1; iBinZ<=denHist->GetNbinsZ(); ++iBinZ) {
	  
	  iBin = (iBinX-1) + denHist->GetNbinsX() * ( (iBinY-1) + denHist->GetNbinsY() * (iBinZ-1) );
	  // iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	  if (iBin<0 || iBin>nBins-1) {cout<<"ERROR, bin index too high: "<<iBin<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") vs. a max of "<<nBins<<endl; return;}
	  if ( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) > AvgRadSq*10 ) continue;

	  distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          // mirrored entries
          distDen[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]-2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]+2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]-2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]+2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()-2,2) );
          distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()+2,2) );

	}

  }

  cout<<"Numerator preparation"<<endl;
  counter=0;
  for (int iEv=0; iEv<numData->sumEntries(); ++iEv) {
    if ( iEv > counter ) {
      cout<<counter*100/numData->sumEntries()<<"%"<<endl;
      counter+=numData->sumEntries()/10;
    }
    
    const RooArgSet *iPoint = numData->get(iEv);
    xVal = iPoint->getRealValue("ctK");
    yVal = iPoint->getRealValue("ctL");
    zVal = iPoint->getRealValue("phi");

    for (iBinX=1; iBinX<=denHist->GetNbinsX(); ++iBinX)
      for (iBinY=1; iBinY<=denHist->GetNbinsY(); ++iBinY)
	for (iBinZ=1; iBinZ<=denHist->GetNbinsZ(); ++iBinZ) {
	  
	  iBin = (iBinX-1) + denHist->GetNbinsX() * ( (iBinY-1) + denHist->GetNbinsY() * (iBinZ-1) );
	  // iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	  if (iBin<0 || iBin>nBins-1) {cout<<"ERROR, bin index too high: "<<iBin<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") vs. a max of "<<nBins<<endl; return;}
	  if ( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) > AvgRadSq*10 ) continue;
	  
	  distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          // mirrored entries
          distNum[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]-2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]+2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]-2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]+2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()-2,2) );
          distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()+2,2) );
	  
	}

  }

  double numInt, denInt;

  for (iBinX=1; iBinX<=denHist->GetNbinsX(); ++iBinX)
    for (iBinY=1; iBinY<=denHist->GetNbinsY(); ++iBinY)
      for (iBinZ=1; iBinZ<=denHist->GetNbinsZ(); ++iBinZ) {

	iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	int binCount = (iBinX-1) + denHist->GetNbinsX() * ( (iBinY-1) + denHist->GetNbinsY() * (iBinZ-1) );

	int maxBin = 1;
	while ( maxBin < distDen[binCount]->GetNbinsX() && (distDen[binCount]->Integral(1,maxBin) < nev || distNum[binCount]->Integral(1,maxBin) < 10) ) ++maxBin;

	// cout<<"("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") range: 1->"<<maxBin<<
	//   " num="<<distNum[binCount]->Integral(1,maxBin)<<" den="<<distDen[binCount]->Integral(1,maxBin)<<endl;

	numInt = distNum[binCount]->Integral(1,maxBin);
	denInt = distDen[binCount]->Integral(1,maxBin);

	if ( denInt >= numInt ) {
	  numHist->SetBinContent(iBin, numInt);
	  denHist->SetBinContent(iBin, denInt);
	} else {
	  cout<<"ERROR, numerator count exceeding denominator count: "<<
	    numInt<<" vs. "<<denInt<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<")"<<endl;
	  numHist->SetBinContent(iBin, denInt);
	  denHist->SetBinContent(iBin, denInt);
	}

	numHist->SetBinError(iBin, sqrt(numHist->GetBinContent(iBin)));
	denHist->SetBinError(iBin, sqrt(denHist->GetBinContent(iBin)));
	  
      }

}
