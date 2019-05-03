#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooPlot.h"

using namespace RooFit ;
using namespace std ;

static const int nBins = 9;

void composeEff_customKDE_uniformBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins, int kernel, int divx, int divy, int divz, int ndiv, bool doMirror);

void fillHists(TH3D* denHist, TH3D* numHist, RooDataSet* data, RooDataSet* numData, int nev, int xbin_min, int xbin_max, int ybin_min, int ybin_max, int zbin_min, int zbin_max, bool doMirror);

void composeEff_customKDE_uniform_parSub(int q2Bin, int tagFlag, int kernel = 1000, int doMirror = 1, int xbins=25, int ybins = 0, int zbins = 0, int ndiv = 0, int divx = 1, int divy = 0, int divz = 0)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively

  if ( ybins<1 ) ybins = xbins;
  if ( zbins<1 ) zbins = xbins;
  if ( xbins<1 ) return;

  if ( divy<1 ) divy = divx;
  if ( divz<1 ) divz = divx;
  if ( divx<1 ) return;

  if ( ndiv<0 || ndiv>=divx*divy*divz ) return;

  if ( doMirror<0 || doMirror>1 ) return;

  if ( q2Bin > -1 && q2Bin < nBins ) {
    if (tagFlag < 2 && tagFlag > -1) {
      cout<<"Computing efficiency for q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
      composeEff_customKDE_uniformBin(q2Bin, (tagFlag==1), xbins, ybins, zbins, kernel, divx, divy, divz, ndiv, (doMirror==1));
    }
    if (tagFlag == 2) {
      cout<<"Computing efficiency for q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
      composeEff_customKDE_uniformBin(q2Bin, true,  xbins, ybins, zbins, kernel, divx, divy, divz, ndiv, (doMirror==1));
      composeEff_customKDE_uniformBin(q2Bin, false, xbins, ybins, zbins, kernel, divx, divy, divz, ndiv, (doMirror==1));
    }
  }
  if (q2Bin == -1) {
    cout<<"Computing efficiency for all q2 bins"<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if (tagFlag < 2 && tagFlag > -1) {
	cout<<endl<<"q2 bin "<<q2Bin<<(tagFlag==1?" correctly":" wrongly")<<" tagged events"<<endl;
	composeEff_customKDE_uniformBin(q2Bin, (tagFlag==1), xbins, ybins, zbins, kernel, divx, divy, divz, ndiv, (doMirror==1));
      }
      if (tagFlag == 2) {
	cout<<endl<<"q2 bin "<<q2Bin<<" correctly and wrongly tagged events"<<endl;
	composeEff_customKDE_uniformBin(q2Bin, true,  xbins, ybins, zbins, kernel, divx, divy, divz, ndiv, (doMirror==1));
	composeEff_customKDE_uniformBin(q2Bin, false, xbins, ybins, zbins, kernel, divx, divy, divz, ndiv, (doMirror==1));
      }
    }
  }
  
}

void composeEff_customKDE_uniformBin(int q2Bin, bool tagFlag, int xbins, int ybins, int zbins, int kernel, int divx, int divy, int divz, int ndiv, bool doMirror)
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
  int xpos = ndiv % divx;
  int ypos = (ndiv/divx) % divy;
  int zpos = (ndiv/divx/divy) % divz;

  int xbin_min = (xpos*xbins)/divx+1;
  int xbin_max = ((xpos+1)*xbins)/divx;
  int ybin_min = (ypos*ybins)/divy+1;
  int ybin_max = ((ypos+1)*ybins)/divy;
  int zbin_min = (zpos*zbins)/divz+1;
  int zbin_max = ((zpos+1)*zbins)/divz;

  TH3D* denHist  = new TH3D(Form("denHist_%i" ,ndiv),Form("denHist_%i" ,ndiv),xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());
  TH3D* numHist  = new TH3D(Form("numHist_%i" ,ndiv),Form("numHist_%i" ,ndiv),xbins,-1,1,ybins,-1,1,zbins,-TMath::Pi(),TMath::Pi());

  fillHists(denHist,numHist,data,numData,kernel,xbin_min,xbin_max,ybin_min,ybin_max,zbin_min,zbin_max,doMirror);
  denHist->Sumw2();
  numHist->Sumw2();

  // save histograms to file
  RooWorkspace *wsp_out = new RooWorkspace("ws","Workspace with efficiency parameterisation");
  wsp_out->import( *denHist );
  wsp_out->import( *numHist );
  wsp_out->import( vars, Silence() );
  wsp_out->writeToFile( ( "effKDEun_"+shortString+Form((doMirror?"_kn%i_%i_%i_%i_%i-frac-%i-%i-%i.root":"_NoMirr_kn%i_%i_%i_%i_%i-frac-%i-%i-%i.root"),kernel,xbins,ybins,zbins,ndiv,divx,divy,divz)).c_str() );
  
}

void fillHists(TH3D* denHist, TH3D* numHist, RooDataSet* data, RooDataSet* numData, int nev, int xbin_min, int xbin_max, int ybin_min, int ybin_max, int zbin_min, int zbin_max, bool doMirror)
{
  // number of grid points to be probed
  int nBins = (xbin_max-xbin_min+1) * (ybin_max-ybin_min+1) * (zbin_max-zbin_min+1);

  // set (the square of) a standard distance unit,
  // as the radius of a sphere containing the needed amount of events if the numerator distribution would be flat
  double AvgRadSq = TMath::Power( 6.0*nev/TMath::Pi()/TMath::Pi()/numData->sumEntries() , 2.0/3 );
  double maxDistSq = AvgRadSq*20;
  int nDistBins = 200;

  // create three arrays with bin coordinates
  double* xBinCenter = new double[denHist->GetNbinsX()];
  double* yBinCenter = new double[denHist->GetNbinsY()];
  double* zBinCenter = new double[denHist->GetNbinsZ()];
  denHist->GetXaxis()->GetCenter(xBinCenter);
  denHist->GetYaxis()->GetCenter(yBinCenter);
  denHist->GetZaxis()->GetCenter(zBinCenter);

  if ( 4*maxDistSq < pow(xBinCenter[xbin_max-1]-xBinCenter[xbin_min-1],2) + pow(yBinCenter[ybin_max-1]-yBinCenter[ybin_min-1],2) + pow(zBinCenter[zbin_max-1]-zBinCenter[zbin_min-1],2) )
    maxDistSq = 0.25 * pow(xBinCenter[xbin_max-1]-xBinCenter[xbin_min-1],2) + pow(yBinCenter[ybin_max-1]-yBinCenter[ybin_min-1],2) + pow(zBinCenter[zbin_max-1]-zBinCenter[zbin_min-1],2);

  // for each point of the grid, a 1D histo is defined for the distance distribution of data points from the bin center
  vector < TH1F* > distNum;
  vector < TH1F* > distDen;

  int iBin;
  for (iBin=0; iBin<nBins; ++iBin) {
    distNum.push_back( new TH1F(Form("distNum%i",iBin),Form("distNum%i",iBin),nDistBins,0,maxDistSq) );
    distDen.push_back( new TH1F(Form("distDen%i",iBin),Form("distDen%i",iBin),nDistBins,0,maxDistSq) );
  }

  int iBinX, iBinY, iBinZ;
  double xVal, yVal, zVal;

  cout<<"Denominator preparation"<<endl;
  int counter=0;
  double distSq;
  // loop over denominator dataset
  for (int iEv=0; iEv<data->sumEntries(); ++iEv) {
    if ( iEv > counter ) {
      cout<<counter*100/data->sumEntries()<<"%"<<endl;
      counter+=data->sumEntries()/10;
    }
    
    // get data-point coordinates
    const RooArgSet *iPoint = data->get(iEv);
    xVal = iPoint->getRealValue("ctK");
    yVal = iPoint->getRealValue("ctL");
    zVal = iPoint->getRealValue("phi");

    if ( (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) ) continue;

    // loop over grid points to be evaluated
    for (iBinX=xbin_min; iBinX<=xbin_max; ++iBinX)
      for (iBinY=ybin_min; iBinY<=ybin_max; ++iBinY)
	for (iBinZ=zbin_min; iBinZ<=zbin_max; ++iBinZ) {
	  
	  iBin = (iBinX-xbin_min) + (xbin_max-xbin_min+1) * ( (iBinY-ybin_min) + (ybin_max-ybin_min+1) * (iBinZ-zbin_min) );
	  if (iBin<0 || iBin>nBins-1) {cout<<"ERROR, bin index too high: "<<iBin<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") vs. a max of "<<nBins<<endl; return;}

	  distSq = pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2);

	  // discard points that would fall in overflow of the distance histogram
	  if ( distSq > maxDistSq ) continue;

	  // fill the distance histogram
	  distDen[iBin]->Fill( distSq );
	  if (doMirror) {
	    // mirrored entries
	    distDen[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]-2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distDen[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]+2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]-2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]+2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()-2,2) );
	    distDen[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()+2,2) );
	  }

	}

  }

  cout<<"Numerator preparation"<<endl;
  counter=0;
  int tmpcnt=0;
  // loop over denominator dataset
  for (int iEv=0; iEv<numData->sumEntries(); ++iEv) {
    if ( iEv > counter ) {
      cout<<counter*100/numData->sumEntries()<<"%"<<endl;
      counter+=numData->sumEntries()/10;
    }
    
    // get data-point coordinates
    const RooArgSet *iPoint = numData->get(iEv);
    xVal = iPoint->getRealValue("ctK");
    yVal = iPoint->getRealValue("ctL");
    zVal = iPoint->getRealValue("phi");

    if ( (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_min-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_min-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_min-1])/TMath::Pi(),2) > maxDistSq ) &&
	 (pow(xVal-xBinCenter[xbin_max-1],2) + pow(yVal-yBinCenter[ybin_max-1],2) + pow((zVal-zBinCenter[zbin_max-1])/TMath::Pi(),2) > maxDistSq ) ) continue;
    ++tmpcnt;
    // loop over grid points to be evaluated
    for (iBinX=xbin_min; iBinX<=xbin_max; ++iBinX)
      for (iBinY=ybin_min; iBinY<=ybin_max; ++iBinY)
	for (iBinZ=zbin_min; iBinZ<=zbin_max; ++iBinZ) {
	  
	  iBin = (iBinX-xbin_min) + (xbin_max-xbin_min+1) * ( (iBinY-ybin_min) + (ybin_max-ybin_min+1) * (iBinZ-zbin_min) );
	  if (iBin<0 || iBin>nBins-1) {cout<<"ERROR, bin index too high: "<<iBin<<" ("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") vs. a max of "<<nBins<<endl; return;}

	  distSq = pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2);

	  // discard points that would fall in overflow of the distance histogram
	  if ( distSq > maxDistSq ) continue;
	  
	  // fill the distance histogram
	  distNum[iBin]->Fill( distSq );
	  if (doMirror) {
	    // mirrored entries
	    distNum[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]-2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distNum[iBin]->Fill( pow(xVal+xBinCenter[iBinX-1]+2,2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]-2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal+yBinCenter[iBinY-1]+2,2) + pow((zVal-zBinCenter[iBinZ-1])/TMath::Pi(),2) );
	    distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()-2,2) );
	    distNum[iBin]->Fill( pow(xVal-xBinCenter[iBinX-1],2) + pow(yVal-yBinCenter[iBinY-1],2) + pow((zVal+zBinCenter[iBinZ-1])/TMath::Pi()+2,2) );
	  }

	}

  }
  cout<<"++"<<tmpcnt<<" "<<sqrt(maxDistSq)<<endl;
  // define minimal kernel radius as the bin emi-width
  int minRadius = distDen[0]->FindBin(pow(denHist->GetXaxis()->GetBinWidth(1)/2.0,2));

  double numInt, denInt;
  // loop over grid points, to evaluate numerator and denominator values for each
  for (iBinX=xbin_min; iBinX<=xbin_max; ++iBinX)
    for (iBinY=ybin_min; iBinY<=ybin_max; ++iBinY)
      for (iBinZ=zbin_min; iBinZ<=zbin_max; ++iBinZ) {

	// define global bin number, and the one 'reduced' to the grid to be evaluated
	iBin = denHist->GetBin(iBinX,iBinY,iBinZ);
	int binCount = (iBinX-xbin_min) + (xbin_max-xbin_min+1) * ( (iBinY-ybin_min) + (ybin_max-ybin_min+1) * (iBinZ-zbin_min) );

	int maxBin = minRadius;
	// determine kernel radius, applying conditions on numerator and denominator counts
	while ( (maxBin < nDistBins/4 || (maxBin < nDistBins && distNum[binCount]->Integral(1,maxBin) < 4)) &&
		(distDen[binCount]->Integral(1,maxBin) < 400 || distNum[binCount]->Integral(1,maxBin) < nev)
		) ++maxBin;

	cout<<"("<<iBinX<<" "<<iBinY<<" "<<iBinZ<<") range: 1->"<<maxBin<<
	  " num="<<distNum[binCount]->Integral(1,maxBin)<<" den="<<distDen[binCount]->Integral(1,maxBin)<<" - "<<distNum[binCount]->GetBinCenter(maxBin)<<endl;

	// fill numerator and denominator description with event counts inside the kernel
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

	// apply poissonian statistics
	numHist->SetBinError(iBin, sqrt(numHist->GetBinContent(iBin)));
	denHist->SetBinError(iBin, sqrt(denHist->GetBinContent(iBin)));
	  
      }

}
