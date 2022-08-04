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


void plotEffBin(int q2Bin, int parity, bool doClosure, int year, int vers)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  string confString = "plotEff_d/effKDEComparison_"+shortString+Form("_v%i",vers);
  gStyle->SetOptStat(0);

  bool doCT = true;
  bool doWT = true;

  int confIndex = nBins*parity + q2Bin;

  const char *varCoord[3];
  varCoord[0] = "X";
  varCoord[1] = "Y";
  varCoord[2] = "Z";

  const char *varNames[3];
  varNames[0] = "ctK";
  varNames[1] = "ctL";
  varNames[2] = "phi";

  // Load variables and dataset
  string filename_data = Form("effDataset_b%i_%i.root",q2Bin,year);
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

  // import KDE efficiency histograms
  vector<RooAbsReal*> effC (0);
  vector<RooAbsReal*> effW (0);
  double intValC = 0;
  double intValW = 0;
  for (int iVers=vers%10; iVers<40; iVers+=10) {
    string filename = Form((parity==0?"files/KDEeff_b%i_ev_%i_v%i.root":"files/KDEeff_b%i_od_%i_v%i.root"),q2Bin, year, iVers);
    TFile* fin = new TFile( filename.c_str(), "READ" );
    if ( !fin || !fin->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      continue;
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
    if (doCT) {
      if (effC.size()==0)
	intValC = effCHist->Integral();
      effCHist->Scale(intValC/effCHist->Integral());
      auto effCData = new RooDataHist("effCData","effCData",vars,effCHist);
      // sara: interpolation order set to 1 -> is it ok? enough?
      effC.push_back( new RooHistFunc("effC","effC",vars,*effCData,1) );
    }
    if (doWT) {
      if (effW.size()==0)
	intValW = effWHist->Integral();
      effWHist->Scale(intValW/effWHist->Integral());
      auto effWData = new RooDataHist("effWData","effWData",vars,effWHist);
      effW.push_back( new RooHistFunc("effW","effW",vars,*effWData,1) );
    }
  }

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
  vector <vector <RooPlot*>> fsC(3,std::vector<RooPlot*>(0) );
  vector <vector <RooPlot*>> fsW(3,std::vector<RooPlot*>(0) );

  // TLegend* leg = new TLegend (0.35,0.8,0.9,0.9);

  // variables to be filled with global efficiency maximum
  double maxEffC[3] = {0,0,0};
  double maxEffW[3] = {0,0,0};

  // loop over slice grid
  for (int i=0; i<5; ++i) for (int j=0; j<5; ++j) {
      cout<<i<<j<<endl;

      std::vector<string> ctCuts, wtCuts;
      // central values and borders of the slices in the hidden variables
      double centA = -0.8 + 1.6*i/4;
      double centB = -0.8 + 1.6*j/4;

      // producing 1D slices of efficiency description 
      for (int ivar = 0; ivar < 3; ivar++){
	cout<<ivar<<endl;
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
  	  fsC[ivar].push_back( vars_vec[ivar]->frame( Name( Form("fs%sC_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()) ),
						      Title( Form("Correct-tag efficiency (q2-bin %i, %s) slice %s=%1.2f %s=%1.2f;%s;Efficiency",q2Bin,(parity==0?"even":"odd"),(ivar==0?"ctL":"ctK"),centA,(ivar==2?"ctL":"phi"),centB*(ivar==2?1:TMath::Pi()),(ivar==0?"cos(#theta_{K})":(ivar==0?"cos(#theta_{L})":"#phi"))) )) );
	  for (uint iEff=0; iEff<effC.size(); ++iEff)
	    effC[iEff]->plotOn( fsC[ivar].back(), LineColor(2+iEff), Name(Form("effC%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()))) ;
        }
        if (doWT) {
  	  fsW[ivar].push_back( vars_vec[ivar]->frame( Name( Form("fs%sW_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()) ),
						      Title( Form("Mistag efficiency (q2-bin %i, %s) slice %s=%1.2f %s=%1.2f;%s;Efficiency",q2Bin,(parity==0?"even":"odd"),(ivar==0?"ctL":"ctK"),centA,(ivar==2?"ctL":"phi"),centB*(ivar==2?1:TMath::Pi()),(ivar==0?"cos(#theta_{K})":(ivar==0?"cos(#theta_{L})":"#phi"))) )) );
	  for (uint iEff=0; iEff<effW.size(); ++iEff)
	    effW[iEff]->plotOn( fsW[ivar].back(), LineColor(2+iEff), Name(Form("effW%s_%i_%i_%s",varCoord[ivar],i,j,shortString.c_str()))) ;
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
          fsC[ivar].back()->Draw();
          // checking maximum values
	  if ( maxEffC[ivar]<fsC[ivar].back()->GetMaximum() )
	    maxEffC[ivar] = fsC[ivar].back()->GetMaximum();
        }
        if (doWT) {
	  csW[ivar][confIndex]->cd(5*j+i+1);
          fsW[ivar].back()->Draw();
          // checking maximum values
	  if ( maxEffW[ivar]<fsW[ivar].back()->GetMaximum() )
	    maxEffW[ivar] = fsW[ivar].back()->GetMaximum();
        }
      }
  } // end of slicing

  //set  uniform y-axis ranges and save
  for (int ivar=0; ivar< 3; ivar++){
    if (doCT) {
      for (vector<RooPlot*>::iterator pad = fsC[ivar].begin(); pad != fsC[ivar].end(); ++pad) {
        (*pad)->SetMaximum(maxEffC[ivar]);
        (*pad)->SetMinimum(0);
      }
      csC[ivar][confIndex]->SaveAs( (confString+Form("_eff-ct_%sslices_comp_%i.pdf",varNames[ivar], year)).c_str() );
    }
    if (doWT) {
      for (vector<RooPlot*>::iterator pad = fsW[ivar].begin(); pad != fsW[ivar].end(); ++pad) {
        (*pad)->SetMaximum(maxEffW[ivar]);
        (*pad)->SetMinimum(0);
      }
      csW[ivar][confIndex]->SaveAs( (confString+Form("_eff-wt_%sslices_comp_%i.pdf",varNames[ivar], year)).c_str() );
    }
  }

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
  int vers   = -1;

  if ( argc > 1 ) q2Bin  = atoi(argv[1]);
  if ( argc > 2 ) parity = atoi(argv[2]);

  if ( q2Bin  < -1 || q2Bin  >= nBins ) return 1;
  if ( parity < -1 || parity > 1      ) return 1;

  bool doClosure = true;
  if ( argc > 3 && atoi(argv[3]) == 0 ) doClosure = false;
  if ( argc > 4 ) year = atoi(argv[4]);
  if ( argc > 5 ) vers = atoi(argv[5]);

  if ( q2Bin > -1 ) {
    if ( parity > -1 ) {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<(parity==1?" - odd":" - even")<<" events"<<endl;
      plotEffBin( q2Bin, parity, doClosure, year, vers );
    } else {
      cout<<"Plotting efficiency for q2 bin "<<q2Bin<<" - both event parities"<<endl;
      plotEffBin( q2Bin, 0, doClosure, year, vers );
      plotEffBin( q2Bin, 1, doClosure, year, vers );
    }
  } else {
    cout<<"Plotting efficiency for all q2 bins - "<<(parity==1?"odd events":(parity==0?"even events":"both event parities"))<<endl;
    for (q2Bin=0; q2Bin<nBins; ++q2Bin) {
      if ( parity > -1 ) plotEffBin( q2Bin, parity, doClosure, year, vers );
      else {
	plotEffBin( q2Bin, 0, doClosure, year, vers );
	plotEffBin( q2Bin, 1, doClosure, year, vers );
      }
    }
  }

  return 0;

}
