using namespace RooFit;
using namespace std;

// want to produce and save plots of distributions?
bool plot = true;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [2*nBins];

RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
RooArgSet vars (ctK, ctL, phi);

void createDataset(int q2Bin, int tagFlag)
{
  // q2-bin format: [0-8] for one bin, [-1] for each bin recursively
  // tag format:    [1] correctly-tagged. [0] mistagged, [2] each tag, recursively
  if ( q2Bin<-1 || q2Bin>=nBins ) {
    cout<<"Invalid bin number!"<<endl;
    return;
  }
  if ( tagFlag<0 || tagFlag>2 ) {
    cout<<"Invalid tag flag!"<<endl;
    return;
  }

  bool runBin [2*nBins];
  string shortString [2*nBins];
  string longString  [2*nBins];
  for (int i=0; i<2*nBins; ++i) {
    runBin [i] = false;
    if ( (q2Bin!=-1 && q2Bin!=i%nBins) || (tagFlag!=2 && tagFlag==i/nBins) ) continue;
    runBin [i] = true;
    shortString [i] = Form(i/nBins==0?"b%ict":"b%iwt",i%nBins);
    longString  [i] = Form(i/nBins==0?"q2 bin %i correct-tag":"q2 bin %i wrong-tag",i%nBins);
  }

  // Load ntuples
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/GEN_BFilter_B0MuMuKstar_p*.root/ntuple");
  t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/NtupleNov21/2016MC_LMNR_bdt0p96.root/ntuple");
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();
  int counter;

  double genCosThetaK, genCosThetaL, genPhi, genDimuMass, genB0pT, genB0eta;
  double recoCosThetaK, recoCosThetaL, recoPhi;
  double recoDimuMass, recoB0pT, recoB0eta, genSignal, tagB0;
  t_den->SetBranchAddress( "cos_theta_k" , &genCosThetaK );
  t_den->SetBranchAddress( "cos_theta_l" , &genCosThetaL );
  t_den->SetBranchAddress( "phi_kst_mumu", &genPhi       );
  t_den->SetBranchAddress( "genq2"       , &genDimuMass  );
  t_den->SetBranchAddress( "genbPt"      , &genB0pT      );
  t_den->SetBranchAddress( "genbEta"     , &genB0eta     );
  t_num->SetBranchAddress( "cos_theta_k" , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l" , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu", &recoPhi       );
  t_num->SetBranchAddress( "mumuMass"    , &recoDimuMass  );
  t_num->SetBranchAddress( "bPt"         , &recoB0pT      );
  t_num->SetBranchAddress( "bEta"        , &recoB0eta     );
  t_num->SetBranchAddress( "genSignal"   , &genSignal     );
  t_num->SetBranchAddress( "tagB0"       , &tagB0         );

  // Define datasets
  RooDataSet* data    [2*nBins];
  RooDataSet* numData [2*nBins];
  for (int i=0; i<2*nBins; ++i) if (runBin[i]) {
      data    [i] = new RooDataSet( ("data_"   +shortString[i]).c_str(), "GEN distribution before GEN-filter", vars );
      numData [i] = new RooDataSet( ("numData_"+shortString[i]).c_str(), "RECO distribution after selections", vars ); 
    }

  // Prepare denominator datasets
  cout<<"Starting denominator dataset filling..."<<endl;
  counter=0;
  int xBin;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
    // find q2 bin of current candidate
    xBin=-1;
    for (int i=0; i<nBins; ++i)
      if ( runBin[i] || runBin[i+nBins] )
	if ( ( pow(genDimuMass,2) < binBorders[i+1] ) &&
	     ( pow(genDimuMass,2) > binBorders[i]   ) ) {
	  xBin = i;
	  break;
	}
    if (xBin<0) continue;
    // status display
    if ( iCand > 1.0*counter*denEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill dataset
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
    if (runBin[xBin])       data[xBin]      ->add(vars);
    if (runBin[xBin+nBins]) data[xBin+nBins]->add(vars);
  }

  // Prepare numerator datasets
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
    // find q2 bin of current candidate
    xBin=-1;
    for (int i=0; i<2*nBins; ++i)
      if ( runBin[i] )
	if ( ( pow(recoDimuMass,2) < binBorders[(i%nBins)+1] ) &&
	     ( pow(recoDimuMass,2) > binBorders[(i%nBins)]   ) )
	  if ( ( (i/nBins==1) && (genSignal == tagB0+1) ) ||
	       ( (i/nBins==0) && (genSignal != tagB0+1) ) ) {
	    xBin = i;
	    break;
	  }
    if (xBin<0) continue;
    // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill dataset
    ctK.setVal(recoCosThetaK);
    ctL.setVal(recoCosThetaL);
    phi.setVal(recoPhi);
    numData[xBin]->add(vars);
  }
  cout<<"Dataset prepared"<<endl;

  // Save datasets
  RooWorkspace *ws [2*nBins];
  for (int i=0; i<2*nBins; ++i) if (runBin[i]) {
      ws [i] = new RooWorkspace(("ws_"+shortString[i]).c_str(),"Workspace with single-bin datasets");
      ws[i]->import( *data   [i], Silence() );
      ws[i]->import( *numData[i], Silence() );
      ws[i]->writeToFile( ( "effDataset_"+shortString[i]+".root" ).c_str() );
    }

  // compute and print average efficiency
  for (int i=0; i<2*nBins; ++i) if (runBin[i]) {
      double avgEff = 1.0 * numData[i]->sumEntries() / data[i]->sumEntries();
      cout<<"Average efficiency bin "<<i%nBins<<(i/nBins==0?" correct-tag = ":" mistag = ")<<avgEff<<endl;
    }

  if (plot) {
    
    // Plot 1D distributions of numerator and denominator datasets
    double rescFac = 0.03;
    TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);

    RooPlot* xframe [2*nBins];
    RooPlot* yframe [2*nBins];
    RooPlot* zframe [2*nBins];

    bool legFilled = false;

    for (int i=0; i<2*nBins; ++i) if (runBin[i]) {
	c [i] = new TCanvas(("c_"+shortString[i]).c_str(),("Num and Den 1D projections - "+longString[i]).c_str(),2000,700);
	xframe [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions").c_str()));
	yframe [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions").c_str()));
	zframe [i] = phi.frame(Title((longString[i]+" #phi distributions").c_str()));
	if (!legFilled) data[i]->plotOn(xframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),Name("plDenDist"));
	else            data[i]->plotOn(xframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
	data[i]->plotOn(yframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
	data[i]->plotOn(zframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
	if (!legFilled) numData[i]->plotOn(xframe[i],Binning(30),Name("plNumDist"));
	else            numData[i]->plotOn(xframe[i],Binning(30));
	numData[i]->plotOn(yframe[i],Binning(30));
	numData[i]->plotOn(zframe[i],Binning(30));
	xframe[i]->GetYaxis()->SetTitleOffset(1.6);
	yframe[i]->GetYaxis()->SetTitleOffset(1.6);
	zframe[i]->GetYaxis()->SetTitleOffset(1.6);
	xframe[i]->SetMaximum(xframe[i]->GetMaximum()*rescFac*1.15);
	yframe[i]->SetMaximum(yframe[i]->GetMaximum()*rescFac*1.15);
	zframe[i]->SetMaximum(zframe[i]->GetMaximum()*rescFac*1.15);
	if (!legFilled) {
	  leg->AddEntry(xframe[i]->findObject("plDenDist"),Form("Generator-level distributions (*%1.2f)",rescFac),"lep");
	  leg->AddEntry(xframe[i]->findObject("plNumDist"),"Post-selection RECO distributions","lep");
	  legFilled = true;
	}
	c[i]->Divide(3,1);
	c[i]->cd(1);
	gPad->SetLeftMargin(0.17); 
	xframe[i]->Draw();
	leg->Draw("same");
	c[i]->cd(2);
	gPad->SetLeftMargin(0.17); 
	yframe[i]->Draw();
	leg->Draw("same");
	c[i]->cd(3);
	gPad->SetLeftMargin(0.17); 
	zframe[i]->Draw();
	leg->Draw("same");

	c[i]->SaveAs( ("dist_GEN_RECO_"+shortString[i]+".pdf").c_str() );
      }
  }

}
