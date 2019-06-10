using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

TCanvas* c [nBins];

void createWeightedDataset(int year, int q2Bin = -1, bool plot = false)
{
  // year format: [6] for 2016
  //              [7] for 2017
  //              [8] for 2018
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( year<6 || year>8 ) return;

  RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
  RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
  RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars (ctK, ctL, phi);

  bool runBin [nBins];
  string shortString [nBins];
  string longString  [nBins];
  for (int i=0; i<nBins; ++i) {
    runBin [i] = false;
    if ( q2Bin!=-1 && q2Bin!=i ) continue;
    runBin [i] = true;
    shortString [i] = Form("b%i",i);
    longString  [i] = Form("q2 bin %i",i);
  }

  // Load GEN-correction weights
  string filename = "correctionWeights.root";
  TFile* fin = TFile::Open( filename.c_str() );
  if ( !fin || !fin->IsOpen() ) {
    cout<<"File not found: "<<filename<<endl;
    return;
  }
  TH3D* hWeight[nBins*2];
  for (int i=0; i<nBins; ++i) if (runBin[i]) {
      hWeight[i] = (TH3D*)fin->Get(Form("hWeight_b%ip0",i));
      if ( !hWeight[i] || hWeight[i]->IsZombie() ) {
	cout<<"Weight histo b"<<i<<"p0 not found in file: "<<filename<<endl;
	cout<<"Using no weight"<<endl;
	hWeight[i] = new TH3D(Form("hWeight_b%ip0",i),Form("hWeight_b%ip0",i),1,-1,1,1,-1,1,1,-TMath::Pi(),TMath::Pi());
	hWeight[i]->SetBinContent(1,1,1,1);
      }
      hWeight[i+nBins] = (TH3D*)fin->Get(Form("hWeight_b%ip1",i));
      if ( !hWeight[i+nBins] || hWeight[i+nBins]->IsZombie() ) {
	cout<<"Weight histo b"<<i<<"p1 not found in file: "<<filename<<endl;
	cout<<"Using no weight"<<endl;
	hWeight[i+nBins] = new TH3D(Form("hWeight_b%ip1",i),Form("hWeight_b%ip1",i),1,-1,1,1,-1,1,1,-TMath::Pi(),TMath::Pi());
	hWeight[i+nBins]->SetBinContent(1,1,1,1);
      }
    }


  // Load ntuples
  TChain* t_gen = new TChain();
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_gen->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/GEN_BFilter_B0MuMuKstar_p*.root/ntuple");
  // 2016
  if ( year==6 ) {
    return; 			// to remove when samples available
    t_den->Add("");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/NtupleNov21/2016MC_LMNR_bdt0p96.root/ntuple");
  // 2017
  } else if ( year==7 ) {
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/2017GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/2017MC_LMNR.root/ntuple");
  } else if ( year==8 ) {
    return; 			// to remove when samples available
    t_den->Add("");
    t_num->Add("");
  }
  int genEntries = t_gen->GetEntries();
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();
  int counter;

  // angular variables
  double genCosThetaK, genCosThetaL, genPhi;
  double recoCosThetaK, recoCosThetaL, recoPhi;
  t_gen->SetBranchAddress( "cos_theta_k"     , &genCosThetaK  );
  t_gen->SetBranchAddress( "cos_theta_l"     , &genCosThetaL  );
  t_gen->SetBranchAddress( "phi_kst_mumu"    , &genPhi        );
  t_den->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK  );
  t_den->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL  );
  t_den->SetBranchAddress( "gen_phi_kst_mumu", &genPhi        );
  t_num->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK  );
  t_num->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL  );
  t_num->SetBranchAddress( "gen_phi_kst_mumu", &genPhi        );
  t_num->SetBranchAddress( "cos_theta_k"     , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l"     , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu"    , &recoPhi       );

  // variables for applying GEN-filter
  double genmupEta, genmumEta, genkstTrkpEta, genkstTrkmEta, genmupPt, genmumPt, genkstTrkpPt, genkstTrkmPt;
  t_gen->SetBranchAddress( "genmupEta", &genmupEta );
  t_gen->SetBranchAddress( "genmumEta", &genmumEta );
  t_gen->SetBranchAddress( "genkstTrkpEta", &genkstTrkpEta );
  t_gen->SetBranchAddress( "genkstTrkmEta", &genkstTrkmEta );
  t_gen->SetBranchAddress( "genmupPt", &genmupPt );
  t_gen->SetBranchAddress( "genmumPt", &genmumPt );
  t_gen->SetBranchAddress( "genkstTrkpPt", &genkstTrkpPt );
  t_gen->SetBranchAddress( "genkstTrkmPt", &genkstTrkmPt );

  // dimuon mass variables
  double genDimuMass2, recoDimuMass;
  t_gen->SetBranchAddress( "genQ2"   , &genDimuMass2 );
  t_den->SetBranchAddress( "genQ2"   , &genDimuMass2 );
  t_num->SetBranchAddress( "mumuMass", &recoDimuMass );

  // B0-kinematic variables
  // double genB0pT, genB0eta;
  // double recoB0pT, recoB0eta;
  // t_gen->SetBranchAddress( "genbPt" , &genB0pT   );
  // t_gen->SetBranchAddress( "genbEta", &genB0eta  );
  // t_den->SetBranchAddress( "genbPt" , &genB0pT   );
  // t_den->SetBranchAddress( "genbEta", &genB0eta  );
  // t_num->SetBranchAddress( "bPt"    , &recoB0pT  );
  // t_num->SetBranchAddress( "bEta"   , &recoB0eta );

  // flavour tagging variables
  double genSignal, tagB0;
  t_num->SetBranchAddress( "genSignal", &genSignal );
  t_num->SetBranchAddress( "tagB0"    , &tagB0     );

  // event number for even/odd splitting
  double eventN_Dou;
  Long64_t eventN;
  t_gen->SetBranchAddress( "eventN", &eventN_Dou );
  t_den->SetBranchAddress( "eventN", &eventN     );
  t_num->SetBranchAddress( "eventN", &eventN     );

  // event pileup weight
  float PUweight = 1;
  t_den->SetBranchAddress( "weight", &PUweight );
  t_num->SetBranchAddress( "weight", &PUweight );

  // Define datasets
  RooDataSet* data_genDen_ev [nBins];
  RooDataSet* data_genDen_od [nBins];
  RooDataSet* data_genNum_ev [nBins];
  RooDataSet* data_genNum_od [nBins];
  RooDataSet* data_den_ev    [nBins];
  RooDataSet* data_den_od    [nBins];
  RooDataSet* data_ctRECO_ev [nBins];
  RooDataSet* data_ctRECO_od [nBins];
  RooDataSet* data_wtRECO_ev [nBins];
  RooDataSet* data_wtRECO_od [nBins];
  for (int i=0; i<nBins; ++i) if (runBin[i]) {
      data_genDen_ev [i] = new RooDataSet( ("data_genDen_ev_"+shortString[i]).c_str(), "GEN distribution before GEN-filter (even)", vars );
      data_genDen_od [i] = new RooDataSet( ("data_genDen_od_"+shortString[i]).c_str(), "GEN distribution before GEN-filter (odd)", vars );
      data_genNum_ev [i] = new RooDataSet( ("data_genNum_ev_"+shortString[i]).c_str(), "GEN distribution after GEN-filter (even)", vars );
      data_genNum_od [i] = new RooDataSet( ("data_genNum_od_"+shortString[i]).c_str(), "GEN distribution after GEN-filter (odd)", vars );
      data_den_ev    [i] = new RooDataSet( ("data_den_ev_"   +shortString[i]).c_str(), "GEN candidates after GEN-filter in full MC sample (even)", vars );
      data_den_od    [i] = new RooDataSet( ("data_den_od_"   +shortString[i]).c_str(), "GEN candidates after GEN-filter in full MC sample (odd)", vars );
      data_ctRECO_ev [i] = new RooDataSet( ("data_ctRECO_ev_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (even)", vars ); 
      data_ctRECO_od [i] = new RooDataSet( ("data_ctRECO_od_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (odd)", vars ); 
      data_wtRECO_ev [i] = new RooDataSet( ("data_wtRECO_ev_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (even)", vars ); 
      data_wtRECO_od [i] = new RooDataSet( ("data_wtRECO_od_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (odd)", vars ); 
    }

  // Prepare GEN-level datasets
  cout<<"Starting GEN datasets filling..."<<endl;
  counter=0;
  int xBin;
  double corrWei=1;
  for (int iCand=0; iCand<genEntries; ++iCand) {
    t_gen->GetEntry(iCand);
    // find q2 bin of current candidate
    xBin=-1;
    for (int i=0; i<nBins; ++i)
      if ( runBin[i] )
	if ( ( genDimuMass2 < binBorders[i+1] ) &&
	     ( genDimuMass2 > binBorders[i]   ) ) {
	  xBin = i;
	  break;
	}
    if (xBin<0) continue;
    // status display
    if ( iCand > 1.0*counter*genEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // compute weight
    // Caveat: weights are saved according to the "efficiency parity", so inverted wrt "dataset parity"
    // thus, using here coherent "weight parity" and "dataset parity", the uncorrelated weights are applied
    if (((int)eventN_Dou)%2==0)
      corrWei=hWeight[xBin]->GetBinContent(hWeight[xBin]->GetXaxis()->FindBin(genCosThetaK),
					   hWeight[xBin]->GetYaxis()->FindBin(genCosThetaL),
					   hWeight[xBin]->GetZaxis()->FindBin(genPhi));
    else
      corrWei=hWeight[xBin+nBins]->GetBinContent(hWeight[xBin+nBins]->GetXaxis()->FindBin(genCosThetaK),
						 hWeight[xBin+nBins]->GetYaxis()->FindBin(genCosThetaL),
						 hWeight[xBin+nBins]->GetZaxis()->FindBin(genPhi));
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
    // fill genDen dataset
    if (((int)eventN_Dou)%2==0) data_genDen_ev[xBin]->add(vars,corrWei);
    else data_genDen_od[xBin]->add(vars,corrWei);
    // fill genNum dataset
    if ( fabs(genmupEta)<2.5 && fabs(genmumEta)<2.5 &&
	 fabs(genkstTrkpEta)<2.5 && fabs(genkstTrkmEta)<2.5 &&
	 genmupPt>2.5 && genmumPt>2.5 &&
	 genkstTrkpPt>0.4 && genkstTrkmPt>0.4) {
      if (((int)eventN_Dou)%2==0) data_genNum_ev[xBin]->add(vars,corrWei);
      else data_genNum_od[xBin]->add(vars,corrWei);
    }
  }
  // Prepare denominator dataset
  cout<<"Starting denominator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
    // find q2 bin of current candidate
    xBin=-1;
    for (int i=0; i<nBins; ++i)
      if ( runBin[i] )
	if ( ( genDimuMass2 < binBorders[i+1] ) &&
	     ( genDimuMass2 > binBorders[i]   ) ) {
	  xBin = i;
	  break;
	}
    if (xBin<0) continue;
    // status display
    if ( iCand > 1.0*counter*denEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // compute weight
    // Caveat: weights are saved according to the "efficiency parity", so inverted wrt "dataset parity"
    // thus, using here coherent "weight parity" and "dataset parity", the uncorrelated weights are applied
    if (eventN%2==0)
      corrWei=hWeight[xBin]->GetBinContent(hWeight[xBin]->GetXaxis()->FindBin(genCosThetaK),
					   hWeight[xBin]->GetYaxis()->FindBin(genCosThetaL),
					   hWeight[xBin]->GetZaxis()->FindBin(genPhi));
    else
      corrWei=hWeight[xBin+nBins]->GetBinContent(hWeight[xBin+nBins]->GetXaxis()->FindBin(genCosThetaK),
						 hWeight[xBin+nBins]->GetYaxis()->FindBin(genCosThetaL),
						 hWeight[xBin+nBins]->GetZaxis()->FindBin(genPhi));
    // fill dataset
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
    if (eventN%2==0) data_den_ev[xBin]->add(vars,PUweight*corrWei);
    else data_den_od[xBin]->add(vars,PUweight*corrWei);
  }

  // Prepare numerator dataset
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
    // find q2 bin of current candidate
    xBin=-1;
    for (int i=0; i<nBins; ++i)
      if ( runBin[i] )
	if ( ( pow(recoDimuMass,2) < binBorders[i+1] ) &&
	     ( pow(recoDimuMass,2) > binBorders[i]   ) ) {
	  xBin = i;
	  break;
	}
    if (xBin<0) continue;
    // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // compute weight
    // Caveat: weights are saved according to the "efficiency parity", so inverted wrt "dataset parity"
    // thus, using here coherent "weight parity" and "dataset parity", the uncorrelated weights are applied
    if (eventN%2==0)
      corrWei=hWeight[xBin]->GetBinContent(hWeight[xBin]->GetXaxis()->FindBin(genCosThetaK),
					   hWeight[xBin]->GetYaxis()->FindBin(genCosThetaL),
					   hWeight[xBin]->GetZaxis()->FindBin(genPhi));
    else
      corrWei=hWeight[xBin+nBins]->GetBinContent(hWeight[xBin+nBins]->GetXaxis()->FindBin(genCosThetaK),
						 hWeight[xBin+nBins]->GetYaxis()->FindBin(genCosThetaL),
						 hWeight[xBin+nBins]->GetZaxis()->FindBin(genPhi));
    // fill dataset
    ctK.setVal(recoCosThetaK);
    ctL.setVal(recoCosThetaL);
    phi.setVal(recoPhi);
    if (genSignal != tagB0+1) {
      if (eventN%2==0) data_ctRECO_ev[xBin]->add(vars,PUweight*corrWei);
      else data_ctRECO_od[xBin]->add(vars,PUweight*corrWei);
    } else {
      if (eventN%2==0) data_wtRECO_ev[xBin]->add(vars,PUweight*corrWei);
      else data_wtRECO_od[xBin]->add(vars,PUweight*corrWei);
    }
  }
  cout<<"Dataset prepared"<<endl;

  // Save datasets
  RooWorkspace *ws [nBins];
  for (int i=0; i<nBins; ++i) if (runBin[i]) {
      ws [i] = new RooWorkspace(("ws_"+shortString[i]).c_str(),"Workspace with single-bin datasets");
      ws[i]->import( *data_genDen_ev[i] );
      ws[i]->import( *data_genDen_od[i] );
      ws[i]->import( *data_genNum_ev[i] );
      ws[i]->import( *data_genNum_od[i] );
      ws[i]->import( *data_den_ev   [i] );
      ws[i]->import( *data_den_od   [i] );
      ws[i]->import( *data_ctRECO_ev[i] );
      ws[i]->import( *data_ctRECO_od[i] );
      ws[i]->import( *data_wtRECO_ev[i] );
      ws[i]->import( *data_wtRECO_od[i] );
      ws[i]->writeToFile( ( "effWeightedDataset_"+shortString[i]+".root" ).c_str() );
    }

  // compute and print average efficiency
  for (int i=0; i<nBins; ++i) if (runBin[i]) {
      double den_ev = data_genNum_ev[i]->sumEntries() / data_genDen_ev[i]->sumEntries() / data_den_ev[i]->sumEntries();
      double den_od = data_genNum_od[i]->sumEntries() / data_genDen_od[i]->sumEntries() / data_den_od[i]->sumEntries();
      double avgEff_ct_ev = data_ctRECO_ev[i]->sumEntries() * den_ev;
      double avgEff_ct_od = data_ctRECO_od[i]->sumEntries() * den_od;
      double avgEff_wt_ev = data_wtRECO_ev[i]->sumEntries() * den_ev;
      double avgEff_wt_od = data_wtRECO_od[i]->sumEntries() * den_od;
      double avgEff_ev = avgEff_ct_ev + avgEff_wt_ev;
      double avgEff_od = avgEff_ct_od + avgEff_wt_od;
      double avgMis_ev = avgEff_wt_ev / avgEff_ev;
      double avgMis_od = avgEff_wt_od / avgEff_od;
      cout<<"Averages bin "<<nBins<<" (even): eps="<<avgEff_ev<<"\tm="<<avgMis_ev<<"\teps_c="<<avgEff_ct_ev<<"\teps_m="<<avgEff_wt_ev<<endl;
      cout<<"Averages bin "<<nBins<<" (odd) : eps="<<avgEff_od<<"\tm="<<avgMis_od<<"\teps_c="<<avgEff_ct_od<<"\teps_m="<<avgEff_wt_od<<endl;
    }

  if (plot) {

    // Plot 1D distributions of datasets
    double rescFac1 = 1.0/12;
    double rescFac2 = 1.0;
    double rescFac3 = 1.0/25;

    TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);

    RooPlot* xframe_ev [nBins];
    RooPlot* yframe_ev [nBins];
    RooPlot* zframe_ev [nBins];
    RooPlot* xframe_od [nBins];
    RooPlot* yframe_od [nBins];
    RooPlot* zframe_od [nBins];

    bool legFilled = false;

    RooDataSet* data_RECO_ev [nBins];
    RooDataSet* data_RECO_od [nBins];

    for (int i=0; i<nBins; ++i) if (runBin[i]) {
	data_RECO_ev [i] = new RooDataSet( ("data_RECO_ev_"+shortString[i]).c_str(), "Reconstructed candidates after selections (even)", data_ctRECO_ev[i], vars );
	data_RECO_od [i] = new RooDataSet( ("data_RECO_od_"+shortString[i]).c_str(), "Reconstructed candidates after selections (odd)" , data_ctRECO_od[i], vars );
	data_RECO_ev[i]->append(*(data_wtRECO_ev[i]));
	data_RECO_od[i]->append(*(data_wtRECO_od[i]));

	c [i] = new TCanvas(("c_"+shortString[i]).c_str(),("Num and Den 1D projections - "+longString[i]).c_str(),2000,1400);
	xframe_ev [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (even)").c_str()));
	yframe_ev [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (even)").c_str()));
	zframe_ev [i] = phi.frame(Title((longString[i]+" #phi distributions (even)").c_str()));
	xframe_od [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (odd)").c_str()));
	yframe_od [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (odd)").c_str()));
	zframe_od [i] = phi.frame(Title((longString[i]+" #phi distributions (odd)").c_str()));

	// if (!legFilled) data[i]->plotOn(xframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30),Name("plDenDist"));
	// else            data[i]->plotOn(xframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
	// data[i]->plotOn(yframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
	// data[i]->plotOn(zframe[i],Rescale(rescFac),MarkerColor(kRed+1),LineColor(kRed+1),Binning(30));
	// if (!legFilled) numData[i]->plotOn(xframe[i],Binning(30),Name("plNumDist"));
	// else            numData[i]->plotOn(xframe[i],Binning(30));
	// numData[i]->plotOn(yframe[i],Binning(30));
	// numData[i]->plotOn(zframe[i],Binning(30));

	data_genDen_ev[i]->plotOn(xframe_ev[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	data_genNum_ev[i]->plotOn(xframe_ev[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	data_den_ev   [i]->plotOn(xframe_ev[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_genDen_od[i]->plotOn(xframe_od[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	data_genNum_od[i]->plotOn(xframe_od[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	data_den_od   [i]->plotOn(xframe_od[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	data_ctRECO_od[i]->plotOn(xframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_od[i]->plotOn(xframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_od  [i]->plotOn(xframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_genDen_ev[i]->plotOn(yframe_ev[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	data_genNum_ev[i]->plotOn(yframe_ev[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	data_den_ev   [i]->plotOn(yframe_ev[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	data_ctRECO_ev[i]->plotOn(yframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_ev[i]->plotOn(yframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_ev  [i]->plotOn(yframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_genDen_od[i]->plotOn(yframe_od[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	data_genNum_od[i]->plotOn(yframe_od[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	data_den_od   [i]->plotOn(yframe_od[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	data_ctRECO_od[i]->plotOn(yframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_od[i]->plotOn(yframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_od  [i]->plotOn(yframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_genDen_ev[i]->plotOn(zframe_ev[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	data_genNum_ev[i]->plotOn(zframe_ev[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	data_den_ev   [i]->plotOn(zframe_ev[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	data_ctRECO_ev[i]->plotOn(zframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_ev[i]->plotOn(zframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_ev  [i]->plotOn(zframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_genDen_od[i]->plotOn(zframe_od[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	data_genNum_od[i]->plotOn(zframe_od[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	data_den_od   [i]->plotOn(zframe_od[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	data_ctRECO_od[i]->plotOn(zframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_od[i]->plotOn(zframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_od  [i]->plotOn(zframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));

	xframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
	yframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
	zframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
	xframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
	yframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
	zframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
	xframe_ev[i]->SetMaximum(xframe_ev[i]->GetMaximum()*rescFac1);
	yframe_ev[i]->SetMaximum(yframe_ev[i]->GetMaximum()*rescFac1);
	zframe_ev[i]->SetMaximum(zframe_ev[i]->GetMaximum()*rescFac1);
	xframe_od[i]->SetMaximum(xframe_od[i]->GetMaximum()*rescFac1);
	yframe_od[i]->SetMaximum(yframe_od[i]->GetMaximum()*rescFac1);
	zframe_od[i]->SetMaximum(zframe_od[i]->GetMaximum()*rescFac1);
	// if (!legFilled) {
	//   leg->AddEntry(xframe[i]->findObject("plDenDist"),Form("Generator-level distributions (*%1.2f)",rescFac),"lep");
	//   leg->AddEntry(xframe[i]->findObject("plNumDist"),"Post-selection RECO distributions","lep");
	//   legFilled = true;
	// }
	c[i]->Divide(3,2);
	c[i]->cd(1);
	gPad->SetLeftMargin(0.17); 
	xframe_ev[i]->Draw();
	// leg->Draw("same");
	c[i]->cd(2);
	gPad->SetLeftMargin(0.17); 
	yframe_ev[i]->Draw();
	// leg->Draw("same");
	c[i]->cd(3);
	gPad->SetLeftMargin(0.17); 
	zframe_ev[i]->Draw();
	// leg->Draw("same");
	c[i]->cd(4);
	gPad->SetLeftMargin(0.17); 
	xframe_od[i]->Draw();
	// leg->Draw("same");
	c[i]->cd(5);
	gPad->SetLeftMargin(0.17); 
	yframe_od[i]->Draw();
	// leg->Draw("same");
	c[i]->cd(6);
	gPad->SetLeftMargin(0.17); 
	zframe_od[i]->Draw();
	// leg->Draw("same");

	c[i]->SaveAs( ("distWeighted_GEN_RECO_"+shortString[i]+".pdf").c_str() );
      }
  }

}
