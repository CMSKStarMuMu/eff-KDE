using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

double PDGB0Mass = 5.27958;
double PDGJpsiMass = 3.096916;
double PDGPsiPrimeMass = 3.686109;
double PDGKstMass = 0.896;

TCanvas* c [nBins];

void createDataset(int year, int q2Bin = -1, bool plot = false)
{
  // year format: [6] for 2016
  //              [7] for 2017
  //              [8] for 2018
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( year<6 || year>8 ) return;

  bool isJpsi = false;
  bool isPsi  = false;
  bool isLMNR = false;
  
  if (q2Bin==4)      isJpsi = true;
  else if (q2Bin==6) isPsi = true;
  else isLMNR = true;
  
  // define if need to split into more files due to workspace memory 
  bool saveDen = true;
  if (q2Bin==4) saveDen = false;

  // define angular variables and variable for PU-reweighting
  RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
  RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
  RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars (ctK, ctL, phi);
  RooRealVar wei ("weight","weight",1);
  RooRealVar rand("rand", "rand", 0,1);
  TRandom rand_gen(1029);
  RooArgSet all_vars (ctK, ctL, phi, rand);

  // flags to mark which q2 bins should be filled
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

  // Load ntuples
  TChain* t_gen = new TChain();
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  //first for 2017
  if (isLMNR){
    t_gen->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0MuMuKstar_p*.root/ntuple");
    if ( year==6 ) {
      // 2016
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016GEN_MC_LMNR.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016MC_LMNR.root/ntuple");
    } else if ( year==7 ) {
      // 2017
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017GEN_MC_LMNR.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017MC_LMNR.root/ntuple");
    } else if ( year==8 ) {
      // 2018
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018GEN_MC_LMNR.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018MC_LMNR.root/ntuple");
    }
  }
  else if (isJpsi){
    t_gen->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0JpsiKstar.root/ntuple");
    if ( year==6 ) {
      // 2016
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016GEN_MC_JPSI.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016MC_JPSI.root/ntuple");
    } else if ( year==7 ) {
      // 2017
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017GEN_MC_JPSI.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017MC_JPSI.root/ntuple");
    } else if ( year==8 ) {
      // 2018
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018GEN_MC_JPSI.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018MC_JPSI.root/ntuple");
    }  
  }
  else if (isPsi){
    t_gen->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0PsiKstar.root/ntuple");
    if ( year==6 ) {
      // 2016
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016GEN_MC_PSI.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016MC_PSI.root/ntuple");
    } else if ( year==7 ) {
      // 2017
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017GEN_MC_PSI.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017MC_PSI.root/ntuple");
    } else if ( year==8 ) {
      // 2018
      t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018GEN_MC_PSI.root/ntuple");
      t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018MC_PSI.root/ntuple");
    }  
  }
  else 
    return;
  int genEntries = t_gen->GetEntries();
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();
  int counter;

  // Import branches from ntuples:
  // angular variables
  double genCosThetaK, genCosThetaL, genPhi;
  double recoCosThetaK, recoCosThetaL, recoPhi;
  t_gen->SetBranchAddress( "cos_theta_k"     , &genCosThetaK  );
  t_gen->SetBranchAddress( "cos_theta_l"     , &genCosThetaL  );
  t_gen->SetBranchAddress( "phi_kst_mumu"    , &genPhi        );
  t_den->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK  );
  t_den->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL  );
  t_den->SetBranchAddress( "gen_phi_kst_mumu", &genPhi        );
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

  // B0 mass variable
  double recoB0Mass;
  bool passB0Psi_lmnr, passB0Psi_jpsi, passB0Psi_psip;
  t_num->SetBranchAddress( "tagged_mass", &recoB0Mass );
  t_num->SetBranchAddress( "passB0Psi_lmnr", &passB0Psi_lmnr );
  t_num->SetBranchAddress( "passB0Psi_jpsi", &passB0Psi_jpsi );
  t_num->SetBranchAddress( "passB0Psi_psip", &passB0Psi_psip );
  
  
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
//   double eventN_Dou;
  Long64_t eventN;
  t_gen->SetBranchAddress( "eventN", &eventN     );
//   t_gen->SetBranchAddress( "eventN", &eventN_Dou );
  t_den->SetBranchAddress( "eventN", &eventN     );
  t_num->SetBranchAddress( "eventN", &eventN     );

  // event pileup weight
  double PUweight = 1;
  double fPUweight = 1;
  
  t_den->SetBranchAddress( "weight", &fPUweight );
  t_num->SetBranchAddress( "weight", &PUweight );

  // final state radiation flag
  double genSignHasFSR;
  t_gen->SetBranchAddress( "genSignHasFSR", &genSignHasFSR );

  
  // cut to remove B+->Psi(2S)K->Jpsi pi pi K
  // will be a boolean in ntuples in the future
  // keep it here for now since not finalized
  // as from https://github.com/CMSKStarMuMu/RooSideBand/blob/master/testSidebandFit.cc#L2592-L2661
  double wt_mass, wt_kstarmass, kaonPt, pionPt, mmpiMass, mmkMass;
  t_num->SetBranchAddress( "wt_mass",      &wt_mass      );
  t_num->SetBranchAddress( "wt_kstarmass", &wt_kstarmass );
  t_num->SetBranchAddress( "kaonPt",       &kaonPt       );
  t_num->SetBranchAddress( "pionPt",       &pionPt       );
  t_num->SetBranchAddress( "mmpiMass",     &mmpiMass     );
  t_num->SetBranchAddress( "mmkMass",      &mmkMass      );

  double x0Cut=-0.4;
  double y0Cut= 0.3;
  double x1Cut= 0.6;
  double y1Cut=-0.1;
  
  double x_0Cut=3;
  double y_0Cut=3.8;
  double x_1Cut=3.6;
  double y_1Cut=4.8;
  
  double CutX1=3.2;
  double CutX2=3.6;
  double CutY1=4.7;
  double CutY2=4.9;


  // Define datasets for five efficiency terms
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
      data_genDen_ev [i] = new RooDataSet( ("data_genDen_ev_"+shortString[i]).c_str(), "GEN distribution before GEN-filter (even)",
					   all_vars );
      data_genDen_od [i] = new RooDataSet( ("data_genDen_od_"+shortString[i]).c_str(), "GEN distribution before GEN-filter (odd)",
					   all_vars );
      data_genNum_ev [i] = new RooDataSet( ("data_genNum_ev_"+shortString[i]).c_str(), "GEN distribution after GEN-filter (even)",
					   vars );
      data_genNum_od [i] = new RooDataSet( ("data_genNum_od_"+shortString[i]).c_str(), "GEN distribution after GEN-filter (odd)",
					   vars );
      data_den_ev    [i] = new RooDataSet( ("data_den_ev_"   +shortString[i]).c_str(), "GEN candidates after GEN-filter in full MC sample (even)",
					   RooArgSet(ctK,ctL,phi,wei), "weight" );
      data_den_od    [i] = new RooDataSet( ("data_den_od_"   +shortString[i]).c_str(), "GEN candidates after GEN-filter in full MC sample (odd)",
					   RooArgSet(ctK,ctL,phi,wei), "weight" );
      data_ctRECO_ev [i] = new RooDataSet( ("data_ctRECO_ev_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (even)",
					   RooArgSet(ctK,ctL,phi,wei), "weight" );
      data_ctRECO_od [i] = new RooDataSet( ("data_ctRECO_od_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (odd)",
					   RooArgSet(ctK,ctL,phi,wei), "weight" );
      data_wtRECO_ev [i] = new RooDataSet( ("data_wtRECO_ev_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (even)",
					   RooArgSet(ctK,ctL,phi,wei), "weight" );
      data_wtRECO_od [i] = new RooDataSet( ("data_wtRECO_od_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (odd)",
					   RooArgSet(ctK,ctL,phi,wei), "weight" );
    }

  // Define counter of total genDen events (w/o FSR veto) for correct efficiency normalisation
  // saved as a two-bin TH1I object
  TH1I* n_genDen [nBins];
  for (int i=0; i<nBins; ++i)
    if (runBin[i])
      n_genDen [i] = new TH1I( ("n_genDen_"+shortString[i]).c_str(), ("n_genDen_"+shortString[i]).c_str(),
			       2, -0.5, 1.5);

  // Prepare GEN-level datasets
  cout<<"Starting GEN datasets filling..."<<endl;
  counter=0;
  int xBin;
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
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
    rand.setVal(rand_gen.Uniform(1));
    // fill genDen dataset
    n_genDen[xBin]->Fill((eventN)%2);
    if ( genSignHasFSR<0.5 ) {
      if ((eventN)%2==0) data_genDen_ev[xBin]->add(all_vars);
      else data_genDen_od[xBin]->add(all_vars);
    }
    // apply same selection as in GEN-filter of recoDen MC sample
    // and fill genNum dataset
    if ( fabs(genmupEta)<2.5 && fabs(genmumEta)<2.5 &&
	 fabs(genkstTrkpEta)<2.5 && fabs(genkstTrkmEta)<2.5 &&
	 genmupPt>2.5 && genmumPt>2.5 &&
	 genkstTrkpPt>0.4 && genkstTrkmPt>0.4) {
      if ((eventN)%2==0) data_genNum_ev[xBin]->add(vars);
      else data_genNum_od[xBin]->add(vars);
    }
  }
  
  delete t_gen;

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
    // fill dataset
    ctK.setVal(genCosThetaK);
    ctL.setVal(genCosThetaL);
    phi.setVal(genPhi);
    if (eventN%2==0) data_den_ev[xBin]->add(vars,fPUweight);
    else data_den_od[xBin]->add(vars,fPUweight);
  }

  delete t_den;

  // Prepare numerator dataset
  cout<<"Starting numerator dataset filling..."<<endl;
  counter=0;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
    // anti-radiation cut
    if (isLMNR && passB0Psi_lmnr == 0) continue;
    else if (isJpsi && passB0Psi_jpsi == 0) continue;
    else if (isPsi  && passB0Psi_psip == 0)  continue;

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

    // apply cut for bin 4 
    bool XCut= (( (PDGB0Mass - wt_mass) - y0Cut ) / (y1Cut-y0Cut)) < (((wt_kstarmass-PDGKstMass)-x0Cut) / (x1Cut-x0Cut)) && \
                  kaonPt > pionPt && \
                  (wt_kstarmass-PDGKstMass)>0 && \
                  (mmpiMass > CutX1 && mmpiMass < CutX2) && \
                  (mmkMass >  CutY1 && mmkMass  < CutY2) && \
                  ((mmkMass - y_0Cut) / (y_1Cut - y_0Cut)) > ((mmpiMass-x_0Cut)/(x_1Cut-x_0Cut));
  
    if (XCut && xBin == 4) continue;

    // status display
    if ( iCand > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }
    // fill dataset
    ctK.setVal(recoCosThetaK);
    ctL.setVal(recoCosThetaL);
    phi.setVal(recoPhi);
    if (genSignal != tagB0+1) { // correctly tagged events
      if (eventN%2==0) data_ctRECO_ev[xBin]->add(vars,PUweight);
      else data_ctRECO_od[xBin]->add(vars,PUweight);
    } else { // wrongly tagged events
      if (eventN%2==0) data_wtRECO_ev[xBin]->add(vars,PUweight);
      else data_wtRECO_od[xBin]->add(vars,PUweight);
    }
  }
  delete t_num;
  cout<<"Dataset prepared"<<endl;

  // Save datasets in workspaces
  RooWorkspace *ws_ev [nBins];
  RooWorkspace *ws_od [nBins];
  for (int i=0; i<nBins; ++i) if (runBin[i]) {
      // Skip the creation of a file when the correct-tag efficiency cannot be computed (empty numerators)
      // which usually means that either you are using a resonant MC, which does not fill signal q2 bins,
      // or using a bin too fine, or out of range
      // If correct-tag numerator is filled and wrong-tag is not, a warning is returned
      if ( data_genNum_ev[i]->numEntries()==0 || data_genNum_od[i]->numEntries()==0 ) {
	cout<<"Error: genNum is empty in q2 bin "<<i<<endl;
	continue;
      }
      if ( data_ctRECO_ev[i]->numEntries()==0 || data_ctRECO_od[i]->numEntries()==0 ) {
	cout<<"Error: ctRECO is empty in q2 bin "<<i<<endl;
	continue;
      }
      if ( data_wtRECO_ev[i]->numEntries()==0 || data_wtRECO_od[i]->numEntries()==0 ) {
	cout<<"Warning: wtRECO is empty in q2 bin "<<i<<endl;
      }
      ws_ev[i] = new RooWorkspace(("ws_"+shortString[i]+"p0").c_str(),"Workspace with single-bin even datasets");
      ws_od[i] = new RooWorkspace(("ws_"+shortString[i]+"p1").c_str(),"Workspace with single-bin odd datasets");
      ws_ev[i]->import( *data_genDen_ev[i] );
      ws_od[i]->import( *data_genDen_od[i] );
      ws_ev[i]->import( *data_genNum_ev[i] );
      ws_od[i]->import( *data_genNum_od[i] );
      if (saveDen){
        ws_ev[i]->import( *data_den_ev[i] );
        ws_od[i]->import( *data_den_od[i] );
      }  
      ws_ev[i]->import( *data_ctRECO_ev[i] );
      ws_od[i]->import( *data_ctRECO_od[i] );
      ws_ev[i]->import( *data_wtRECO_ev[i] );
      ws_od[i]->import( *data_wtRECO_od[i] );
      TFile* fout = new TFile( ( "effDataset_"+shortString[i]+Form("_201%i.root",year) ).c_str(), "RECREATE" );
      ws_ev[i]->Write();
      ws_od[i]->Write();
      n_genDen[i]->Write();
      fout->Close();

      if (!saveDen){
        TFile* fout = new TFile( ( "effDataset_"+shortString[i]+Form("_201%i_den.root",year) ).c_str(), "RECREATE" );
        ws_ev[i] = new RooWorkspace(("ws_"+shortString[i]+"p0").c_str(),"Workspace with single-bin even datasets");
        ws_od[i] = new RooWorkspace(("ws_"+shortString[i]+"p1").c_str(),"Workspace with single-bin odd datasets");
        ws_ev[i]->import( *data_den_ev[i] );
        ws_od[i]->import( *data_den_od[i] );
        ws_ev[i]->Write();
        ws_od[i]->Write();
        n_genDen[i]->Write();
        fout->Close();
      }
    }

  // compute and print average efficiency (merged and individual tag configurations) and mistag fraction
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
      cout<<"Averages bin "<<i<<" (even): eps="<<avgEff_ev<<"\tm="<<avgMis_ev<<"\teps_c="<<avgEff_ct_ev<<"\teps_m="<<avgEff_wt_ev<<endl;
      cout<<"Averages bin "<<i<<" (odd) : eps="<<avgEff_od<<"\tm="<<avgMis_od<<"\teps_c="<<avgEff_ct_od<<"\teps_m="<<avgEff_wt_od<<endl;
    }

  // Plot 1D distributions of datasets
  if (plot) {

    // to keep all distributions visible in the same plot, the ones with higher stats (tipically denominators) need to be rescaled
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
	// Create dataset containing both correct-tag and wrong-tag events
	data_RECO_ev [i] = new RooDataSet( ("data_RECO_ev_"+shortString[i]).c_str(), "Reconstructed candidates after selections (even)", data_ctRECO_ev[i], vars );
	data_RECO_od [i] = new RooDataSet( ("data_RECO_od_"+shortString[i]).c_str(), "Reconstructed candidates after selections (odd)" , data_ctRECO_od[i], vars );
	data_RECO_ev[i]->append(*(data_wtRECO_ev[i]));
	data_RECO_od[i]->append(*(data_wtRECO_od[i]));

	// create frames (one for each bin/parity/variable, but all the six efficiency terms are plotted together)
	c [i] = new TCanvas(("c_"+shortString[i]).c_str(),("Num and Den 1D projections - "+longString[i]).c_str(),2000,1400);
	xframe_ev [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (even)").c_str()));
	yframe_ev [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (even)").c_str()));
	zframe_ev [i] = phi.frame(Title((longString[i]+" #phi distributions (even)").c_str()));
	xframe_od [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (odd)").c_str()));
	yframe_od [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (odd)").c_str()));
	zframe_od [i] = phi.frame(Title((longString[i]+" #phi distributions (odd)").c_str()));

	// plot datasets on frames
	if (!legFilled) { // the first time assign names to tag them in the legend
	  data_genDen_ev[i]->plotOn(xframe_ev[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1),Name("plGenDen"));
	  data_genNum_ev[i]->plotOn(xframe_ev[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2),Name("plGenNum"));
	  data_den_ev   [i]->plotOn(xframe_ev[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3),Name("plRecoDen"));
	  data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40),Name("plCTrecoNum"));
	  data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40),Name("plWTrecoNum"));
	  data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40),Name("plRecoNum"));
	} else {
	  data_genDen_ev[i]->plotOn(xframe_ev[i],MarkerColor(kRed+1)   ,LineColor(kRed+1)   ,Binning(40),Rescale(rescFac1));
	  data_genNum_ev[i]->plotOn(xframe_ev[i],MarkerColor(kBlue)    ,LineColor(kBlue)    ,Binning(40),Rescale(rescFac2));
	  data_den_ev   [i]->plotOn(xframe_ev[i],MarkerColor(kGreen+2) ,LineColor(kGreen+2) ,Binning(40),Rescale(rescFac3));
	  data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	  data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	  data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	}
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

	// Fill legend (only one time)
	if (!legFilled) {
	  string strRescFac1 = (rescFac1<1?Form(" (*%1.2f)",rescFac1):"");
	  string strRescFac2 = (rescFac2<1?Form(" (*%1.2f)",rescFac2):"");
	  string strRescFac3 = (rescFac3<1?Form(" (*%1.2f)",rescFac3):"");
          leg->AddEntry(xframe_ev[i]->findObject("plGenDen"   ),("Generator-level distribution"+strRescFac1).c_str(),"lep");
          leg->AddEntry(xframe_ev[i]->findObject("plGenNum"   ),("Post-GEN-filter distribution of special MC sample"+strRescFac2).c_str(),"lep");
          leg->AddEntry(xframe_ev[i]->findObject("plRecoDen"  ),("Post-GEN-filter distribution of full MC sample"+strRescFac3).c_str(),"lep");
	  leg->AddEntry(xframe_ev[i]->findObject("plRecoNum"  ),"Post-selection distribution","lep");
	  leg->AddEntry(xframe_ev[i]->findObject("plCTrecoNum"),"Post-selection correct-tag distribution","lep");
	  leg->AddEntry(xframe_ev[i]->findObject("plWTrecoNum"),"Post-selection wrong-tag distribution","lep");
	  legFilled = true;
	}

	// Plot even distributions in the top row and odd ones in the bottom row
	c[i]->Divide(3,2);
	c[i]->cd(1);
	gPad->SetLeftMargin(0.17); 
	xframe_ev[i]->Draw();
	leg->Draw("same");
	c[i]->cd(2);
	gPad->SetLeftMargin(0.17); 
	yframe_ev[i]->Draw();
	leg->Draw("same");
	c[i]->cd(3);
	gPad->SetLeftMargin(0.17); 
	zframe_ev[i]->Draw();
	leg->Draw("same");
	c[i]->cd(4);
	gPad->SetLeftMargin(0.17); 
	xframe_od[i]->Draw();
	leg->Draw("same");
	c[i]->cd(5);
	gPad->SetLeftMargin(0.17); 
	yframe_od[i]->Draw();
	leg->Draw("same");
	c[i]->cd(6);
	gPad->SetLeftMargin(0.17); 
	zframe_od[i]->Draw();
	leg->Draw("same");

	c[i]->SaveAs( ("plotDist_d/dist_GEN_RECO_"+shortString[i]+".pdf").c_str() );
      }
  }

}
