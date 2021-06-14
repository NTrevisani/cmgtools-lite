/*

root -l -b -q 'my_tmva_reader_2.C("2016")'
root -l -b -q 'my_tmva_reader_2.C("2016_VFP")'

root -l -b -q 'my_tmva_reader_2.C("2016_TTG")'
root -l -b -q 'my_tmva_reader_2.C("2016_TTG_VFP")'

*/


#include <algorithm>

using namespace std;

TH1F* h_sig_ttW = new TH1F("h_sig_ttW", "h_sig_ttW", 1000, -1., 1.);
TH1F* h_bkg_ttW = new TH1F("h_bkg_ttW", "h_bkg_ttW", 1000, -1., 1.);

TH1F* h_sig_ttH = new TH1F("h_sig_ttH", "h_sig_ttH", 1000, -1., 1.);
TH1F* h_bkg_ttH = new TH1F("h_bkg_ttH", "h_bkg_ttH", 1000, -1., 1.);

TH1F* h_sig_legacy = new TH1F("h_sig_legacy", "h_sig_legacy", 1000, -1., 1.);
TH1F* h_bkg_legacy = new TH1F("h_bkg_legacy", "h_bkg_legacy", 1000, -1., 1.);

TH1F* h_sig_ttW_alt3 = new TH1F("h_sig_ttW_alt3", "h_sig_ttW_alt3", 1000, -1., 1.);
TH1F* h_bkg_ttW_alt3 = new TH1F("h_bkg_ttW_alt3", "h_bkg_ttW_alt3", 1000, -1., 1.);

TH1F* h_sig_ttH_alt3 = new TH1F("h_sig_ttH_alt3", "h_sig_ttH_alt3", 1000, -1., 1.);
TH1F* h_bkg_ttH_alt3 = new TH1F("h_bkg_ttH_alt3", "h_bkg_ttH_alt3", 1000, -1., 1.);


TGraph ROC_curve(TH1F* h_sig, TH1F* h_bkg, TString name){

  float integral_sig = 0.;
  float integral_bkg = 0.;
  
  integral_sig = h_sig -> Integral();
  integral_bkg = h_bkg -> Integral();
  
  cout << h_sig -> Integral() << endl;
  cout << h_bkg -> Integral() << endl;

  float eff_sig = 0.;
  float eff_bkg = 0.;
  
  // ROC Curve
  TGraph ROC;
  int nbins = 1000;
  for (int b = 0; b < nbins; ++b){
    float S = h_sig -> Integral(b, nbins);
    float B = h_bkg -> Integral(b, nbins);
    eff_sig = S / integral_sig;
    eff_bkg = B / integral_bkg;
    ROC.SetPoint(b, 1 - eff_bkg, eff_sig);
  }
  ROC.SetName(name);
  return ROC;
}


void my_tmva_reader_mu(TString file){

  TFile* f = new TFile(file);

  TTree* t = (TTree*) f -> Get("Events");

  TMVA::Reader* reader = new TMVA::Reader();
  
  float mvaTTH;
  float miniPFRelIso_all;
  float genPartFlav;
  float isGlobal;
  float isTracker;
  float isPFcand;
  float looseId;
  float mediumId;
  float dxy;
  float dz;

  float pt;
  float eta;
  float pfRelIso03_all;
  float miniPFRelIso_chg;
  float miniRelIsoNeutral;
  float jetNDauCharged;
  float jetPtRelv2;
  float jetPtRatio;
  float jetBTagDeepFlavB;
  float sip3d;
  float log_dxy;
  float log_dz;
  float segmentComp;
  
  // Spectators
  reader->AddSpectator("Muon_mvaTTH",                   &mvaTTH);
  reader->AddSpectator("Muon_miniPFRelIso_all",         &miniPFRelIso_all);
  reader->AddSpectator("Muon_looseId",                  &looseId);
  reader->AddSpectator("Muon_genPartFlav",              &genPartFlav);
  reader->AddSpectator("Muon_isGlobal",                 &isGlobal);
  reader->AddSpectator("Muon_isTracker",                &isTracker);
  reader->AddSpectator("Muon_isPFcand",                 &isPFcand);
  reader->AddSpectator("Muon_mediumId",                 &mediumId);
  reader->AddSpectator("Muon_looseId",                  &looseId);
  reader->AddSpectator("Muon_dxy",                      &dxy);
  reader->AddSpectator("Muon_dz",                       &dz);
  // Training variables
  reader->AddVariable("Muon_pt",                                                                        &pt);
  reader->AddVariable("Muon_eta",                                                                       &eta);
  reader->AddVariable("Muon_pfRelIso03_all",                                                            &pfRelIso03_all);
  reader->AddVariable("Muon_miniPFRelIso_chg",                                                          &miniPFRelIso_chg);
  reader->AddVariable("Muon_miniRelIsoNeutral := Muon_miniPFRelIso_all - Muon_miniPFRelIso_chg",        &miniRelIsoNeutral);
  reader->AddVariable("Muon_jetNDauCharged",                                                            &jetNDauCharged);
  reader->AddVariable("Muon_jetPtRelv2",                                                                &jetPtRelv2);
  reader->AddVariable("Muon_jetPtRatio := min(1 / (1 + Muon_jetRelIso), 1.5)",                          &jetPtRatio);
  reader->AddVariable("Muon_jetBTagDeepFlavB := Muon_jetIdx > -1 ? Jet_btagDeepFlavB[Muon_jetIdx] : 0", &jetBTagDeepFlavB);
  reader->AddVariable("Muon_sip3d",                                                                     &sip3d);
  reader->AddVariable("Muon_log_dxy := log(abs(Muon_dxy))",                                             &log_dxy);
  reader->AddVariable("Muon_log_dz  := log(abs(Muon_dz))",                                              &log_dz);
  reader->AddVariable("Muon_segmentComp",                                                               &segmentComp);
  
  // cout << "Variables added to the reader" << endl;

  TString readers_names[2];
  if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2")){
    readers_names[0] = "dataset/weights/UL20_mu_TTH-like_2016_BDTG.weights.xml";
    readers_names[1] = "dataset/weights/UL20_mu_default_2016_BDTG.weights.xml";
  }
  else if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20")){
    readers_names[0] = "dataset/weights/UL20_mu_TTH-like_2016_VFP_BDTG.weights.xml";
    readers_names[1] = "dataset/weights/UL20_mu_default_2016_VFP_BDTG.weights.xml";
  }

  reader->BookMVA("BDTG_ttH-like",      readers_names[0]);
  reader->BookMVA("BDTG_default",       readers_names[1]);

  // cout << "MVAs booked by the reader" << endl;
  
  float mvaTTH_[40];
  float miniPFRelIso_all_[40];
  char  genPartFlav_[40];
  bool  isGlobal_[40];
  bool  isTracker_[40];
  bool  isPFcand_[40];
  bool  looseId_[40];
  bool  mediumId_[40];
  float dxy_[40];
  float dz_[40];

  float pt_[40];
  float eta_[40];
  float pfRelIso03_all_[40];
  float miniPFRelIso_chg_[40];
  float miniRelIsoAll_[40];
  char  jetNDauCharged_[40];
  float jetPtRelv2_[40];
  float jetPtRatio_[40];
  float jetBTagDeepFlavB_[40];
  float sip3d_[40];
  float segmentComp_[40];
  float jetRelIso_[40];
  int   jetIdx_[40];

  uint  nMuon_;

  t->SetBranchAddress("Muon_mvaTTH",           &mvaTTH_);
  t->SetBranchAddress("Muon_miniPFRelIso_all", &miniPFRelIso_all_);
  t->SetBranchAddress("Muon_looseId",          &looseId_);
  t->SetBranchAddress("Muon_genPartFlav",      &genPartFlav_);
  t->SetBranchAddress("Muon_isGlobal",         &isGlobal_);
  t->SetBranchAddress("Muon_isTracker",        &isTracker_);
  t->SetBranchAddress("Muon_isPFcand",         &isPFcand_);
  t->SetBranchAddress("Muon_mediumId",         &mediumId_);
  t->SetBranchAddress("Muon_looseId",          &looseId_);
  t->SetBranchAddress("Muon_dxy",              &dxy_);
  t->SetBranchAddress("Muon_dz",               &dz_);

  t->SetBranchAddress("Muon_pt",               &pt_);
  t->SetBranchAddress("Muon_eta",              &eta_);
  t->SetBranchAddress("Muon_pfRelIso03_all",   &pfRelIso03_all_);
  t->SetBranchAddress("Muon_miniPFRelIso_chg", &miniPFRelIso_chg_);
  t->SetBranchAddress("Muon_miniPFRelIso_all", &miniRelIsoAll_);
  t->SetBranchAddress("Muon_jetNDauCharged",   &jetNDauCharged_);
  t->SetBranchAddress("Muon_jetPtRelv2",       &jetPtRelv2_);
  t->SetBranchAddress("Muon_jetRelIso",        &jetRelIso_);
  t->SetBranchAddress("Muon_jetIdx",           &jetIdx_);
  t->SetBranchAddress("Jet_btagDeepFlavB",     &jetBTagDeepFlavB_);
  t->SetBranchAddress("Muon_sip3d",            &sip3d_);
  t->SetBranchAddress("Muon_segmentComp",      &segmentComp_);

  t->SetBranchAddress("nMuon",                 &nMuon_);

  // cout << "Branches set, starting loop" << endl;

  for (Long64_t ievt = 0; ievt < t->GetEntries(); ievt++) {

    if (ievt % 100000 == 0) cout << "--- ... Processing event: " << ievt <<std::endl;
    
    t->GetEntry(ievt);

    // cout << "Entry read" << endl;

    for (int i = 0; i < nMuon_; ++i){

      mvaTTH = mvaTTH_[i];
      miniPFRelIso_all = miniPFRelIso_all_[i];
      genPartFlav = genPartFlav_[i];
      isGlobal = isGlobal_[i];
      isTracker = isTracker_[i];
      isPFcand = isPFcand_[i];
      looseId = looseId_[i];
      mediumId = mediumId_[i];
      dxy = dxy_[i];
      dz = dz_[i];
      
      // cout << "Spectator variables assigned" << endl;

      pt = pt_[i];
      // cout << "pT assigned" << endl;
      eta = eta_[i];
      // cout << "eta assigned" << endl;
      pfRelIso03_all = pfRelIso03_all_[i];
      // cout << "pfrelisoall assigned" << endl;
      miniPFRelIso_chg = miniPFRelIso_chg_[i];
      // cout << "minipfrelisochg assigned" << endl;
      miniRelIsoNeutral = pfRelIso03_all_[i] - miniPFRelIso_chg_[i];
      // cout << "minirelisoneutral assigned" << endl;
      jetNDauCharged = jetNDauCharged_[i];
      // cout << "jetndau assigned" << endl;
      jetPtRelv2 = jetPtRelv2_[i];
      // cout << "jetptrel assigned" << endl;
      jetPtRatio = min(1. / (1. + jetRelIso_[i]), 1.5);
      // cout << "jetptratio assigned" << endl;
      jetBTagDeepFlavB = jetIdx_[i] > -1 ? jetBTagDeepFlavB_[jetIdx_[i]] : 0;
      // cout << "jetbtag assigned" << endl;
      sip3d = sip3d_[i];
      // cout << "sip3d assigned" << endl;
      log_dxy = log(abs(dxy));
      // cout << "logdxy assigned" << endl;
      log_dz = log(abs(dz));
      // cout << "logdz assigned" << endl;
      segmentComp = segmentComp_[i];
      
      // cout << "Variables assigned" << endl;
      
      // Preselections
      if (pt < 5)                          continue;
      if (miniPFRelIso_all > 0.4)          continue;
      if (sip3d > 8)                       continue;
      if (fabs(eta) > 2.4)                 continue;
      if (isGlobal != 1 && isTracker != 1) continue;
      if (isPFcand != 1)                   continue;
      if (fabs(dxy) > 0.05)                continue;
      if (fabs(dz) > 0.1)                  continue; 
      if (looseId != 1)                    continue;

      // cout << "Cuts passed" << endl;
      
      float BDTG_ttH;
      float BDTG_ttW;
      
      BDTG_ttH = reader->EvaluateMVA("BDTG_ttH-like");
      BDTG_ttW = reader->EvaluateMVA("BDTG_default");
      
      // if (BDTG_ttH == -999 || BDTG_ttW == -999)
      //   // cout << jetRelIso_ << " - " << jetPtRatio << " - " << BDTG_ttH << " - " << BDTG_ttW << " - " << mvaTTH << endl;
      
      // cout << "Filling histograms" << endl;
      
      // Fill signal histograms
      if (genPartFlav == 1 || genPartFlav == 15){
	h_sig_ttW    -> Fill((-1)*(mediumId == 0) + (mediumId == 1)*BDTG_ttW);
	h_sig_ttH    -> Fill(BDTG_ttH);
	h_sig_legacy -> Fill(mvaTTH);   
      }
      
      // Fill background histograms
      if (genPartFlav != 1 && genPartFlav != 15){
	h_bkg_ttW    -> Fill((-1)*(mediumId == 0) + (mediumId == 1)*BDTG_ttW);
	h_bkg_ttH    -> Fill(BDTG_ttH);
	h_bkg_legacy -> Fill(mvaTTH);  
      }
    }
  }

  // cout << "Preparing ROC curves" << endl;

  // Now prepare ROC curves
  TGraph ROC_ttW    = ROC_curve(h_sig_ttW,    h_bkg_ttW,    "ttW-like");
  TGraph ROC_ttH    = ROC_curve(h_sig_ttH,    h_bkg_ttH,    "ttH-like");
  TGraph ROC_legacy = ROC_curve(h_sig_legacy, h_bkg_legacy, "legacy");
  
  // cout << "Printing" << endl;

  // Cosmetics
  ROC_legacy.SetLineColor(kBlack);
  ROC_ttW.SetLineColor(kGreen+1);
  ROC_ttH.SetLineColor(kRed+1);

  ROC_legacy.SetLineWidth(2);
  ROC_ttW.SetLineWidth(2);
  ROC_ttH.SetLineWidth(2);
  

  TLegend* leg = new TLegend(0.15, 0.55, 0.75, 0.70);
  leg->SetLineColor(0);
  leg->AddEntry(&ROC_legacy, "Legacy mvaTTH", "l");  
  leg->AddEntry(&ROC_ttW,  "ttW-like presel", "l");  
  leg->AddEntry(&ROC_ttH,  "ttH-like presel", "l");  

  TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
  c2 -> cd();
  c2 -> SetLogx();
  ROC_ttW.SetTitle("ROC Curve for Muon MVA");
  ROC_ttW.GetXaxis()->SetRangeUser(0.90, 1.0);
  ROC_ttW.GetYaxis()->SetRangeUser(0.50, 1.0);
  ROC_ttW.GetXaxis()->SetTitle("Background Rejection");
  ROC_ttW.GetYaxis()->SetTitle("Signal Efficiency");
  ROC_ttW.Draw("AL");
  ROC_ttH.Draw("same");
  ROC_legacy.Draw("same");
  leg->Draw();
  
  TString canvas_name;
  if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2"))
    canvas_name = "ROC_mu_2016.png";
  else if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20"))
    canvas_name = "ROC_mu_2016_VFP.png";

  c2 -> Print(canvas_name);

}

void my_tmva_reader_el(TString file){

  TFile* f = new TFile(file);
  // cout << "File opened" << endl;

  TTree* t = (TTree*) f -> Get("Events");
  // cout << "Tree opened" << endl;

  TMVA::Reader* reader = new TMVA::Reader();
  
  float mvaTTH;
  float miniPFRelIso_all;
  float mvaFall17V2noIso_WPL;
  float lostHits;
  float genPartFlav;
  float genPartFlavAlt;
  float genPartFlavAlt2;
  float genPartFlavAlt3;
  float genPartFlavAlt4;
  float dxy;
  float dz;

  float pt;
  float eta;
  float pfRelIso03_all;
  float miniPFRelIso_chg;
  float miniRelIsoNeutral;
  float jetNDauCharged;
  float jetPtRelv2;
  float jetPtRatio;
  float jetBTagDeepFlavB;
  float sip3d;
  float log_dxy;
  float log_dz;
  float mvaFall17V2noIso;

  // Spectators
  reader->AddSpectator("Electron_mvaTTH",                   &mvaTTH);
  reader->AddSpectator("Electron_miniPFRelIso_all",         &miniPFRelIso_all);
  reader->AddSpectator("Electron_mvaFall17V2noIso_WPL",     &mvaFall17V2noIso_WPL);
  reader->AddSpectator("Electron_lostHits",                 &lostHits);
  reader->AddSpectator("Electron_genPartFlav",              &genPartFlav);
  reader->AddSpectator("Electron_genPartFlavAlt",           &genPartFlavAlt);
  reader->AddSpectator("Electron_genPartFlavAlt2",          &genPartFlavAlt2);
  reader->AddSpectator("Electron_genPartFlavAlt3",          &genPartFlavAlt3);
  reader->AddSpectator("Electron_genPartFlavAlt4",          &genPartFlavAlt4);
  reader->AddSpectator("Electron_dxy",                      &dxy);
  reader->AddSpectator("Electron_dz",                       &dz);
  // Training variables
  reader->AddVariable("Electron_pt",                                                                                &pt);
  reader->AddVariable("Electron_eta",                                                                               &eta);
  reader->AddVariable("Electron_pfRelIso03_all",                                                                    &pfRelIso03_all);
  reader->AddVariable("Electron_miniPFRelIso_chg",                                                                  &miniPFRelIso_chg);
  reader->AddVariable("Electron_miniRelIsoNeutral := Electron_miniPFRelIso_all - Electron_miniPFRelIso_chg",        &miniRelIsoNeutral);
  reader->AddVariable("Electron_jetNDauCharged",                                                                    &jetNDauCharged);
  reader->AddVariable("Electron_jetPtRelv2",                                                                        &jetPtRelv2);
  reader->AddVariable("Electron_jetPtRatio := min(1 / (1 + Electron_jetRelIso), 1.5)",                              &jetPtRatio);
  reader->AddVariable("Electron_jetBTagDeepFlavB := Electron_jetIdx > -1 ? Jet_btagDeepFlavB[Electron_jetIdx] : 0", &jetBTagDeepFlavB);
  reader->AddVariable("Electron_sip3d",                                                                             &sip3d);
  reader->AddVariable("Electron_log_dxy := log(abs(Electron_dxy))",                                                 &log_dxy);
  reader->AddVariable("Electron_log_dz  := log(abs(Electron_dz))",                                                  &log_dz);
  reader->AddVariable("Electron_mvaFall17V2noIso",                                                                  &mvaFall17V2noIso);
  
  // cout << "Variables added to the reader" << endl;

  TString readers_names[4];
  if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2")){
    readers_names[0] = "dataset/weights/UL20_el_TTH-like_2016_BDTG.weights.xml";
    readers_names[1] = "dataset/weights/UL20_el_default_2016_BDTG.weights.xml";
    readers_names[2] = "dataset/weights/UL20_el_TTH-like_2016_matchingAlt3_BDTG.weights.xml";
    readers_names[3] = "dataset/weights/UL20_el_default_2016_matchingAlt3_BDTG.weights.xml";
  }
  else if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20")){
    readers_names[0] = "dataset/weights/UL20_el_TTH-like_2016_VFP_BDTG.weights.xml";
    readers_names[1] = "dataset/weights/UL20_el_default_2016_VFP_BDTG.weights.xml";
    readers_names[2] = "dataset/weights/UL20_el_TTH-like_2016_VFP_matchingAlt3_BDTG.weights.xml";
    readers_names[3] = "dataset/weights/UL20_el_default_2016_VFP_matchingAlt3_BDTG.weights.xml";
  }

  reader->BookMVA("BDTG_ttH-like",      readers_names[0]);
  reader->BookMVA("BDTG_default",       readers_names[1]);
  reader->BookMVA("BDTG_ttH-like_alt3", readers_names[2]);
  reader->BookMVA("BDTG_default_alt3",  readers_names[3]);

  // cout << "MVAs booked by the reader" << endl;
  
  float mvaTTH_[40];
  float miniPFRelIso_all_[40];
  bool  mvaFall17V2noIso_WPL_[40];
  char  lostHits_[40];
  char  genPartFlav_[40];
  int   genPartFlavAlt_[40];
  int   genPartFlavAlt2_[40];
  int   genPartFlavAlt3_[40];
  int   genPartFlavAlt4_[40];
  float dxy_[40];
  float dz_[40];

  float pt_[40];
  float eta_[40];
  float pfRelIso03_all_[40];
  float miniPFRelIso_chg_[40];
  float miniRelIsoAll_[40];
  char  jetNDauCharged_[40];
  float jetPtRelv2_[40];
  float jetPtRatio_[40];
  float jetBTagDeepFlavB_[40];
  float sip3d_[40];
  float mvaFall17V2noIso_[40];
  float jetRelIso_[40];
  int   jetIdx_[40];

  uint  nElectron_;

  t->SetBranchAddress("Electron_mvaTTH",               &mvaTTH_);
  t->SetBranchAddress("Electron_miniPFRelIso_all",     &miniPFRelIso_all_);
  t->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", &mvaFall17V2noIso_WPL_);
  t->SetBranchAddress("Electron_lostHits",             &lostHits_);
  t->SetBranchAddress("Electron_genPartFlav",          &genPartFlav_);
  t->SetBranchAddress("Electron_genPartFlavAlt",       &genPartFlavAlt_);
  t->SetBranchAddress("Electron_genPartFlavAlt2",      &genPartFlavAlt2_);
  t->SetBranchAddress("Electron_genPartFlavAlt3",      &genPartFlavAlt3_);
  t->SetBranchAddress("Electron_genPartFlavAlt4",      &genPartFlavAlt4_);
  t->SetBranchAddress("Electron_dxy",                  &dxy_);
  t->SetBranchAddress("Electron_dz",                   &dz_);

  t->SetBranchAddress("Electron_pt",                   &pt_);
  t->SetBranchAddress("Electron_eta",                  &eta_);
  t->SetBranchAddress("Electron_pfRelIso03_all",       &pfRelIso03_all_);
  t->SetBranchAddress("Electron_miniPFRelIso_chg",     &miniPFRelIso_chg_);
  t->SetBranchAddress("Electron_miniPFRelIso_all",     &miniRelIsoAll_);
  t->SetBranchAddress("Electron_jetNDauCharged",       &jetNDauCharged_);
  t->SetBranchAddress("Electron_jetPtRelv2",           &jetPtRelv2_);
  t->SetBranchAddress("Electron_jetRelIso",            &jetRelIso_);
  t->SetBranchAddress("Electron_jetIdx",               &jetIdx_);
  t->SetBranchAddress("Jet_btagDeepFlavB",             &jetBTagDeepFlavB_);
  t->SetBranchAddress("Electron_sip3d",                &sip3d_);
  t->SetBranchAddress("Electron_mvaFall17V2noIso",     &mvaFall17V2noIso_);

  t->SetBranchAddress("nElectron",                     &nElectron_);

  // cout << "Branches set, starting loop" << endl;

  for (Long64_t ievt = 0; ievt < t->GetEntries(); ievt++) {

    //cout << "Event " << ievt << endl;
    
    if (ievt % 100000 == 0) cout << "--- ... Processing event: " << ievt <<std::endl;
    
    t->GetEntry(ievt);

    for (int i = 0; i < nElectron_; ++i){

      mvaTTH = mvaTTH_[i];
      miniPFRelIso_all = miniPFRelIso_all_[i];
      mvaFall17V2noIso_WPL = mvaFall17V2noIso_WPL_[i];
      lostHits = lostHits_[i];
      genPartFlav     = genPartFlav_[i];
      genPartFlavAlt  = genPartFlavAlt_[i];
      genPartFlavAlt2 = genPartFlavAlt2_[i];
      genPartFlavAlt3 = genPartFlavAlt3_[i];
      genPartFlavAlt4 = genPartFlavAlt4_[i];
      dxy = dxy_[i];
      dz = dz_[i];
    
      // cout << "Spectator variables assigned" << endl;
    
      pt = pt_[i];
      // cout << "pT assigned: " << pt << endl;
      eta = eta_[i];
      // cout << "eta assigned" << endl;
      pfRelIso03_all = pfRelIso03_all_[i];
      // cout << "pfrelisoall assigned" << endl;
      miniPFRelIso_chg = miniPFRelIso_chg_[i];
      // cout << "minipfrelisochg assigned" << endl;
      miniRelIsoNeutral = pfRelIso03_all_[i] - miniPFRelIso_chg_[i];
      // cout << "minirelisoneutral assigned" << endl;
      jetNDauCharged = jetNDauCharged_[i];
      // cout << "jetndau assigned" << endl;
      jetPtRelv2 = jetPtRelv2_[i];
      // cout << "jetptrel assigned" << endl;
      jetPtRatio = min(1. / (1. + jetRelIso_[i]), 1.5);
      // cout << "jetptratio assigned" << endl;
      jetBTagDeepFlavB = jetIdx_[i] > -1 ? jetBTagDeepFlavB_[jetIdx_[i]] : 0;
      // cout << "jetbtag assigned" << endl;
      sip3d = sip3d_[i];
      // cout << "sip3d assigned" << endl;
      log_dxy = log(abs(dxy));
      // cout << "logdxy assigned" << endl;
      log_dz = log(abs(dz));
      // cout << "logdz assigned" << endl;
      mvaFall17V2noIso = mvaFall17V2noIso_[i];
    
      // cout << "Variables assigned" << endl;
    
      // Preselections
      if (pt <= 5)                         continue;
      if (miniPFRelIso_all >= 0.4)         continue;
      if (sip3d >= 8)                      continue;
      if (fabs(eta) >= 2.5)                continue;
      if (lostHits >= 2)                    continue;
      if (fabs(dxy) >= 0.05)               continue;
      if (fabs(dz) >= 0.1)                 continue; 
      
      // cout << "Cuts passed" << endl;

      // cout << genPartFlavAlt3 << endl;

      float BDTG_ttH;
      float BDTG_ttW;
      float BDTG_ttH_alt3;
      float BDTG_ttW_alt3;

      BDTG_ttH = reader->EvaluateMVA("BDTG_ttH-like");
      BDTG_ttW = reader->EvaluateMVA("BDTG_default");
      BDTG_ttH_alt3 = reader->EvaluateMVA("BDTG_ttH-like_alt3");
      BDTG_ttW_alt3 = reader->EvaluateMVA("BDTG_default_alt3");
      
      // if (BDTG_ttH == -999 || BDTG_ttW == -999)
      //   // cout << jetRelIso_ << " - " << jetPtRatio << " - " << BDTG_ttH << " - " << BDTG_ttW << " - " << mvaTTH << endl;
      
      // cout << "Filling histograms" << endl;
      
      // Fill signal histograms
      if (genPartFlavAlt3 == 1 || genPartFlavAlt3 == 15){
	h_sig_ttW      -> Fill(BDTG_ttW);
	h_sig_ttH      -> Fill( (-1.)*(mvaFall17V2noIso_WPL == 0) + BDTG_ttH*(mvaFall17V2noIso_WPL == 1) );
	h_sig_ttW_alt3 -> Fill(BDTG_ttW_alt3);
	h_sig_ttH_alt3 -> Fill( (-1.)*(mvaFall17V2noIso_WPL == 0) + BDTG_ttH_alt3*(mvaFall17V2noIso_WPL == 1) );
	h_sig_legacy   -> Fill( (-1.)*(mvaFall17V2noIso_WPL == 0) + mvaTTH*(mvaFall17V2noIso_WPL == 1) );
      }

      // Fill background histograms
      if (genPartFlavAlt3 != 1 && genPartFlavAlt3 != 15){
	h_bkg_ttW      -> Fill(BDTG_ttW);
	h_bkg_ttH      -> Fill( (-1.)*(mvaFall17V2noIso_WPL == 0) + BDTG_ttH*(mvaFall17V2noIso_WPL == 1) ); // BDTG_ttH);
	h_bkg_ttW_alt3 -> Fill(BDTG_ttW_alt3);
	h_bkg_ttH_alt3 -> Fill( (-1.)*(mvaFall17V2noIso_WPL == 0) + BDTG_ttH_alt3*(mvaFall17V2noIso_WPL == 1) ); // BDTG_ttH_alt3);
	h_bkg_legacy   -> Fill( (-1.)*(mvaFall17V2noIso_WPL == 0) + mvaTTH*(mvaFall17V2noIso_WPL == 1) ); // mvaTTH);  
      }
    }
  } // end nElectron loop

  // cout << "Preparing ROC curves" << endl;

  // Now prepare ROC curves
  TGraph ROC_ttW      = ROC_curve(h_sig_ttW,      h_bkg_ttW,      "ttW-like");
  TGraph ROC_ttH      = ROC_curve(h_sig_ttH,      h_bkg_ttH,      "ttH-like");
  TGraph ROC_ttW_alt3 = ROC_curve(h_sig_ttW_alt3, h_bkg_ttW_alt3, "ttW-like alt3");
  TGraph ROC_ttH_alt3 = ROC_curve(h_sig_ttH_alt3, h_bkg_ttH_alt3, "ttH-like alt3");
  TGraph ROC_legacy   = ROC_curve(h_sig_legacy,   h_bkg_legacy,   "legacy");
  
  // cout << "Printing" << endl;

  // Cosmetics
  ROC_legacy.SetLineColor(kBlack);
  ROC_ttW.SetLineColor(kGreen+1);
  ROC_ttH.SetLineColor(kRed+1);
  ROC_ttW_alt3.SetLineColor(kBlue+1);
  ROC_ttH_alt3.SetLineColor(kYellow+1);

  ROC_legacy.SetLineWidth(2);
  ROC_ttW.SetLineWidth(2);
  ROC_ttH.SetLineWidth(2);
  ROC_ttW_alt3.SetLineWidth(2);
  ROC_ttH_alt3.SetLineWidth(2);
  

  TLegend* leg = new TLegend(0.15, 0.45, 0.75, 0.70);
  leg->SetLineColor(0);
  leg->AddEntry(&ROC_legacy,   "Legacy mvaTTH",   "l");  
  leg->AddEntry(&ROC_ttW,      "ttW-like presel", "l");  
  leg->AddEntry(&ROC_ttH,      "ttH-like presel", "l");  
  leg->AddEntry(&ROC_ttW_alt3, "ttW-like presel alt3", "l");  
  leg->AddEntry(&ROC_ttH_alt3, "ttH-like presel alt3", "l");  

  TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
  c2 -> cd();
  c2 -> SetLogx();
  ROC_ttW.SetTitle("ROC Curve for Electron MVA");
  ROC_ttW.GetXaxis()->SetRangeUser(0.95, 1.0);
  ROC_ttW.GetXaxis()->SetNdivisions(100, 200, 400);
  ROC_ttW.GetYaxis()->SetRangeUser(0.50, 1.0);
  ROC_ttW.GetXaxis()->SetTitle("Background Rejection");
  ROC_ttW.GetYaxis()->SetTitle("Signal Efficiency");
  ROC_ttW.Draw("AL");
  ROC_ttH.Draw("same");
  ROC_ttW_alt3.Draw("same");
  ROC_ttH_alt3.Draw("same");
  ROC_legacy.Draw("same");
  leg->Draw();
  
  TString canvas_name;
  if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2"))
    canvas_name = "ROC_el_2016.png";
  else if (file.Contains("TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20"))
    canvas_name = "ROC_el_2016_VFP.png";

  c2 -> Print(canvas_name);

}


void my_tmva_reader_2(TString year = "2016"){

  TString input_file;
  if (year == "2016")
    input_file = "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/09EA2721-660A-0F49-80CC-E8AB8E038BFA_Skim.root";
  else if (year == "2016_VFP")
    input_file = "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/E6F97EA6-EE6A-CA47-9116-E465AD9E8E2A_Skim.root";
  
  my_tmva_reader_mu(input_file);
  my_tmva_reader_el(input_file);

}
