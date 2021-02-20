/* 

root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL_mu_default_2016", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016.root", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016.root", "")'

root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL_el_default_2016", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016.root", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016.root", "")'

*/

#include <assert.h>

void trainLeptonID_UL(TString folder,
		      TString name, 
		      TString sig1file, 
		      TString sig2file, 
		      TString bkg1file, 
		      TString bkg2file, 
		      bool    doMultiClass    = false, 
		      TString file_for_sigW_1 = "", 
		      TString file_for_sigW_2 = "", 
		      TString file_for_bkgW_1 = "", 
		      TString file_for_bkgW_2 = "", 
		      double  int1s = 0, 
		      double  int2s = 0, 
		      double  int1b = 0,
		      double  int2b = 0) {
  
  
  TFile *_f_s1 = TFile::Open(sig1file.Data(), "read");
  TFile *_f_s2 = (sig2file == "") ? NULL : TFile::Open(sig2file.Data(), "read");
  TFile *_f_b1 = TFile::Open(bkg1file.Data(), "read");
  TFile *_f_b2 = (bkg2file == "") ? NULL : TFile::Open(bkg2file.Data(), "read");
 
  TTree *dSig1 = (TTree*) _f_s1->Get("Events");
  TTree *dSig2 = (_f_s2) ? ((TTree*) _f_s2->Get("Events")) : NULL;

  if (file_for_sigW_1 != "") 
    dSig1->AddFriend("wtree", file_for_sigW_1.Data());
  if (file_for_sigW_2 != "") 
    dSig2->AddFriend("wtree", file_for_sigW_2.Data());

  TTree *dBg1 = (TTree*) _f_b1->Get("Events");
  TTree *dBg2 = (_f_b2) ? ((TTree*) _f_b2->Get("Events")) : NULL;

  if (file_for_bkgW_1!="") 
    dBg1->AddFriend("wtree",file_for_bkgW_1.Data());
  if (file_for_bkgW_2!="") 
    dBg2->AddFriend("wtree",file_for_bkgW_2.Data());
  
  TString mkdir = "mkdir -p " + folder;
  gSystem -> Exec(mkdir);

  TFile *fOut = new TFile(folder + "/" + name + ".root", "RECREATE");

  TString factory_conf = (!doMultiClass) ? "!V:!Color:Transformations=I" : "!V:!Color:Transformations=I:AnalysisType=Multiclass";
  
 
  TMVA::Factory *factory       = new TMVA::Factory(name, fOut, factory_conf.Data());
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");


  TCut lepton = "1";
    
  if (name.Contains("UL")) {
    if (name.Contains("_mu")) {
      dataloader->AddVariable("Muon_pt",                        'D'); // pt
      dataloader->AddVariable("Muon_eta",                       'D'); // eta
      dataloader->AddVariable("Muon_pfRelIso03_all",            'D'); // relative isolation
      dataloader->AddVariable("Muon_miniPFRelIso_chg",          'D'); // mini relative isolation - charged component
      dataloader->AddVariable("Muon_miniRelIsoNeutral := Muon_miniPFRelIso_all - Muon_miniPFRelIso_chg", 'D'); // mini relative isolation - neutral component
      dataloader->AddVariable("Muon_jetNDauCharged",            'D'); // number of charged daughters of the closest jet
      dataloader->AddVariable("Muon_jetPtRelv2",                'D'); // relative momentum of the lepton with respect to the closest jet after subtracting the lepton
      dataloader->AddVariable("Muon_jetPtRatio := min(1 / (1 + Muon_jetRelIso), 1.5)", 'D'); // ratio between the lepton and the closest jet pTs (1+(PFIsoAll04/pt) if no matching jet)
      dataloader->AddVariable("Muon_jetBTagDeepFlavB := Muon_jetIdx > -1 ? Jet_btagDeepFlavB[Muon_jetIdx] : 0", 'D'); // closest jet b-tagging probability 
      dataloader->AddVariable("Muon_sip3d",                     'D'); // 3D impact parameter significance wrt first PV
      dataloader->AddVariable("Muon_dxy := log(abs(Muon_dxy))", 'D'); // transverse impact parameter wrt PV
      dataloader->AddVariable("Muon_dz  := log(abs(Muon_dz))",  'D'); // longitudinal impact parameter wrt PV
      dataloader->AddVariable("Muon_segmentComp",               'D'); // segment compatibility
      lepton += "Muon_miniPFRelIso_all < 0.4 && Muon_sip3d < 8";
      lepton += "Muon_pt > 10 && abs(Muon_eta) < 2.4 && (Muon_isGlobal == 1 || Muon_isTracker == 1) && Muon_isPFcand == 1 && abs(Muon_dxy) < 0.05 && abs(Muon_dz) < 0.1 && Muon_mediumId == 1";
    }
    else if (name.Contains("_el")) {
      dataloader->AddVariable("Electron_pt",                            'D'); // pt
      dataloader->AddVariable("Electron_eta",                           'D'); // eta
      dataloader->AddVariable("Electron_pfRelIso03_all",                'D'); // relative isolation
      dataloader->AddVariable("Electron_miniPFRelIso_chg",              'D'); // mini relative isolation - charged component
      dataloader->AddVariable("Electron_miniRelIsoNeutral := Electron_miniPFRelIso_all - Electron_miniPFRelIso_chg", 'D'); // mini relative isolation - neutral component
      dataloader->AddVariable("Electron_jetNDauCharged",                'D'); // number of charged daughters of the closest jet
      dataloader->AddVariable("Electron_jetPtRelv2",                    'D'); // relative momentum of the lepton with respect to the closest jet after subtracting the lepton
      dataloader->AddVariable("Electron_jetPtRatio := min(1 / (1 + Electron_jetRelIso), 1.5)", 'D'); // ratio between the lepton and the closest jet pTs (1+(PFIsoAll04/pt) if no matching jet)
      dataloader->AddVariable("Electron_jetBTagDeepFlavB := Electron_jetIdx > -1 ? Jet_btagDeepFlavB[Electron_jetIdx] : 0", 'D'); // closest jet b-tagging probability 
      dataloader->AddVariable("Electron_sip3d",                         'D'); // 3D impact parameter significance wrt first PV
      dataloader->AddVariable("Electron_dxy := log(abs(Electron_dxy))", 'D'); // transverse impact parameter wrt PV
      dataloader->AddVariable("Electron_dz  := log(abs(Electron_dz))",  'D'); // longitudinal impact parameter wrt PV
      dataloader->AddVariable("Electron_mvaFall17V2noIso",              'D'); // electron MVA ID value
      lepton += "Electron_miniPFRelIso_all < 0.4 && Electron_sip3d < 8";
      lepton += "Electron_pt > 10 && abs(Electron_eta) < 2.5 && Electron_lostHits < 2 && abs(Electron_dxy) < 0.05 && abs(Electron_dz) < 0.1";
    }
    else { 
      std::cerr << "ERROR: must either be electron or muon." << std::endl; 
      return; 
    }
  }


  // if (name.Contains("mu")) {
  //   lepton += "abs(LepGood_pdgId) == 13";
  // } 
  // else if (name.Contains("el")) {
  //   lepton += "abs(LepGood_pdgId) == 11";
  // }
  
  double wSig = 1.0, wBkg = 1.0;
  dataloader->AddSignalTree(dSig1, wSig);
  dataloader->AddBackgroundTree(dBg1, wBkg);


  // if (!doMultiClass){
  //   if (!dSig2) 
  //     dataloader->AddSignalTree(dSig1, wSig);
  //   else {
  //     std::cout << "Adding signal tree 1 and 2 with respective weights: " << wSig/int1s/2. << " " << wSig/int2s/2. << std::endl;
  //     dataloader->AddSignalTree(dSig1, wSig/int1s/2.);
  //     dataloader->AddSignalTree(dSig2, wSig/int2s/2.);
  //   }
  //   if (!dBg2) 
  //     dataloader->AddBackgroundTree(dBg1, wBkg);
  //   else {
  //     std::cout << "Adding background tree 1 and 2 with respective weights: " << wBkg/int1b/2. << " " << wBkg/int2b/2. << std::endl;
  //     dataloader->AddBackgroundTree(dBg1, wBkg/int1b/2.);
  //     dataloader->AddBackgroundTree(dBg2, wBkg/int2b/2.);
  //   }
  // }
  // else {
  //   if (!dSig2) 
  //     dataloader->AddTree(dSig1,"signal",wSig, "LepGood_mcMatchId!=0");
  //   else {
  //     dataloader->AddTree(dSig1, "signal", wSig/int1s/2., "LepGood_mcMatchId!=0");
  //     dataloader->AddTree(dSig2, "signal", wSig/int2s/2., "LepGood_mcMatchId!=0");
  //   }
  //   if (!dBg2) {
  //     dataloader->AddTree(dBg1, "bfake", wBkg, "LepGood_mcMatchId==0 && (abs(LepGood_mcMatchAny)==4 || abs(LepGood_mcMatchAny)==5)");
  //     dataloader->AddTree(dBg1, "light", wBkg, "LepGood_mcMatchId==0 && (abs(LepGood_mcMatchAny)<4 || abs(LepGood_mcMatchAny)>5)");
  //   }
  //   else {
  //     dataloader->AddTree(dBg1, "bfake", wBkg/int1b/2., "LepGood_mcMatchId==0 && (abs(LepGood_mcMatchAny)==4 || abs(LepGood_mcMatchAny)==5)");
  //     dataloader->AddTree(dBg1, "light", wBkg/int1b/2., "LepGood_mcMatchId==0 && (abs(LepGood_mcMatchAny)<4 || abs(LepGood_mcMatchAny)>5)");
  //     dataloader->AddTree(dBg2, "bfake", wBkg/int2b/2., "LepGood_mcMatchId==0 && (abs(LepGood_mcMatchAny)==4 || abs(LepGood_mcMatchAny)==5)");
  //     dataloader->AddTree(dBg2, "light", wBkg/int2b/2., "LepGood_mcMatchId==0 && (abs(LepGood_mcMatchAny)<4 || abs(LepGood_mcMatchAny)>5)");
  //   }
  // }
  

  // Do not use weights, for the moment (we are using the same sample as signal and background!)
  // if (file_for_sigW_1 != "" || file_for_sigW_2 != "")
  //   dataloader->SetSignalWeightExpression("addW*xsec*genWeight");
  // // else
  // //   dataloader->SetSignalWeightExpression("xsec*genWeight");
  // if (file_for_bkgW_1!="" || file_for_bkgW_2!="") 
  //   dataloader->SetBackgroundWeightExpression("addW*xsec*genWeight");
  // //  else dataloader->SetBackgroundWeightExpression("xsec*genWeight");
    
  //  if (!doMultiClass)

  // Flavour of genParticle for MC matching (for signal and bkg)
  if (name.Contains("_mu"))
    dataloader->PrepareTrainingAndTestTree(lepton+" Muon_genPartFlav == 1", lepton+" Muon_genPartFlav != 1", "" );
  else if (name.Contains("_el"))
    dataloader->PrepareTrainingAndTestTree(lepton+" Electron_genPartFlav == 1", lepton+" Electron_genPartFlav != 1", "" );
  else{
    std::cerr << "ERROR: must either be electron or muon." << std::endl; 
    return;
  }
  
  // else 
  //   dataloader->PrepareTrainingAndTestTree(lepton,"SplitMode=Random:NormMode=NumEvents:!V");

  //    if (!doMultiClass) factory->BookMethod( TMVA::Types::kLD, "LD", "!H:!V:VarTransform=None" );
    
  // Boosted Decision Trees with gradient boosting
  //TString BDTGopt = "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=9:UseNvars=9:MaxDepth=8";
  TString BDTGopt = "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:UseNvars=9:MaxDepth=8";

  // alternative options
  //TString BDTGopt = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:MaxDepth=3";
  //TString BDTGopt = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000";
    

  if (!doMultiClass) 
    BDTGopt += ":CreateMVAPdfs"; // Create Rarity distribution
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG", BDTGopt);
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  fOut->Close();
}
