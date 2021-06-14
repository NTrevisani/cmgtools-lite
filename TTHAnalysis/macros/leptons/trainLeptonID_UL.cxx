/* 

Latest directory with Alternatie Matching:
/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/

----------
-- 2016 --
----------

## Using selections as from tW talk

---------------
--   Muons   --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_mu_default_2016", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "2016")'
---------------
-- Electrons --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_default_2016", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "2016")'
---------------------------------------
-- Alternative electron gen-matching --
---------------------------------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_default_2016_matchingAlt3", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "2016")'


## Using selections as from ttH leptons module

---------------
--   Muons   --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_mu_TTH-like_2016", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "2016")'
---------------
-- Electrons --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_TTH-like_2016", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "2016")'
---------------------------------------
-- Alternative electron gen-matching --
---------------------------------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_TTH-like_2016_matchingAlt3", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_UL20_part2/", "", "2016")'


--------------
-- 2016_VFP --
--------------

## Using selections as from tW talk

---------------
--   Muons   --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_mu_default_2016_VFP", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "2016_VFP")'
---------------
-- Electrons --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_default_2016_VFP", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "2016_VFP")'
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_default_2016_VFP_matchingAlt3", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "2016_VFP")'

## Using selections as from ttH leptons module

---------------
--   Muons   --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_mu_TTH-like_2016_VFP", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "2016_VFP")'
---------------
-- Electrons --
---------------
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_TTH-like_2016_VFP", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "2016_VFP")'
root -l -b -q 'trainLeptonID_UL.cxx("training_nanoAOD", "UL20_el_TTH-like_2016_VFP_matchingAlt3", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "/eos/user/n/ntrevisa/ttH/leptonMVA/rootFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_2016_VFP_UL20/", "", "2016_VFP")'

*/

#include <assert.h>

void trainLeptonID_UL(TString folder,
		      TString name, 
		      TString sig1file, 
		      TString sig2file, 
		      TString bkg1file, 
		      TString bkg2file, 
		      TString year,
		      bool    doMultiClass    = false, 
		      TString file_for_sigW_1 = "", 
		      TString file_for_sigW_2 = "", 
		      TString file_for_bkgW_1 = "", 
		      TString file_for_bkgW_2 = "", 
		      double  int1s = 0, 
		      double  int2s = 0, 
		      double  int1b = 0,
		      double  int2b = 0) {
  
  
  std::vector<TString> files_names_strings;

  if (year == "2016"){
    files_names_strings.push_back("09EA2721-660A-0F49-80CC-E8AB8E038BFA_Skim.root");
    files_names_strings.push_back("3432CF59-5280-054C-82EC-01132DD8E196_Skim.root");
    files_names_strings.push_back("3A50DB1B-526E-F84B-B423-C25886D8B9E2_Skim.root");
    files_names_strings.push_back("45F2F0C0-8DA3-494D-B4F9-F0CA15746068_Skim.root");
    files_names_strings.push_back("51452850-C2FF-0640-815B-4F11464A5750_Skim.root");
    files_names_strings.push_back("529F2EDE-2C63-A846-A6F0-7804D4E4AFF9_Skim.root");
    files_names_strings.push_back("5C750601-4EDC-6F46-9B2A-D3B60F4BE9DF_Skim.root");
    files_names_strings.push_back("5F5FBB5C-170A-5E4E-A370-4EA224E273F3_Skim.root");
    files_names_strings.push_back("6546823D-A1FF-ED4F-BD37-17BA3BB84942_Skim.root");
    files_names_strings.push_back("67CECA18-E87E-524C-B356-EC9151A056AD_Skim.root");
    files_names_strings.push_back("861F48B6-DCE9-0D47-8A67-0FB0CCCF4337_Skim.root");
    files_names_strings.push_back("9255B0C9-7BE4-304C-9EBC-5416E9C6196E_Skim.root");
    files_names_strings.push_back("98F148C6-FD40-6A4C-A704-5E3398A792C6_Skim.root");
    files_names_strings.push_back("DC985944-2EBF-5344-A1D7-AA5B26156DAD_Skim.root");
    files_names_strings.push_back("E42CF18F-07E9-E444-AB5E-A5D3DAA83425_Skim.root");
    files_names_strings.push_back("E5EABD17-62DC-3746-A9C4-1F0FFA7606EB_Skim.root");
  }
  else if (year == "2016_VFP"){
    files_names_strings.push_back("293A49D5-E719-7C42-9EA0-82DCC0370101_Skim.root");
    files_names_strings.push_back("2E471CD5-EB2E-C34A-A520-295D8C8D2724_Skim.root");
    files_names_strings.push_back("5856A3BC-B4E7-6549-8619-56F0D4556DFE_Skim.root");
    // files_names_strings.push_back("5EADD416-E3D5-7E49-866C-B3D50B7C4D98_Skim.root");
    files_names_strings.push_back("7D1A92EC-C079-C24A-9F1D-796B96A5D168_Skim.root");
    files_names_strings.push_back("9D0F5653-79E6-8745-A43F-13A3DF09DD54_Skim.root");
    files_names_strings.push_back("A4659DA6-FE45-524A-BF9C-800BE082FCD7_Skim.root");
    files_names_strings.push_back("A5B1B5C6-E28B-9D42-96EB-A82D64C67186_Skim.root");
    files_names_strings.push_back("C0E66DF4-622B-404B-A6AF-8F5A57D30C54_Skim.root");
    files_names_strings.push_back("C737A51F-00A9-2E4A-90C6-13715F8BEA8D_Skim.root");
    files_names_strings.push_back("DD1CC8F8-D09D-5149-82C2-07F2776CF378_Skim.root");
    files_names_strings.push_back("E6F97EA6-EE6A-CA47-9116-E465AD9E8E2A_Skim.root");
    files_names_strings.push_back("EC87D37D-E22E-1545-B9D8-B4D00AFA0CB0_Skim.root");
    files_names_strings.push_back("EE4D1744-7A21-B94E-9387-203F249E3180_Skim.root");
    files_names_strings.push_back("F34167E6-D611-4A4A-8D7E-6FBACA1CC4A7_Skim.root");
  }

  TChain *dSig1 = new TChain("Events");
  TChain *dBg1  = new TChain("Events");
  for (int i = 0; i < files_names_strings.size(); ++i){
    dSig1 -> AddFile(sig1file + files_names_strings.at(i));
    dBg1  -> AddFile(bkg1file + files_names_strings.at(i));
  }
  
  TString mkdir = "mkdir -p " + folder;
  gSystem -> Exec(mkdir);

  TFile *fOut = new TFile(folder + "/" + name + ".root", "RECREATE");

  TString factory_conf = (!doMultiClass) ? "!V:!Color:Transformations=I" : "!V:!Color:Transformations=I:AnalysisType=Multiclass";
  
 
  TMVA::Factory *factory       = new TMVA::Factory(name, fOut, factory_conf.Data());
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");


  TCut lepton = "1";
    
  if (name.Contains("UL")) {
    if (name.Contains("_mu")) {
      // Spectators
      dataloader->AddSpectator("Muon_mvaTTH",                   'D');
      dataloader->AddSpectator("Muon_miniPFRelIso_all",         'D'); // relative isolation
      dataloader->AddSpectator("Muon_looseId",                  'D');
      dataloader->AddSpectator("Muon_genPartFlav",              'D');
      dataloader->AddSpectator("Muon_isGlobal",                 'D');
      dataloader->AddSpectator("Muon_isTracker",                'D');
      dataloader->AddSpectator("Muon_isPFcand",                 'D');
      dataloader->AddSpectator("Muon_mediumId",                 'D');
      dataloader->AddSpectator("Muon_looseId",                  'D');
      dataloader->AddSpectator("Muon_dxy",                      'D');
      dataloader->AddSpectator("Muon_dz",                       'D');
      // Training variables
      dataloader->AddVariable("Muon_pt",                        'D'); // pt
      dataloader->AddVariable("Muon_eta",                       'D'); // eta
      dataloader->AddVariable("Muon_pfRelIso03_all",            'D'); // PF relative isolation
      dataloader->AddVariable("Muon_miniPFRelIso_chg",          'D'); // mini relative isolation - charged component
      dataloader->AddVariable("Muon_miniRelIsoNeutral := Muon_miniPFRelIso_all - Muon_miniPFRelIso_chg", 'D'); // mini relative isolation - neutral component
      dataloader->AddVariable("Muon_jetNDauCharged",            'D'); // number of charged daughters of the closest jet
      dataloader->AddVariable("Muon_jetPtRelv2",                'D'); // relative momentum of the lepton with respect to the closest jet after subtracting the lepton
      dataloader->AddVariable("Muon_jetPtRatio := min(1 / (1 + Muon_jetRelIso), 1.5)", 'D'); // ratio between the lepton and the closest jet pTs (1+(PFIsoAll04/pt) if no matching jet)
      dataloader->AddVariable("Muon_jetBTagDeepFlavB := Muon_jetIdx > -1 ? Jet_btagDeepFlavB[Muon_jetIdx] : 0", 'D'); // Jet b-tagging probability of the associated jet 
      dataloader->AddVariable("Muon_sip3d",                         'D'); // 3D impact parameter significance wrt first PV
      dataloader->AddVariable("Muon_log_dxy := log(abs(Muon_dxy))", 'D'); // transverse impact parameter wrt PV
      dataloader->AddVariable("Muon_log_dz  := log(abs(Muon_dz))",  'D'); // longitudinal impact parameter wrt PV
      dataloader->AddVariable("Muon_segmentComp",                   'D'); // segment compatibility
      lepton += "Muon_miniPFRelIso_all < 0.4 && Muon_sip3d < 8";
      // From tW talk
      if (name.Contains("default")){
	lepton += "Muon_pt > 10 && abs(Muon_eta) < 2.4 && (Muon_isGlobal == 1 || Muon_isTracker == 1) && Muon_isPFcand == 1 && abs(Muon_dxy) < 0.05 && abs(Muon_dz) < 0.1 && Muon_mediumId == 1";
      }
      else if (name.Contains("TTH-like")){
	// From https://github.com/sscruz/cmgtools-lite/blob/104X_dev_nano/TTHAnalysis/python/tools/nanoAOD/ttH_modules.py
	lepton += "abs(Muon_eta) < 2.4 && Muon_pt > 5 && Muon_miniPFRelIso_all < 0.4 && Muon_sip3d < 8 && abs(Muon_dxy) < 0.05 && abs(Muon_dz) < 0.1"; // line 18
	lepton += "Muon_looseId == 1"; // line 97
      }
    }
    else if (name.Contains("_el")) {
      // Spectators
      dataloader->AddSpectator("Electron_mvaTTH",                       'D');
      dataloader->AddSpectator("Electron_miniPFRelIso_all",             'D'); // relative isolation
      dataloader->AddSpectator("Electron_mvaFall17V2noIso_WPL",         'D');
      dataloader->AddSpectator("Electron_lostHits",                     'D');
      dataloader->AddSpectator("Electron_genPartFlav",                  'D');
      dataloader->AddSpectator("Electron_genPartFlavAlt",               'D');
      dataloader->AddSpectator("Electron_genPartFlavAlt2",              'D');
      dataloader->AddSpectator("Electron_genPartFlavAlt3",              'D');
      dataloader->AddSpectator("Electron_genPartFlavAlt4",              'D');
      dataloader->AddSpectator("Electron_dxy",                          'D');
      dataloader->AddSpectator("Electron_dz",                           'D');
      // Training variables
      dataloader->AddVariable("Electron_pt",                            'D'); // pt
      dataloader->AddVariable("Electron_eta",                           'D'); // eta
      dataloader->AddVariable("Electron_pfRelIso03_all",                'D'); // PF relative isolation
      dataloader->AddVariable("Electron_miniPFRelIso_chg",              'D'); // mini relative isolation - charged component
      dataloader->AddVariable("Electron_miniRelIsoNeutral := Electron_miniPFRelIso_all - Electron_miniPFRelIso_chg", 'D'); // mini relative isolation - neutral component
      dataloader->AddVariable("Electron_jetNDauCharged",                'D'); // number of charged daughters of the closest jet
      dataloader->AddVariable("Electron_jetPtRelv2",                    'D'); // relative momentum of the lepton with respect to the closest jet after subtracting the lepton
      dataloader->AddVariable("Electron_jetPtRatio := min(1 / (1 + Electron_jetRelIso), 1.5)", 'D'); // ratio between the lepton and the closest jet pTs (1+(PFIsoAll04/pt) if no matching jet)
      dataloader->AddVariable("Electron_jetBTagDeepFlavB := Electron_jetIdx > -1 ? Jet_btagDeepFlavB[Electron_jetIdx] : 0", 'D'); // closest jet b-tagging probability 
      dataloader->AddVariable("Electron_sip3d",                             'D'); // 3D impact parameter significance wrt first PV
      dataloader->AddVariable("Electron_log_dxy := log(abs(Electron_dxy))", 'D'); // transverse impact parameter wrt PV
      dataloader->AddVariable("Electron_log_dz  := log(abs(Electron_dz))",  'D'); // longitudinal impact parameter wrt PV
      dataloader->AddVariable("Electron_mvaFall17V2noIso",                  'D'); // electron MVA ID value
      lepton += "Electron_miniPFRelIso_all < 0.4 && Electron_sip3d < 8";
      // From tW talk
      if (name.Contains("default")){
	lepton += "Electron_pt > 10 && abs(Electron_eta) < 2.5 && Electron_lostHits < 2 && abs(Electron_dxy) < 0.05 && abs(Electron_dz) < 0.1";
      }
      else if (name.Contains("TTH-like")){
	// From https://github.com/sscruz/cmgtools-lite/blob/104X_dev_nano/TTHAnalysis/python/tools/nanoAOD/ttH_modules.py
	lepton += "abs(Electron_eta) < 2.5 && Electron_pt > 5 && Electron_miniPFRelIso_all < 0.4 && Electron_sip3d < 8 && abs(Electron_dxy) < 0.05 && abs(Electron_dz) < 0.1 && Electron_mvaFall17V2noIso_WPL == 1"; // line 19
	lepton += "Electron_lostHits <= 1"; // line 97
      }
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
    dataloader->PrepareTrainingAndTestTree(lepton+" Muon_genPartFlav == 1 || Muon_genPartFlav == 15", lepton+" Muon_genPartFlav != 1 && Muon_genPartFlav != 15", "");
  else if (name.Contains("_el")){
    if (name.Contains("Alt")){
      if (name.Contains("Alt2"))
	// Alt2 --> matching with gen particles before FSR  
	dataloader->PrepareTrainingAndTestTree(lepton+" Electron_genPartFlavAlt2 == 1 || Electron_genPartFlavAlt2 == 15", lepton+" Electron_genPartFlavAlt2 != 1 && Electron_genPartFlavAlt2 != 15", "");
      else if (name.Contains("Alt3"))
	// Alt3 --> matching with electrons or photons with status 1 or 23 -- USE THIS ONE!
	dataloader->PrepareTrainingAndTestTree(lepton+" Electron_genPartFlavAlt3 == 1 || Electron_genPartFlavAlt3 == 15", lepton+" Electron_genPartFlavAlt3 != 1 && Electron_genPartFlavAlt3 != 15", "");
      else if (name.Contains("Alt4"))
	// Alt4 --> matching with electrons with status 1 or 23  
	dataloader->PrepareTrainingAndTestTree(lepton+" Electron_genPartFlavAlt4 == 1 || Electron_genPartFlavAlt4 == 15", lepton+" Electron_genPartFlavAlt4 != 1 && Electron_genPartFlavAlt4 != 15", "");
      else // Just "Alt" in name
	// Alt --> dressed lepton matching  
     	dataloader->PrepareTrainingAndTestTree(lepton+" Electron_genPartFlavAlt == 1 || Electron_genPartFlavAlt == 15", lepton+" Electron_genPartFlavAlt != 1 && Electron_genPartFlavAlt != 15", "");
    }
    else // Standard nanoAOD matching
      dataloader->PrepareTrainingAndTestTree(lepton+" Electron_genPartFlav == 1 || Electron_genPartFlav == 15", lepton+" Electron_genPartFlav != 1 && Electron_genPartFlav != 15", "");
  }
  else{
    std::cerr << "ERROR: must either be electron or muon." << std::endl; 
    return;
  }
  
  // else 
  //   dataloader->PrepareTrainingAndTestTree(lepton,"SplitMode=Random:NormMode=NumEvents:!V");

  //    if (!doMultiClass) factory->BookMethod( TMVA::Types::kLD, "LD", "!H:!V:VarTransform=None" );
    
  // Boosted Decision Trees with gradient boosting
  //TString BDTGopt = "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=9:UseNvars=9:MaxDepth=8";

  TString BDTGopt  = "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:UseNvars=9:MaxDepth=8";
  TString BDTGopt2 = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=1000:UseNvars=9:MaxDepth=4:MinNodeSize=0.05";

  // alternative options
  //TString BDTGopt = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:MaxDepth=3";
  //TString BDTGopt = "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000";
    

  // if (!doMultiClass) 
  BDTGopt += ":CreateMVAPdfs"; // Create Rarity distribution
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG", BDTGopt);
  // factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG2", BDTGopt2);
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  
  fOut->Close();
}
