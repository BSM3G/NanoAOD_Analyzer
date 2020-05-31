#include "Analyzer.h"
#include "Compression.h"
#include <regex>
#include <sstream>
#include <cmath>
#include <map>
//// Used to convert Enums to integers
#define ival(x) static_cast<int>(x)
//// BIG_NUM = sqrt(sizeof(int)) so can use diparticle convention of
//// index = BIG_NUM * i1 + i2
//// This insures easy way to extract indices
//// Needs to be changed if go to size_t instead (if want to play safe
#define BIG_NUM 46340

///// Macros defined to shorten code.  Made since lines used A LOT and repeative.  May change to inlines
///// if tests show no loss in speed
#define histAddVal2(val1, val2, name) ihisto.addVal(val1, val2, group, max, name, wgt);
#define histAddVal(val, name) ihisto.addVal(val, group, max, name, wgt);
#define SetBranch(name, variable) BOOM->SetBranchStatus(name, 1);  BOOM->SetBranchAddress(name, &variable);

typedef std::vector<int>::iterator vec_iter;



//////////////////////////////////////////////////////////////////
///////////////////CONSTANTS DEFINITONS///////////////////////////
//////////////////////////////////////////////////////////////////

//Filespace that has all of the .in files
const std::string PUSPACE = "Pileup/";


//////////PUBLIC FUNCTIONS////////////////////

const std::vector<CUTS> Analyzer::genCuts = {
  CUTS::eGTau, CUTS::eGNuTau, CUTS::eGTop,
  CUTS::eGElec, CUTS::eGMuon, CUTS::eGZ,
  CUTS::eGW, CUTS::eGHiggs, CUTS::eGJet, CUTS::eGBJet, CUTS::eGHadTau, CUTS::eGMatchedHadTau
};

const std::vector<CUTS> Analyzer::jetCuts = {
  CUTS::eRJet1,  CUTS::eRJet2,   CUTS::eRCenJet,
  CUTS::eR1stJet, CUTS::eR2ndJet, CUTS::eRBJet
};

const std::vector<CUTS> Analyzer::nonParticleCuts = {
  CUTS::eRVertex,CUTS::eRTrig1, CUTS::eRTrig2,
};
//01.16.19
const std::unordered_map<std::string, CUTS> Analyzer::cut_num = {
  {"NGenTau", CUTS::eGTau},                             {"NGenTop", CUTS::eGTop},
  {"NGenElectron", CUTS::eGElec},                       {"NGenMuon", CUTS::eGMuon},
  {"NGenZ", CUTS::eGZ},                                 {"NGenW", CUTS::eGW},
  {"NGenHiggs", CUTS::eGHiggs},                         {"NGenJet", CUTS::eGJet},
  {"NGenBJet", CUTS::eGBJet},                           {"NGenHadTau", CUTS::eGHadTau},
  {"NMatchedGenHadTau", CUTS::eGMatchedHadTau},
  {"NRecoMuon1", CUTS::eRMuon1},                        {"NRecoMuon2", CUTS::eRMuon2},
  {"NRecoElectron1", CUTS::eRElec1},                    {"NRecoElectron2",CUTS::eRElec2},
  {"NRecoTau1", CUTS::eRTau1},                          {"NRecoTau2", CUTS::eRTau2},
  {"NRecoJet1", CUTS::eRJet1},                          {"NRecoJet2", CUTS::eRJet2},
  {"NRecoCentralJet", CUTS::eRCenJet},                  {"NRecoBJet", CUTS::eRBJet},
  {"NRecoTriggers1", CUTS::eRTrig1},                    {"NRecoTriggers2", CUTS::eRTrig2},
  {"NRecoFirstLeadingJet", CUTS::eR1stJet},             {"NRecoSecondLeadingJet", CUTS::eR2ndJet},
  {"NDiMuonCombinations", CUTS::eDiMuon},               {"NDiElectronCombinations", CUTS::eDiElec},
  {"NDiTauCombinations", CUTS::eDiTau},                 {"NDiJetCombinations", CUTS::eDiJet},
  {"NMuon1Tau1Combinations", CUTS::eMuon1Tau1},         {"NMuon1Tau2Combinations", CUTS::eMuon1Tau2},
  {"NMuon2Tau1Combinations", CUTS::eMuon2Tau1},         {"NMuon2Tau2Combinations", CUTS::eMuon2Tau2},
  {"NElectron1Tau1Combinations", CUTS::eElec1Tau1},     {"NElectron1Tau2Combinations", CUTS::eElec1Tau2},
  {"NElectron2Tau1Combinations", CUTS::eElec2Tau1},     {"NElectron2Tau2Combinations", CUTS::eElec2Tau2},
  {"NMuon1Electron1Combinations", CUTS::eMuon1Elec1},   {"NMuon1Electron2Combinations", CUTS::eMuon1Elec2},
  {"NMuon2Electron1Combinations", CUTS::eMuon2Elec1},   {"NMuon2Electron2Combinations", CUTS::eMuon2Elec2},
  {"NElectron1Jet1Combinations", CUTS::eElec1Jet1},     {"NElectron1Jet2Combinations", CUTS::eElec1Jet2},
  {"NElectron2Jet1Combinations", CUTS::eElec2Jet1},     {"NElectron2Jet2Combinations", CUTS::eElec2Jet2},
  {"NLeadJetCombinations", CUTS::eSusyCom},             {"METCut", CUTS::eMET},
  {"NRecoWJet", CUTS::eRWjet},                          {"NRecoVertex", CUTS::eRVertex}
};



//////////////////////////////////////////////////////
//////////////////PUBLIC FUNCTIONS////////////////////
//////////////////////////////////////////////////////

///Constructor
Analyzer::Analyzer(std::vector<std::string> infiles, std::string outfile, bool setCR, std::string configFolder, std::string year) : goodParts(getArray()), genName_regex(".*([A-Z][^[:space:]]+)"){
  std::cout << "setup start" << std::endl;
  
  routfile = new TFile(outfile.c_str(), "RECREATE", outfile.c_str(), ROOT::CompressionSettings(ROOT::kLZMA, 9));
  
  add_metadata(infiles);

  BOOM= new TChain("Events");
  infoFile=0;

  for( std::string infile: infiles){
    BOOM->AddFile(infile.c_str());
  }


  nentries = (int) BOOM->GetEntries();
  BOOM->SetBranchStatus("*", 0);
  std::cout << "TOTAL EVENTS: " << nentries << std::endl;

  srand(0);
  
  filespace=configFolder;//"PartDet";
  filespace+="/";

  setupGeneral(year);
  
  CalculatePUSystematics = distats["Run"].bfind("CalculatePUSystematics");
  // New variable to do special PU weight calculation (2017)
  specialPUcalculation = distats["Run"].bfind("SpecialMCPUCalculation");
  
  // If data, always set specialcalculation to false
  if(isData){ 
    specialPUcalculation = false;
    std::cout << "This is Data!! Setting SpecialMCPUCalculation to False." << std::endl;
  }

  // New! Initialize pileup information and take care of exceptions if histograms not found.
  initializePileupInfo(specialPUcalculation, outfile);

  
  syst_names.push_back("orig");
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> tmp;
  syst_parts.push_back(tmp);
  if(!isData && distats["Systematics"].bfind("useSystematics")) {
    for(auto systname : distats["Systematics"].bset) {
      if( systname == "useSystematics")
        doSystematics= true;
      else {
        syst_names.push_back(systname);
        syst_parts.push_back(getArray());
      }
    }
  }else {
    doSystematics=false;
  }

  _Electron = new Electron(BOOM, filespace + "Electron_info.in", syst_names, year);

  _Muon     = new Muon(BOOM, filespace + "Muon_info.in", syst_names, year);

  _Tau      = new Taus(BOOM, filespace + "Tau_info.in", syst_names, year);

  _Jet      = new Jet(BOOM, filespace + "Jet_info.in", syst_names, year);

  _FatJet   = new FatJet(BOOM, filespace + "FatJet_info.in", syst_names, year);

  _MET      = new Met(BOOM, "MET" , syst_names, distats["Run"].dmap.at("MT2Mass"));

	
  // B-tagging scale factor stuff
  setupBJetSFInfo(_Jet->pstats["BJet"], year);  
  // Here the calibration module will be defined, which is then needed to define the readers below.
/*
  btagsfreader.load(calib, BTagEntry::FLAV_B, "comb");
  btagsfreaderup.load(calib, BTagEntry::FLAV_B, "comb");
  btagsfreaderdown.load(calib, BTagEntry::FLAV_B, "comb");
*/

  if(!isData) {
    std::cout<<"This is MC if not, change the flag!"<<std::endl;
    _Gen = new Generated(BOOM, filespace + "Gen_info.in", syst_names);
    _GenHadTau = new GenHadronicTaus(BOOM, filespace + "Gen_info.in", syst_names);
    allParticles= {_Gen,_GenHadTau,_Electron,_Muon,_Tau,_Jet,_FatJet};
  } else {
    std::cout<<"This is Data if not, change the flag!"<<std::endl;
    allParticles= {_Electron,_Muon,_Tau,_Jet,_FatJet};
  }

  particleCutMap[CUTS::eGElec]=_Electron;
  particleCutMap[CUTS::eGMuon]=_Muon;
  particleCutMap[CUTS::eGTau]=_Tau;

  std::vector<std::string> cr_variables;
  if(setCR) {
    char buf[64];
    read_info(filespace + "Control_Regions.in");
    crbins = pow(2.0, distats["Control_Region"].dmap.size());
    for(auto maper: distats["Control_Region"].dmap) {
      cr_variables.push_back(maper.first);
      sprintf(buf, "%.*G", 16, maper.second);
      cr_variables.push_back(buf);
    }
    if(isData) {
      if(distats["Control_Region"].smap.find("SR") == distats["Control_Region"].smap.end()) {
        std::cout << "Using Control Regions with data, but no signal region specified can lead to accidentially unblinding a study  before it should be.  Please specify a SR in the file PartDet/Control_Region.in" << std::endl;
        exit(1);
      } else if(distats["Control_Region"].smap.at("SR").length() != distats["Control_Region"].dmap.size()) {
        std::cout << "Signal Region specified incorrectly: check signal region variable to make sure the number of variables matches the number of signs in SR" << std::endl;
        exit(1);
      }
      int factor = 1;
      SignalRegion = 0;
      for(auto gtltSign: distats["Control_Region"].smap["SR"]) {
        if(gtltSign == '>') SignalRegion += factor;
        factor *= 2;
      }
      if(distats["Control_Region"].smap.find("Unblind") != distats["Control_Region"].smap.end()) {

        blinded = distats["Control_Region"].smap["Unblind"] == "false";
        std::cout << "we have " << blinded << std::endl;
      }
    }
  }
  //we update the root file if it exist so now we have to delete it:
  //std::remove(outfile.c_str()); // delete file
  histo = Histogramer(1, filespace+"Hist_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables);
  if(doSystematics)
    syst_histo=Histogramer(1, filespace+"Hist_syst_entries.in", filespace+"Cuts.in", outfile, isData, cr_variables,syst_names);
  
  systematics = Systematics(distats);
  
  jetScaleRes = JetScaleResolution("Pileup/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt", "",  "Pileup/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt", "Pileup/Spring16_25nsV6_MC_SF_AK4PFchs.txt");



  ///this can be done nicer
  //put the variables that you use here:
  zBoostTree["tau1_pt"] =0;
  zBoostTree["tau1_eta"]=0;
  zBoostTree["tau1_phi"]=0;
  zBoostTree["tau2_pt"] =0;
  zBoostTree["tau2_eta"]=0;
  zBoostTree["tau2_phi"]=0;
  zBoostTree["met"]     =0;
  zBoostTree["mt_tau1"] =0;
  zBoostTree["mt_tau2"] =0;
  zBoostTree["mt2"]     =0;
  zBoostTree["cosDphi1"]=0;
  zBoostTree["cosDphi2"]=0;
  zBoostTree["jet1_pt"] =0;
  zBoostTree["jet1_eta"]=0;
  zBoostTree["jet1_phi"]=0;
  zBoostTree["jet2_pt"] =0;
  zBoostTree["jet2_eta"]=0;
  zBoostTree["jet2_phi"]=0;
  zBoostTree["jet_mass"]=0;


  histo.createTree(&zBoostTree,"TauTauTree");


  if(setCR) {
    cuts_per.resize(histo.get_folders()->size());
    cuts_cumul.resize(histo.get_folders()->size());
  } else {
    cuts_per.resize(histo.get_cuts()->size());
    cuts_cumul.resize(histo.get_cuts()->size());
  }

  create_fillInfo();
  for(auto maper: distats["Control_Region"].dmap) {

    setupCR(maper.first, maper.second);
  }
  // check if we need to make gen level cuts to cross clean the samples:
  for(auto iselect : gen_selection){
    if(iselect.second){
      std::cout<<"Waning: The selection "<< iselect.first<< " is active!"<<std::endl;
    }
  }

  if(distats["Run"].bfind("DiscrByGenDileptonMass")){
    // If the gen-dilepton mass filter is on, check that the sample is DY. Otherwise, the filter won't be applied.
    isZsample = infiles[0].find("DY") != std::string::npos; 
    // if(isZsample){std::cout << "This is a DY sample " << std::endl;}
    // else{std::cout << "Not a DY sample " << std::endl;}
  }

  if(distats["Run"].bfind("InitializeMCSelection")){
     // std::cout << "MC selection initialized." << std::endl;
     initializeMCSelection(infiles);
  }


  initializeWkfactor(infiles);
  setCutNeeds();
  
  

  std::cout << "setup complete" << std::endl << std::endl;
  start = std::chrono::system_clock::now();
}

void Analyzer::add_metadata(std::vector<std::string> infiles){
  std::cout << "------------------------------------------------------------ " << std::endl;
  std::cout << "Copying minimal original trees from input files:"<<std::endl;
  
  // Define all the variables needed for this function
  TFile* rfile;
  std::string keyname;
  TTree* keytree;

  // Loop over the list of input files.
  for( std::string infile: infiles){
    std::cout << "File: " << infile << std::endl;

    // Open the input file
    rfile = TFile::Open(infile.c_str());
    routfile->cd();
    // Loop over all key stored in the current input file
    std::cout << "Processing key: " << std::endl;
    for(const auto&& inkey : *rfile->GetListOfKeys()){
      keyname = inkey->GetName();
      std::cout << "\t" << keyname << std::endl;

      if(keyname == "Events"){
        keytree= ((TTree*) rfile->Get(keyname.c_str()));     // Get the tree from file
        keytree->SetBranchStatus("*",0);                     // Disable all branches
        keytree->SetBranchStatus("run",1);                   // Enable only the branch named run
        originalTrees[keyname] = keytree->CopyTree("1","",1); // Add this tree to the original trees map. No selection nor option (1st and 2nd arg.) are applied, only 1 event is stored (3rd arg.)
      }else if(keyname == "MetaData" or keyname == "ParameterSets" or keyname == "Runs"){ // All branches from these trees are included but only 1 event is stored.
        originalTrees[keyname] = ((TTree*) rfile->Get(keyname.c_str()))->CopyTree("1","",1); 
      }else if(keyname == "LuminosityBlocks"){
        originalTrees[keyname] = ((TTree*) rfile->Get(keyname.c_str()))->CopyTree("1");  // All events for this tree are stored since they are useful when comparing with lumi filtering (JSON)
      }else if( std::string(inkey->ClassName()) == "TTree"){
        std::cout << "Not copying unknown tree " << inkey->GetName() << std::endl;
      }
    }
    routfile->cd();
    for(auto tree : originalTrees){
      tree.second->Write();
    }
    rfile->Close();
    delete rfile;
  }
  
  std::cout << "Finished copying minimal original trees." << std::endl;
  std::cout << "------------------------------------------------------------ " << std::endl;
}

std::unordered_map<CUTS, std::vector<int>*, EnumHash> Analyzer::getArray() {
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> rmap;
  for(auto e: Enum<CUTS>()) {
    rmap[e] = new std::vector<int>();
  }
  return rmap;
}



void Analyzer::create_fillInfo() {

  fillInfo["FillLeadingJet"] = new FillVals(CUTS::eSusyCom, FILLER::Dipart, _Jet, _Jet);
  fillInfo["FillGen"] =        new FillVals(CUTS::eGen, FILLER::Single, _Gen);
  fillInfo["FillTau1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTau2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuon1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuon2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectron1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectron2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);

  fillInfo["FillJet1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJet2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);
  fillInfo["FillBJet"] =       new FillVals(CUTS::eRBJet, FILLER::Single, _Jet);
  fillInfo["FillCentralJet"] = new FillVals(CUTS::eRCenJet, FILLER::Single, _Jet);
  fillInfo["FillWJet"] =       new FillVals(CUTS::eRWjet, FILLER::Single, _FatJet);

  fillInfo["FillDiElectron"] = new FillVals(CUTS::eDiElec, FILLER::Dipart, _Electron, _Electron);
  fillInfo["FillDiMuon"] =     new FillVals(CUTS::eDiMuon, FILLER::Dipart, _Muon, _Muon);
  fillInfo["FillDiTau"] =      new FillVals(CUTS::eDiTau, FILLER::Dipart, _Tau, _Tau);
  fillInfo["FillMetCuts"] =    new FillVals();
  fillInfo["FillDiJet"] =      new FillVals(CUTS::eDiJet, FILLER::Dipart, _Jet, _Jet);

  fillInfo["FillMuon1Tau1"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon1Tau2"] =       new FillVals(CUTS::eMuon1Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau1"] =       new FillVals(CUTS::eMuon2Tau1, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillMuon2Tau2"] =       new FillVals(CUTS::eMuon2Tau2, FILLER::Dipart, _Muon, _Tau);
  fillInfo["FillElectron1Tau1"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron1Tau2"] =   new FillVals(CUTS::eElec1Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau1"] =   new FillVals(CUTS::eElec2Tau1, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillElectron2Tau2"] =   new FillVals(CUTS::eElec2Tau2, FILLER::Dipart, _Electron, _Tau);
  fillInfo["FillMuon1Electron1"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon1Electron2"] =  new FillVals(CUTS::eMuon1Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron1"] =  new FillVals(CUTS::eMuon2Elec1, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillMuon2Electron2"] =  new FillVals(CUTS::eMuon2Elec2, FILLER::Dipart, _Muon, _Electron);
  fillInfo["FillElectron1Jet1"] =   new FillVals(CUTS::eElec1Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron1Jet2"] =   new FillVals(CUTS::eElec1Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron2Jet1"] =   new FillVals(CUTS::eElec2Jet1, FILLER::Dilepjet, _Electron, _Jet);
  fillInfo["FillElectron2Jet2"] =   new FillVals(CUTS::eElec2Jet2, FILLER::Dilepjet, _Electron, _Jet);

  //////I hate this solution so much.  Its terrible
  fillInfo["FillElectron1Electron2"] =     new FillVals(CUTS::eDiElec, FILLER::Single, _Electron, _Electron);
  fillInfo["FillMuon1Muon2"] =             new FillVals(CUTS::eDiMuon, FILLER::Single, _Muon, _Muon);
  fillInfo["FillTau1Tau2"] =               new FillVals(CUTS::eDiTau, FILLER::Single, _Tau, _Tau);

  //efficiency plots
  //In principal the efficiency plots should only be used, when also the object is used, but hey nobody knows!
  fillInfo["FillTauEfficiency1"] =       new FillVals(CUTS::eRTau1, FILLER::Single, _Tau);
  fillInfo["FillTauEfficiency2"] =       new FillVals(CUTS::eRTau2, FILLER::Single, _Tau);
  fillInfo["FillMuonEfficiency1"] =      new FillVals(CUTS::eRMuon1, FILLER::Single, _Muon);
  fillInfo["FillMuonEfficiency2"] =      new FillVals(CUTS::eRMuon2, FILLER::Single, _Muon);
  fillInfo["FillElectronEfficiency1"] =  new FillVals(CUTS::eRElec1, FILLER::Single, _Electron);
  fillInfo["FillElectronEfficiency2"] =  new FillVals(CUTS::eRElec2, FILLER::Single, _Electron);
  fillInfo["FillJetEfficiency1"] =       new FillVals(CUTS::eRJet1, FILLER::Single, _Jet);
  fillInfo["FillJetEfficiency2"] =       new FillVals(CUTS::eRJet2, FILLER::Single, _Jet);



  for(auto it: *histo.get_groups()) {
    if(fillInfo[it] == nullptr) fillInfo[it] = new FillVals();
  }

}

void Analyzer::setupCR(std::string var, double val) {
  std::smatch m;
  std::regex part ("^(.+)_(.+)$");
  if(std::regex_match(var, m, part)) {
    std::string name = m[1];
    std::string cut = "Fill" + name;
    if(fillInfo.find(cut) == fillInfo.end()) {
      std::cout << cut << " not found, put into fillInfo" << std::endl;
      exit(1);
    }
    std::cout << cut << " " << m[2] << " " << val << " " << name << std::endl;
    testVec.push_back(new CRTester(fillInfo.at(cut), m[2], val, name));
  } else {
    std::cout << "Could not process line: " << var << std::endl;
    exit(1);
  }

}





////destructor
Analyzer::~Analyzer() {
  clear_values();
  delete BOOM;
  delete _Electron;
  delete _Muon;
  delete _Tau;
  delete _Jet;
  if(!isData){
    delete _Gen;
    delete _GenHadTau;
  }

  for(auto fpair: fillInfo) {
    delete fpair.second;
    fpair.second=nullptr;
  }

  for(auto e: Enum<CUTS>()) {
    delete goodParts[e];
    goodParts[e]=nullptr;
  }
  //for(auto &it: syst_parts) {
    //for(auto e: Enum<CUTS>()) {
      //if( it[e] != nullptr) {
      //if(it.find(e) != it.end()){
        //delete it[e];
        //it[e]=nullptr;
      //}
      //}
    //}
  //}
  for(auto it: testVec){
    delete it;
    it=nullptr;
  }

}


///resets values so analysis can start
void Analyzer::clear_values() {

  for(auto e: Enum<CUTS>()) {
    goodParts[e]->clear();
  }
  //faster!!
  for(auto &it: syst_parts) {
    if (it.size() == 0) continue;
    for(auto e: Enum<CUTS>()) {
      it[e]->clear();
    }
  }
  if(infoFile!=BOOM->GetFile()){
    std::cout<<"New file!"<<std::endl;
    infoFile=BOOM->GetFile();
  }

  leadIndex=-1;
  maxCut = 0;
}

// New function: this sets up parameters that only need to be called once per event.
void Analyzer::setupEventGeneral(int nevent){

  // This class is an intermediate step called from preprocess that will set up all the variables that are common to the event and not particle specific.
  // We want to set those branches first here and then call BOOM->GetEntry(nevent) so that the variables change properly for each event.

  // For MC samples, set number of true pileup interactions, gen-HT and gen-weights.
  if(!isData){ 
    SetBranch("Pileup_nTrueInt",nTruePU);
    SetBranch("genWeight",gen_weight);
    if (BOOM->FindBranch("LHE_HT") != 0){
    	SetBranch("LHE_HT",generatorht);
    }
  }
  // Get the number of primary vertices, applies to both data and MC
  SetBranch("PV_npvs", bestVertices);

  // Finally, call get entry so all the branches assigned here are filled with the proper values for each event.
  BOOM->GetEntry(nevent);

  // Check that the sample does not have crazy values of nTruePU
  if(nTruePU < 100.0){
       // std::cout << "pileupntrueint = " << pileupntrueint << std::endl;
  }
  else{
    // std::cout << "event with abnormal pileup = " << pileupntrueint << std::endl;
    clear_values();
    return;
  }

  // Calculate the pu_weight for this event.
  pu_weight = (!isData && CalculatePUSystematics) ? hPU[(int)(nTruePU+1)] : 1.0;

  // Get the trigger decision vector.
  triggerDecision = false; // Reset the decision flag for each event.

  for(std::string triggname : triggerBranchesList){
    // std::cout << "Trigger name: " << triggname << std::endl;
    
    TBranch *triggerBranch = BOOM->GetBranch(triggname.c_str());
    triggerBranch->SetStatus(1);
    triggerBranch->SetAddress(&triggerDecision);

    // SetBranch(triggname.c_str(), triggerDecision);
    BOOM->GetEntry(nevent);

    // std::cout << "Decision = " << triggerDecision << std::endl;
    triggernamedecisions.push_back(triggerDecision);
    triggerBranch->ResetAddress();
  }

  /*
  for(size_t i=0; i < triggernamedecisions.size(); i++){
    std::cout << "Trigger decision #" << i << " = " << triggernamedecisions.at(i) << std::endl; 
  }
  */

}

bool Analyzer::passGenHTFilter(float genhtvalue){

  if(genhtvalue >= distats["Run"].dmap.at("LowerGenHtCut") && genhtvalue <= distats["Run"].dmap.at("UpperGenHtCut")){
    //std::cout << "genhtvalue = " << genhtvalue << ", passed genht filter " << std::endl;
    return true;
  }
  else{
    //std::cout << "genhtvalue = " << genhtvalue << ", failed genht filter " << std::endl;
    return false;
  }

}

bool Analyzer::passGenMassFilterZ(float mass_lowbound, float mass_upbound){

  if(!isZsample) return true;

   std::vector<uint> genleptonindices;
   int genpart_id = 0, genmotherpart_idx = 0, genmotherpart_id = 0;

   // std::cout << "---------------------------" << std::endl;
   // Loop over all generator level particles to look for the leptons that come from a Z:
   for(size_t idx = 0; idx < _Gen->size(); idx++) {
     // Get the particle PDG ID
     genpart_id = abs(_Gen->pdg_id[idx]);
     // Find the index of the mother particle
     genmotherpart_idx = _Gen->genPartIdxMother[idx];

     // Find the id of the mother particle using the index we just retrieved
     genmotherpart_id = _Gen->pdg_id[genmotherpart_idx];

     // std::cout << "part_idx = " << idx << ", genpart_id = " << genpart_id << ", status = " << _Gen->status[idx] << ", mother part_idx = " << genmotherpart_idx << ", mother genpart_id = " << genmotherpart_id << std::endl;

     // Only select those particles that are electrons (11), muons (13) or taus (15)
     if(! ( ( (genpart_id == 11 || genpart_id == 13) && _Gen->status[idx] == 1) || (genpart_id == 15 && _Gen->status[idx] == 2) ) )  continue;
     // std::cout << "This is a lepton" << std::endl;

     // Check that the mother particle is a Z boson (23)
     if(genmotherpart_id != 23) continue;
     // std::cout << "The mother is a Z boson" << std::endl;
     // Add the 4-momentum of this particle to the vector used later for dilepton mass calculation.
     genleptonindices.push_back(idx);

   }

   // Check that this vector only contains the information from two leptons.
   if(genleptonindices.size() != 2){
     // std::cout << "Not enough or too much leptons coming from the Z, returning false" << std::endl;
     return false;
   }
   // Check that they have the same mother (same mother index)
   else if(_Gen->genPartIdxMother[genleptonindices.at(0)] != _Gen->genPartIdxMother[genleptonindices.at(1)]){
     // std::cout << "The mothers of these two leptons are different, returning false" << std::endl;
     return false;
   }

   // std::cout << "We got a real Z -> ll event" << std::endl;
   // Get the dilepton mass using the indices from the leptons stored
   TLorentzVector lep1 = _Gen->p4(genleptonindices.at(0));
   TLorentzVector lep2 = _Gen->p4(genleptonindices.at(1));

   gendilepmass = (lep1 + lep2).M();

   if( gendilepmass >= mass_lowbound && gendilepmass <= mass_upbound){
     // std::cout << "Dilepton mass = " << (lep1 + lep2).M() << " between " << mass_lowbound << " and " << mass_upbound << " Passed the genMassfilter! " << std::endl;
     return true;
   }
   else{
     // std::cout << "Dilepton mass = " << (lep1 + lep2).M() << " not between " << mass_lowbound << " and " << mass_upbound << " Failed the genMassfilter! " << std::endl;
     return false;
   }

 }


bool Analyzer::checkGoodRunsAndLumis(int event){

    UInt_t run_num;  //NEW:  create run_num to store the run number.
    UInt_t luminosityBlock_num;  //NEW:  create luminosityBlock_num to store the number of the luminosity block.

    SetBranch("run",run_num);  //NEW:  define the branch for the run number.
    SetBranch("luminosityBlock",luminosityBlock_num);  //NEW:  define the branch for the luminosity block.

    BOOM->GetEntry(event);  //NEW:  get event.

    int key = run_num;
    int element = luminosityBlock_num;
    
    auto search = jsonlinedict.find(key);  //NEW:  see if the run number is in the json dictionary.
    
    if(search != jsonlinedict.end()){  //NEW:  this means that the run is in there.

		std::vector<int> lumivector;  //NEW:  going to make a container to store the lumis for the run.

		for (auto itr = jsonlinedict.begin(); itr != jsonlinedict.end(); itr++){ //NEW:  go through all the pairs of run, lumibound in the dictionary.
			if (itr->first != key) continue; //{ //NEW:  look for the run number.
			lumivector.push_back(itr->second);  //NEW:  grab all of the lumibounds corresponding to the run number.
		}

		// NEW:  going to go through pairs.  good lumi sections defined by [lumibound1, lumibound2], [lumibound3, lumibound4], etc.  
		// This is why we have to step through the check in twos.
		for(size_t lowbound = 0; lowbound < lumivector.size(); lowbound = lowbound+2){
			int upbound = lowbound+1; //NEW:  checking bounds that are side-by-side in the vector.
			if(!(element >= lumivector[lowbound] && element <= lumivector[upbound])) continue; //NEW:  if the lumisection is not within the bounds continue, check all lumi pairs.
			// if the lumisection falls within one of the lumipairs, you will make it here. Then, return true
			return true;
		}    
    }
    else{

	   return false;
    }

    return false;

}

bool Analyzer::passHEMveto2018(){

  bool hasJetinHEM = false;
  // Loop over all jets in the reco Jet collection before applying any selections.
  for(size_t index = 0; index < _Jet->size(); index++){
    // Get the 4-momentum of the jet:
    TLorentzVector jetP4 = _Jet->p4(index);

    //std::cout << "Jet #" << index << ": pt = " << jetP4.Pt() << ", eta = " << jetP4.Eta() << ", phi = " << jetP4.Phi() << std::endl;

    // Check if jet is in the HEM region
    if(jetP4.Pt() >= 30.0 && jetP4.Eta() >= -3.0 && jetP4.Eta() <= -1.3 && jetP4.Phi() >= -1.57 && jetP4.Phi() <= -0.87){
      //std::cout << "....... This jet is in the HEM region" << std::endl;
      hasJetinHEM = true;
      break;
    }
  }

  if(hasJetinHEM == true) return false;

  return true;

}

///Function that does most of the work.  Calculates the number of each particle
void Analyzer::preprocess(int event, std::string year){ // This function no longer needs to get the JSON dictionary as input.

  int test= BOOM->GetEntry(event);
  if(test<0){
    std::cout << "Could not read the event from the following file: "<<BOOM->GetFile()->GetNewUrl().Data() << std::endl;
  }

  for(Particle* ipart: allParticles){
    ipart->init();
  }
  _MET->init();

  active_part = &goodParts;
  if(!select_mc_background()){
    //we will put nothing in good particles
    clear_values();
    return;
  }

  // Call the new function setupEventGeneral: this will set generatorht, pu weight and genweight
  setupEventGeneral(event);

  
  if(!isData){ // Do everything that corresponds only to MC

    // Initialize the lists of generator-level particles.
    _Gen->setOrigReco();
    _GenHadTau->setOrigReco();

    getGoodGen(_Gen->pstats["Gen"]);
    getGoodGenHadronicTaus(_GenHadTau->pstats["Gen"]);
    getGoodGenHadronicTauNeutrinos(_Gen->pstats["Gen"]);
    getGoodGenBJet(); //01.16.19

    //--- filtering inclusive HT-binned samples: must be done after setupEventGeneral --- // 

    if(distats["Run"].bfind("DiscrByGenHT")){ 
      //std::cout << "generatorht = " << generatorht << std::endl;
      if(passGenHTFilter(generatorht) == false){
        clear_values();
        return;
      }
    }

    //--- filtering inclusive HT-binned samples: must be done after  _Gen->setOrigReco() --- // 
     if(distats["Run"].bfind("DiscrByGenDileptonMass")){
       if(passGenMassFilterZ(distats["Run"].pmap.at("GenDilepMassRange").first, distats["Run"].pmap.at("GenDilepMassRange").second) == false){
         clear_values();
         return;
       }
     }

  }
  else if(isData){

  	// If you want to filter data by good run and lumi sections, turn on this option in Run_info.in
  	// SingleMuon data sets need to be filtered. SingleElectron does not (?).
  	if(distats["Run"].bfind("FilterDataByGoldenJSON")){
	    if(checkGoodRunsAndLumis(event) == false){
        //std::cout << "Thrown away by JSON filter" << std::endl;
	    	clear_values();
	    	return;
	    }
	}
  }
	
  // Call the new function passMetFilters
  // Apply here the MET filters, in case the option is turned on. It applies to both data and MC
  applymetfilters = distats["Run"].bfind("ApplyMetFilters");
  if(applymetfilters){
  	passedmetfilters = passMetFilters(year, event);	
  	if(!passedmetfilters){
      // std::cout << "Thrown away by MET filter" << std::endl;
  		clear_values();
  		return;
  	}
  }

  // Apply HEM veto for 2018 if the flag is on.
  bool checkHEM = distats["Run"].bfind("ApplyHEMVeto2018");
  if(checkHEM == 1 || checkHEM == true){
    if(passHEMveto2018() == false){
      // std::cout << "This event didn't pass the HEM veto... throwing it away." << std::endl;
      clear_values();
      return;
    }
  }

  // std::cout << "------------" << std::endl;

  // ------- Number of primary vertices requirement -------- // 
  active_part->at(CUTS::eRVertex)->resize(bestVertices);

  // ---------------- Trigger requirement ------------------ //
  TriggerCuts(CUTS::eRTrig1);

  ////check update met is ok
  for(size_t i=0; i < syst_names.size(); i++) {
     //////Smearing
    smearLepton(*_Electron, CUTS::eGElec, _Electron->pstats["Smear"], distats["Electron_systematics"], i);
    smearLepton(*_Muon, CUTS::eGMuon, _Muon->pstats["Smear"], distats["Muon_systematics"], i);
    smearLepton(*_Tau, CUTS::eGTau, _Tau->pstats["Smear"], distats["Tau_systematics"], i);

    smearJet(*_Jet,CUTS::eGJet,_Jet->pstats["Smear"], i);
    smearJet(*_FatJet,CUTS::eGJet,_FatJet->pstats["Smear"], i);
    updateMet(i);

  }

  for(size_t i=0; i < syst_names.size(); i++) {
    std::string systname = syst_names.at(i);
    for( auto part: allParticles) part->setCurrentP(i);
    _MET->setCurrentP(i);
    getGoodParticles(i);
  }
  active_part = &goodParts;
  
  if( event < 10 || ( event < 100 && event % 10 == 0 ) ||
  ( event < 1000 && event % 100 == 0 ) ||
  ( event < 10000 && event % 1000 == 0 ) ||
  ( event >= 10000 && event % 10000 == 0 ) ) {
    std::cout << std::setprecision(2)<<event << " Events analyzed "<< static_cast<double>(event)/nentries*100. <<"% done"<<std::endl;
    std::cout << std::setprecision(5);
  }
}


void Analyzer::getGoodParticles(int syst){

  std::string systname=syst_names.at(syst);
  if(syst == 0) active_part = &goodParts;
  else active_part=&syst_parts.at(syst);
    //    syst=syst_names[syst];



  // // SET NUMBER OF RECO PARTICLES
  // // MUST BE IN ORDER: Muon/Electron, Tau, Jet
  getGoodRecoLeptons(*_Electron, CUTS::eRElec1, CUTS::eGElec, _Electron->pstats["Elec1"],syst);
  getGoodRecoLeptons(*_Electron, CUTS::eRElec2, CUTS::eGElec, _Electron->pstats["Elec2"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon1, CUTS::eGMuon, _Muon->pstats["Muon1"],syst);
  getGoodRecoLeptons(*_Muon, CUTS::eRMuon2, CUTS::eGMuon, _Muon->pstats["Muon2"],syst);
  //getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGTau, _Tau->pstats["Tau1"],syst);
  //getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGTau, _Tau->pstats["Tau2"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau1, CUTS::eGHadTau, _Tau->pstats["Tau1"],syst);
  getGoodRecoLeptons(*_Tau, CUTS::eRTau2, CUTS::eGHadTau, _Tau->pstats["Tau2"],syst);
  getGoodRecoBJets(CUTS::eRBJet, _Jet->pstats["BJet"],syst); //01.16.19
  //getGoodRecoJets(CUTS::eRBJet, _Jet->pstats["BJet"],syst);
  getGoodRecoJets(CUTS::eRJet1, _Jet->pstats["Jet1"],syst);
  getGoodRecoJets(CUTS::eRJet2, _Jet->pstats["Jet2"],syst);
  getGoodRecoJets(CUTS::eRCenJet, _Jet->pstats["CentralJet"],syst);
  getGoodRecoJets(CUTS::eR1stJet, _Jet->pstats["FirstLeadingJet"],syst);
  getGoodRecoJets(CUTS::eR2ndJet, _Jet->pstats["SecondLeadingJet"],syst);

  getGoodRecoFatJets(CUTS::eRWjet, _FatJet->pstats["Wjet"],syst);
  //  treatMuons_Met(systname);

  ///VBF Susy cut on leadin jets
  VBFTopologyCut(distats["VBFSUSY"],syst);

  /////lepton lepton topology cuts
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1,CUTS::eRTau1, CUTS::eElec1Tau1, distats["Electron1Tau1"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau1, CUTS::eElec2Tau1, distats["Electron2Tau1"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec1, CUTS::eRTau2, CUTS::eElec1Tau2, distats["Electron1Tau2"],syst);
  getGoodLeptonCombos(*_Electron, *_Tau, CUTS::eRElec2, CUTS::eRTau2, CUTS::eElec2Tau2, distats["Electron2Tau2"],syst);

  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau1, CUTS::eMuon1Tau1, distats["Muon1Tau1"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon1, CUTS::eRTau2, CUTS::eMuon1Tau2, distats["Muon1Tau2"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau1, CUTS::eMuon2Tau1, distats["Muon2Tau1"],syst);
  getGoodLeptonCombos(*_Muon, *_Tau, CUTS::eRMuon2, CUTS::eRTau2, CUTS::eMuon2Tau2, distats["Muon2Tau2"],syst);

  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec1, CUTS::eMuon1Elec1, distats["Muon1Electron1"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon1, CUTS::eRElec2, CUTS::eMuon1Elec2, distats["Muon1Electron2"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec1, CUTS::eMuon2Elec1, distats["Muon2Electron1"],syst);
  getGoodLeptonCombos(*_Muon, *_Electron, CUTS::eRMuon2, CUTS::eRElec2, CUTS::eMuon2Elec2, distats["Muon2Electron2"],syst);

  ////DIlepton topology cuts
  getGoodLeptonCombos(*_Tau, *_Tau, CUTS::eRTau1, CUTS::eRTau2, CUTS::eDiTau, distats["DiTau"],syst);
  getGoodLeptonCombos(*_Electron, *_Electron, CUTS::eRElec1, CUTS::eRElec2, CUTS::eDiElec, distats["DiElectron"],syst);
  getGoodLeptonCombos(*_Muon, *_Muon, CUTS::eRMuon1, CUTS::eRMuon2, CUTS::eDiMuon, distats["DiMuon"],syst);

  //
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec1, CUTS::eRJet1, CUTS::eElec1Jet1, distats["Electron1Jet1"],syst);
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec1, CUTS::eRJet2, CUTS::eElec1Jet2, distats["Electron1Jet2"],syst);
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec2, CUTS::eRJet1, CUTS::eElec2Jet1, distats["Electron2Jet1"],syst);
  getGoodLeptonJetCombos(*_Electron, *_Jet, CUTS::eRElec2, CUTS::eRJet2, CUTS::eElec2Jet2, distats["Electron2Jet2"],syst);

  ////Dijet cuts
  getGoodDiJets(distats["DiJet"],syst);

}


void Analyzer::fill_efficiency() {
  //cut efficiency
  const std::vector<CUTS> goodGenLep={CUTS::eGElec,CUTS::eGMuon,CUTS::eGTau};
  //just the lepton 1 for now
  const std::vector<CUTS> goodRecoLep={CUTS::eRElec1,CUTS::eRMuon1,CUTS::eRTau1};



  for(size_t igen=0;igen<goodGenLep.size();igen++){
    Particle* part =particleCutMap.at(goodGenLep[igen]);
    CUTS cut=goodRecoLep[igen];
    std::smatch mGen;
    std::string tmps=part->getName();
    std::regex_match(tmps, mGen, genName_regex);
    //loop over all gen leptons
    for(int iigen : *active_part->at(goodGenLep[igen])){


      int foundReco=-1;
      for(size_t ireco=0; ireco<part->size(); ireco++){
        if(part->p4(ireco).DeltaR(_Gen->p4(iigen))<0.3){
          foundReco=ireco;
        }
      }
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Pt", _Gen->pt(iigen), foundReco>=0,0);
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Eta",_Gen->eta(iigen),foundReco>=0,0);
      histo.addEffiency("eff_Reco_"+std::string(mGen[1])+"Phi",_Gen->phi(iigen),foundReco>=0,0);
      if(foundReco>=0){
        bool id_particle= (find(active_part->at(cut)->begin(),active_part->at(cut)->end(),foundReco)!=active_part->at(cut)->end());
        histo.addEffiency("eff_"+std::string(mGen[1])+"Pt", _Gen->pt(iigen), id_particle,0);
        histo.addEffiency("eff_"+std::string(mGen[1])+"Eta",_Gen->eta(iigen),id_particle,0);
        histo.addEffiency("eff_"+std::string(mGen[1])+"Phi",_Gen->phi(iigen),id_particle,0);
      }
    }
  }

  //for(Particle* part : allParticles){





    //regex genName_regex(".*([A-Z][^[:space:]]+)");
    //smatch mGen;
    //std::string tmps=part->getName();
    //regex_match(tmps, mGen, genName_regex);
    ////no efficiency for gen particles
    //if(part->getName().find("Gen") != std::string::npos)
      //continue;
    ////we don't want to make met efficiency plots
    //if(particleCutMap.find(part) == particleCutMap.end())
      //continue;
    //if(part->cutMap.find(part->type) == part->cutMap.end())
      //continue;
    //for(size_t i=0; i < part->size(); i++){
      ////make match to gen
      //if(matchLeptonToGen(part->p4(i), part->pstats.at("Smear") ,part->cutMap.at(part->type)) == TLorentzVector(0,0,0,0)) continue;
      ////check if the particle is part of the reco
      //for(CUTS cut:  particleCutMap.at(part).first){
        //bool id_particle= (find(active_part->at(cut)->begin(),active_part->at(cut)->end(),i)!=active_part->at(cut)->end());
        //histo.addEffiency("eff_"+std::string(mGen[1])+"Pt",part->pt(i),id_particle,0);
        //histo.addEffiency("eff_"+std::string(mGen[1])+"Eta",part->eta(i),id_particle,0);
        //histo.addEffiency("eff_"+std::string(mGen[1])+"Phi",part->phi(i),id_particle,0);
      //}
    //}
  //}
}


////Reads cuts from Cuts.in file and see if the event has enough particles
bool Analyzer::fillCuts(bool fillCounter) {
  const std::unordered_map<std::string,std::pair<int,int> >* cut_info = histo.get_cuts();
  const std::vector<std::string>* cut_order = histo.get_cutorder();

  bool prevTrue = true;

  maxCut=0;
  //  std::cout << active_part << std::endl;;

  for(size_t i = 0; i < cut_order->size(); i++) {
    std::string cut = cut_order->at(i);
    if(isData && cut.find("Gen") != std::string::npos){
      maxCut += 1;
      continue;
    }
    int min= cut_info->at(cut).first;
    int max= cut_info->at(cut).second;
    int nparticles = active_part->at(cut_num.at(cut))->size();
    //if(!fillCounter) std::cout << cut << ": " << nparticles << " (" << min << ", " << max << ")" <<std::endl;
    if( (nparticles >= min) && (nparticles <= max || max == -1)) {
      if((cut_num.at(cut) == CUTS::eR1stJet || cut_num.at(cut) == CUTS::eR2ndJet) && active_part->at(cut_num.at(cut))->at(0) == -1 ) {
        //cout<<"here   "<<std::endl;
        prevTrue = false;
        continue;  ////dirty dirty hack
      }
      if(fillCounter && crbins == 1) {
        cuts_per[i]++;
        cuts_cumul[i] += (prevTrue) ? 1 : 0;
        maxCut += (prevTrue) ? 1 : 0;
      }else{
        maxCut += (prevTrue) ? 1 : 0;
      }
    }else {
      //cout<<"here 2  "<<std::endl;
      prevTrue = false;
    }
  }

  if(crbins != 1) {
    if(!prevTrue) {
      maxCut = -1;
      return prevTrue;
    }

    int factor = crbins;
    for(auto tester: testVec) {
      factor /= 2;
      /////get variable value from maper.first.
      if(tester->test(this)) { ///pass cut
        maxCut += factor;
      }
    }
    if(isData && blinded && maxCut == SignalRegion) return false;
    cuts_per[maxCut]++;
  }


  return prevTrue;
}



///Prints the number of events that passed each cut per event and cumulatively
//done at the end of the analysis
void Analyzer::printCuts() {
  std::vector<std::string> cut_order;
  if(crbins > 1) cut_order = *(histo.get_folders());
  else cut_order = *(histo.get_cutorder());
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  double run_time_real=elapsed_seconds.count();


  std::cout.setf(std::ios::floatfield,std::ios::fixed);
  std::cout<<std::setprecision(3);
  std::cout << "\n";
  std::cout << "Selection Efficiency " << "\n";
  std::cout << "Total events: " << nentries << "\n";
  std::cout << "\n";
  std::cout << "Run Time (real): " <<run_time_real <<" s\n";
  std::cout << "Time per 1k Events (real): " << run_time_real/(nentries/1000) <<" s\n";
  std::cout << "Events/s: " << static_cast<double>(nentries)/(run_time_real) <<" 1/s (real) \n";
  std::cout << "                        Name                  Indiv.";
  if(crbins == 1) std::cout << "            Cumulative";
  std::cout << std::endl << "---------------------------------------------------------------------------\n";
  for(size_t i = 0; i < cut_order.size(); i++) {
    std::cout << std::setw(28) << cut_order.at(i) << "    ";
    if(isData && cut_order.at(i).find("Gen") != std::string::npos) std::cout << "Skipped" << std::endl;
    else if(crbins != 1 && blinded && i == (size_t)SignalRegion) std::cout << "Blinded Signal Region" << std::endl;
    else {
      std::cout << std::setw(10) << cuts_per.at(i) << "  ( " << std::setw(5) << ((float)cuts_per.at(i)) / nentries << ") ";
      if(crbins == 1) std::cout << std::setw(12) << cuts_cumul.at(i) << "  ( " << std::setw(5) << ((float)cuts_cumul.at(i)) / nentries << ") ";

      std::cout << std::endl;
    }
  }
  std::cout <<std::setprecision(5);
  std::cout << "---------------------------------------------------------------------------\n";

  //write all the histograms
  //attention this is not the fill_histogram method from the Analyser
  histo.fill_histogram(routfile);
  if(doSystematics)
    syst_histo.fill_histogram(routfile);

}

/////////////PRIVATE FUNCTIONS////////////////


bool Analyzer::select_mc_background(){
  //will return true if Z* mass is smaller than 200GeV
  if(_Gen == nullptr){
    return true;
  }
  if(gen_selection["DY_noMass_gt_200"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id[i])==11 or abs(_Gen->pdg_id[i])==13 or abs(_Gen->pdg_id[i])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<200;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //will return true if Z* mass is smaller than 200GeV
  if(gen_selection["DY_noMass_gt_100"]){
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t i=0; i<_Gen->size(); i++){
      if(abs(_Gen->pdg_id[i])==11 or abs(_Gen->pdg_id[i])==13 or abs(_Gen->pdg_id[i])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(i);
          return (lep1+lep2).M()<100;
        }else{
          lep1= _Gen->p4(i);
        }
      }
    }
  }
  //cout<<"Something is rotten in the state of Denmark."<<std::endl;
  //cout<<"could not find gen selection particle"<<std::endl;
  return true;
}


double Analyzer::getTauDataMCScaleFactor(int updown){
  double sf=1.;
  //for(size_t i=0; i<_Tau->size();i++){
  for(auto i : *active_part->at(CUTS::eRTau1)){
    if(matchTauToGen(_Tau->p4(i),0.4)!=TLorentzVector()){

      if(updown==-1) sf*=  _Tau->pstats["Smear"].dmap.at("TauSF") * (1.-(0.35*_Tau->pt(i)/1000.0));
      else if(updown==0) sf*=  _Tau->pstats["Smear"].dmap.at("TauSF");
      else if(updown==1) sf*=  _Tau->pstats["Smear"].dmap.at("TauSF") * (1.+(0.05*_Tau->pt(i)/1000.0));
    }
  }
  return sf;
}

///Calculates met from values from each file plus smearing and treating muons as neutrinos
void Analyzer::updateMet(int syst) {
  _MET->update(distats["Run"], *_Jet,  syst);

  /////MET CUTS

  if(!passCutRange(_MET->pt(), distats["Run"].pmap.at("MetCut"))) return;
  if(distats["Run"].bfind("DiscrByHT") && _MET->HT() < distats["Run"].dmap.at("HtCut")) return;

  if(syst==0){
    active_part->at(CUTS::eMET)->push_back(1);
  }else{
    syst_parts.at(syst).at(CUTS::eMET)->push_back(1);
  }
}

bool Analyzer::passMetFilters(std::string year, int ievent){

  // Set the branches accordingly:
  // good vertices filter
  SetBranch("Flag_goodVertices", primaryvertexfilter);
  // beam halo filter
  SetBranch("Flag_globalSuperTightHalo2016Filter", beamhalofilter);
  // HBHE noise filter
  SetBranch("Flag_HBHENoiseFilter", hbhenoisefilter);
  // ECAL trigger primitives filter
  SetBranch("Flag_EcalDeadCellTriggerPrimitiveFilter", ecaltpfilter);
  // Bad PF muon filter
  SetBranch("Flag_BadPFMuonFilter", badpfmuonfilter);
  // Bad charged hadron filter - not recommended.
  // SetBranch("Flag_BadChargedCandidateFilter", badchargedhadronfilter);
  if(year.compare("2016") == 0){
    // in 2016, this filter is not recommended... therefore we set it always to true.
    ecalbadcalibrationfilter = true;
  }
  else{
    // ECAL bad calibration filter (2017 + 2018).
    // SetBranch("Flag_ecalBadCalibFilter", ecalbadcalibrationfilter);
    SetBranch("Flag_ecalBadCalibFilterV2", ecalbadcalibrationfilter);
  }

  // Call get entry so all the branches assigned here are filled with the proper values for each event.
  BOOM->GetEntry(ievent);

  // Check if the current event passed all the flags
  allmetfilters = primaryvertexfilter && beamhalofilter && hbhenoisefilter && ecaltpfilter && badpfmuonfilter && ecalbadcalibrationfilter;

  return allmetfilters;

}

///////////////////////////////////////////////
////////removed for teh time being/////////////
///////////////////////////////////////////////

// void Analyzer::treatMuons_Met(std::string syst) {

//   //syst not implemented for muon as tau or neutrino yet
//   if( syst!="orig" or !( distats["Run"].bfind("TreatMuonsAsNeutrinos") || distats["Run"].bfind("TreatMuonsAsTaus")) ){
//     return;
//   }

//   //  Neutrino update before calculation
//   _MET->addP4Syst(_MET->p4(),"muMET");
//   _MET->systdeltaMEx["muMET"]=0;
//   _MET->systdeltaMEy["muMET"]=0;

//   if(distats["Run"].bfind("TreatMuonsAsNeutrinos")) {
//     for(auto it : *active_part->at(CUTS::eRMuon1)) {
//       if(find(active_part->at(CUTS::eRMuon2)->begin(), active_part->at(CUTS::eRMuon2)->end(), it) != active_part->at(CUTS::eRMuon2)->end() ) continue;
//       _MET->systdeltaMEx["muMET"] += _Muon->p4(it).Px();
//       _MET->systdeltaMEy["muMET"] += _Muon->p4(it).Py();
//     }
//     for(auto it : *active_part->at(CUTS::eRMuon2)) {
//       _MET->systdeltaMEx["muMET"] += _Muon->p4(it).Px();
//       _MET->systdeltaMEy["muMET"] += _Muon->p4(it).Py();
//     }
//   }
//   // else if(distats["Run"].bmap.at("TreatMuonsAsTaus")) {

//   //   if(active_part->at(CUTS::eRMuon1)->size() == 1) {

//   //     int muon = (int)active_part->at(CUTS::eRMuon1)->at(0);

//   //     double rand1 = 1;//Tau_HFrac->GetRandom();
//   //     double rand2 = 0;//Tau_Resol->GetRandom();

//   //     double ETau_Pt = _Muon->p4(muon).Pt()*rand1*(rand2+1.0);
//   //     double ETau_Eta = _Muon->p4(muon).Eta();
//   //     double ETau_Phi=normPhi(_Muon->p4(muon).Phi());//+DeltaNu_Phi->GetRandom());
//   //     double ETau_Energy = 0.;


//   //     // double theta = 2.0*TMath::ATan2(1.0,TMath::Exp(_Muon->p4(muon).Eta()));
//   //     // double sin_theta = TMath::Sin(theta);
//   //     // double P_tau = ETau_Pt/sin_theta;

//   //     // //ETau_Energy = sqrt(pow(P_tau, 2) + pow(1.77699, 2));
//   //     // ETau_Energy = sqrt( pow(1.77699, 2) + pow(ETau_Pt, 2) + pow(_Muon->p4(muon).Pz(), 2));

//   //     /*if(ETau_Pt <= 15.0){
//   //       while(ETau_Pt<=15.0){
//   //       rand1 = Tau_HFrac->GetRandom();
//   //       rand2 = Tau_Resol->GetRandom();
//   //       ETau_Pt = _Muon->p4(muon).Pt()*rand1*(rand2+1.0);
//   //       ENu_Pt = _Muon->p4(muon).Pt()-ETau_Pt;
//   //       }
//   //     }
//   //     */

//   //     TLorentzVector Emu_Tau;
//   //     Emu_Tau.SetPtEtaPhiE(ETau_Pt, ETau_Eta, ETau_Phi, ETau_Energy);
//   //     _Muon->cur_P->clear();

//   //     if (ETau_Pt >= _Muon->pstats["Muon1"].pmap.at("PtCut").first ){
//   //       _Muon->cur_P->push_back(Emu_Tau);
//   //       _MET->systdeltaMEy["muMET"] += (_Muon->p4(muon).Px()-Emu_Tau.Px());
//   //       _MET->systdeltaMEy["muMET"] += (_Muon->p4(muon).Py()-Emu_Tau.Py());

//   //     }
//   //   }
//   //}
//   // recalculate MET
//   //  _MET->update("muMET");

//   /////MET CUTS
//   active_part->at(CUTS::eMET)->clear();

//   if(passCutRange(_MET->pt(), distats["Run"].pmap.at("MetCut"))) {
//     active_part->at(CUTS::eMET)->push_back(1);
//   }
// }

/////Check if a given branch is not found in the file

void Analyzer::branchException(std::string branch){
  if(BOOM->FindBranch(branch.c_str()) == 0 ){
     throw "Branch not found in the current sample. Check the config files associated to this branch.";
  }
}

void Analyzer::getTriggerBranchesList(std::string trigger, bool usewildcard){

 if(usewildcard){
  for(int i=0; i < BOOM->GetListOfBranches()->GetSize(); i++){

    std::string branch_name = BOOM->GetListOfBranches()->At(i)->GetName();
    if(branch_name.find(trigger) == std::string::npos) continue;
    //std::cout << "The branch: " << branch_name << " is a selected trigger branch." << std::endl;
    triggerBranchesList.push_back(branch_name);
  }
 }
 else{
    for(int i=0; i < BOOM->GetListOfBranches()->GetSize(); i++){
       std::string branch_name = BOOM->GetListOfBranches()->At(i)->GetName();
       // Look for branches that match exactly the trigger name
       if(branch_name.compare(trigger.c_str()) != 0) continue;
       // std::cout << "The branch: " << branch_name << " is a selected trigger branch." << std::endl;
       triggerBranchesList.push_back(branch_name);
 }
}

  if(triggerBranchesList.size() == 0) throw "no branches matching this name were found. Check the trigger name requirement.";
}

/////sets up other values needed for analysis that aren't particle specific
void Analyzer::setupGeneral(std::string year) {

  genMaper = {
    {5, new GenFill(2, CUTS::eGJet)},     {6,  new GenFill(2, CUTS::eGTop)},
    {11, new GenFill(1, CUTS::eGElec)},   {13, new GenFill(1, CUTS::eGMuon)},
    {15, new GenFill(2, CUTS::eGTau)},    {23, new GenFill(62, CUTS::eGZ)},
    {24, new GenFill(62, CUTS::eGW)},      {25, new GenFill(2, CUTS::eGHiggs)}, 
    {5, new GenFill(52, CUTS::eGBJet)}
  };
  
  isData=true;
  if(BOOM->FindBranch("Pileup_nTrueInt")!=0){
    isData=false;
  }

  if(BOOM->FindBranch("Electron_mvaFall17Iso")!=0){
    std::cout<<"This file needs the new version of the analyzer"<<std::endl;
  } 
  
  read_info(filespace + "ElectronTau_info.in");
  read_info(filespace + "MuonTau_info.in");
  read_info(filespace + "MuonElectron_info.in");
  read_info(filespace + "DiParticle_info.in");
  read_info(filespace + "ElectronJet_info.in");
  read_info(filespace + "VBFCuts_info.in");
  read_info(filespace + "Run_info.in");
  read_info(filespace + "Systematics_info.in");
	
  // Call readinJSON and save this in the jsonlinedict we declared in Analyzer.h
  jsonlinedict = readinJSON(year);


  for(std::string trigger : inputTriggerNames){

    try{
      getTriggerBranchesList(trigger, distats["Run"].bfind("UseTriggerWildcard"));
    }
    catch(const char* msg){
      std::cout << "ERROR! Trigger " << trigger << ": " << msg << std::endl;
      continue;
    }
  
  }

  // Check that there are no elements on the trigger list that refer to the same trigger 
  // to speed up the process.
  removeDuplicates(triggerBranchesList);

  std::cout << " ---------------------------------------------------------------------- " << std::endl;
  std::cout << "Full list of triggers to be probed: " << std::endl;
  for(std::string name : triggerBranchesList){
    std::cout << name << std::endl;
  }
  std::cout << " ---------------------------------------------------------------------- " << std::endl;
}


void Analyzer::initializeMCSelection(std::vector<std::string> infiles) {
    // check if we need to make gen level cuts to cross clean the samples:

  isVSample = false;
  if(infiles[0].find("DY") != std::string::npos){
    isVSample = true;
    if(infiles[0].find("DYJetsToLL_M-50_HT-") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
      //gen_selection["DY_noMass_gt_200"]=true;
    //get the DY1Jet DY2Jet ...
    }else if(infiles[0].find("JetsToLL_TuneCUETP8M1_13TeV") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
      gen_selection["DY_noMass_gt_200"]=false;
    }

    if(infiles[0].find("DYJetsToLL_M-50_TuneCUETP8M1_13TeV") != std::string::npos){
      gen_selection["DY_noMass_gt_100"]=true;
    }else{
      //set it to false!!
      gen_selection["DY_noMass_gt_100"]=false;
      gen_selection["DY_noMass_gt_200"]=false;
    }
  }else{
    //set it to false!!
    gen_selection["DY_noMass_gt_200"]=false;
    gen_selection["DY_noMass_gt_100"]=false;
  }

  if(infiles[0].find("WJets") != std::string::npos){
    isVSample = true;
  }
}


///parsing method that gets info on diparts and basic run info
//put in std::map called "distats"
void Analyzer::read_info(std::string filename) {
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  std::ifstream info_file(filename);
  boost::char_separator<char> sep(", \t");

  if(!info_file) {
    std::cout << "could not open file " << filename <<std::endl;
    exit(1);
  }


  std::string group, line;
  while(getline(info_file, line)) {
    tokenizer tokens(line, sep);
    std::vector<std::string> stemp;

    for(tokenizer::iterator iter = tokens.begin();iter != tokens.end(); iter++) {
      if( ((*iter)[0] == '/' && (*iter)[0] == '/') || ((*iter)[0] == '#') ) break;
      stemp.push_back(*iter);
    }
    if(stemp.size() == 0) continue;
    else if(stemp.size() == 1) {
      group = stemp[0];
      continue;
    } else if(group == "") {
      std::cout << "error in " << filename << "; no groups specified for data" << std::endl;
      exit(1);
    } else if(stemp.size() == 2) {
      char* p;
      strtod(stemp[1].c_str(), &p);
      if(group.compare("Control_Region") !=0 ){
        if(stemp[1] == "1" || stemp[1] == "true"){
          distats[group].bset.push_back(stemp[0]);
          if(stemp[1] == "1" ){
            distats[group].dmap[stemp[0]]=std::stod(stemp[1]);
          }
        }
        else if(*p) distats[group].smap[stemp[0]] = stemp[1];
        else  distats[group].dmap[stemp[0]]=std::stod(stemp[1]);
      }else{
        if(*p) distats[group].smap[stemp[0]] = stemp[1];
        else  distats[group].dmap[stemp[0]]=std::stod(stemp[1]);
      }
      if(stemp.at(0).find("Trigger") != std::string::npos) {
        for(auto trigger : stemp){
          if(trigger.find("Trigger")== std::string::npos and "="!=trigger ){
            inputTriggerNames.push_back(trigger);
          }
        }
        continue;
      }
    } else if(stemp.size() == 3 and stemp.at(0).find("Trigger") == std::string::npos){
      distats[group].pmap[stemp[0]] = std::make_pair(std::stod(stemp[1]), std::stod(stemp[2]));
    } else{
      if(stemp.at(0).find("Trigger") != std::string::npos) {
        for(auto trigger : stemp){
          if(trigger.find("Trigger")== std::string::npos and "="!=trigger ){
            inputTriggerNames.push_back(trigger);
          }
        }
        continue;
      }
    }
  }
  info_file.close();
  //for( std::pair<const std::basic_string<char>, PartStats> group : distats){
    //for( std::pair<const std::basic_string<char>, double> i : group.second.dmap){
      //if(group.first.find("Tau")!=string::npos){
        //cout<<group.first<<" "<< i.first<< "  "<<i.second<<std::endl;
      //}
    //}
  //}
}


// This code works pretty much (at least in my tests), but dagnabit, its ugly.  They all can't be winners, at least now...
void Analyzer::setCutNeeds() {


  for(auto it: *histo.get_groups()) {
    if(fillInfo[it]->type == FILLER::None) continue;
    neededCuts.loadCuts(fillInfo[it]->ePos);
  }
  for(auto it : *histo.get_cutorder()) {
    try{
      neededCuts.loadCuts(cut_num.at(it));
    }catch(...){
      std::cout<<"The following cut is strange: "<<it<<std::endl;
      exit(2);
    }
  }
  for(auto it: testVec) {
    neededCuts.loadCuts(it->info->ePos);
  }

  if(!isData and distats["Run"].bfind("ApplyZBoostSF") and isVSample){
    neededCuts.loadCuts(CUTS::eGen);
    neededCuts.loadCuts(CUTS::eGZ);
    neededCuts.loadCuts(CUTS::eGW);
  }

  neededCuts.loadCuts(CUTS::eGHadTau);
  neededCuts.loadCuts(CUTS::eGBJet); //01.16.19

  neededCuts.loadCuts(_Jet->findExtraCuts());
  if(doSystematics) {
    neededCuts.loadCuts(CUTS::eGen);
  }

  for(auto it: jetCuts) {
    if(!neededCuts.isPresent(it)) continue;
    neededCuts.loadCuts(_Jet->overlapCuts(it));
  }

  if(neededCuts.isPresent(CUTS::eRWjet)) {
    neededCuts.loadCuts(_FatJet->findExtraCuts());
    neededCuts.loadCuts(_FatJet->overlapCuts(CUTS::eRWjet));
  } else {
    std::cout<<"WJets not needed. They will be deactivated!"<<std::endl;
    _FatJet->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRTau1) || neededCuts.isPresent(CUTS::eRTau2) ) {
    neededCuts.loadCuts(_Tau->findExtraCuts());
  } else {
    std::cout<<"Taus not needed. They will be deactivated!"<<std::endl;
    _Tau->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRElec1) || neededCuts.isPresent(CUTS::eRElec2) ) {
    neededCuts.loadCuts(_Electron->findExtraCuts());
  } else {
    std::cout<<"Electrons not needed. They will be deactivated!"<<std::endl;
    _Electron->unBranch();
  }

  if( neededCuts.isPresent(CUTS::eRMuon1) || neededCuts.isPresent(CUTS::eRMuon2) ) {
    neededCuts.loadCuts(_Muon->findExtraCuts());
  } else {
    std::cout<<"Muons not needed. They will be deactivated!"<<std::endl;
    _Muon->unBranch();
  }

  if( !neededCuts.isPresent(CUTS::eGen) and !isData) {
    std::cout<<"Gen not needed. They will be deactivated!"<<std::endl;
    _Gen->unBranch();

  }

  std::cout << "Cuts being filled: " << std::endl;
  for(auto cut : neededCuts.getCuts()) {
    std::cout << enumNames.at(static_cast<CUTS>(cut)) << "   ";
  }
  std::cout << std::endl;
}


///Smears lepton only if specified and not a data file.  Otherwise, just filles up lorentz std::vectors
//of the data into the std::vector container smearP with is in each lepton object.
void Analyzer::smearLepton(Lepton& lep, CUTS eGenPos, const PartStats& stats, const PartStats& syst_stats, int syst) {
  if( isData) {
    lep.setOrigReco();
    return;
  }

  std::string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) return;

  if(systname=="orig" && !stats.bfind("SmearTheParticle")){
    lep.setOrigReco();
  } else {
    systematics.loadScaleRes(stats, syst_stats, systname);
    for(size_t i = 0; i < lep.size(); i++) {
      TLorentzVector lepReco = lep.RecoP4(i);
      TLorentzVector genVec =  matchLeptonToGen(lepReco, lep.pstats["Smear"],eGenPos);
      systematics.shiftLepton(lep, lepReco, genVec, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst], syst);
    }
  }
}

///Same as smearlepton, just jet specific
void Analyzer::smearJet(Particle& jet, const CUTS eGenPos, const PartStats& stats, int syst) {
  //at the moment
  if(isData || jet.type != PType::Jet ){
    //|| !stats.bfind("SmearTheJet")
    jet.setOrigReco();
    return;
  }
  if(!jet.needSyst(syst)){
    return;
  }
  //add energy scale uncertainty


  std::string systname = syst_names.at(syst);

  for(size_t i=0; i< jet.size(); i++) {
    TLorentzVector jetReco = jet.RecoP4(i);
    if(JetMatchesLepton(*_Muon, jetReco, stats.dmap.at("MuonMatchingDeltaR"), CUTS::eGMuon) ||
       JetMatchesLepton(*_Tau, jetReco, stats.dmap.at("TauMatchingDeltaR"), CUTS::eGTau) ||
       JetMatchesLepton(*_Electron, jetReco,stats.dmap.at("ElectronMatchingDeltaR"), CUTS::eGElec)){
      jet.addP4Syst(jetReco,syst);
      continue;
    }

    double sf=1.;
    //only apply corrections for jets not for FatJets

    TLorentzVector genJet=matchJetToGen(jetReco, jet.pstats["Smear"],eGenPos);
    if(systname=="orig" && stats.bfind("SmearTheJet")){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, 0);
    }else if(systname=="Jet_Res_Up"){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, 1);
    }else if(systname=="Jet_Res_Down"){
      sf=jetScaleRes.GetRes(jetReco,genJet, rho, -1);
    }else if(systname=="Jet_Scale_Up"){
      sf = jetScaleRes.GetScale(jetReco, false, +1.);
    }else if(systname=="Jet_Scale_Down"){
      sf = jetScaleRes.GetScale(jetReco, false, -1) ;
    }
    //cout<<systname<<"  "<<sf<<"  "<<jetReco.Pt()<<"  "<<genJet.Pt()<<std::endl;
    systematics.shiftParticle(jet, jetReco, sf, _MET->systdeltaMEx[syst], _MET->systdeltaMEy[syst], syst);
  }
}


/////checks if jet is close to a lepton and the lepton is a gen particle, then the jet is a lepton object, so
//this jet isn't smeared
bool Analyzer::JetMatchesLepton(const Lepton& lepton, const TLorentzVector& jetV, double partDeltaR, CUTS eGenPos) {
  for(size_t j = 0; j < lepton.size(); j++) {
    if(jetV.DeltaR(lepton.RecoP4(j)) < partDeltaR && matchLeptonToGen(lepton.RecoP4(j), lepton.pstats.at("Smear"), eGenPos) != TLorentzVector(0,0,0,0)) return true;
  }
  return false;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchLeptonToGen(const TLorentzVector& recoLepton4Vector, const PartStats& stats, CUTS ePos) {
  if(ePos == CUTS::eGTau) {
    return matchTauToGen(recoLepton4Vector, stats.dmap.at("GenMatchingDeltaR"));
  }
  if(ePos == CUTS::eGHadTau){
    return matchHadTauToGen(recoLepton4Vector, stats.dmap.at("GenMatchingDeltaR"));
  }
  for(auto it : *active_part->at(ePos)) {
    if(recoLepton4Vector.DeltaR(_Gen->p4(it)) <= stats.dmap.at("GenMatchingDeltaR")) {
      if(stats.bfind("UseMotherID") && abs(_Gen->pdg_id[_Gen->genPartIdxMother[it]]) != stats.dmap.at("MotherID")) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}

///Tau specific matching function.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchHadTauToGen(const TLorentzVector& recoTau4Vector, double recogenDeltaR) {
  //TLorentzVector genTau4Vector(0,0,0,0);
  //int i = 0; // This will give us the position of the element at CUTS::eGHadTau we are looking at
  for(vec_iter genhadtau_it = active_part->at(CUTS::eGHadTau)->begin(); genhadtau_it != active_part->at(CUTS::eGHadTau)->end(); genhadtau_it++){ // (genhadtau_it) is the index of the gen-level hadronic tau in the gen-hadtau vector.
    //genTau4Vector = _GenHadTau->p4(*genhadtau_it);
    // Compare the separation between the reco and gen hadronic tau candidates. If it's greather than the requirement, continue with the next gen-tau_h candidate.
    if(recoTau4Vector.DeltaR(_GenHadTau->p4(*genhadtau_it)) > recogenDeltaR) continue;
    // If the requirement DeltaR <= recogenDeltaR, the code will get to this point. Then we want to fill a vector that only stores the matched gen taus
    active_part->at(CUTS::eGMatchedHadTau)->push_back(*genhadtau_it); 
    // Remove this gen-tau from CUTS::eGHadTau to avoid getting matched again with any other tau present in the event.
    // active_part->at(CUTS::eGHadTau)->erase(genhadtau_it);
    // And we also return the gen-tau p4 vector, which we are interested in.
    //return genTau4Vector;
    return _GenHadTau->p4(*genhadtau_it);
  }
  //return genTau4Vector;
  return TLorentzVector(0,0,0,0);
}


///Tau specific matching function.  Works by seeing if a tau doesn't decay into a muon/electron and has
//a matching tau neutrino showing that the tau decayed and decayed hadronically
TLorentzVector Analyzer::matchTauToGen(const TLorentzVector& lvec, double lDeltaR) {
  TLorentzVector genVec(0,0,0,0);
  int i = 0;
  for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it !=active_part->at(CUTS::eGTau)->end();it++, i++) {
    int nu = active_part->at(CUTS::eGNuTau)->at(i);
    if(nu == -1) continue;

    genVec = _Gen->p4(*it) - _Gen->p4(nu);
    if(lvec.DeltaR(genVec) <= lDeltaR) {
      return genVec;
    }
  }
  return genVec;
}


////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
TLorentzVector Analyzer::matchJetToGen(const TLorentzVector& lvec, const PartStats& stats, CUTS ePos) {
  //for the future store gen jets
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(_Gen->p4(it)) <= stats.dmap.at("GenMatchingDeltaR")) {
      //nothing more than b quark or gluon
      if( !( (abs(_Gen->pdg_id[it])<5) || (abs(_Gen->pdg_id[it])==9) ||  (abs(_Gen->pdg_id[it])==21) ) ) continue;
      return _Gen->p4(it);
    }
  }
  return TLorentzVector(0,0,0,0);
}



////checks if reco object matchs a gen object.  If so, then reco object is for sure a correctly identified particle
int Analyzer::matchToGenPdg(const TLorentzVector& lvec, double minDR) {
  double _minDR=minDR;
  int found=-1;
  for(size_t i=0; i< _Gen->size(); i++) {

    if(lvec.DeltaR(_Gen->p4(i)) <=_minDR) {
      //only hard interaction
      if( _Gen->status[i]<10){
        found=i;
        _minDR=lvec.DeltaR(_Gen->p4(i));
      }
    }
  }
  if (found>=0){
    return _Gen->pdg_id[found];
  }
  return 0;
}


////Calculates the number of gen particles.  Based on id number and status of each particle
// The file that corresponds to PartStats& stats is Gen_info.in
void Analyzer::getGoodGen(const PartStats& stats) {

  if(! neededCuts.isPresent(CUTS::eGen)) return;

  int particle_id = 0;

  for(size_t j = 0; j < _Gen->size(); j++) {

    particle_id = abs(_Gen->pdg_id[j]);

    if( (particle_id < 5 || particle_id == 9 || particle_id == 21) && genMaper.find(particle_id) != genMaper.end() && _Gen->status[j] == genMaper.at(5)->status){
      active_part->at(genMaper.at(5)->ePos)->push_back(j);
    }
    else if(genMaper.find(particle_id) != genMaper.end() && _Gen->status[j] == genMaper.at(particle_id)->status) {
      // Cuts on gen-taus (before decaying)
      if(stats.bfind("DiscrTauByPtAndEta")){
        if(particle_id == 15 && (_Gen->pt(j) < stats.pmap.at("TauPtCut").first || _Gen->pt(j) > stats.pmap.at("TauPtCut").second || abs(_Gen->eta(j)) > stats.dmap.at("TauEtaCut"))) continue;
      }
      // Cuts on light lepton mother IDs
      else if(stats.bfind("DiscrLightLepByMotherID")){
       if( (particle_id == 11 || particle_id == 13) && (abs(_Gen->pdg_id[_Gen->genPartIdxMother[j]]) != stats.pmap.at("LightLepMotherIDs").first || abs(_Gen->pdg_id[_Gen->genPartIdxMother[j]]) != stats.pmap.at("LightLepMotherIDs").second)) continue;
      }

      active_part->at(genMaper.at(particle_id)->ePos)->push_back(j);    
    }

  }
}

// --- Function that applies selections to hadronic taus at gen-level (stored in the GenVisTau list) --- //
void Analyzer::getGoodGenHadronicTaus(const PartStats& stats){

  // Loop over all gen-level hadronic taus stored in the corresponding list to apply certain selections
  for(size_t i=0; i < _GenHadTau->size(); i++){

    if(_GenHadTau->pt(i) < stats.pmap.at("HadTauPtCut").first || _GenHadTau->pt(i) > stats.pmap.at("HadTauPtCut").second || abs(_GenHadTau->eta(i)) > stats.dmap.at("HadTauEtaCut")) continue;
    else if( stats.bfind("DiscrByTauDecayMode") && (_GenHadTau->decayMode[i] < stats.pmap.at("TauDecayModes").first || _GenHadTau->decayMode[i] > stats.pmap.at("TauDecayModes").second)) continue;

    active_part->at(CUTS::eGHadTau)->push_back(i);
  } 
}

// --- Function that gets the Lorentz vector of taus that decayed hadronically using the tagging method in Analyzer::getGenHadronicTauNeutrinos() --- //
TLorentzVector Analyzer::getGenVisibleTau4Vector(int gentau_idx, int gentaunu_idx){
  TLorentzVector visTau4Vector(0,0,0,0);

  if(gentaunu_idx != -1 && gentau_idx != -1){
    visTau4Vector = _Gen->p4(gentau_idx) - _Gen->p4(gentaunu_idx);
  } 
  return visTau4Vector;
}

// --- Function that looks for tau neutrinos coming from hadronic tau decays in the full gen-level particle list --- //
void Analyzer::getGoodGenHadronicTauNeutrinos(const PartStats& stats){
  int genTauNeutrino_idx = -1, genHadTauNeutrino_idx = -1;   // integers for the particle indices in the gen-particle vector.
  TLorentzVector genVisHadTau4Vector(0,0,0,0);

  // Loop over the gen-level taus that satisfied the conditions imposed in getGoodGen, which are stored at CUTS::eGTau.
  for(auto gentau_it : *active_part->at(CUTS::eGTau)){ // (gentau_it) is the index of the gen-level tau in the gen-particles vector.
    // For each iteration, reset the tau neutrino index and the leptonic flag.
    genTauNeutrino_idx = -1, genHadTauNeutrino_idx = -1;
    // isTauLeptonicDecay = false;

    // Make sure that the status code of the current gen-tau is 2, meaning, it's an unstable particle about to decay
    if(_Gen->status[(gentau_it)] != 2) continue;

    // Loop over all gen-level particles to find those that come from the decay of the current gen-tau.
    for(size_t genpart_idx = 0; genpart_idx < _Gen->size(); genpart_idx++) {

      // Check that the mother particle index matches the index of the current gen-tau and that the particle is any neutrino:
      if( abs(_Gen->genPartIdxMother[genpart_idx]) != (gentau_it)) continue; // || abs(_Gen->pdg_id[genpart_idx]) != 12 || abs(_Gen->pdg_id[genpart_idx]) != 14 || abs(_Gen->pdg_id[genpart_idx]) != 16) continue;

      // Look specifically for electron, muon and tau neutrinos. They are enough to tell if a decay was hadronic or leptonic.
      // If in the particle decays there is an electron or muon neutrino, the decay is leptonic and we break this for loop and continue with the next gen-tau in CUTS::eGTau
      if( ( abs(_Gen->pdg_id[genpart_idx]) == 12 || ( abs(_Gen->pdg_id[genpart_idx]) == 14 ) ) ){
        genTauNeutrino_idx = -1;
        break;
      }
      // Since we reject all taus that have e/mu neutrinos, we will only get hadronic tau decays and we want to store the index of the neutrino for these events.
      // As soon as it's found, then break the loop to not go over the entire list of gen particles for this events.
      else if( abs(_Gen->pdg_id[genpart_idx]) == 16){ // genTauNeutrino_idx = genpart_idx;
        genTauNeutrino_idx = genpart_idx;
      }
    }

    // Check if the current tau decayed hadronically or leptonically and assign the index accordingly
    genHadTauNeutrino_idx = genTauNeutrino_idx;

    // Apply kinematic cuts on the visible hadronic tau vector
    genVisHadTau4Vector = getGenVisibleTau4Vector(gentau_it, genHadTauNeutrino_idx);
    if(genVisHadTau4Vector.Pt() < stats.pmap.at("HadTauPtCut").first || genVisHadTau4Vector.Pt() > stats.pmap.at("HadTauPtCut").second || abs(genVisHadTau4Vector.Eta()) > stats.dmap.at("HadTauEtaCut")){
      genHadTauNeutrino_idx = -1;
    }

    active_part->at(CUTS::eGNuTau)->push_back(genHadTauNeutrino_idx);
  }

}

void Analyzer::getGoodGenBJet() { //01.16.19
  for (size_t j=0; j < _Gen->size(); j++){
    int id = abs(_Gen->pdg_id[j]);
    int motherid = abs(_Gen->pdg_id[_Gen->genPartIdxMother[j]]);
    int motherind = abs(_Gen->genPartIdxMother[j]);
    if(id == 5 && motherid == 6)
      {for(size_t k=0; k < _Gen->size(); k++){
    if(abs(_Gen->pdg_id[k]) == 24 && _Gen->genPartIdxMother[k] == motherind)
      {active_part->at(CUTS::eGBJet)->push_back(j);
      }
  }
      }
  }
}

double Analyzer::getTopBoostWeight(){ //01.15.19
  double topPt; //initialize a value to hold the p_T of the top
  double topBarPt; //initialize a value to hold the p_T of the tbar
  double SFtop = 1; //initialize a value to hold the top SF
  double SFtopBar = 1; //initialize a value to hold the tbar SF
  double SFttbar = 1; //holds SF for both
  
  for(size_t j = 0; j < _Gen->size(); j++) {
    int id = _Gen->pdg_id[j]; //grab the particle ID
    int daught = _Gen->numDaught[j]; //grab the number of daughter particles
    if(id == 6 && daught == 2){ //check that a top has two daughters 
      topPt = _Gen->pt(j); //grab its p_T
      SFtop = exp(0.0615 - (0.0005 * topPt));} //calculate the top SF
    if(id == -6 && daught == 2){ //check that a tbar has two daughters
      topBarPt = _Gen->pt(j); //grab its p_T
      SFtopBar = exp(0.0615 - (0.0005 * topBarPt));} //calculate the tbar SF
    SFttbar = sqrt(SFtop * SFtopBar);} //calculate the total SF
  
  return SFttbar;
}

///Function used to find the number of reco leptons that pass the various cuts.
///Divided into if blocks for the different lepton requirements.
void Analyzer::getGoodRecoLeptons(const Lepton& lep, const CUTS ePos, const CUTS eGenPos, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!lep.needSyst(syst)) {
    active_part->at(ePos) = goodParts[ePos];
    return;
  }

  int i = 0;
  for(auto lvec: lep) {
    bool passCuts = true;
    if (fabs(lvec.Eta()) > stats.dmap.at("EtaCut")) passCuts = passCuts && false;
    else if (lvec.Pt() < stats.pmap.at("PtCut").first || lvec.Pt() > stats.pmap.at("PtCut").second) passCuts = passCuts && false;

    if((lep.pstats.at("Smear").bfind("MatchToGen")) && (!isData)) {   /////check
      if(matchLeptonToGen(lvec, lep.pstats.at("Smear") ,eGenPos) == TLorentzVector(0,0,0,0)) continue;
    }

    for( auto cut: stats.bset) {
      if(!passCuts) break;
      else if(cut == "DoDiscrByIsolation") {
        double firstIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").first : ival(ePos) - ival(CUTS::eRTau1) + 1;
	//if (lep.type == PType::Electron){std::cout << "firstIso: " << firstIso << std::endl;}
        double secondIso = (stats.pmap.find("IsoSumPtCutValue") != stats.pmap.end()) ? stats.pmap.at("IsoSumPtCutValue").second : stats.bfind("FlipIsolationRequirement");
	//if (lep.type == PType::Electron){std::cout << "secondIso: " << secondIso << std::endl;}
	passCuts = passCuts && lep.get_Iso(i, firstIso, secondIso);
	//if (lep.type == PType::Electron){std::cout << "pass cuts: " << passCuts << std::endl;}
      }
      else if(cut == "DiscrIfIsZdecay" && lep.type != PType::Tau ) passCuts = passCuts && isZdecay(lvec, lep);
      else if(cut == "DiscrByMetDphi") passCuts = passCuts && passCutRange(absnormPhi(lvec.Phi() - _MET->phi()), stats.pmap.at("MetDphiCut"));
      else if(cut == "DiscrByMetMt") passCuts = passCuts && passCutRange(calculateLeptonMetMt(lvec), stats.pmap.at("MetMtCut"));
      /////muon cuts
      else if(lep.type == PType::Muon){
        if(cut == "DoDiscrByTightID") passCuts = passCuts && _Muon->tight[i];
        else if(cut == "DoDiscrBySoftID") passCuts = passCuts && _Muon->soft[i];
      }
      ////electron cuts
      else if(lep.type == PType::Electron){
	//std::cout << "got an electron to check..." << std::endl;
        if(cut == "DoDiscrByHLTID"){
          std::bitset<8> idvariable(_Electron->cutBased_HLTPreSel[i]);
          if(ival(ePos) - ival(CUTS::eRElec2)){ //test if it is electron1
            passCuts = passCuts && (_Electron->cbHLTIDele1&idvariable).count();
	    //std::cout << "it's an elec1..." << endl;
          }else{
            passCuts = passCuts && (_Electron->cbHLTIDele2&idvariable).count();
          }
        }
        if(cut == "DoDiscrByCBID"){
	  //std::cout << "cutBased[i]: " << (_Electron->cutBased[i]) << std::endl;
          std::bitset<8> idvariable(_Electron->cutBased[i]);
	  //std::cout << "ival(ePos): " << ival(ePos) << std::endl;
	  //std::cout << "ival(CUTS::eRElec2): " << ival(CUTS::eRElec2) << std::endl;
          if(ival(ePos) - ival(CUTS::eRElec2)){ //test if it is electron1
	    //std::cout << "cbIDele1: " << (_Electron->cbIDele1) << std::endl;
	    //std::cout << "idvariable: " << (idvariable) << std::endl;
	    //std::cout << "the AND: " << (_Electron->cbIDele1&idvariable) << std::endl;
            passCuts = passCuts && (_Electron->cbIDele1&idvariable).count();
	    //std::cout << "it's an elec1..." << std::endl;
	    //std::cout << "passCuts value2: " << passCuts << std::endl;
          }else{
            passCuts = passCuts && (_Electron->cbIDele2&idvariable).count();
          }
        }
        else if(cut == "DoDiscrByHEEPID")
         passCuts = passCuts && _Electron->isPassHEEPId[i];
      }
      else if(lep.type == PType::Tau){
        if(cut == "DoDiscrByCrackCut") passCuts = passCuts && !isInTheCracks(lvec.Eta());
        /////tau cuts
        else if(cut == "DoDzCut") passCuts = passCuts && (abs(_Tau->dz[i]) <= stats.dmap.at("DzCutThreshold"));
        else if(cut == "DoDiscrByLeadTrack") passCuts = passCuts && (_Tau->leadTkPtOverTauPt[i]*_Tau->pt(i) >= stats.dmap.at("LeadTrackThreshold"));
             // ----Electron and Muon vetos
        else if (cut == "DoDiscrAgainstElectron") passCuts = passCuts && _Tau->pass_against_Elec(ePos, i);
        else if (cut == "SelectTausThatAreElectrons") passCuts = passCuts && !_Tau->pass_against_Elec(ePos, i);

        else if (cut == "DoDiscrAgainstMuon") passCuts = passCuts && _Tau->pass_against_Muon(ePos, i);
        else if (cut == "SelectTausThatAreMuons") passCuts = passCuts &&  !_Tau->pass_against_Muon(ePos, i);

        else if(cut == "DiscrByProngType") {
          passCuts = passCuts && (stats.smap.at("ProngType").find("hps") == std::string::npos || _Tau->DecayModeNewDMs[i] != 0);
          // passCuts = passCuts && passProng(stats.smap.at("ProngType"), _Tau->decayMode[i]); //original.
          passCuts = passCuts && passProng(stats.smap.at("ProngType"), _Tau->decayModeInt[i]);
        }
        else if(cut == "decayModeFindingNewDMs") passCuts = passCuts && _Tau->DecayModeNewDMs[i] != 0;
        // else if(cut == "decayModeFinding") passCuts = passCuts && _Tau->DecayMode[i] != 0; // original
        else if(cut == "decayModeFinding") passCuts = passCuts && _Tau->DecayModeOldDMs[i] != 0;
              // ----anti-overlap requirements
        else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
        else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      }
      else std::cout << "cut: " << cut << " not listed" << std::endl;
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }

  return;
}

////Jet specific function for finding the number of jets that pass the cuts.
//used to find the nubmer of good jet1, jet2, central jet, 1st and 2nd leading jets and bjet.
void Analyzer::getGoodRecoJets(CUTS ePos, const PartStats& stats, const int syst) {

  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_Jet) {
    if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet){
      break;
    }
    bool passCuts = true;
    double dphi1rjets = normPhi(lvec.Phi() - _MET->phi());
    if( ePos == CUTS::eRCenJet) passCuts = passCuts && (fabs(lvec.Eta()) < 2.5);
    else  passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    for( auto cut: stats.bset) {
      if(!passCuts) break;

    /// BJet specific
      // else if(cut == "ApplyJetBTagging") passCuts = passCuts && (_Jet->bDiscriminator[i] > stats.dmap.at("JetBTaggingCut")); // original	
      else if(cut == "ApplyJetBTaggingCSVv2") passCuts = passCuts && (_Jet->bDiscriminatorCSVv2[i] > stats.dmap.at("JetBTaggingCut")); 
      else if(cut == "ApplyJetBTaggingDeepCSV") passCuts = passCuts && (_Jet->bDiscriminatorDeepCSV[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyJetBTaggingDeepFlav") passCuts = passCuts && (_Jet->bDiscriminatorDeepFlav[i] > stats.dmap.at("JetBTaggingCut"));
      //else if(cut == "MatchBToGen") passCuts = passCuts && (isData ||  abs(_Jet->partonFlavour->at(i)) == 5);
      else if(cut == "ApplyLooseID") passCuts = passCuts && _Jet->passedLooseJetID(i);
      else if(cut == "ApplyTightID") passCuts = passCuts && _Jet->passedTightJetID(i);
      else if(cut == "RemoveOverlapWithBs") passCuts = passCuts && !isOverlapingB(lvec, *_Jet, CUTS::eRBJet, stats.dmap.at("BJMatchingDeltaR"));
    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));
      else if(cut == "DiscrByDphi1") {//cout << "dphi1rjets in function: " << dphi1rjets << endl;
	passCuts = passCuts && passCutRange(fabs(dphi1rjets), stats.pmap.at("Dphi1CutMet"));} //01.17.19
      else if(cut == "UseBtagSF") {
        //double bjet_SF = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, lvec.Eta(), lvec.Pt());
        //passCuts = passCuts && (isData || ((double) rand()/(RAND_MAX)) <  bjet_SF);
      }
    }
    if(_Jet->pstats["BJet"].bfind("RemoveBJetsFromJets") and ePos!=CUTS::eRBJet){
      passCuts = passCuts && find(active_part->at(CUTS::eRBJet)->begin(), active_part->at(CUTS::eRBJet)->end(), i) == active_part->at(CUTS::eRBJet)->end();
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;

  }

  //clean up for first and second jet
  //note the leading jet has to be selected fist!
  if(ePos == CUTS::eR1stJet || ePos == CUTS::eR2ndJet) {

    std::vector<std::pair<double, int> > ptIndexVector;
    for(auto it : *active_part->at(CUTS::eRJet1)) {
      ptIndexVector.push_back(std::make_pair(_Jet->pt(it),it));
    }
    sort(ptIndexVector.begin(),ptIndexVector.end());
    if(ePos == CUTS::eR1stJet && ptIndexVector.size()>0){
      active_part->at(ePos)->push_back(ptIndexVector.back().second);
    }
    else if(ePos == CUTS::eR2ndJet && ptIndexVector.size()>1){
      active_part->at(ePos)->push_back(ptIndexVector.at(ptIndexVector.size()-2).second);
    }
    
  }

}

void Analyzer::getGoodRecoBJets(CUTS ePos, const PartStats& stats, const int syst) { //01.16.19   
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_Jet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;
  for(auto lvec: *_Jet) {
    bool passCuts = true;
    if( ePos == CUTS::eRCenJet) passCuts = passCuts && (fabs(lvec.Eta()) < 2.5);
    else  passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    for( auto cut: stats.bset) {
      if(!passCuts) break;

    /// BJet specific
      //else if(cut == "ApplyJetBTagging") passCuts = passCuts && (_Jet->bDiscriminator[i] > stats.dmap.at("JetBTaggingCut")); //original
      else if(cut == "ApplyJetBTaggingCSVv2") passCuts = passCuts && (_Jet->bDiscriminatorCSVv2[i] > stats.dmap.at("JetBTaggingCut")); 
      else if(cut == "ApplyJetBTaggingDeepCSV") passCuts = passCuts && (_Jet->bDiscriminatorDeepCSV[i] > stats.dmap.at("JetBTaggingCut"));
      else if(cut == "ApplyJetBTaggingDeepFlav") passCuts = passCuts && (_Jet->bDiscriminatorDeepFlav[i] > stats.dmap.at("JetBTaggingCut"));
      //else if(cut == "MatchBToGen") passCuts = passCuts && (isData ||  abs(_Jet->partonFlavour->at(i)) == 5);
      else if(cut == "ApplyLooseID") passCuts = passCuts && _Jet->passedLooseJetID(i);
      else if(cut == "ApplyTightID") passCuts = passCuts && _Jet->passedTightJetID(i);

    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));

      else if(cut == "UseBtagSF") {
        //double bjet_SF = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, lvec.Eta(), lvec.Pt());
        //passCuts = passCuts && (isData || ((double) rand()/(RAND_MAX)) <  bjet_SF);
      }
    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }
}

// The function below sets up the information from the right CSV file in the Pileup folder
// to obtain the functions needed to apply b-tagging SF in an automatic way.
void Analyzer::setupBJetSFInfo(const PartStats& stats, std::string year){

   // Always check that the filenames are up to date!!
   // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
   static std::map<std::string, std::string> csvv2btagsfcsvfiles = {
    {"2016", "CSVv2_Moriond17_B_H.csv"},
    {"2017", "CSVv2_94XSF_WP_V2_B_F.csv"},
    {"2018", "none.csv"} // No CSVv2 available for 2018
   }; 

   static std::map<std::string, std::string> deepcsvbtagsfcsvfiles = {
    {"2016", "DeepCSV_2016LegacySF_WP_V1.csv"},
    {"2017", "DeepCSV_94XSF_WP_V4_B_F.csv"},
    {"2018", "DeepCSV_102XSF_WP_V1.csv"}
   };

  static std::map<std::string, std::string> deepflavbtagsfcsvfiles = {
    {"2016", "DeepJet_2016LegacySF_WP_V1.csv"},
    {"2017", "DeepFlavour_94XSF_WP_V3_B_F.csv"},
    {"2018", "DeepJet_102XSF_WP_V1.csv"}
  };

   // This will be the default.
  std::string btagalgoname = "DeepCSV";
  std::string btagsffilename = deepcsvbtagsfcsvfiles.begin()->second;
  
  // std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
  // std::cout << "Setting up the b-jet scale factors... " << std::endl;
  // Get the b-tagging algorithm to use to initialize the appropriate csv file:
  if(stats.bfind("ApplyJetBTaggingCSVv2")){
    btagalgoname = "CSVv2"; 
    btagsffilename = csvv2btagsfcsvfiles[year]; 
  }
  else if(stats.bfind("ApplyJetBTaggingDeepCSV")){
    btagalgoname = "DeepCSV";
    btagsffilename = deepcsvbtagsfcsvfiles[year];

  }
  else if(stats.bfind("ApplyJetBTaggingDeepFlav")){
    btagalgoname = "DeepFlav";
    btagsffilename = deepflavbtagsfcsvfiles[year];
  }

  try{
    btagcalib = BTagCalibration(btagalgoname, (PUSPACE+"BJetDatabase/"+btagsffilename).c_str());

    if(stats.bfind("UseBtagSF")){  
      std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
      std::cout << "Setting up the b-jet scale factors... " << std::endl;
      std::cout << "B-jet ID algorithm selected: " << btagalgoname << std::endl;
      std::cout << "Selected working point: " << stats.smap.at("JetBTaggingWP") << std::endl;
      std::cout << "B-tagging SF csv file: " << btagsffilename << std::endl;
      std::cout << "-----------------------------------------------------------------------------------------" << std::endl;
    }
  }
  catch(std::exception& e){
    std::cerr << "ERROR in setupBJetSFInfo: " << e.what() << std::endl;
    std::cout << "\t Setting dummy b-tagging from " << (PUSPACE+"btagging.csv").c_str() << std::endl;
    btagcalib = BTagCalibration("DeepCSV", (PUSPACE+"btagging.csv").c_str());   
  }

 // Check BTagCalibrationStandalone.h for more information

   static std::map<std::string, int> bjetflavors = {
    {"bjet", 0}, {"cjet", 1}, {"lightjet", 2},
   };

   bjetflavor = (BTagEntry::JetFlavor) bjetflavors["bjet"];

   static std::map<std::string, int> btagoperatingpoints = {
    {"loose", 0}, {"medium", 1}, {"tight", 2}, {"reshaping", 3}
   };

   b_workingpoint = (BTagEntry::OperatingPoint) btagoperatingpoints[stats.smap.at("JetBTaggingWP").c_str()];

}

double Analyzer::getBJetSF(CUTS ePos, const PartStats& stats) {
  double bjetSFall = 1.0, bjetSFtemp = 1.0; 

  if(!neededCuts.isPresent(ePos)) return bjetSFall;

  // Load the info from the btaggin reader
  btagsfreader = BTagCalibrationReader(b_workingpoint, "central");
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  if(!stats.bfind("UseBtagSF")){
    bjetSFtemp = 1.0;
  }
  else{

    if(active_part->at(CUTS::eRBJet)->size() == 0){
      bjetSFtemp = 1.0;
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 1){
      bjetSFtemp = btagsfreader.eval_auto_bounds("central", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 2){
      // Get the SF for the first b-jet
      bjetSFtemp = btagsfreader.eval_auto_bounds("central", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      // Now multiply by the SF of the second b-jet
      bjetSFtemp = bjetSFtemp * btagsfreader.eval_auto_bounds("central", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Pt());
    }

  }

  // Calculate the full SF
  bjetSFall = bjetSFall * bjetSFtemp;

  return bjetSFall;
}

double Analyzer::getBJetSFResUp(CUTS ePos, const PartStats& stats) {
  double bjetSFall = 1.0, bjetSFtemp = 1.0; 

  if(!neededCuts.isPresent(ePos)) return bjetSFall;

  // Load the info from the btaggin reader
  btagsfreader = BTagCalibrationReader(b_workingpoint, "up");
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  if(!stats.bfind("UseBtagSF")){
    bjetSFtemp = 1.0;
  }
  else{

    if(active_part->at(CUTS::eRBJet)->size() == 0){
      bjetSFtemp = 1.0;
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 1){
      bjetSFtemp = btagsfreader.eval_auto_bounds("up", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 2){
      // Get the SF for the first b-jet
      bjetSFtemp = btagsfreader.eval_auto_bounds("up", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      // Now multiply by the SF of the second b-jet
      bjetSFtemp = bjetSFtemp * btagsfreader.eval_auto_bounds("up", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Pt());
    }

  }

  // Calculate the full SF
  bjetSFall = bjetSFall * bjetSFtemp;

  return bjetSFall;
}

double Analyzer::getBJetSFResDown(CUTS ePos, const PartStats& stats) { 
  double bjetSFall = 1.0, bjetSFtemp = 1.0; 

  if(!neededCuts.isPresent(ePos)) return bjetSFall;

  btagsfreader = BTagCalibrationReader(b_workingpoint, "down");
  btagsfreader.load(btagcalib, bjetflavor, "comb");

  if(!stats.bfind("UseBtagSF")){
    bjetSFtemp = 1.0;
  }
  else{
    if(active_part->at(CUTS::eRBJet)->size() == 0){
      bjetSFtemp = 1.0;
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 1){
      bjetSFtemp = btagsfreader.eval_auto_bounds("down", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
    }
    else if(active_part->at(CUTS::eRBJet)->size() == 2){
      // Get the SF for the first b-jet
      bjetSFtemp = btagsfreader.eval_auto_bounds("down", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(0))).Pt());
      // Now multiply by the SF of the second b-jet
      bjetSFtemp = bjetSFtemp * btagsfreader.eval_auto_bounds("down", bjetflavor, (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Eta(), (_Jet->p4(active_part->at(CUTS::eRBJet)->at(1))).Pt());
    }
  }

  // Calculate the full SF
  bjetSFall = bjetSFall * bjetSFtemp;

  return bjetSFall;
}

////FatJet specific function for finding the number of V-jets that pass the cuts.
void Analyzer::getGoodRecoFatJets(CUTS ePos, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePos)) return;

  std::string systname = syst_names.at(syst);
  if(!_FatJet->needSyst(syst)) {
    active_part->at(ePos)=goodParts[ePos];
    return;
  }

  int i=0;

  for(auto lvec: *_FatJet) {
    bool passCuts = true;
    passCuts = passCuts && passCutRange(fabs(lvec.Eta()), stats.pmap.at("EtaCut"));
    passCuts = passCuts && (lvec.Pt() > stats.dmap.at("PtCut")) ;

    ///if else loop for central jet requirements
    for( auto cut: stats.bset) {
      if(!passCuts) break;

      else if(cut == "ApplyJetWTagging") passCuts = passCuts && (passCutRange(_FatJet->tau2[i]/_FatJet->tau1[i], stats.pmap.at("JetTau2Tau1Ratio")) &&
                                                     passCutRange(_FatJet->PrunedMass[i], stats.pmap.at("JetWmassCut")));
    // ----anti-overlap requirements
      else if(cut == "RemoveOverlapWithMuon1s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon1, stats.dmap.at("Muon1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithMuon2s") passCuts = passCuts && !isOverlaping(lvec, *_Muon, CUTS::eRMuon2, stats.dmap.at("Muon2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron1s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec1, stats.dmap.at("Electron1MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithElectron2s") passCuts = passCuts && !isOverlaping(lvec, *_Electron, CUTS::eRElec2, stats.dmap.at("Electron2MatchingDeltaR"));
      else if(cut == "RemoveOverlapWithTau1s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau1, stats.dmap.at("Tau1MatchingDeltaR"));
      else if (cut =="RemoveOverlapWithTau2s") passCuts = passCuts && !isOverlaping(lvec, *_Tau, CUTS::eRTau2, stats.dmap.at("Tau2MatchingDeltaR"));

    }
    if(passCuts) active_part->at(ePos)->push_back(i);
    i++;
  }
}

///function to see if a lepton is overlapping with another particle.  Used to tell if jet or tau
//came ro decayed into those leptons
bool Analyzer::isOverlaping(const TLorentzVector& lvec, Lepton& overlapper, CUTS ePos, double MatchingDeltaR) {
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
  }
  return false;
}

bool Analyzer::isOverlapingB(const TLorentzVector& lvec, Jet& overlapper, CUTS ePos, double MatchingDeltaR) { //01.17.19
  for(auto it : *active_part->at(ePos)) {
    if(lvec.DeltaR(overlapper.p4(it)) < MatchingDeltaR) return true;
  }
  return false;
}

///Tests if tau decays into the specified number of jet prongs.
bool Analyzer::passProng(std::string prong, int value) {
  return ( (prong.find("1") != std::string::npos &&  (value<5)) ||
  (prong.find("2") != std::string::npos &&  (value>=5 && value<10)) ||
  (prong.find("3") != std::string::npos && (value>=10 && value<12)) );
}


////Tests if tau is within the cracks of the detector (the specified eta ranges)
bool Analyzer::isInTheCracks(float etaValue){
  return (fabs(etaValue) < 0.018 ||
  (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
  (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
  (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
  (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}


///sees if the event passed one of the two cuts provided
void Analyzer::TriggerCuts(CUTS ePos) {
	if(! neededCuts.isPresent(ePos)) return;

	// Loop over all elements of the trigger decisions vector
	for(bool decision : triggernamedecisions){
		if(decision){
			// If one element is true (1), then store back the event in the triggers vector
			active_part->at(ePos)->push_back(0);
			// Clean up the trigger decisions vector to reduce memory usage and have an empty vector for the next event
			triggernamedecisions.clear();
			// End of the function
			return;
		}
	}
	// If all the elements of the trigger decisions vector are false, then just clean up the trigger decisions vector to reduce memory usage.
	triggernamedecisions.clear();
}


////VBF specific cuts dealing with the leading jets.
void Analyzer::VBFTopologyCut(const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eSusyCom)) return;
  std::string systname = syst_names.at(syst);


  if(systname!="orig"){
    //only jet stuff is affected
    //save time to not rerun stuff
    if( systname.find("Jet")==std::string::npos){
      active_part->at(CUTS::eSusyCom)=goodParts[CUTS::eSusyCom];
      return;
    }
  }

  if(active_part->at(CUTS::eR1stJet)->size()==0 || active_part->at(CUTS::eR2ndJet)->size()==0) return;

  TLorentzVector ljet1 = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0));
  TLorentzVector ljet2 = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));
  TLorentzVector dijet = ljet1 + ljet2;
  double dphi1 = normPhi(ljet1.Phi() - _MET->phi());
  double dphi2 = normPhi(ljet2.Phi() - _MET->phi());

  bool passCuts = true;
  for(auto cut: stats.bset) {
    if(!passCuts) break;
    else if(cut == "DiscrByMass") passCuts = passCuts && passCutRange(dijet.M(), stats.pmap.at("MassCut"));
    else if(cut == "DiscrByPt") passCuts = passCuts && passCutRange(dijet.Pt(), stats.pmap.at("PtCut"));
    else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(ljet1.Eta() - ljet2.Eta()), stats.pmap.at("DeltaEtaCut"));
    else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(ljet1.Phi() - ljet2.Phi()), stats.pmap.at("DeltaPhiCut"));
    else if(cut == "DiscrByOSEta") passCuts = passCuts && (ljet1.Eta() * ljet2.Eta() < 0);
    else if(cut == "DiscrByR1") passCuts = passCuts && passCutRange(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0)), stats.pmap.at("R1Cut"));
    else if(cut == "DiscrByR2") passCuts = passCuts && passCutRange(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), stats.pmap.at("R2Cut"));
    else if(cut == "DiscrByAlpha") {
      double alpha = (dijet.M() > 0) ? ljet2.Pt() / dijet.M() : -1;
      passCuts = passCuts && passCutRange(alpha, stats.pmap.at("AlphaCut"));
    }
    else if(cut == "DiscrByDphi1") passCuts = passCuts && passCutRange(abs(dphi1), stats.pmap.at("Dphi1Cut"));
    else if(cut == "DiscrByDphi2") passCuts = passCuts && passCutRange(abs(dphi2), stats.pmap.at("Dphi2Cut"));

    else std::cout << "cut: " << cut << " not listed" << std::endl;
  }

  if(passCuts)  active_part->at(CUTS::eSusyCom)->push_back(0);
  return;
}

bool Analyzer::passCutRange(double value, const std::pair<double, double>& cuts) {
  return (value > cuts.first && value < cuts.second);
}

bool Analyzer::passCutRangeAbs(double value, const std::pair<double, double>& cuts) {
  return ((value > cuts.first && value < cuts.second) || (value > (cuts.second * -1) && value < (cuts.first * -1)));
}

//-----Calculate lepton+met transverse mass
double Analyzer::calculateLeptonMetMt(const TLorentzVector& Tobj) {
  double px = Tobj.Px() + _MET->px();
  double py = Tobj.Py() + _MET->py();
  double et = Tobj.Et() + _MET->energy();
  double mt2 = et*et - (px*px + py*py);
  return (mt2 >= 0) ? sqrt(mt2) : -1;
}


/////Calculate the diparticle mass based on how to calculate it
///can use Collinear Approximation, which can fail (number failed available in a histogram)
///can use VectorSumOfVisProductAndMet which is sum of particles and met
///Other which is adding without met
double Analyzer::diParticleMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string howCalc) {
  bool ratioNotInRange = false;
  TLorentzVector The_LorentzVect;

  if(howCalc == "InvariantMass") {
    return (Tobj1 + Tobj2).M();
  }


  //////check this equation/////
  if(howCalc == "CollinearApprox") {
    double denominator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1 = (Tobj2.Py()*_MET->px() - Tobj2.Px()*_MET->py())/denominator;
    double x2 = (Tobj1.Px()*_MET->py() - Tobj1.Py()*_MET->px())/denominator;
    ratioNotInRange=!((x1 < 0.) && (x2 < 0.));
    if (!ratioNotInRange) {
      The_LorentzVect.SetPxPyPzE( (Tobj1.Px()*(1 + x1) + Tobj2.Px()*(1+x2)), (Tobj1.Py()*(1+x1) + Tobj2.Py()*(1+x2)), (Tobj1.Pz()*(1+x1) + Tobj2.Pz()*(1+x2)), (Tobj1.Energy()*(1+x1) + Tobj2.Energy()*(1+x2)) );
      return The_LorentzVect.M();
    }
  }

  if(howCalc == "VectorSumOfVisProductsAndMet" || ratioNotInRange) {
    return (Tobj1 + Tobj2 + _MET->p4()).M();
  }

  return (Tobj1 + Tobj2).M();
}

////Tests if the CollinearApproximation works for finding the mass of teh particles
bool Analyzer::passDiParticleApprox(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string howCalc) {
  if(howCalc == "CollinearApprox") {
    double x1_numerator = (Tobj1.Px() * Tobj2.Py()) - (Tobj2.Px() * Tobj1.Py());
    double x1_denominator = (Tobj2.Py() * (Tobj1.Px() + _MET->px())) - (Tobj2.Px() * (Tobj1.Py() + _MET->py()));
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
    double x2_numerator = x1_numerator;
    double x2_denominator = (Tobj1.Px() * (Tobj2.Py() + _MET->py())) - (Tobj1.Py() * (Tobj2.Px() + _MET->px()));
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
    return (x1 > 0. && x1 < 1.) && (x2 > 0. && x2 < 1.);
  } else {
    return true;
  }
}


/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonCombos(Lepton& lep1, Lepton& lep2, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  std::string systname = syst_names.at(syst);

  if(!lep1.needSyst(syst) && !lep2.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  bool sameParticle = (&lep1 == &lep2);
  TLorentzVector part1, part2, llep;

  for(auto i1 : *active_part->at(ePos1)) {
    part1 = lep1.p4(i1);
    for(auto i2 : *active_part->at(ePos2)) {
      if(sameParticle && i2 <= i1) continue;
      part2 = lep2.p4(i2);
      bool passCuts = true;
      for(auto cut : stats.bset) {
        if(!passCuts) break;
        else if (cut == "DiscrByDeltaR") passCuts = passCuts && (part1.DeltaR(part2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(part1.Phi() - part2.Phi())), stats.pmap.at("CosDphiCut"));
        else if(cut == "DiscrByDeltaPt") passCuts = passCuts && passCutRange(part1.Pt() - part2.Pt(), stats.pmap.at("DeltaPtCutValue"));
        else if(cut == "DiscrByCDFzeta2D") {
          double CDFzeta = stats.dmap.at("PZetaCutCoefficient") * getPZeta(part1, part2).first
            + stats.dmap.at("PZetaVisCutCoefficient") * getPZeta(part1, part2).second;
          passCuts = passCuts && passCutRange(CDFzeta, stats.pmap.at("CDFzeta2DCutValue"));
        }
        else if(cut == "DiscrByDeltaPtDivSumPt") {
          double ptDiv = (part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt());
          passCuts = passCuts && passCutRange(ptDiv, stats.pmap.at("DeltaPtDivSumPtCutValue"));
        }
        else if (cut == "DiscrByMassReco") {
          double diMass = diParticleMass(part1,part2, stats.smap.at("HowCalculateMassReco"));
          passCuts = passCuts && passCutRange(diMass, stats.pmap.at("MassCut"));
        }
        else if(cut == "DiscrByCosDphiPtAndMet"){
          double CosDPhi1 = cos(absnormPhi(part1.Phi() - _MET->phi()));
          passCuts = passCuts && passCutRange(CosDPhi1, stats.pmap.at("CosDphiPtAndMetCut"));
        }
	// ---------- New: DY Z' team ---------- //
	else if(cut == "DiscrByCosDphiLeadPtAndMet"){
	  if(part1.Pt() > part2.Pt()){
	    llep = part1;
	  }
	  else{
	    llep = part2;
	  }
	  double CosDPhiLead = cos(absnormPhi(llep.Phi() - _MET->phi()));
	  passCuts = passCuts && passCutRange(CosDPhiLead, stats.pmap.at("CosDphiLeadPtAndMetCut"));
	}
  else if (cut == "DiscrByAbsCosDphiLeadPtandMet"){
    if (part1.Pt() > part2.Pt()){
      llep = part1;
    }
    else{
      llep = part2;
    }
    double CosDPhiLead_forabs = cos(absnormPhi(llep.Phi() - _MET->phi()));
    passCuts = passCuts && passCutRangeAbs(CosDPhiLead_forabs, stats.pmap.at("AbsCosDphiLeadPtAndMetCut"));
  }
	else if(cut == "DiscrByMtLeadPtAndMet"){
	   if(part1.Pt() > part2.Pt()){
	     llep = part1;
	  }
	  else{
	    llep = part2;
	  }
	  double mTlead = calculateLeptonMetMt(llep);
	  passCuts = passCuts && passCutRange(mTlead, stats.pmap.at("MtLeadPtAndMetCut"));
	}
	else if(cut == "DiscrByDiLepMassDeltaPt"){
	  double dilepmass = CalculateDiLepMassDeltaPt(part1, part2);
	  passCuts = passCuts && passCutRange(dilepmass, stats.pmap.at("DiLeadMassDeltaPtCut"));
	}
	// ---------------------------------------- //

        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      if (stats.bfind("DiscrByOSLSType")){
        //   if it is 1 or 0 it will end up in the bool std::map!!
        if(stats.bfind("DiscrByOSLSType") && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
      }else if (stats.dmap.find("DiscrByOSLSType") != stats.dmap.end() ){
        if(lep1.charge(i1) * lep2.charge(i2) > 0) continue;
      }else if (stats.smap.find("DiscrByOSLSType") != stats.smap.end() ){
        if(stats.smap.at("DiscrByOSLSType") == "LS" && (lep1.charge(i1) * lep2.charge(i2) <= 0)) continue;
        else if(stats.smap.at("DiscrByOSLSType") == "OS" && (lep1.charge(i1) * lep2.charge(i2) >= 0)) continue;
      }

      ///Particles that lead to good combo are nGen * part1 + part2
      /// final / nGen = part1 (make sure is integer)
      /// final % nGen = part2
      if(passCuts)
        active_part->at(ePosFin)->push_back(i1*BIG_NUM + i2);
    }
  }
}

// ---------- New function: DY Z' team ---------- //	     
double Analyzer::CalculateDiLepMassDeltaPt(const TLorentzVector& part1, const TLorentzVector& part2){
      double pt1  = part1.Pt();
      double eta1 = part1.Eta();
      double phi1 = part1.Phi();
      double mass1= 0.000511; // electron mass in GeV                                                           
      double pt2  = part2.Pt();
      double eta2 = part2.Eta();
      double phi2 = part2.Phi();
      double mass2= 1.777;    // tau mass in GeV
      double mass3= 0.0;

      double px1 = pt1*cos(phi1);
      double py1 = pt1*sin(phi1);
      double pz1 = pt1*sinh(eta1);
      double E1  = sqrt( pow(pt1*cosh(eta1),2) + pow(mass1,2) );
      double px2 = pt2*cos(phi2);
      double py2 = pt2*sin(phi2);
      double pz2 = pt2*sinh(eta2);
      double E2  = sqrt( pow(pt2*cosh(eta2),2) + pow(mass2,2) );
      double px3 = -(px1 + px2);
      double py3 = -(py1 + py2);
      double pz3 = -0.0;
      double E3  = sqrt( pow(px3,2) + pow(py3,2) + pow(pz3,2) + pow(mass3,2) );

      double diMass= sqrt(pow(E1+E2+E3,2)-pow(px1+px2+px3,2)-pow(py1+py2+py3,2)-pow(pz1+pz2+pz3,2));

      return diMass;
}	     
// ---------------------------------------- //
	     
/////abs for values
///Find the number of lepton combos that pass the dilepton cuts
void Analyzer::getGoodLeptonJetCombos(Lepton& lep1, Jet& jet1, CUTS ePos1, CUTS ePos2, CUTS ePosFin, const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(ePosFin)) return;
  std::string systname = syst_names.at(syst);
  if(!lep1.needSyst(syst) && !jet1.needSyst(syst)) {
    active_part->at(ePosFin)=goodParts[ePosFin];
    return;
  }

  TLorentzVector llep1, ljet1;
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(ePos1)) {
    llep1 = lep1.p4(ij2);
    for(auto ij1 : *active_part->at(ePos2)) {
      ljet1 = _Jet->p4(ij1);

      bool passCuts = true;
      for(auto cut : stats.bset) {
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = (ljet1.DeltaR(llep1) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByDeltaEta") passCuts = passCutRange(abs(ljet1.Eta() - llep1.Eta()), stats.pmap.at("DeltaEtaCut"));
        else if(cut == "DiscrByDeltaPhi") passCuts = passCutRange(absnormPhi(ljet1.Phi() - llep1.Phi()), stats.pmap.at("DeltaPhiCut"));
        else if(cut == "DiscrByOSEta") passCuts = (ljet1.Eta() * llep1.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCutRange((ljet1+llep1).M(), stats.pmap.at("MassCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCutRange(cos(absnormPhi(ljet1.Phi() - llep1.Phi())), stats.pmap.at("CosDphiCut"));
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(ePosFin)->push_back(ij1*_Jet->size() + ij2);
    }
  }
}


//////////////LOOK INTO DIJET PICKING
///////HOW TO GET RID OF REDUNCENCIES??

/////Same as gooddilepton, just jet specific
void Analyzer::getGoodDiJets(const PartStats& stats, const int syst) {
  if(! neededCuts.isPresent(CUTS::eDiJet)) return;
  std::string systname = syst_names.at(syst);
  if(systname!="orig"){
    //save time to not rerun stuff
    if( systname.find("Jet")==std::string::npos){
      active_part->at(CUTS::eDiJet)=goodParts[CUTS::eDiJet];
      return;
    }
  }
  TLorentzVector jet1, jet2;
  // ----Separation cut between jets (remove overlaps)
  for(auto ij2 : *active_part->at(CUTS::eRJet2)) {
    jet2 = _Jet->p4(ij2);
    for(auto ij1 : *active_part->at(CUTS::eRJet1)) {
      if(ij1 == ij2) continue;
      jet1 = _Jet->p4(ij1);

      bool passCuts = true;
      //cout<<"---------------------"<<std::endl;
      for(auto cut : stats.bset) {
        //cout<<cut<<"    "<<passCuts<<std::endl;;
        if(!passCuts) break;
        else if(cut == "DiscrByDeltaR") passCuts = passCuts && (jet1.DeltaR(jet2) >= stats.dmap.at("DeltaRCut"));
        else if(cut == "DiscrByDeltaEta") passCuts = passCuts && passCutRange(abs(jet1.Eta() - jet2.Eta()), stats.pmap.at("DeltaEtaCut"));
        else if(cut == "DiscrByDeltaPhi") passCuts = passCuts && passCutRange(absnormPhi(jet1.Phi() - jet2.Phi()), stats.pmap.at("DeltaPhiCut"));
        else if(cut == "DiscrByOSEta") passCuts = passCuts && (jet1.Eta() * jet2.Eta() < 0);
        else if(cut == "DiscrByMassReco") passCuts = passCuts && passCutRange((jet1+jet2).M(), stats.pmap.at("MassCut"));
        else if(cut == "DiscrByCosDphi") passCuts = passCuts && passCutRange(cos(absnormPhi(jet1.Phi() - jet2.Phi())), stats.pmap.at("CosDphiCut"));
        else std::cout << "cut: " << cut << " not listed" << std::endl;
      }
      ///Particlesp that lead to good combo are totjet * part1 + part2
      /// final / totjet = part1 (make sure is integer)
      /// final % totjet = part2
      if(passCuts) active_part->at(CUTS::eDiJet)->push_back(ij1*_Jet->size() + ij2);
    }
  }
}

///////Only tested for if is Zdecay, can include massptasymmpair later?
/////Tests to see if a light lepton came form a zdecay
bool Analyzer::isZdecay(const TLorentzVector& theObject, const Lepton& lep) {
  bool eventIsZdecay = false;
  const float zMass = 90.1876;
  const float zWidth = 2.4952;
  float zmmPtAsymmetry = -10.;

  // if mass is within 3 sigmas of z or pt asymmetry is small set to true.
  for(std::vector<TLorentzVector>::const_iterator lepit= lep.begin(); lepit != lep.end(); lepit++) {
    if(theObject.DeltaR(*lepit) < 0.3) continue;
    if(theObject == (*lepit)) continue;

    TLorentzVector The_LorentzVect = theObject + (*lepit);
    zmmPtAsymmetry = (theObject.Pt() - lepit->Pt()) / (theObject.Pt() + lepit->Pt());

    if( (abs(The_LorentzVect.M() - zMass) < 3.*zWidth) || (fabs(zmmPtAsymmetry) < 0.20) ) {
      eventIsZdecay = true;
      break;
    }
  }

  return eventIsZdecay;
}


///Calculates the Pzeta value
std::pair<double, double> Analyzer::getPZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2) {
  double zetaX = cos(Tobj1.Phi()) + cos(Tobj2.Phi());
  double zetaY = sin(Tobj1.Phi()) + sin(Tobj2.Phi());
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) { zetaX /= zetaR; zetaY /= zetaR; }
  double visPx = Tobj1.Px() + Tobj2.Px();
  double visPy = Tobj1.Py() + Tobj2.Py();
  double px = visPx + _MET->px();
  double py = visPy + _MET->py();
  return std::make_pair(px*zetaX + py*zetaY, visPx*zetaX + visPy*zetaY);
}

void Analyzer::checkParticleDecayList(){
  std::fstream file;
  file.open("particle_decay_list.txt", std::fstream::in | std::fstream::out);
  int s;
  if (file.is_open()){
    std::cout << "Warning, file already exists.";
    std::cout << "Do you wish to clear the file? 1 for yes; 0 for no.";
    std::cin >> s;
    if (s == 1)
      {file.open("particle_decay_list.txt", std::ios::out | std::ios::trunc);
	file.close();
	std::cout << "You have cleared the file.";}
  }
}

void Analyzer::writeParticleDecayList(int event){  //01.16.19
  BOOM->GetEntry(event);
  std::fstream file("particle_decay_list.txt", std::fstream::in | std::fstream::out | std::fstream::app);
  if (file.is_open()){
  for (unsigned p=0; p < _Gen->size(); p++){  
    file << std::setw(2) << p << std::setw(2) << " " << std::setw(8) << "pdg_id: " << std::setw(4) << abs(_Gen->pdg_id[p]) << std::setw(2) << "  " << std::setw(5) << "p_T: " << std::setw(10) << _Gen->pt(p) << std::setw(2) << "  " << std::setw(5) << "phi: " << std::setw(10) << _Gen->phi(p) << std::setw(10) << "status: " << std::setw(3) << _Gen->status[p] << std::setw(2) << "  " << std::setw(7) << "mind:" << std::setw(1) << " " << std::setw(2) << _Gen->genPartIdxMother[p] << "\n";
  }
  file << "----------" << "\n";
  file.close();}
  else std::cout << "Unable to open file." << std::endl;
  return;
}

std::multimap<int,int> Analyzer::readinJSON(std::string year){ //NEW:  function for reading in the JSON file in the format that I established.  This happens once before events start entering preprocess.
  // Newer: this function will take the year as an input to set the JSON file to be read.
  // The naming convention is json201X.txt, make sure your JSON files match this convention.
  std::fstream fs((PUSPACE+"json"+year+".txt").c_str(), std::fstream::in); // Take the info of the JSON file as input.	
  std::string line; //NEW:  need this since we'll be analyzing line by line.
  std::multimap<int,int> json_line_dict; //NEW:  we'll work with the information as a multimap (C++ equivalent of a python dictionary).
  
  while(!fs.eof()){ //NEW:  while the file is open...
      getline(fs, line); //NEW:  grab each line.
      if(line.size() == 0) continue;  //NEW:  if there's nothing in the line, skip it.
      std::vector<std::string> vals = string_split(line, {","});  //NEW:  put the contents of the line into a vector.
      int run_num_line = (stringtotype<int>(vals[0])); //NEW:  the first element of the line is the run_number.
      for(size_t p=1; p < vals.size(); p=p+1){  //NEW:  step through the remaining elements of the line.
	json_line_dict.insert (std::make_pair(run_num_line,stringtotype<int>(vals[p]))); //NEW:  push the remaining lumisection bounds as each an element corresponding to the run_number key.
      }
  }
  fs.close(); // Close the file once the loop has gone over all lines in the file.

  return json_line_dict;  //NEW:  returns the multimap that you need for proper filtering.
}

double Analyzer::getZBoostWeight(){
  double boostweigth=1.;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGZ)->size() ==1 || active_part->at(CUTS::eGW)->size() ==1)){
    //cout<<" Z or W " <<std::endl;
    double boostz = 0;
    if(active_part->at(CUTS::eGZ)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGZ)->at(0));
    }
    if(active_part->at(CUTS::eGW)->size() ==1){
      boostz = _Gen->pt(active_part->at(CUTS::eGW)->at(0));
    }
    if(boostz > 0 && boostz <= 50) {boostweigth = 1.1192;}
    else if (boostz > 50 && boostz <= 100) {boostweigth = 1.1034;}
    else if (boostz > 100 && boostz <= 150) {boostweigth = 1.0675;}
    else if (boostz > 150 && boostz <= 200) {boostweigth = 1.0637;}
    else if (boostz > 200 && boostz <= 300) {boostweigth = 1.0242;}
    else if (boostz > 300 && boostz <= 400) {boostweigth = 0.9453;}
    else if (boostz > 400 && boostz <= 600) {boostweigth = 0.8579;}
    else if (boostz >= 600) {boostweigth = 0.7822;}
    else {boostweigth = 1;}
  }
  return boostweigth;
}


double Analyzer::getWkfactor(){
  double kfactor=1.;
  if(!isWSample)
    return kfactor;
  if((active_part->at(CUTS::eGElec)->size() + active_part->at(CUTS::eGTau)->size() + active_part->at(CUTS::eGMuon)->size()) >=1 && (active_part->at(CUTS::eGW)->size() ==1)){
    //this k-factor is not computed for the low mass W!
    double wmass=_Gen->p4(active_part->at(CUTS::eGW)->at(0)).M();
    if(wmass<100){
      return 1.;
    }
    if(active_part->at(CUTS::eGTau)->size()){
      kfactor=k_ele_h->GetBinContent(k_ele_h->FindBin(wmass));
    }
    else if(active_part->at(CUTS::eGMuon)->size()){
      kfactor=k_mu_h->GetBinContent(k_mu_h->FindBin(wmass));
    }
    else if(active_part->at(CUTS::eGElec)->size()){
      kfactor=k_tau_h->GetBinContent(k_tau_h->FindBin(wmass));
    }
  }
  return kfactor;
}


////Grabs a list of the groups of histograms to be filled and asked Fill_folder to fill up the histograms
void Analyzer::fill_histogram() {
  
  if(!isData && distats["Run"].bfind("ApplyGenWeight") && gen_weight == 0.0) return;

  if(isData && blinded && maxCut == SignalRegion) return;

  const std::vector<std::string>* groups = histo.get_groups();
  if(!isData){
    wgt = 1.;
    //wgt *= getTopBoostWeight(); //01.15.19
    if(distats["Run"].bfind("UsePileUpWeight")) wgt*= pu_weight;
    if(distats["Run"].bfind("ApplyGenWeight")) wgt *= (gen_weight > 0) ? 1.0 : -1.0;
    //add weight here
    if(distats["Run"].bfind("ApplyTauIDSF")) wgt *= getTauDataMCScaleFactor(0);

    if(distats["Run"].bfind("ApplyZBoostSF") && isVSample){
      wgt *= getZBoostWeight();
    }
    if(distats["Run"].bfind("ApplyWKfactor")){
      wgt *= getWkfactor();
    }
    wgt *= getBJetSF(CUTS::eRBJet, _Jet->pstats["BJet"]); //01.16.19
  }else  wgt=1.;
  //backup current weight
  backup_wgt=wgt;

  for(size_t i = 0; i < syst_names.size(); i++) {
    for(Particle* ipart: allParticles) ipart->setCurrentP(i);
    _MET->setCurrentP(i);
    active_part =&syst_parts.at(i);
    //////i == 0 is orig or no syst case
    if(i == 0) {
      active_part = &goodParts;
      fillCuts(true);
      for(auto it: *groups) {
        fill_Folder(it, maxCut, histo, false);
      }
      if(!fillCuts(false)) {
        fill_Tree();
      }
    }else{
      wgt=backup_wgt;
      if(syst_names[i].find("weight")!=std::string::npos){
        if(syst_names[i]=="Tau_weight_Up"){
          if(distats["Run"].bfind("ApplyTauIDSF")) {
            wgt/=getTauDataMCScaleFactor(0);
            wgt *= getTauDataMCScaleFactor(1);
          }
        }else if(syst_names[i]=="Tau_weight_Down"){
          if(distats["Run"].bfind("ApplyTauIDSF")) {
            wgt/=getTauDataMCScaleFactor(0);
            wgt *= getTauDataMCScaleFactor(-1);
          }
        }
        if(syst_names[i]=="Pileup_weight_Up"){
          if(distats["Run"].bfind("UsePileUpWeight")) {
            wgt/=   pu_weight;
            wgt *=  hPU_up[(int)(nTruePU+1)];
          }
        }else if(syst_names[i]=="Pileup_weight_Down"){
          if(distats["Run"].bfind("UsePileUpWeight")) {
            wgt/=   pu_weight;
            wgt *=  hPU_down[(int)(nTruePU+1)];
          }
        }
      }

      if(syst_names[i].find("Btag")!=std::string::npos){ //01.16.19
        if(syst_names[i]=="Btag_Up"){
	  wgt/=getBJetSF(CUTS::eRBJet, _Jet->pstats["BJet"]);
	  wgt*=getBJetSFResUp(CUTS::eRBJet, _Jet->pstats["BJet"]);
        }else if(syst_names[i]=="Btag_Down"){
	  wgt/=getBJetSF(CUTS::eRBJet, _Jet->pstats["BJet"]);
          wgt*=getBJetSFResDown(CUTS::eRBJet, _Jet->pstats["BJet"]);
        }
      }

      //get the non particle conditions:
      for(auto itCut : nonParticleCuts){
        active_part->at(itCut)=goodParts.at(itCut);
      }
      if(!fillCuts(false)) continue;
      for(auto it: *syst_histo.get_groups()) {
        fill_Folder(it, i, syst_histo, true);
      }
      wgt=backup_wgt;
    }
  }
  for(Particle* ipart: allParticles) ipart->setCurrentP(0);
  _MET->setCurrentP(0);
  active_part = &goodParts;
}

///Function that fills up the histograms
void Analyzer::fill_Folder(std::string group, const int max, Histogramer &ihisto, bool issyst) {
  /*be aware in this function
   * the following definition is used:
   * histAddVal(val, name) histo.addVal(val, group, max, name, wgt)
   * so each histogram knows the group, max and weight!
   */
  if(group == "FillRun" && (&ihisto==&histo)) {
    if(crbins != 1) {
      for(int i = 0; i < crbins; i++) {
        ihisto.addVal(false, group, i, "Events", 1);
        if(distats["Run"].bfind("ApplyGenWeight")) {
          //put the weighted events in bin 3
          ihisto.addVal(2, group,i, "Events", (gen_weight > 0) ? 1.0 : -1.0);
        }
        ihisto.addVal(wgt, group, i, "Weight", 1);
      }
    }
    else{
      ihisto.addVal(false, group,ihisto.get_maxfolder(), "Events", 1);
      if(distats["Run"].bfind("ApplyGenWeight")) {
        //put the weighted events in bin 3
        ihisto.addVal(2, group,ihisto.get_maxfolder(), "Events", (gen_weight > 0) ? 1.0 : -1.0);
      }
      ihisto.addVal(wgt, group, ihisto.get_maxfolder(), "Weight", 1);
      ihisto.addVal(nTruePU, group, ihisto.get_maxfolder(), "PUWeight", 1);
    }
    histAddVal(true, "Events");
    histAddVal(bestVertices, "NVertices");
  } else if(!isData && group == "FillGen") {

    int nhadtau = 0;
    TLorentzVector genVec;
    int i = 0;
    for(vec_iter it=active_part->at(CUTS::eGTau)->begin(); it!=active_part->at(CUTS::eGTau)->end(); it++, i++) {
      int nu = active_part->at(CUTS::eGNuTau)->at(i);
      if(nu != -1) {
        genVec = _Gen->p4(*it) - _Gen->p4(nu);
        histAddVal(genVec.Pt(), "HadTauPt");
        histAddVal(genVec.Eta(), "HadTauEta");
        nhadtau++;
      }
      histAddVal(_Gen->energy(*it), "TauEnergy");
      histAddVal(_Gen->pt(*it), "TauPt");
      histAddVal(_Gen->eta(*it), "TauEta");
      histAddVal(_Gen->phi(*it), "TauPhi");
      for(vec_iter it2=it+1; it2!=active_part->at(CUTS::eGTau)->end(); it2++) {
        histAddVal(diParticleMass(_Gen->p4(*it),_Gen->p4(*it2), "none"), "DiTauMass");
      }
    }

    for(vec_iter genhadtau_it = active_part->at(CUTS::eGHadTau)->begin(); genhadtau_it != active_part->at(CUTS::eGHadTau)->end(); genhadtau_it++){
        histAddVal(_GenHadTau->pt(*genhadtau_it), "VisHadTauPt");
        histAddVal(_GenHadTau->eta(*genhadtau_it), "VisHadTauEta");
        histAddVal(_GenHadTau->phi(*genhadtau_it), "VisHadTauPhi");
        histAddVal(_GenHadTau->p4(*genhadtau_it).M(), "VisHadTauMass");
        histAddVal(_GenHadTau->decayMode[*genhadtau_it], "VisHadTauDecayMode");
    }
    histAddVal(active_part->at(CUTS::eGHadTau)->size(), "NVisHadTau");

    for(vec_iter matchedgentauh_it = active_part->at(CUTS::eGMatchedHadTau)->begin(); matchedgentauh_it != active_part->at(CUTS::eGMatchedHadTau)->end(); matchedgentauh_it++){
        histAddVal(_GenHadTau->pt(*matchedgentauh_it), "MatchedVisTauHPt");
        histAddVal(_GenHadTau->eta(*matchedgentauh_it), "MatchedVisTauHEta");
        histAddVal(_GenHadTau->phi(*matchedgentauh_it), "MatchedVisTauHPhi");
        histAddVal(_GenHadTau->p4(*matchedgentauh_it).M(), "MatchedVisTauHMass");
        histAddVal(_GenHadTau->decayMode[*matchedgentauh_it], "MatchedVisTauHDecayMode");
    }
    histAddVal(active_part->at(CUTS::eGHadTau)->size(), "NMatchedVisTauH");

    int grbj = 0;
    TLorentzVector genVec2;
    for(auto it : *active_part->at(CUTS::eGBJet)){
    histAddVal(_Gen->pt(it), "BJPt");
    histAddVal(_Gen->eta(it), "BJEta");
    grbj = grbj + 1;
    }

    histAddVal(active_part->at(CUTS::eGTau)->size(), "NTau");
    histAddVal(nhadtau, "NHadTau");

    histAddVal(active_part->at(CUTS::eGBJet)->size(), "NBJ");

    // Add these to the default histograms to perform quick checks if necessary.
    histAddVal(nTruePU, "PUNTrueInt");
    histAddVal(generatorht, "HT");
    histAddVal(gen_weight, "Weight");

    for(auto it : *active_part->at(CUTS::eGZ)) {
      histAddVal(_Gen->pt(it), "ZPt");
      histAddVal(_Gen->eta(it), "ZEta");
      histAddVal(_Gen->p4(it).M(), "ZMass");
    }
    histAddVal(active_part->at(CUTS::eGZ)->size(), "NZ");

    for(auto it : *active_part->at(CUTS::eGW)) {
      histAddVal(_Gen->pt(it), "WPt");
      histAddVal(_Gen->eta(it), "WEta");
      histAddVal(_Gen->p4(it).M(), "WMass");
    }
    histAddVal(active_part->at(CUTS::eGW)->size(), "NW");


    for(auto it : *active_part->at(CUTS::eGMuon)) {
      histAddVal(_Gen->energy(it), "MuonEnergy");
      histAddVal(_Gen->pt(it), "MuonPt");
      histAddVal(_Gen->eta(it), "MuonEta");
      histAddVal(_Gen->phi(it), "MuonPhi");
    }
    histAddVal(active_part->at(CUTS::eGMuon)->size(), "NMuon");

    histAddVal(gendilepmass, "ZDiLepMass");

    double mass=0;
    TLorentzVector lep1;
    TLorentzVector lep2;
    for(size_t igen=0; igen<_Gen->size(); igen++){
      //if a Z boson is explicitly there
      if(abs(_Gen->pdg_id[igen])==11 or abs(_Gen->pdg_id[igen])==13 or abs(_Gen->pdg_id[igen])==15){
        if(lep1!=TLorentzVector(0,0,0,0)){
          lep2= _Gen->p4(igen);
          mass=(lep1+lep2).M();
          //cout<<"mass  leptons "<<mass<<std::endl;
          break;
        }else{
          //cout<<_Gen->size()<<"   "<<igen<<std::endl;
          //if(_Gen->size()>_Gen->cur_P.size()){
           //_Gen->init();
          //}
          lep1= _Gen->RecoP4(igen);
        }
      }
    }
    histAddVal(mass, "LeptonMass");
  } else if(fillInfo[group]->type == FILLER::Single) {
    Particle* part = fillInfo[group]->part;
    CUTS ePos = fillInfo[group]->ePos;

    for(auto it : *active_part->at(ePos)) {
      histAddVal(part->p4(it).Energy(), "Energy");
      histAddVal(part->p4(it).Pt(), "Pt");
      histAddVal(part->p4(it).Eta(), "Eta");
      histAddVal(part->p4(it).Phi(), "Phi");
      histAddVal(part->p4(it).DeltaPhi(_MET->p4()), "MetDphi");
      if(part->type == PType::Tau) {
        //if(_Tau->nProngs->at(it) == 1){
          //histAddVal(part->pt(it), "Pt_1prong");
        //}else if(_Tau->nProngs->at(it) == 3){
          //histAddVal(part->pt(it), "Pt_3prong");
        //}
        //histAddVal(_Tau->nProngs->at(it), "NumSignalTracks");
        histAddVal(_Tau->charge(it), "Charge");
        histAddVal(_Tau->againstElectron[it], "againstElectron");
        histAddVal(_Tau->againstMuon[it], "againstMuon");
        //histAddVal(_Tau->DecayMode[it], "DecayMode"); // original
        histAddVal(_Tau->DecayModeOldDMs[it], "DecayMode");
        histAddVal(_Tau->DecayModeNewDMs[it], "DecayModeNewDMs");
        // histAddVal(_Tau->MVAoldDM[it], "MVAoldDM"); //original
        histAddVal(_Tau->TauIdDiscr[it], "MVAoldDM");
        //histAddVal(_Tau->decayMode[it], "decayMode"); // original
        histAddVal(_Tau->decayModeInt[it], "decayMode"); 
        histAddVal(_Tau->leadTkDeltaEta[it], "leadTkDeltaEta");
        histAddVal(_Tau->leadTkDeltaPhi[it], "leadTkDeltaPhi");
        histAddVal(_Tau->leadTkPtOverTauPt[it], "leadTkPtOverTauPt");
        histAddVal(_Tau->dz[it], "dz");
        //histAddVal(_Tau->leadChargedCandPt->at(it), "SeedTrackPt");
        //histAddVal(_Tau->leadChargedCandDz_pv->at(it), "leadChargedCandDz");
      }
      if(part->type != PType::Jet) {
        histAddVal(calculateLeptonMetMt(part->p4(it)), "MetMt");
      }
      if(part->type == PType::FatJet ) {
        histAddVal(_FatJet->PrunedMass[it], "PrunedMass");
        histAddVal(_FatJet->SoftDropMass[it], "SoftDropMass");
        histAddVal(_FatJet->tau1[it], "tau1");
        histAddVal(_FatJet->tau2[it], "tau2");
        histAddVal(_FatJet->tau2[it]/_FatJet->tau1[it], "tau2Overtau1");
      }
    }

    if((part->type != PType::Jet ) && active_part->at(ePos)->size() > 0) {
      std::vector<std::pair<double, int> > ptIndexVector;
      for(auto it : *active_part->at(ePos)) {
        ptIndexVector.push_back(std::make_pair(part->pt(it),it));
      }
      sort(ptIndexVector.begin(),ptIndexVector.end());
      if(ptIndexVector.size()>0){
        histAddVal(part->pt(ptIndexVector.back().second), "FirstLeadingPt");
        histAddVal(part->eta(ptIndexVector.back().second), "FirstLeadingEta");
      }
      if(ptIndexVector.size()>1){
        histAddVal(part->pt(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingPt");
        histAddVal(part->eta(ptIndexVector.at(ptIndexVector.size()-2).second), "SecondLeadingEta");
      }
    }

    histAddVal(active_part->at(ePos)->size(), "N");


  } else if(group == "FillMetCuts") {
    histAddVal(_MET->MHT(), "MHT");
    histAddVal(_MET->HT(), "HT");
    histAddVal(_MET->HT() + _MET->MHT(), "Meff");
    histAddVal(_MET->pt(), "Met");
    histAddVal(_MET->phi(), "MetPhi");

  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() == 0) {

    if(active_part->at(CUTS::eR1stJet)->size()>1) { //01.17.19
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Pt(), "FirstPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Eta(), "FirstEta");
      Double_t dphi1new = normPhi(_Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size()-1)).Phi() - _MET->phi());
      histAddVal(dphi1new,"Dphi1");
    }
    if(active_part->at(CUTS::eR2ndJet)->size()>2) {
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size()-1)).Pt(), "SecondPt");
      histAddVal(_Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size()-1)).Eta(), "SecondEta");
    }


  } else if(group == "FillLeadingJet" && active_part->at(CUTS::eSusyCom)->size() != 0) {

    TLorentzVector first = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(active_part->at(CUTS::eR1stJet)->size() - 1));
    //cout << "first and SUS: " << first.Pt() << endl;
    TLorentzVector second = _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(active_part->at(CUTS::eR2ndJet)->size() - 1));
    //cout << "second and SUS: " << second.Pt() << endl;

    histAddVal(first.Pt(), "FirstPt");
    histAddVal(second.Pt(), "SecondPt");

    histAddVal(first.Eta(), "FirstEta");
    histAddVal(second.Eta(), "SecondEta");

    TLorentzVector LeadDiJet = first + second;

    histAddVal(LeadDiJet.M(), "Mass");
    histAddVal(LeadDiJet.Pt(), "Pt");
    histAddVal(fabs(first.Eta() - second.Eta()), "DeltaEta");
    histAddVal(first.DeltaR(second), "DeltaR");

    double dphiDijets = absnormPhi(first.Phi() - second.Phi());
    double dphi1 = normPhi(first.Phi() - _MET->phi());
    double dphi2 = normPhi(second.Phi() - _MET->phi());
    double alpha = (LeadDiJet.M() > 0) ? second.Pt() / LeadDiJet.M() : 999999999.0;

    histAddVal(dphiDijets, "LeadSublDijetDphi");
    histAddVal2(_MET->pt(),dphiDijets, "MetVsDiJetDeltaPhiLeadSubl");
    histAddVal2(fabs(first.Eta()-second.Eta()), dphiDijets, "DeltaEtaVsDeltaPhiLeadSubl");

    histAddVal(absnormPhi(_MET->phi() - LeadDiJet.Phi()), "MetDeltaPhi");


    histAddVal(sqrt( pow(dphi1,2.0) + pow((TMath::Pi() - dphi2),2.0) ), "R1");
    histAddVal(sqrt( pow(dphi2,2.0) + pow((TMath::Pi() - dphi1),2.0)), "R2");
    histAddVal(normPhi(first.Phi() - _MET->MHTphi()), "Dphi1MHT");
    histAddVal(normPhi(second.Phi() - _MET->MHTphi()), "Dphi2MHT");
    histAddVal(dphi1, "Dphi1");
    histAddVal(dphi2, "Dphi2");
    histAddVal2(dphi1,dphi2, "Dphi1VsDphi2");
    histAddVal(alpha, "Alpha");


    //dijet info
  } else if(group == "FillDiJet") {
    double leaddijetmass = 0;
    double leaddijetpt = 0;
    double leaddijetdeltaR = 0;
    double leaddijetdeltaEta = 0;
    double etaproduct = 0;
    for(auto it : *active_part->at(CUTS::eDiJet)) {
      int p1 = (it) / _Jet->size();
      int p2 = (it) % _Jet->size();
      TLorentzVector jet1 = _Jet->p4(p1);
      TLorentzVector jet2 = _Jet->p4(p2);
      TLorentzVector DiJet = jet1 + jet2;

      if(DiJet.M() > leaddijetmass) {
        leaddijetmass = DiJet.M();
        etaproduct = (jet1.Eta() * jet2.Eta() > 0) ? 1 : -1;
      }
      if(DiJet.Pt() > leaddijetpt) leaddijetpt = DiJet.Pt();
      if(fabs(jet1.Eta() - jet2.Eta()) > leaddijetdeltaEta) leaddijetdeltaEta = fabs(jet1.Eta() - jet2.Eta());
      if(jet1.DeltaR(jet2) > leaddijetdeltaR) leaddijetdeltaR = jet1.DeltaR(jet2);

      histAddVal(DiJet.M(), "Mass");
      histAddVal(DiJet.Pt(), "Pt");
      histAddVal(fabs(jet1.Eta() - jet2.Eta()), "DeltaEta");
      histAddVal(absnormPhi(jet1.Phi() - jet2.Phi()), "DeltaPhi");
      histAddVal(jet1.DeltaR(jet2), "DeltaR");
    }


    histAddVal(leaddijetmass, "LargestMass");
    histAddVal(leaddijetpt, "LargestPt");
    histAddVal(leaddijetdeltaEta, "LargestDeltaEta");
    histAddVal(leaddijetdeltaR, "LargestDeltaR");
    histAddVal(etaproduct, "LargestMassEtaProduct");

    for(auto index : *(active_part->at(CUTS::eRTau1)) ) {
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetmass, "mTvsLeadingMass");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaEta, "mTvsLeadingDeltaEta");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetdeltaR, "mTvsLeadingDeltaR");
      histAddVal2(calculateLeptonMetMt(_Tau->p4(index)), leaddijetpt, "mTvsLeadingPt");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetmass, "MetDphiVSLeadingMass");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaEta, "MetDphiVSLeadingDeltaEta");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetdeltaR, "MetDphiVSLeadingDeltaR");
      histAddVal2((absnormPhi(_Tau->p4(index).Phi()-_MET->phi())), leaddijetpt, "MetDphiVSLeadingPt");
    }



    ////diparticle stuff

  } else if(fillInfo[group]->type == FILLER::Dilepjet) {
    Jet* jet = static_cast<Jet*>(fillInfo[group]->part);
    Lepton* lep = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    std::string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1;
    TLorentzVector part2;

    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / _Jet->size();;
      int p2= (it) % _Jet->size();;

      part1 = lep->p4(p1);
      part2 = jet->p4(p2);

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != std::string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");
      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }
      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt");
      histAddVal(PZeta, "PZeta");
      histAddVal(PZetaVis, "PZetaVis");
      histAddVal2(PZetaVis,PZeta, "Zeta2D");
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

        histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
        histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
        histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass");
      }
    }
  } else if(fillInfo[group]->type == FILLER::Dipart) {
    Lepton* lep1 = static_cast<Lepton*>(fillInfo[group]->part);
    Lepton* lep2 = static_cast<Lepton*>(fillInfo[group]->part2);
    CUTS ePos = fillInfo[group]->ePos;
    std::string digroup = group;
    digroup.erase(0,4);

    TLorentzVector part1;
    TLorentzVector part2;
    TLorentzVector llep;

    for(auto it : *active_part->at(ePos)) {

      int p1= (it) / BIG_NUM;
      int p2= (it) % BIG_NUM;

      part1 = lep1->p4(p1);
      part2 = lep2->p4(p2);

      if(part1.Pt() > part2.Pt()) llep = lep1->p4(p1);
      else if(part1.Pt() < part2.Pt()) llep = lep2->p4(p2);

      histAddVal2(part1.Pt(),part2.Pt(), "Part1PtVsPart2Pt");
      histAddVal(part1.DeltaR(part2), "DeltaR");
      if(group.find("Di") != std::string::npos) {
        histAddVal((part1.Pt() - part2.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part1.Pt() - part2.Pt(), "DeltaPt");
      } else {
        histAddVal((part2.Pt() - part1.Pt()) / (part1.Pt() + part2.Pt()), "DeltaPtDivSumPt");
        histAddVal(part2.Pt() - part1.Pt(), "DeltaPt");
      }
      histAddVal(cos(absnormPhi(part2.Phi() - part1.Phi())), "CosDphi");

      histAddVal(cos(absnormPhi(part1.Phi() - _MET->phi())), "Part1CosDphiPtandMet");
      histAddVal(cos(absnormPhi(part2.Phi() - _MET->phi())), "Part2CosDphiPtandMet");


      histAddVal(absnormPhi(part1.Phi() - _MET->phi()), "Part1MetDeltaPhi");
      histAddVal2(absnormPhi(part1.Phi() - _MET->phi()), cos(absnormPhi(part2.Phi() - part1.Phi())), "Part1MetDeltaPhiVsCosDphi");
      histAddVal(absnormPhi(part2.Phi() - _MET->phi()), "Part2MetDeltaPhi");
      histAddVal(cos(absnormPhi(atan2(part1.Py() - part2.Py(), part1.Px() - part2.Px()) - _MET->phi())), "CosDphi_DeltaPtAndMet");

      double diMass = diParticleMass(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"));
      if(passDiParticleApprox(part1,part2, distats[digroup].smap.at("HowCalculateMassReco"))) {
        histAddVal(diMass, "ReconstructableMass");
      } else {
        histAddVal(diMass, "NotReconstructableMass");
      }

      double diMass1 = CalculateDiLepMassDeltaPt(part1, part2);
      histAddVal(diMass1, "ReconstructableMassDeltaPt");
	    
      double InvMass = diParticleMass(part1,part2, "InvariantMass");
      histAddVal(InvMass, "InvariantMass");

      double AbscosDphiLeadLepMet = cos(absnormPhi(llep.Phi() - _MET->phi()));                                                                                                                                    
      histAddVal(AbscosDphiLeadLepMet, "RecoAbsCosDphiPtLeadLepandMet"); 

      double ptSum = part1.Pt() + part2.Pt();
      histAddVal(ptSum, "SumOfPt");

      double PZeta = getPZeta(part1,part2).first;
      double PZetaVis = getPZeta(part1,part2).second;
      histAddVal(calculateLeptonMetMt(part1), "Part1MetMt");
      histAddVal(calculateLeptonMetMt(part2), "Part2MetMt");
      histAddVal(lep2->charge(p2) * lep1->charge(p1), "OSLS");
      histAddVal(PZeta, "PZeta");
      histAddVal(PZetaVis, "PZetaVis");
      histAddVal2(PZetaVis,PZeta, "Zeta2D");
      histAddVal((distats.at(digroup).dmap.at("PZetaCutCoefficient") * PZeta) + (distats.at(digroup).dmap.at("PZetaVisCutCoefficient") * PZetaVis), "Zeta1D");

      if ((active_part->at(CUTS::eR1stJet)->size()>0 && active_part->at(CUTS::eR1stJet)->at(0) != -1) && (active_part->at(CUTS::eR2ndJet)->size()>0 && active_part->at(CUTS::eR2ndJet)->at(0) != -1)) {
        TLorentzVector TheLeadDiJetVect = _Jet->p4(active_part->at(CUTS::eR1stJet)->at(0)) + _Jet->p4(active_part->at(CUTS::eR2ndJet)->at(0));

        histAddVal(absnormPhi(part1.Phi() - TheLeadDiJetVect.Phi()), "Part1DiJetDeltaPhi");
        histAddVal(absnormPhi(part2.Phi() - TheLeadDiJetVect.Phi()), "Part2DiJetDeltaPhi");
        histAddVal(diParticleMass(TheLeadDiJetVect, part1+part2, "VectorSumOfVisProductsAndMet"), "DiJetReconstructableMass");
      }

      if(lep1->type != PType::Tau) {
        histAddVal(isZdecay(part1, *lep1), "Part1IsZdecay");
      }
      if(lep2->type != PType::Tau){
        histAddVal(isZdecay(part2, *lep2), "Part2IsZdecay");
      }


      //electron tau stuff:
      if(lep1->type == PType::Electron && lep2->type == PType::Electron){
        //loop over taus to find a match in the unisolated taus:
        int matchedTauInd=-1;
        TLorentzVector matchedEle;
        TLorentzVector unmatchedEle;
        for( size_t itau =0; itau< _Tau->size(); itau++){
          if(part2.DeltaR(_Tau->p4(itau))<0.3){
            //we are sure that part1 passes the tight id
            matchedTauInd=itau;
            matchedEle=part2;
            unmatchedEle=part1;
          }else if(part1.DeltaR(_Tau->p4(itau))<0.3){
            //check if part2 passes the tight id:
            if(find(active_part->at(CUTS::eRElec1)->begin(),active_part->at(CUTS::eRElec1)->end(),p2)!=active_part->at(CUTS::eRElec1)->end()){
              matchedTauInd=itau;
              matchedEle=part1;
              unmatchedEle=part2;
            }
          }
        }
        if(matchedTauInd>=0){
          if(find(active_part->at(CUTS::eRTau1)->begin(),active_part->at(CUTS::eRTau1)->end(),matchedTauInd)!=active_part->at(CUTS::eRTau1)->end()){
            histAddVal(_Tau->p4(matchedTauInd).Pt(), "DiEleGoodTauMatchPt");
            histAddVal(_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt(), "DiEleGoodTauMatchDeltaPt");
            histAddVal((_Tau->p4(matchedTauInd)+unmatchedEle).M(), "DiEleGoodTauMatchMass");
            histAddVal((matchedEle+unmatchedEle).M(), "DiEleEleGoodMatchMass");
            histAddVal(matchedEle.Pt(), "DiEleEleGoodMatchPt");
            //histAddVal(_Tau->leadChargedCandPtError->at(matchedTauInd),"DiEleleadChargedCandPtErrorGoodMatched");
            //histAddVal(_Tau->leadChargedCandValidHits->at(matchedTauInd),"DiEleleadChargedCandValidHitGoodMatched");
            histAddVal2( matchedEle.Pt(),   (_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt())/matchedEle.Pt(), "DiEleTauGoodMatchPt_vs_DeltaPt");
            histAddVal2( matchedEle.Pt(),   matchedEle.Eta(), "DiEleTauGoodMatchPt_vs_eta");
            // histAddVal2( _Tau->pt(matchedTauInd),   _Tau->decayMode[matchedTauInd], "DiEleTauGoodMatchPt_vs_Decay"); //original
            histAddVal2( _Tau->pt(matchedTauInd),   _Tau->decayModeInt[matchedTauInd], "DiEleTauGoodMatchPt_vs_Decay");
          }else{
            histAddVal(_Tau->p4(matchedTauInd).Pt(), "DiEleTauMatchPt");
            histAddVal(_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt(), "DiEleTauMatchDeltaPt");
            histAddVal((_Tau->p4(matchedTauInd)+unmatchedEle).M(), "DiEleTauMatchMass");
            histAddVal((matchedEle+unmatchedEle).M(), "DiEleEleMatchMass");
            histAddVal(matchedEle.Pt(), "DiEleEleMatchPt");
            //histAddVal(_Tau->leadChargedCandPtError->at(matchedTauInd),"DiEleleadChargedCandPtErrorMatched");
            //histAddVal(_Tau->leadChargedCandValidHits->at(matchedTauInd),"DiEleleadChargedCandValidHitsMatched");
            histAddVal2( matchedEle.Pt(),   (_Tau->p4(matchedTauInd).Pt()-matchedEle.Pt())/matchedEle.Pt(), "DiEleTauMatchPt_vs_DeltaPt");
            histAddVal2( matchedEle.Pt(),   matchedEle.Eta(), "DiEleTauMatchPt_vs_eta");
            // histAddVal2( _Tau->pt(matchedTauInd),   _Tau->decayMode[matchedTauInd], "DiEleTauMatchPt_vs_Decay"); //original
            histAddVal2( _Tau->pt(matchedTauInd),   _Tau->decayModeInt[matchedTauInd], "DiEleTauMatchPt_vs_Decay");
          }
        }else{
          histAddVal((part1+part2).M(), "DiEleEleUnMatchMass");
          histAddVal(part2.Pt(), "DiEleEleUnMatchPt");
          histAddVal2( part2.Pt(),   part2.Eta(), "DiEleUnMatchPt_vs_eta");
          if(!isData){
            histAddVal(part2.Pt(), "DiEleEleUnMatchPt_gen_"+std::to_string(abs(matchToGenPdg(part2,0.3))));
          }
          int found=-1;
          for(size_t i=0; i< _Jet->size(); i++) {
            if(part2.DeltaR(_Jet->p4(i)) <=0.4) {
              found=i;
            }
          }
          if (found>=0){
            //histAddVal(_Jet->chargedMultiplicity[found], "DiEleEleUnMatchJetMultiplicity");
          }else{
            histAddVal(-1, "DiEleEleUnMatchJetMultiplicity");
          }

        }
      }
    }
  }
}

void Analyzer::fill_Tree(){

  if(0){
    //do our dirty tree stuff here:
    int p1=-1;
    int p2=-1;
    if(active_part->at(CUTS::eDiTau)->size()==1){
      p1= active_part->at(CUTS::eDiTau)->at(0) / BIG_NUM;
      p2= active_part->at(CUTS::eDiTau)->at(0) % BIG_NUM;
    } else{
      return;
    }
    int j1=-1;
    int j2=-1;
    double mass=0;
    for(auto it : *active_part->at(CUTS::eDiJet)) {
      int j1tmp= (it) / _Jet->size();
      int j2tmp= (it) % _Jet->size();
      if(diParticleMass(_Jet->p4(j1tmp),_Jet->p4(j2tmp),"")>mass){
        j1=j1tmp;
        j2=j2tmp;
        mass=diParticleMass(_Jet->p4(j1tmp),_Jet->p4(j2tmp),"");
      }
    }
    if(p1<0 or p2<0 or j1<0 or j2 <0)
      return;
    zBoostTree["tau1_pt"]   = _Tau->pt(p1);
    zBoostTree["tau1_eta"]  = _Tau->eta(p1);
    zBoostTree["tau1_phi"]  = _Tau->phi(p1);
    zBoostTree["tau2_pt"]   = _Tau->pt(p2);
    zBoostTree["tau2_eta"]  = _Tau->eta(p2);
    zBoostTree["tau2_phi"]  = _Tau->phi(p2);
    zBoostTree["tau_mass"]  = diParticleMass(_Tau->p4(p1),_Tau->p4(p2),"");
    zBoostTree["met"]       = _MET->pt();
    zBoostTree["mt_tau1"]   = calculateLeptonMetMt(_Tau->p4(p1));
    zBoostTree["mt_tau2"]   = calculateLeptonMetMt(_Tau->p4(p2));
    zBoostTree["mt2"]       = _MET->MT2(_Tau->p4(p1),_Tau->p4(p2));
    zBoostTree["cosDphi1"]  = absnormPhi(_Tau->phi(p1) - _MET->phi());
    zBoostTree["cosDphi2"]  = absnormPhi(_Tau->phi(p2) - _MET->phi());
    zBoostTree["jet1_pt"]   = _Jet->pt(j1);
    zBoostTree["jet1_eta"]  = _Jet->eta(j1);
    zBoostTree["jet1_phi"]  = _Jet->phi(j1);
    zBoostTree["jet2_pt"]   = _Jet->pt(j2);
    zBoostTree["jet2_eta"]  = _Jet->eta(j2);
    zBoostTree["jet2_phi"]  = _Jet->phi(j2);
    zBoostTree["jet_mass"]  = mass;
    zBoostTree["weight"]    = wgt;

    //put it accidentally in the tree
    histo.fillTree("TauTauTree");
  }
}

void Analyzer::initializePileupInfo(const bool& specialPU, std::string outfilename){
   // If specialPUcalculation is true, then take the name of the output file (when submitting jobs) 
   // to retrieve the right MC nTruePU distribution to calculate the PU weights, otherwise, it will do
   // the calculation as usual, having a single MC nTruePU histogram (MCHistos option).
   // If you need to run the Analyzer interactively and use the specialPUcalculation option, make sure that you 
   // include the name of the sample to analyze. For example, you can look at the names given to output files in past
   // runs sumbitted as jobs and use those names instead. This will be improved later.

   // try-catch block for initialize pileup info
   try{

     if(!specialPU){// No special PU calculation.
       initializePileupWeights(distats["Run"].smap.at("MCHistos"),distats["Run"].smap.at("DataHistos"),distats["Run"].smap.at("DataPUHistName"),distats["Run"].smap.at("MCPUHistName"));
     }
     else{ // Special PU calculation
    std::string outputname = outfilename;

    std::string delimitertune = "_Tune";
    std::string delimiterenergy = "_13TeV";

    bool istuneinname = outputname.find(delimitertune.c_str()) != std::string::npos;

    std::string samplename;

    if((samplename.length() == 0) && istuneinname){
      unsigned int pos = outputname.find(delimitertune.c_str());
      samplename = outputname.erase(pos, (outputname.substr(pos).length()));
    }
    else if((samplename.length() == 0) && (!istuneinname)){
      unsigned int pos = outputname.find(delimiterenergy.c_str());
      samplename = outputname.erase(pos, (outputname.substr(pos).length()));
    }

         initializePileupWeights(distats["Run"].smap.at("SpecialMCPUHistos"),distats["Run"].smap.at("DataHistos"),distats["Run"].smap.at("DataPUHistName"), samplename);
     }
   }// end of try block
   catch(std::runtime_error& err){
    std::cerr << "ERROR in initializePileupInfo! " << std::endl;
    std::cerr << "\t" << err.what()  << std::endl;
    std::cout << "\tMake sure you are using the right file and the pileup histogram name is correct in Run_info.in." << std::endl;
    if(specialPUcalculation){
      std::cerr << "\tWARNING: You are using SpecialPUCalculation. Check that the output filename (" << outfilename << ") contains the name of the analyzed MC sample." << std::endl;
    }
    std::cerr << "\tAborting Analyzer..." << std::endl;
    std::abort();
   }
    catch(std::out_of_range& err){
    std::cerr << "ERROR in initializePileupInfo! " << std::endl;
    std::cerr << "\t" << err.what()  << std::endl;
    if(specialPU){
      std::cerr << "\tYou are using SpecialPUCalculation. Name of the MC sample was not found in the output filename (" << outfilename << ")." << std::endl;
      std::cerr << "\tTry with a different name that contains the name of your MC sample followed by _Tune or _13TeV (e.g. DYJetsToLL_M-50_TuneCP5.root)." << std::endl;
    }
    std::cerr << "\tAborting Analyzer..." << std::endl;
    std::abort();
   }
   catch(...){
    std::cerr << "ERROR in initializePileupInfo! Unknown exception." << std::endl; 
    std::cerr << "\tAborting Analyzer..." << std::endl;
    std::abort();
   } // end of catch blocks

 }


void Analyzer::initializePileupWeights(std::string MCHisto, std::string DataHisto, std::string DataHistoName, std::string MCHistoName) {

  TFile *file1 = new TFile((PUSPACE+MCHisto).c_str());
  TH1D* histmc = (TH1D*)file1->FindObjectAny(MCHistoName.c_str());
  if(!histmc) throw std::runtime_error(("Failed to extract histogram "+MCHistoName+" from "+PUSPACE+MCHisto+"!").c_str());

  TFile* file2 = new TFile((PUSPACE+DataHisto).c_str());
  TH1D* histdata = (TH1D*)file2->FindObjectAny(DataHistoName.c_str());
  if(!histdata) throw std::runtime_error(("Failed to extract histogram "+DataHistoName+" from "+PUSPACE+DataHisto+"!").c_str());

  TH1D* histdata_up = (TH1D*)file2->FindObjectAny((DataHistoName+"Up").c_str());
  TH1D* histdata_down = (TH1D*)file2->FindObjectAny((DataHistoName+"Down").c_str());


  histmc->Scale(1./histmc->Integral());
  histdata->Scale(1./histdata->Integral());
  if(histdata_up){
    histdata_up->Scale(1./histdata_up->Integral());
    histdata_down->Scale(1./histdata_down->Integral());
  }

  //double factor = histmc->Integral() / histdata->Integral();
  double value,valueUp,valueDown;
  // The bin corresponding to nTruePU = 0 will be bin = 1 and, the calculated weight will be stored 
  // in hPU[bin]. That's why when we calculate the pu_weight in setupEventGeneral, we require hPU[(int)nTruePU+1].
  for(int bin=0; bin <= histmc->GetNbinsX(); bin++) {
    if(histmc->GetBinContent(bin) == 0){
      value = 1;
      valueUp = 1;
      valueDown = 1;
    }else{
      value = histdata->GetBinContent(bin) / histmc->GetBinContent(bin);
      if(histdata_up){
        valueUp = histdata->GetBinContent(bin) / histmc->GetBinContent(bin);
        valueDown = histdata->GetBinContent(bin) / histmc->GetBinContent(bin);
      }
    }
    hPU[bin]      = value;
    if(histdata_up){
      hPU_up[bin]   = valueUp;
      hPU_down[bin] = valueDown;
    }else{
      hPU_up[bin]   = value;
      hPU_down[bin] = value;
    }
  }

  file1->Close();
  file2->Close();

}

void Analyzer::initializeWkfactor(std::vector<std::string> infiles) {
  if(infiles[0].find("WJets") != std::string::npos){
    isWSample = true;
  }else{
    isWSample=false;
    return;
  }
  //W-jet k-factor Histograms:
  TFile k_ele("Pileup/k_faktors_ele.root");
  TFile k_mu("Pileup/k_faktors_mu.root");
  TFile k_tau("Pileup/k_faktors_tau.root");

  k_ele_h =dynamic_cast<TH1D*>(k_ele.FindObjectAny("k_fac_m"));
  k_mu_h  =dynamic_cast<TH1D*>(k_mu.FindObjectAny("k_fac_m"));
  k_tau_h =dynamic_cast<TH1D*>(k_tau.FindObjectAny("k_fac_m"));

  k_ele.Close();
  k_mu.Close();
  k_tau.Close();

}

///Normalizes phi to be between -PI and PI
double normPhi(double phi) {
  static double const TWO_PI = TMath::Pi() * 2;
  while ( phi <= -TMath::Pi() ) phi += TWO_PI;
  while ( phi >  TMath::Pi() ) phi -= TWO_PI;
  return phi;
}


///Takes the absolute value of of normPhi (made because constant use)
double absnormPhi(double phi) {
  return abs(normPhi(phi));
}