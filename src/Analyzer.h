#ifndef Analyzer_h
#define Analyzer_h

struct CRTester;

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>

#include <TDirectory.h>
#include <TEnv.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1.h>

#include "Particle.h"
#include "MET.h"
#include "Histo.h"

/////fix
#include "./btagging/BTagCalibrationStandalone.h"
#include "TauIDSFTool.h"

#include "Cut_enum.h"
#include "FillInfo.h"
#include "CRTest.h"
#include "Systematics.h"
#include "DepGraph.h"
#include "JetScaleResolution.h"
#include "JetRecalibrator.h"
#include "L1ECALPrefiringWgtProd.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"

double normPhi(double phi);
double absnormPhi(double phi);

//#define const
//using namespace std;

template <typename T>
void removeDuplicates(std::vector<T>& vec)
{
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

static const int nTrigReq = 2;

class Analyzer {
  friend class CRTester;
public:
  Analyzer(std::vector<std::string>, std::string, bool setCR = false, std::string configFolder="PartDet", std::string year="2016");
  ~Analyzer();
  void add_metadata(std::vector<std::string> infiles);
  void clear_values();
  void preprocess(int, std::string); //05.28.19
  bool fillCuts(bool);
  void printCuts();
  void writeout();
  int nentries;
  void fill_efficiency();
  void fill_histogram(std::string);
  // void fill_Tree();
  void setControlRegions() { histo.setControlRegions();}
  void checkParticleDecayList(); //01.16.19
  /*----------- ReadinJSON class variables -------------*/
  std::multimap<int,int> readinJSON(std::string);
  std::multimap<int,int> jsonlinedict;
  /*----------------------------------------------------*/
  void writeParticleDecayList(int); //01.16.19
  void getGoodGenBJet(); //01.16.19
  std::vector<int>* getList(CUTS ePos) {return goodParts[ePos];}
  double getMet() {return _MET->pt();}
  double getHT() {return _MET->HT();}
  double getMHT() {return _MET->MHT();}
  double getMass(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string partName) {
    return diParticleMass(Tobj1, Tobj2, distats[partName].smap.at("HowCalculateMassReco"));
  }
  double getZeta(const TLorentzVector& Tobj1, const TLorentzVector& Tobj2, std::string partName) {
    return distats[partName].dmap.at("PZetaCutCoefficient") * getPZeta(Tobj1, Tobj2).first;

  }
// private:
  std::string runyear;
  void CRfillCuts();
  ///// Functions /////
  //void fill_Folder(std::string, const int, std::string syst="");
  void fill_Folder(std::string, const int, Histogramer& ihisto, bool issyst);

  void getInputs();
  void setupJob(std::string);
  bool specialPUcalculation = false;
  void initializePileupInfo(const bool&, std::string);
  void initializePileupWeights(std::string, std::string, std::string, std::string);
  void initializeMCSelection(std::vector<std::string> infiles);
  void initializeWkfactor(std::vector<std::string> infiles);

  void read_info(std::string);
  void setupGeneral(std::string);
  void setupEventGeneral(int);
  void getTriggerBranchesList(CUTS, std::string, bool);
  bool passGenHTFilter(float);
  bool passGenMassFilterZ(float mass_lowbound, float mass_upbound);
  bool checkGoodRunsAndLumis(int);
  void branchException(std::string);
  void initializeTrigger();
  void setCutNeeds();

  bool passHEMveto2018();
  bool passJetVetoEEnoise2017(int);

  bool skimSignalMC(int);
  std::string inputSignalModel="", inputSignalMassParam="";
  bool finalInputSignal = false;

  void smearLepton(Lepton&, CUTS, const PartStats&, const PartStats&, int syst=0);
  void smearTaus(Lepton&, const PartStats&, const PartStats&, int syst=0);
  //void smearJet(Particle&, CUTS, const PartStats&, int syst=0);
  void smearJetRes(Particle&, CUTS, const PartStats&, int syst=0);

  bool JetMatchesLepton(const Lepton&, const TLorentzVector&, double, CUTS);
  TLorentzVector matchLeptonToGen(const TLorentzVector&, const PartStats&, CUTS);
  // TLorentzVector matchLeptonToGen(int, const PartStats& stats, CUTS ePos);
  TLorentzVector matchTauToGen(const TLorentzVector&, double);
  TLorentzVector matchHadTauToGen(const TLorentzVector&, double);
  TLorentzVector matchJetToGen(const TLorentzVector&, const double&, CUTS, bool);

  int matchToGenPdg(const TLorentzVector& lvec, double minDR);


  void getGoodParticles(int);
  void getGoodGenHadronicTaus(const PartStats&);
  void getGoodGenJets(const PartStats&);
  void getGoodGenBJets(const PartStats&);
  void getGoodGenHadronicTauNeutrinos(const PartStats&);
  TLorentzVector getGenVisibleTau4Vector(int, int);
  void getGoodGen(const PartStats&);
  void getGoodRecoLeptons(const Lepton&, const CUTS, const CUTS, const PartStats&, const int);
  void getGoodRecoJets(CUTS, const PartStats&, const int);
  void getGoodRecoLeadJets(CUTS, const PartStats&, const int);
  void getGoodRecoBJets(CUTS, const CUTS, const PartStats&, const int); //01.16.19
  void getGoodRecoFatJets(CUTS, const PartStats&, const int);

  bool discrBySUSYSoftElectrons(int, bool, bool);
  bool discrBySUSYSoftMuons(int, bool, bool);

  void getGoodLeptonCombos(Lepton&, Lepton&, CUTS, CUTS, CUTS, const PartStats&, const int);
  double CalculateDiLepMassDeltaPt(const TLorentzVector&, const TLorentzVector&, const float, const float);
  void getGoodLeptonJetCombos(Lepton&, Jet&, CUTS, CUTS, CUTS, const PartStats&, const int);
  void getGoodDiJets(const PartStats&, const int);

  void VBFTopologyCut(const PartStats&, const int);
  void fastTriggerCuts(CUTS);
  void TriggerCuts(CUTS);

  inline bool passCutRangeAbs(std::string, double, const PartStats&);
  bool passCutRangeAbs(double, const std::pair<double, double>&);

  double calculateLeptonMetMt(const TLorentzVector&);
  double diParticleMass(const TLorentzVector&, const TLorentzVector&, std::string);
  bool passDiParticleApprox(const TLorentzVector&, const TLorentzVector&, std::string);
  bool isZdecay(const TLorentzVector&, const Lepton&);

  bool isOverlaping(const TLorentzVector&, Lepton&, CUTS, double);
  bool isOverlapingB(const TLorentzVector&, Jet&, CUTS, double); //01.17.19
  bool isOverlapingC(const TLorentzVector&, FatJet&, CUTS, double);
  bool passProng(std::string, int);
  bool isInTheCracks(float);
  bool passedLooseJetID(int);
  bool select_mc_background();
  double getTauIdSFs(bool, bool, bool, bool, std::string);
  double getWkfactor();
  double getZBoostWeight();
  double getZpTWeight();
  double getZpTWeight_vbfSusy(std::string);
  double getZBoostWeightSyst(int ud); // 06.02.20
  double getTopBoostWeight(); //01.15.19
  void setupBJetSFInfo(const PartStats&, std::string); // new function that sets up the b-tagging SF info
  double getBJetSF(CUTS, const PartStats&); //01.16.19
  double getBJetSFResUp(CUTS, const PartStats&); //01.16.19
  double getBJetSFResDown(CUTS, const PartStats&); //01.16.19
  std::pair<double, double> getPZeta(const TLorentzVector&, const TLorentzVector&);
  void create_fillInfo();

  void setupJetCorrections(std::string, std::string);
  // void getJetEnergyResSFs(Particle& jet, const CUTS eGenPos, bool);
  void applyJetEnergyCorrections(Particle&, CUTS, const PartStats&, std::string, int syst=0);

  void setupTauIDSFsInfo(std::string, std::string, bool, bool);
  void setupTauResSFsInfo(bool);

  inline bool passCutRange(std::string, double, const PartStats&);
  bool passCutRange(double, const std::pair<double, double>&);
  bool findCut(const std::vector<std::string>&, std::string);

  void updateMet(int syst=0);
  void selectMet(int syst=0);
  bool passMetFilters(std::string, int);
  void treatMuonsAsMet(int);
  double getPileupWeight(float);
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> getArray();

  double getCRVal(std::string);
  void setupCR(std::string, double);

  ///// values /////

  TChain* BOOM;
  TTree* BAAM;
  TFile* infoFile;
  TFile* routfile = new TFile;
  std::string filespace = "";
  float hPU[200] = { };        // initialize this array to zero.
  float hPU_up[200] = { };     // initialize this array to zero.
  float hPU_down[200] = { };   // initialize this array to zero.
  int version=0;
  std::map<std::string,TTree* > originalTrees;
  //std::map<std::string,*TObject> otherObjects;


  Generated* _Gen;
  GenHadronicTaus* _GenHadTau;
  GenJets* _GenJet;
  Electron* _Electron;
  Muon* _Muon;
  Taus* _Tau;
  Jet* _Jet;
  FatJet* _FatJet;
  Met* _MET;
  Photon* _Photon;
  Histogramer histo;
  Histogramer syst_histo;
  std::unordered_map<CUTS, std::vector<int>*, EnumHash>* active_part;
  static const std::unordered_map<std::string, CUTS> cut_num;

  Systematics systematics;
  JetRecalibrator jetRecalib;
  JetRecalibrator jetRecalibL1;
  JetScaleResolution jetScaleRes;
  PartStats genStat;

  std::unordered_map<std::string, PartStats> distats;
  std::unordered_map<std::string, FillVals*> fillInfo;
  std::unordered_map<std::string, double> genMap;
  std::unordered_map<CUTS, std::vector<int>*, EnumHash> goodParts;
  std::vector<std::unordered_map<CUTS, std::vector<int>*, EnumHash>> syst_parts;
  std::unordered_map<CUTS, bool, EnumHash> need_cut;

  std::unordered_map<std::string,bool> gen_selection;
  std::regex genName_regex;

  TH1D* k_ele_h;
  TH1D* k_mu_h;
  TH1D* k_tau_h;

  bool isVSample;
  bool isZsample;
  bool isWSample;
  double boosters[3] = { }; //06.02.20

  std::string btagalgoname;

  std::vector<Particle*> allParticles;
  std::vector<std::string> syst_names;
  std::map<CUTS, Particle* >  particleCutMap;
  DepGraph neededCuts;

  static const std::unordered_map<CUTS, std::vector<CUTS>, EnumHash> adjList;

  std::vector<int>* trigPlace[nTrigReq];
  bool setTrigger = false;
  std::vector<std::string> trigger1BranchesList, trigger2BranchesList;
  bool triggerDecision = false;
  std::vector<std::string> inputTrigger1Names, inputTrigger2Names; // Brenda: This will take the triggers from the configuration file Run_info.in
  std::vector<bool> trigger1namedecisions, trigger2namedecisions; // Brenda
  std::vector<int> cuts_per, cuts_cumul;

  std::vector<std::pair<double, int> > jetPtIndexVector;

  // std::unordered_map< std::string,float > zBoostTree;

  double maxIso, minIso;
  int leadIndex, maxCut, crbins=1;
  bool isData, isSignalMC, CalculatePUSystematics, doSystematics;

  float nTruePU = 0;
  int bestVertices = 0;
  float gen_weight = 0;
  float generatorht = 0;
  float gendilepmass = 0;

  // Met filters' variables
  bool applymetfilters = false;
  bool primaryvertexfilter = false;
  bool beamhalofilter = false;
  bool hbhenoisefilter = false;
  bool ecaltpfilter = false;
  bool badpfmuonfilter = false;
  bool badchargedhadronfilter = false;
  bool ecalbadcalibrationfilter = false;
  bool allmetfilters = false;
  bool passedmetfilters = false;

  // B-tagging scale factors - calibration + readers
  BTagCalibration btagcalib;
  BTagCalibrationReader btagsfreader;
  BTagEntry::JetFlavor bjetflavor;
  BTagEntry::OperatingPoint b_workingpoint;

  // Tau ID SFs
  TauIDSFTool tau1idSFs, tau2idSFs;
  TauIDSFTool tau1id_antiEleSFs, tau2id_antiEleSFs;
  TauIDSFTool tau1id_antiMuSFs, tau2id_antiMuSFs;
  TauESTool tauesSFs;
  TauFESTool taufesSFs;
  std::string tauidyear;
  std::map<int, std::string> tauidwpsmap;
  std::map<int, std::string> antielewpsmap;
  std::map<int, std::string> antimuwpsmap;
  std::string tauid_algo, antiele_algo, antimu_algo;
  int tauidwp, antielewp, antimuwp;
  bool failtau1iso = false, failtau2iso = false;

  // Prefiring weights
  // L1ECALPrefiringWgtProd prefiringwgtprod; // commented by Brenda FE, Dec 2020
  float l1prefiringwgt = 1.0, l1prefiringwgt_up = 1.0, l1prefiringwgt_dn = 1.0;

  Float_t jec_rho =20.;
  std::string runera;
  // std::vector< std::vector<float> > jets_jer_sfs;
  std::map<int, std::vector<float> > jets_jer_sfs;
  std::map<int, float> jets_jetptres;
  std::map<int, float> jets_datatoMCSFs_jer;
  std::vector<int> genMatchedJets;
  std::vector<int> genUnmatchedJets;


  const static std::vector<CUTS> genCuts;
  const static std::vector<CUTS> jetCuts;
  const static std::vector<CUTS> nonParticleCuts;
  double pu_weight, wgt, backup_wgt;
  std::unordered_map<int, GenFill*> genMaper;

  std::vector<CRTester*> testVec;
  int SignalRegion = -1;
  bool blinded = true;
  clock_t start_time;
  std::chrono::time_point<std::chrono::system_clock> start;

};



#endif
