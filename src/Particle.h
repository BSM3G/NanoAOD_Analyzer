#ifndef Particle_h
#define Particle_h

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <regex>

#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>

#include "tokenizer.hpp"
#include "Cut_enum.h"

//using namespace std;
typedef unsigned int uint;

struct PartStats {
  std::unordered_map<std::string,double> dmap;
  std::unordered_map<std::string,std::string> smap;
  std::unordered_map<std::string,std::pair<double,double> > pmap;
  std::unordered_map<std::string,std::vector<int> > vmap;
  //  std::unordered_map<std::string,bool> bmap;
  std::vector<std::string> bset;

  bool bfind(std::string cut) const {
    return find(bset.begin(), bset.end(), cut) != bset.end();
  }

};


enum class PType { Electron, Muon, Tau, Jet, FatJet, Photon, None};

const Int_t MAXINDEX = 1000;

class Particle {

public:
  Particle();
  Particle(TTree*, std::string, std::string, std::vector<std::string>);
  virtual ~Particle() {}

  virtual std::vector<CUTS> findExtraCuts() {return std::vector<CUTS>();}
  void init();
  void unBranch();
  double pt(uint) const;
  double eta(uint) const;
  double phi(uint) const;
  double energy(uint) const;
  double mass(uint) const;
  virtual double charge(uint) const;
  TLorentzVector p4(uint) const;
  TLorentzVector& p4(uint);
  TLorentzVector RecoP4(uint) const;
  TLorentzVector& RecoP4(uint);

  uint size() const;
  std::vector<TLorentzVector>::iterator begin();
  std::vector<TLorentzVector>::iterator end();
  std::vector<TLorentzVector>::const_iterator begin() const;
  std::vector<TLorentzVector>::const_iterator end() const;

  bool needSyst(int) const;

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setOrigReco();
  void setCurrentP(int);
  std::string getName() const {return GenName;};

  bool findCut(const std::vector<std::string>&, std::string);

  PType type;
  std::unordered_map<std::string, PartStats> pstats;
  const std::map<PType,CUTS> cutMap = {{PType::Electron, CUTS::eGElec}, {PType::Muon, CUTS::eGMuon},
				  {PType::Tau, CUTS::eGTau}};


protected:
  void getPartStats(std::string);
  TTree* BOOM;
  std::string GenName;
  std::unordered_map<CUTS, std::string, EnumHash> jetNameMap = {
    {CUTS::eRJet1, "Jet1"},               {CUTS::eRJet2, "Jet2"},
    {CUTS::eRCenJet, "CentralJet"},      {CUTS::eRBJet, "BJet"},
    {CUTS::eR1stJet, "FirstLeadingJet"},  {CUTS::eR2ndJet, "SecondLeadingJet"},
    {CUTS::eRWjet, "WJet"}
  };

 private:
  //vector<double>* mpt = 0;
  //vector<double>* meta = 0;
  //vector<double>* mphi = 0;
  //vector<double>* menergy = 0;
  uint  m_n;
  float m_pt[MAXINDEX];
  float m_phi[MAXINDEX];
  float m_eta[MAXINDEX];
  float m_mass[MAXINDEX];

  std::vector<TLorentzVector> Reco;
  std::vector<TLorentzVector> *cur_P;
  std::vector<std::string> syst_names;
  std::vector<std::vector<TLorentzVector>* > systVec;

  std::string activeSystematic;
};

class Photon : public Particle {
public:
  Photon();
  Photon(TTree*, std::string, std::vector<std::string>, std::string);

  bool get_Iso(int, double, double) const;
  std::bitset<8> cbIDphot1;
  std::bitset<8> cbIDphot2;

  float hoverE[MAXINDEX];
  float phoR[MAXINDEX];
  float sigmaIEtaIEta[MAXINDEX];
  float pfIso_all[MAXINDEX];
  float pfIso_chg[MAXINDEX];
  bool eleVeto[MAXINDEX];
  bool hasPixelSeed[MAXINDEX];
  Int_t cutBasedID[MAXINDEX];
  bool mvaID_WP80[MAXINDEX];
  bool mvaID_WP90[MAXINDEX];
};


/////////////////////////////////////////////////////////////////
class Generated : public Particle {

public:
  Generated();
  Generated(TTree*, std::string, std::vector<std::string>);

  int  genPartIdxMother[MAXINDEX];
  int  pdg_id[MAXINDEX];
  int  status[MAXINDEX];
  int  statusFlags[MAXINDEX];
  int  numDaught[MAXINDEX]; //01.15.19

};

/////////////////////////////////////////////////////////////////////////
class GenHadronicTaus : public Particle {

public:
  GenHadronicTaus();
  GenHadronicTaus(TTree*, std::string, std::vector<std::string>);

  int  genPartIdxMother[MAXINDEX];
  int  decayMode[MAXINDEX];

};

/////////////////////////////////////////////////////////////////////////
class GenJets : public Particle {

public:
  GenJets();
  GenJets(TTree*, std::string, std::vector<std::string>);

UChar_t  genHadronFlavor[MAXINDEX];
int      genPartonFlavor[MAXINDEX];

};

/////////////////////////////////////////////////////////////////////////
class Jet : public Particle {

public:
  Jet(TTree*, std::string, std::vector<std::string>, std::string);

  std::vector<CUTS> findExtraCuts();
  std::vector<CUTS> overlapCuts(CUTS);
  bool passedLooseJetID(int);
  bool passedTightJetID(int);
  bool passedTightLepVetoJetID(int);
  bool getPileupJetID(int, int);

  float area[MAXINDEX];
  // float bDiscriminator[MAXINDEX];
  float bDiscriminatorCSVv2[MAXINDEX];
  float bDiscriminatorDeepCSV[MAXINDEX];
  float bDiscriminatorDeepFlav[MAXINDEX];

  float chargedEmEnergyFraction[MAXINDEX];
  float chargedHadronEnergyFraction[MAXINDEX];
  float neutralEmEnergyFraction[MAXINDEX];
  float neutralHadEnergyFraction[MAXINDEX];
  int jetId[MAXINDEX];
  int nMuons[MAXINDEX];
  int numberOfConstituents[MAXINDEX];
  int puID[MAXINDEX];
  int partonFlavour[MAXINDEX];
  float rawFactor[MAXINDEX];
  int matchingMuonIdx1[MAXINDEX];
  int matchingMuonIdx2[MAXINDEX];
  int genJetIdx[MAXINDEX];

 protected:

};

class FatJet : public Particle {

public:
  FatJet(TTree*, std::string, std::vector<std::string>, std::string);

  std::vector<CUTS> findExtraCuts();
  std::vector<CUTS> overlapCuts(CUTS);

  float tau1[MAXINDEX];
  float tau2[MAXINDEX];
  float tau3[MAXINDEX];
  float tau4[MAXINDEX];
  float PrunedMass[MAXINDEX];
  float SoftDropMass[MAXINDEX];

};

class Lepton : public Particle {

public:
  Lepton(TTree*, std::string, std::string, std::vector<std::string>);

  std::vector<CUTS> findExtraCuts();

  double charge(uint)const;
  int _charge[MAXINDEX];

  virtual bool get_Iso(int, double, double) const {return false;}
};

class Electron : public Lepton {

public:
  Electron(TTree*, std::string, std::vector<std::string>, std::string);

  bool get_Iso(int, double, double) const;

  std::bitset<8> cbIDele1;
  std::bitset<8> cbIDele2;
  std::bitset<8> cbHLTIDele1;
  std::bitset<8> cbHLTIDele2;

  // Electron MVA ID fall 2017
  float miniPFRelIso_all[MAXINDEX];
  float miniPFRelIso_chg[MAXINDEX];
  float mvaFall17V2Iso[MAXINDEX];
  float mvaFall17V2Iso_WP80[MAXINDEX];
  float mvaFall17V2Iso_WP90[MAXINDEX];
  float mvaFall17V2Iso_WPL[MAXINDEX];
  float mvaFall17V2noIso[MAXINDEX];
  float mvaFall17V2noIso_WP80[MAXINDEX];
  float mvaFall17V2noIso_WP90[MAXINDEX];
  float mvaFall17V2noIso_WPL[MAXINDEX];
  bool mvaTTH[MAXINDEX];
  float pfRelIso03_all[MAXINDEX];
  float pfRelIso03_chg[MAXINDEX];
  int cutBased[MAXINDEX];
  bool isPassHEEPId[MAXINDEX];

};



class Muon : public Lepton {

public:
  Muon(TTree*, std::string, std::vector<std::string>, std::string);

  bool get_Iso(int, double, double) const;

  bool soft[MAXINDEX];
  bool loose[MAXINDEX];
  bool medium[MAXINDEX];
  bool tight[MAXINDEX];
  float miniPFRelIso_all[MAXINDEX];
  float miniPFRelIso_chg[MAXINDEX];
  float pfRelIso03_all[MAXINDEX];
  float pfRelIso03_chg[MAXINDEX];
  float pfRelIso04_all[MAXINDEX];
  bool isGlobal[MAXINDEX];
  bool isPFCand[MAXINDEX];
  bool isTracker[MAXINDEX];
};

class Taus : public Lepton {

public:
  Taus(TTree*, std::string, std::vector<std::string>, std::string);

  std::vector<CUTS> findExtraCuts();
  bool get_Iso(int, double, double) const;
  bool pass_against_Elec(CUTS, int);
  bool pass_against_Muon(CUTS, int);

  std::bitset<8> tau1minIso;
  std::bitset<8> tau1maxIso;

  std::bitset<8> tau2minIso;
  std::bitset<8> tau2maxIso;

  std::bitset<8> tau1ele;
  std::bitset<8> tau1mu;

  std::bitset<8> tau2ele;
  std::bitset<8> tau2mu;

   // --------  Anti-particle discriminators ------- //
  UChar_t againstElectron[MAXINDEX];
  UChar_t againstMuon[MAXINDEX];

  // --------  Decay modes ------- //
  bool DecayModeOldDMs[MAXINDEX];
  bool DecayModeNewDMs[MAXINDEX];
  int decayModeInt[MAXINDEX];

  // -------- Tau isolation ----------- //
  UChar_t TauIdDiscr[MAXINDEX];

  // ------- Tau-related quantities ---------- //
  float leadTkDeltaEta[MAXINDEX];
  float leadTkDeltaPhi[MAXINDEX];
  float leadTkPtOverTauPt[MAXINDEX];
  float dz[MAXINDEX];
  float dxy[MAXINDEX];
  float chargedIsoPtSum[MAXINDEX];
  float neutralIso[MAXINDEX];
  float puCorr[MAXINDEX];

  // ----- Tau gen-matching for ID SFs ----- //
  UChar_t genPartFlav[MAXINDEX];
  Int_t genPartIdx[MAXINDEX];

};



#endif
