#ifndef Met_h
#define Met_h

// system include files
#include <memory>

// user include files
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

#include <TTree.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include "Particle.h"
#include "mt2/mt2_bisect.hh"


#include "tokenizer.hpp"
#include "Cut_enum.h"

//using namespace std;
typedef unsigned int uint;


class Met {

public:
  Met(){};
  Met(TTree*, std::string, std::string, std::vector<std::string>, std::string){};
  Met(TTree*, std::string, std::vector<std::string>, double, std::string);
  virtual ~Met() {}

  virtual std::vector<CUTS> findExtraCuts(){return std::vector<CUTS>();}
  void init();
  void unBranch();
  double pt() const;
  double px() const;
  double py() const;
  double eta() const;
  double phi() const;
  double energy() const;
  double HT() const;
  double MHT() const;
  double MHTphi() const;
  double MT2(TLorentzVector&, TLorentzVector&);
  TLorentzVector p4() const;
  TLorentzVector& p4();

  
  void propagateJER(TLorentzVector recoJet, double const& jer_sf_nom, double const& jer_sf_shift, int syst);
  void propagateJES(TLorentzVector recoJet, double const& jer_sf_nom, double const& jes_delta, double const& jes_sigma, int syst);
  void propagateJetEnergyCorr(TLorentzVector recoJet, double const& jer_sf_nom, double const& jec_param, std::string& systname, int syst);
  

  void addPtEtaPhiESyst(double, double, double, double, int);
  void addP4Syst(TLorentzVector, int);
  void setMT2Mass(double);
  void setCurrentP(int);
  std::string getName() {return GenName;};
  void update(PartStats&, Jet&, int);

  TLorentzVector RecoMet;
  TLorentzVector *cur_P;

  std::vector<TLorentzVector* > systVec;
  std::vector<double> systdeltaMEx;
  std::vector<double> systdeltaMEy;
  std::vector<double> syst_HT;
  std::vector<double> syst_MHT;
  std::vector<double> syst_MHTphi;


  int activeSystematic;

protected:
  TTree* BOOM;
  std::string GenName;
  std::vector<std::string> syst_names;
  double MT2mass;
  
  float met_pt;
  float met_phi;
  float met_px;
  float met_py;
  float met_px_nom;
  float met_py_nom;
  float met_pt_nom;
  float met_phi_nom;
  float met_px_jerShifted;
  float met_py_jerShifted;
  float met_px_jesShifted;
  float met_py_jesShifted;
  float met_px_shifted;
  float met_py_shifted;
  float finalmet_px_shifted;
  float finalmet_py_shifted;

  double jet_unclEnThreshold = 15.0; // energy threshold below which jets are considered as "unclustered energy" (cf. JetMETCorrections/Type1MET/python/correctionTermsPfMetType1Type2_cff.py )
  //note this is only for pt and phi
  double MetUnclUp[2] = {0, 0};
  double MetUnclDown[2] = {0, 0};
  int Unclup=-1;
  int Uncldown=-1;
  mt2_bisect::mt2 mt2_event;
  

};



#endif
