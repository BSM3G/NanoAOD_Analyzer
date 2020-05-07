#ifndef Systematics_hh
#define Systematics_hh

#include <string>
#include <iostream>
#include <functional>
#include <unordered_map>
#include "Particle.h"
//#include <boost/unordered_map.hpp>

// we will put stuff from the main analyser here once we know what


class TRandom3;

class Systematics {

public:
  Systematics();
  Systematics(std::unordered_map<std::string, PartStats> const &distats);
  ~Systematics();

  void init();

  void shiftParticle(Particle& jet, double const& corrJetPt, double& corrJetMass, std::string& syst_name, int syst);
  void shiftLepton(Lepton& lepton, TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy, int syst);
  void loadScaleRes(const PartStats& smear, const PartStats& syst, std::string syst_name);

  // void shiftParticle(Particle& jet, TLorentzVector recJet, double const& ratio, double& dPx, double& dPy, int syst);
  // void shiftParticleRes(Particle& jet, TLorentzVector recJet, double const& shiftSF, int syst);
  // void shiftParticleScale(Particle& jet, TLorentzVector recoJet, double const& jer_sf_nom, double const& jes_delta, double const& sigma, int syst);

private:
  double scale;
  double resolution;
};
#endif /*Systematics_hh*/
