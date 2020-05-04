#include "Systematics.h"


Systematics::Systematics(){}

Systematics::Systematics(std::unordered_map<std::string, PartStats> const &distats){

}
Systematics::~Systematics(){}

void Systematics::init(){

}


void Systematics::shiftParticleRes(Particle& jet, TLorentzVector recoJet, double const& shiftSF, int syst){

  // Recalculate pt and mass values:
  double jet_pt_res_shifted = shiftSF * recoJet.Pt();
  double jet_mass_res_shifted = shiftSF * recoJet.M();

  // std::cout << "jet_pt_res_shifted = " << jet_pt_res_shifted << std::endl;
  
  // If nominal calculation, check the sign of these quantities
  if(syst == 0){
    if(jet_pt_res_shifted < 0.0) jet_pt_res_shifted *= -1.0;
    if(jet_mass_res_shifted < 0.0) jet_mass_res_shifted *= -1.0;
  }
  // Set the new components of the 4-momentum
  TLorentzVector shiftedRecoJet;
  shiftedRecoJet.SetPtEtaPhiM(jet_pt_res_shifted, recoJet.Eta(), recoJet.Phi(), jet_mass_res_shifted); // Eta and phi components shouldn't change
  // Change the particle content for the particle
  jet.addP4Syst(shiftedRecoJet, syst);
  return;
}



void Systematics::shiftParticleScale(Particle& jet, TLorentzVector recoJet, double const& jer_sf_nom, double const& jes_delta, double const& sigma, int syst){
  // Recalculate pt and mass values from nominal JER correction:
  double jet_pt_nom = jer_sf_nom * recoJet.Pt();
  double jet_mass_nom = jer_sf_nom * recoJet.M();

  // Make sure they are positive
  if(jet_pt_nom < 0.0) jet_pt_nom *= -1.0;
  if(jet_mass_nom < 0.0) jet_mass_nom *= -1.0;

  double jet_pt_jes_shifted = jet_pt_nom * (1.0 + sigma*jes_delta);
  double jet_mass_jes_shifted = jet_mass_nom * (1.0 + sigma*jes_delta);

  // Set the new components of the 4-momentum
  TLorentzVector shiftedRecoJet; 
  shiftedRecoJet.SetPtEtaPhiM(jet_pt_jes_shifted, recoJet.Eta(), recoJet.Phi(), jet_mass_jes_shifted);

  jet.addP4Syst(shiftedRecoJet, syst);
  return;

}

void Systematics::shiftParticle(Particle& jet, TLorentzVector recJet, double const& ratio, double& dPx, double& dPy, int syst){

   //add the shifted part up
   dPx+=recJet.Px()*(ratio-1);
   dPy+=recJet.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recJet*=ratio;
   jet.addP4Syst(recJet, syst);
   return;
}

void Systematics::shiftLepton(Lepton& lepton, TLorentzVector recoLep, TLorentzVector genLep, double& dPx, double& dPy, int syst){
  if (genLep == TLorentzVector(0,0,0,0)) {
    lepton.addP4Syst(recoLep, syst);
    return;
  }
  double ratio = ((genLep.Pt()*scale) + (recoLep.Pt() - genLep.Pt())*resolution)/recoLep.Pt();
  //cout<<"ratio  "<<ratio<<"  "<<scale<<"  "<<resolution    <<std::endl;
   //add the shifted part up
   dPx+=recoLep.Px()*(ratio-1);
   dPy+=recoLep.Py()*(ratio-1);
   //WARNING change the particle content for the particle
   recoLep*=ratio;
   lepton.addP4Syst(recoLep, syst);
   return;
}


void Systematics::loadScaleRes(const PartStats& smear, const PartStats& syst, std::string syst_name) {
  scale = 1;
  resolution = 1;
  if(smear.bfind("SmearTheParticle")) {
    scale = smear.dmap.at("PtScaleOffset");
    resolution = smear.dmap.at("PtResolutionOffset");
  } 
  if(syst_name.find("_Res_")) {
    resolution = syst_name.find("_Up") ? 1 + syst.dmap.at("res") : 1 - syst.dmap.at("res");
    scale=1;
  } else if(syst_name.find("_Scale_")) {
    scale = syst_name.find("_Up") ? 1+syst.dmap.at("scale") : 1- syst.dmap.at("scale");
    resolution=1;
  }
}

