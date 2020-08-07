#include "MET.h"
#include <algorithm>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, true);  BOOM->SetBranchAddress(name, &variable);

//particle is a objet that stores multiple versions of the particle candidates
Met::Met(TTree* _BOOM, std::string _GenName,  std::vector<std::string> _syst_names, double _MT2mass) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names), MT2mass(_MT2mass)  {

  //std::cout << "Met branch = " << GenName << std::endl;
  SetBranch((GenName+"_pt").c_str(), met_pt);
  SetBranch((GenName+"_phi").c_str(), met_phi);

  // Then we get the default met
  if(GenName.compare("MET") != 0){
    samedeft1met = false; 
    SetBranch("MET_pt", def_met_pt);
    SetBranch("MET_phi", def_met_phi);
  }
  // Then we get the raw met for Type-I corrections
  SetBranch("RawMET_pt", raw_met_pt);
  SetBranch("RawMET_phi", raw_met_phi);

  // Initialize the systRawMet vector   
  for(auto name: syst_names){
    if(name == "orig") 
      systRawMetVec.push_back(new TLorentzVector);
    else if(name.find("Met")!=std::string::npos){
      systRawMetVec.push_back(new TLorentzVector);
    }
    else if(name.find("weight")!=std::string::npos){
      systRawMetVec.push_back(nullptr);
    }else if(name.find("Tau_qcd")!=std::string::npos){
      systRawMetVec.push_back(nullptr);
    }else
      systRawMetVec.push_back(new TLorentzVector);
  }
  
  // For TreatMuonsAsNeutrinos
  systdeltaMEx.resize(syst_names.size());
  systdeltaMEy.resize(syst_names.size());
  
  // For HT, MHT
  syst_HT.resize(syst_names.size());
  syst_MHT.resize(syst_names.size());
  syst_MHTphi.resize(syst_names.size());
 
  if( std::find(syst_names.begin(), syst_names.end(), "MetUncl_Up") != syst_names.end() && _BOOM->GetListOfBranches()->FindObject((GenName+"_MetUnclustEnUpDeltaX").c_str()) !=0){
    SetBranch((GenName+"_MetUnclustEnUpDeltaX").c_str(), MetUnclUp[0]);
    SetBranch((GenName+"_MetUnclustEnUpDeltaY").c_str(), MetUnclUp[1]);
    Unclup = std::find(syst_names.begin(), syst_names.end(), "MetUncl_Up") -syst_names.begin();
  }
  if( std::find(syst_names.begin(), syst_names.end(), "MetUncl_Down") != syst_names.end() && _BOOM->GetListOfBranches()->FindObject((GenName+"_MetUnclustEnUpDeltaY").c_str()) !=0){
    SetBranch((GenName+"_MetUnclustEnUpDeltaX").c_str(), MetUnclDown[0]);
    SetBranch((GenName+"_MetUnclustEnUpDeltaY").c_str(), MetUnclDown[1]);
    Uncldown = std::find(syst_names.begin(), syst_names.end(), "MetUncl_Down")- syst_names.begin();
  }
  activeSystematic=0;
}

bool Met::needSys(int syst) const {
  return systRawMetVec.at(syst) != nullptr;
}

void Met::addPtEtaPhiESyst(double ipt,double ieta, double iphi, double ienergy, int syst){
  systRawMetVec.at(syst)->SetPtEtaPhiE(ipt,ieta,iphi,ienergy);
}


void Met::addP4Syst(TLorentzVector mp4, int syst){
  systRawMetVec.at(syst)->SetPtEtaPhiE(mp4.Pt(),mp4.Eta(),mp4.Phi(),mp4.E());
}


void Met::init(){
  //cleanup of the particles
  //keep this if there is any ever some need for a unchanged met
  Reco.SetPtEtaPhiM(met_pt,0,met_phi,met_pt);
  // Define the four vectors for the various kinds of Met that will be needed.
  T1Met.SetPtEtaPhiM(met_pt,0,met_phi,met_pt);
  if(samedeft1met){
    DefMet.SetPtEtaPhiM(met_pt,0,met_phi, met_pt);
  }
  else{
    DefMet.SetPtEtaPhiM(def_met_pt,0,def_met_phi,def_met_pt);
  }
  RawMet.SetPtEtaPhiM(raw_met_pt,0,raw_met_phi,raw_met_pt);

  for(int i=0; i < (int) syst_names.size(); i++) {
     // initialize the Raw MET vector
     if(systRawMetVec.at(i) != nullptr) addP4Syst(RawMet, i);
     // initialize the deltaMuMet vectors to zero
     fill(systdeltaMEx.begin(), systdeltaMEx.end(), 0);
     fill(systdeltaMEy.begin(), systdeltaMEy.end(), 0);
  } 

  cur_P=&Reco;

}

void Met::propagateJetEnergyCorr(TLorentzVector recoJet, double const& jet_pt_up, double const& jet_pt_down, std::string& systname, int syst){

  if(systRawMetVec.at(syst) == nullptr) return;

  float jet_cosPhi = cos(recoJet.Phi());
  float jet_sinPhi = sin(recoJet.Phi());

  // Update the nominal values:
  float met_px_shifted = systRawMetVec.at(syst)->Px();
  float met_py_shifted = systRawMetVec.at(syst)->Py();

  // update shifted values
  met_px_shifted = met_px_shifted - (jet_pt_up - jet_pt_down) * jet_cosPhi;
  met_py_shifted = met_py_shifted - (jet_pt_up - jet_pt_down) * jet_sinPhi;

  // Change the Raw Met vector:
  systRawMetVec.at(syst)->SetPxPyPzE(met_px_shifted, met_py_shifted, systRawMetVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_shifted,2) + pow(met_py_shifted,2)));

}

// This is only used for 2017 data/MC for the EE noise 
void Met::propagateUnclEnergyUnctyEE(double const& delta_x_T1Jet, double const& delta_y_T1Jet, double const& delta_x_rawJet, double const& delta_y_rawJet, std::string& systname, int syst){

  if(systRawMetVec.at(syst) == nullptr) return;

  // Remove the L1L2L3 - L1 corrected jets in the EE region from the default MET branch
  float new_def_met_px = DefMet.Px() + delta_x_T1Jet;
  float new_def_met_py = DefMet.Py() + delta_y_T1Jet;

  DefMet.SetPxPyPzE(new_def_met_px, new_def_met_py, DefMet.Pz(), sqrt(pow(new_def_met_px,2) + pow(new_def_met_py, 2)));

  // Get the unclustered energy part that is removed in the v2 recipe
  float met_unclEE_x = DefMet.Px() - T1Met.Px(); //t1met_px;
  float met_unclEE_y = DefMet.Py() - T1Met.Py(); //t1met_py;

  // Finalize the v2 recipe for the rawMET by removing the unclustered part in the EE region
  float met_px_unclshift = systRawMetVec.at(syst)->Px();
  float met_py_unclshift = systRawMetVec.at(syst)->Py();

  met_px_unclshift = met_px_unclshift + delta_x_rawJet - met_unclEE_x;
  met_py_unclshift = met_py_unclshift + delta_y_rawJet - met_unclEE_y;

  // Add this to the systematics vector
  systRawMetVec.at(syst)->SetPxPyPzE(met_px_unclshift, met_py_unclshift, systRawMetVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_unclshift,2) + pow(met_py_unclshift,2)));

}

void Met::propagateUnclEnergyUncty(std::string& systname, int syst){

  if(systRawMetVec.at(syst) == nullptr) return;

  if(systname.find("MetUncl") == std::string::npos) return;

  // This will only apply for the MetUncl uncertainty in MC
  float met_px_unclEnshift = systRawMetVec.at(0)->Px(); // 0 refers to the nominal value, which at this point should already have all corrections applied
  float met_py_unclEnshift = systRawMetVec.at(0)->Py();
  
  if(systname.find("_Up") != std::string::npos){

    met_px_unclEnshift = met_px_unclEnshift + MetUnclUp[0];
    met_py_unclEnshift = met_py_unclEnshift + MetUnclUp[1];
    }
    else if(systname.find("_Down") != std::string::npos){

    met_px_unclEnshift = met_px_unclEnshift + MetUnclDown[0];
    met_py_unclEnshift = met_py_unclEnshift + MetUnclDown[1];
  }

  systRawMetVec.at(syst)->SetPxPyPzE(met_px_unclEnshift, met_py_unclEnshift, systRawMetVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_unclEnshift,2) + pow(met_py_unclEnshift,2)));

}

void Met::update(int syst=0){

  if(systRawMetVec.at(syst) == nullptr) return;

  // Treat muons as neutrinos. This is done on the systRawMetVec which is the vector that has all the JERC propagated. 
  systRawMetVec.at(syst)->SetPxPyPzE(systRawMetVec.at(syst)->Px()+systdeltaMEx[syst], 
                               systRawMetVec.at(syst)->Py()+systdeltaMEy[syst], 
                               systRawMetVec.at(syst)->Pz(), 
                               TMath::Sqrt(pow(systRawMetVec.at(syst)->Px()+systdeltaMEx[syst],2) + pow(systRawMetVec.at(syst)->Py()+systdeltaMEy[syst],2)));

}

void Met::calculateHtAndMHt(PartStats& stats, Jet& jet, int syst=0){

  if(systRawMetVec.at(syst) == nullptr) return;

  float sumpxForMht=0;
  float sumpyForMht=0;
  float sumptForHt=0;

  // Calculates HT and MHT.
  int i=0;
  for(auto jetVec: jet){
    bool add = true;
    if( (jetVec.Pt() < stats.dmap.at("JetPtForMhtAndHt")) || 
    	(abs(jetVec.Eta()) > stats.dmap.at("JetEtaForMhtAndHt")) || 
    	(stats.bfind("ApplyJetLooseIDforMhtAndHt") && !jet.passedLooseJetID(i)) ||
    	(stats.bfind("ApplyJetTightIDforMhtAndHt") && !jet.passedTightJetID(i)) ) add = false;
    
    if(add) {
      sumpxForMht -= jetVec.Px();
      sumpyForMht -= jetVec.Py();
      sumptForHt  += jetVec.Pt();
    }

    i++;
  }
  syst_HT.at(syst)=sumptForHt;
  syst_MHT.at(syst)= sqrt( pow(sumpxForMht,2.0) + pow(sumpyForMht,2.0) );
  syst_MHTphi.at(syst)=atan2(sumpyForMht,sumpxForMht);

}

double Met::pt()const         {return cur_P->Pt();}
double Met::px()const         {return cur_P->Px();}
double Met::py()const         {return cur_P->Py();}
double Met::eta()const        {return cur_P->Eta();}
double Met::phi()const        {return cur_P->Phi();}
double Met::energy()const     {return cur_P->E();}


TLorentzVector Met::p4()const {return (*cur_P);}
TLorentzVector& Met::p4() {return *cur_P;}

double Met::HT() const {return syst_HT.at(activeSystematic);};
double Met::MHT() const {return syst_MHT.at(activeSystematic);};
double Met::MHTphi() const {return syst_MHTphi.at(activeSystematic);};

double Met::MT2(TLorentzVector& pa, TLorentzVector& pb){
  mt2_event.set_momenta( pa, pb, px(), py());
  mt2_event.set_mn( MT2mass );
  return mt2_event.get_mt2();
}

void Met::setMT2Mass(double mass){
  MT2mass=mass;
}

TLorentzVector Met::getNominalP(){
	TLorentzVector nominalP4 = *systRawMetVec.at(0); 
	// nominalP4.SetPxPyPzE(systRawMetVec.at(0)->Px(), systRawMetVec.at(0)->Py(), systRawMetVec.at(0)->Pz(), systRawMetVec.at(0)->E());
	return nominalP4;
}

void Met::setCurrentP(int syst){
  // Use the information from the raw met vector instead, since this is the vector all jet corrections are propagated to.
  if(systRawMetVec.at(syst) == nullptr) {
    cur_P = systRawMetVec.at(0);
    activeSystematic = 0;
  } else {
    cur_P = systRawMetVec[syst];
    activeSystematic=syst;
  }

}

void Met::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
}

