#include "MET.h"
#include <algorithm>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, true);  BOOM->SetBranchAddress(name, &variable);

//particle is a objet that stores multiple versions of the particle candidates
Met::Met(TTree* _BOOM, std::string _GenName,  std::vector<std::string> _syst_names, double _MT2mass) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names), MT2mass(_MT2mass)  {

  std::cout << "Met branch = " << GenName << std::endl;

  SetBranch((GenName+"_pt").c_str(), met_pt);
  SetBranch((GenName+"_phi").c_str(), met_phi);

  // Then we get the default met
  if(GenName.compare("MET") != 0){
    samedeft1met = false;
    // std::cout << "Not here! " << std::endl;
    SetBranch("MET_pt", def_met_pt);
    SetBranch("MET_phi", def_met_phi);
  }
  // Then we get the raw met for Type-I corrections
  SetBranch("RawMET_pt", raw_met_pt);
  SetBranch("RawMET_phi", raw_met_phi);

  // Initialize the systRawMet vector 
  for(auto name: syst_names){
    systRawMetVec.push_back(new TLorentzVector);
  }

  systdeltaMEx.resize(syst_names.size());
  systdeltaMEy.resize(syst_names.size());
  syst_HT.resize(syst_names.size());
  syst_MHT.resize(syst_names.size());
  syst_MHTphi.resize(syst_names.size());
  

  for( auto name : syst_names) {
    if(name == "orig") 
      systVec.push_back(new TLorentzVector);
    else if(name.find("Met")!=std::string::npos){
      systVec.push_back(new TLorentzVector);
    }
    else if(name.find("weight")!=std::string::npos){
      systVec.push_back(nullptr);
    }else if(name.find("Tau_qcd")!=std::string::npos){
      systVec.push_back(nullptr);
    }else
      systVec.push_back(new TLorentzVector);
  }
  
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


void Met::addPtEtaPhiESyst(double ipt,double ieta, double iphi, double ienergy, int syst){
  systVec.at(syst)->SetPtEtaPhiE(ipt,ieta,iphi,ienergy);
}


void Met::addP4Syst(TLorentzVector mp4, int syst){
  systVec.at(syst)->SetPtEtaPhiE(mp4.Pt(),mp4.Eta(),mp4.Phi(),mp4.E());
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
  
  // Get the x and y components of the raw MET
  met_px = RawMet.Px();
  met_py = RawMet.Py();
  def_met_px = DefMet.Px();
  def_met_py = DefMet.Py();
  t1met_px = Reco.Px();
  t1met_py = Reco.Py();

  for(int i=0; i < (int) syst_names.size(); i++) {
    systRawMetVec.at(i)->SetPxPyPzE(RawMet.Px(), RawMet.Py(), RawMet.Pz(), RawMet.E());
  } 

  for(int i=0; i < (int) syst_names.size(); i++) {
    
    if(i == Unclup) systVec.at(i)->SetPxPyPzE(MetUnclUp[0]+Reco.Px(),MetUnclUp[1]+Reco.Py(),0,sqrt(pow(MetUnclUp[0]+Reco.Px(),2)+pow(MetUnclUp[1]+Reco.Py(),2)));
    else if(i == Uncldown) systVec.at(i)->SetPxPyPzE(MetUnclDown[0]+Reco.Px(),MetUnclDown[1]+Reco.Py(),0,sqrt(pow(MetUnclDown[0]+Reco.Px(),2)+pow(MetUnclDown[1]+Reco.Py(),2)));
    else if(systVec.at(i) != nullptr) addP4Syst(Reco, i);
    // else if(systVec.at(i) != nullptr) addP4Syst(RawMet, i);
    
    fill(systdeltaMEx.begin(), systdeltaMEx.end(), 0);
    fill(systdeltaMEy.begin(), systdeltaMEy.end(), 0);
  }
   cur_P=&Reco;
   // cur_P=&RawMet;

  // syst_HT[activeSystematic]=0.;
  // syst_MHT[activeSystematic]=0.;
  // syst_MHTphi[activeSystematic]=99.;
}


void Met::update(PartStats& stats, Jet& jet, int syst=0){
  ///Calculates met from values from each file plus smearing and treating muons as neutrinos
  if(systVec.at(syst) == nullptr) return;
  double sumpxForMht=0;
  double sumpyForMht=0;
  double sumptForHt=0;

  int i=0;
  for(auto jetVec: jet) {
    bool add = true;
    if( (jetVec.Pt() < stats.dmap.at("JetPtForMhtAndHt")) ||
        (abs(jetVec.Eta()) > stats.dmap.at("JetEtaForMhtAndHt")) ||
        ( stats.bfind("ApplyJetLooseIDforMhtAndHt") &&
          !jet.passedLooseJetID(i) ) ) add = false;
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

  systVec.at(syst)->SetPxPyPzE(systVec.at(syst)->Px()+systdeltaMEx[syst], 
                               systVec.at(syst)->Py()+systdeltaMEy[syst], 
                               systVec.at(syst)->Pz(), 
                               TMath::Sqrt(pow(systVec.at(syst)->Px()+systdeltaMEx[syst],2) + pow(systVec.at(syst)->Py()+systdeltaMEy[syst],2)));

}

void Met::propagateJetEnergyCorr(TLorentzVector recoJet, double const& jet_pt_up, double const& jet_pt_down, std::string& systname, int syst){

  if(systRawMetVec.at(syst) == nullptr) return;

  double jet_cosPhi = cos(recoJet.Phi());
  double jet_sinPhi = sin(recoJet.Phi());

  // std::cout << "Before: RawMet.px() = " << systRawMetVec.at(syst)->Px() << ", RawMet.py() = " << systRawMetVec.at(syst)->Py() << std::endl;
  // Update the nominal values:
  met_px_shifted = systRawMetVec.at(syst)->Px();
  met_py_shifted = systRawMetVec.at(syst)->Py();
  // met_px_shifted = RawMet.Px();
  // met_py_shifted = RawMet.Py();

  // std::cout << "met_px_shifted (b) = " << met_px_shifted << ", met_py_shifted (b) = " << met_py_shifted << std::endl;

  // update shifted values
  met_px_shifted = met_px_shifted - (jet_pt_up - jet_pt_down) * jet_cosPhi;
  met_py_shifted = met_py_shifted - (jet_pt_up - jet_pt_down) * jet_sinPhi;

  // std::cout << "met_px_shifted (a) = " << met_px_shifted << ", met_py_shifted (a) = " << met_py_shifted << std::endl;

  // Change the Raw Met vector:
  systRawMetVec.at(syst)->SetPxPyPzE(met_px_shifted, met_py_shifted, systRawMetVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_shifted,2) + pow(met_py_shifted,2)));

  // std::cout << "After: RawMet.px() (=- (jet_pt_L1L2L3 - jet_pt_L1)*jet_cosPhi) = " << systRawMetVec.at(syst)->Px() << ", RawMet.py() (=- (jet_pt_L1L2L3 - jet_pt_L1)*jet_sinPhi) = " << systRawMetVec.at(syst)->Py() << std::endl; 

  // Add this to the systematics vector
  // systVec.at(syst)->SetPxPyPzE(met_px_shifted, met_py_shifted, systVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_shifted,2) + pow(met_py_shifted,2)));

}

// This is only used for 2017 data/MC for the EE noise 
void Met::propagateUnclEnergyUnctyEE(double const& delta_x_T1Jet, double const& delta_y_T1Jet, double const& delta_x_rawJet, double const& delta_y_rawJet, std::string& systname, int syst){

  if(systRawMetVec.at(syst) == nullptr) return;

  // std::cout << std::endl <<  "Before EE: RawMet.px() = " << systRawMetVec.at(syst)->Px() << ", RawMet.py() = " << systRawMetVec.at(syst)->Py() << std::endl;
  // std::cout << "def_met_px (b) = " << DefMet.Px() << ", def_met_py (b) = " << DefMet.Py() << std::endl;
  // std::cout << "t1_met_px (b) = " << T1Met.Px() << ", t1_met_py (b) = " << T1Met.Py() << std::endl;

  // Remove the L1L2L3 - L1 corrected jets in the EE region from the default MET branch
  double new_def_met_px = DefMet.Px() + delta_x_T1Jet;
  double new_def_met_py = DefMet.Py() + delta_y_T1Jet;

  // std::cout << "new_def_met_px (=+ delta_x_T1Jet) = " << new_def_met_px << ", new_def_met_py (=+ delta_y_T1Jet) = " << new_def_met_py << std::endl;

  DefMet.SetPxPyPzE(new_def_met_px, new_def_met_py, DefMet.Pz(), sqrt(pow(new_def_met_px,2) + pow(new_def_met_py, 2)));

  // std::cout << "def_met_px (a) = " << DefMet.Px() << ", def_met_py (a) = " << DefMet.Py() << std::endl;

  // Get the unclustered energy part that is removed in the v2 recipe
  double met_unclEE_x = DefMet.Px() - T1Met.Px(); //t1met_px;
  double met_unclEE_y = DefMet.Py() - T1Met.Py(); //t1met_py;

  // std::cout << "met_unclEE_x (=def_met_px - t1met_px) = " <<  met_unclEE_x << ", met_unclEE_y (=def_met_py - t1met_py) = " <<  met_unclEE_y << std::endl;

  // Finalize the v2 recipe for the rawMET by removing the unclustered part in the EE region
  double met_px_unclshift = systRawMetVec.at(syst)->Px();
  double met_py_unclshift = systRawMetVec.at(syst)->Py();

  //std::cout << "met_px_unclshift (RawMetVec) = " << met_px_unclshift << ", met_py_unclshift (RawMetVec) = " << met_py_unclshift << std::endl;

  met_px_unclshift = met_px_unclshift + delta_x_rawJet - met_unclEE_x;
  met_py_unclshift = met_py_unclshift + delta_y_rawJet - met_unclEE_y;

  // std::cout << "met_px_unclshift (met_px_unclshift +delta_x_rawJet - met_unclEE_x) = " << met_px_unclshift << ", met_py_unclshift (met_py_unclshift +delta_y_rawJet - met_unclEE_y) = " << met_py_unclshift << std::endl;

  // Add this to the systematics vector
  systRawMetVec.at(syst)->SetPxPyPzE(met_px_unclshift, met_py_unclshift, systRawMetVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_unclshift,2) + pow(met_py_unclshift,2)));

  // std::cout << std::endl << "After EE: RawMet.px() = " << systRawMetVec.at(syst)->Px() << ", RawMet.py() = " << systRawMetVec.at(syst)->Py() << std::endl; 

}

void Met::propagateUnclEnergyUncty(std::string& systname, int syst){

  if(systRawMetVec.at(syst) == nullptr) return;

  if(systname.find("MetUncl") == std::string::npos) return;

  // This will only apply for the MetUncl uncertainty in MC
  double met_px_unclEnshift = systRawMetVec.at(0)->Px(); // 0 refers to the nominal value, which at this point should already have all corrections applied
  double met_py_unclEnshift = systRawMetVec.at(0)->Py();
  
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

void Met::setCurrentP(int syst){
  // Use the information from the raw met vector instead, since this is the vector all jet corrections
  // are propagated to.
  if(systRawMetVec.at(syst) == nullptr) {
    cur_P = systRawMetVec.at(0);
    activeSystematic = 0;
  } else {
    cur_P = systRawMetVec[syst];
    activeSystematic=syst;
  }

/*
  if( systVec.at(syst) == nullptr) {
    cur_P = systVec.at(0);
    activeSystematic = 0;
  } else {
    cur_P = systVec[syst];
    activeSystematic=syst;
  }
*/
}


void Met::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
}

