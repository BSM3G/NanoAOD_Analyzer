#include "MET.h"
#include <algorithm>

#define SetBranch(name, variable) BOOM->SetBranchStatus(name, true);  BOOM->SetBranchAddress(name, &variable);

//particle is a objet that stores multiple versions of the particle candidates
Met::Met(TTree* _BOOM, std::string _GenName,  std::vector<std::string> _syst_names, double _MT2mass, std::string year) : BOOM(_BOOM), GenName(_GenName), syst_names(_syst_names), MT2mass(_MT2mass)  {

  if(year.compare("2017") == 0) GenName = (GenName+"FixEE2017").c_str();

  SetBranch((GenName+"_pt").c_str(), met_pt);
  SetBranch((GenName+"_phi").c_str(), met_phi);
  std::cout << "Name of the branch loaded for MET: " << (GenName+"_pt").c_str() << std::endl;

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

// Seems to be unused
void Met::addPtEtaPhiESyst(double ipt,double ieta, double iphi, double ienergy, int syst){
  systVec.at(syst)->SetPtEtaPhiE(ipt,ieta,iphi,ienergy);
}

// Seems to be unused
void Met::addP4Syst(TLorentzVector mp4, int syst){
  systVec.at(syst)->SetPtEtaPhiE(mp4.Pt(),mp4.Eta(),mp4.Phi(),mp4.E());
}

// Initialize MET
void Met::init(){
  //cleanup of the particles
  //keep this if there is any ever some need for a unchanged met
  RecoMet.SetPtEtaPhiM(met_pt,0,met_phi,met_pt);
  
  // Get the x and y components of the MET
  met_px = RecoMet.Px();
  met_py = RecoMet.Py();

  // std::cout << "met_px init (1) = " << met_px << ". met_py init (1) = " << met_py << std::endl;

  for(int i=0; i < (int) syst_names.size(); i++) {
    
    if(i == Unclup) systVec.at(i)->SetPxPyPzE(MetUnclUp[0]+RecoMet.Px(),MetUnclUp[1]+RecoMet.Py(),0,sqrt(pow(MetUnclUp[0]+RecoMet.Px(),2)+pow(MetUnclUp[1]+RecoMet.Py(),2)));
    else if(i == Uncldown) systVec.at(i)->SetPxPyPzE(MetUnclDown[0]+RecoMet.Px(),MetUnclDown[1]+RecoMet.Py(),0,sqrt(pow(MetUnclDown[0]+RecoMet.Px(),2)+pow(MetUnclDown[1]+RecoMet.Py(),2)));
    else if(systVec.at(i) != nullptr) addP4Syst(RecoMet, i);
    
    fill(systdeltaMEx.begin(), systdeltaMEx.end(), 0);
    fill(systdeltaMEy.begin(), systdeltaMEy.end(), 0);
  }
  cur_P=&RecoMet;

  // std::cout << "met_px init (2) = " << RecoMet.Px() << ". met_py init (2) = " << RecoMet.Py() << std::endl;

  // syst_HT[activeSystematic]=0.;
  // syst_MHT[activeSystematic]=0.;
  // syst_MHTphi[activeSystematic]=99.;
}


void Met::propagateJetEnergyCorr(TLorentzVector recoJet, double const& jer_sf_nom, double const& jec_param, std::string& systname, int syst){

  if(systVec.at(syst) == nullptr) return;

  // update shifted values
  met_px_shifted = systVec.at(syst)->Px();
  met_py_shifted = systVec.at(syst)->Py();

  double jet_pt_nom = jer_sf_nom * recoJet.Pt();
  // make sure this is positive
  if(jet_pt_nom < 0.0) jet_pt_nom *= -1.0;

  // Check that the jet pt is greater than the thrshold for unclustered energy in order to do the propagation
  if(jet_pt_nom <= jet_unclEnThreshold) return;

    // update those values we want for the systematic calculated, this also includes the nominal value.
  if(syst == 0){
    met_px_shifted = met_px_shifted - (jet_pt_nom - recoJet.Pt()) * cos(recoJet.Phi());
    met_py_shifted = met_py_shifted - (jet_pt_nom - recoJet.Pt()) * sin(recoJet.Phi());
  }
  else if(systname.find("_Res_") != std::string::npos){ // in this case jec_param = jer_sf_shift
  	met_px_shifted = met_px_shifted - (jec_param * recoJet.Pt() - jet_pt_nom) * cos(recoJet.Phi());
    met_py_shifted = met_py_shifted - (jec_param * recoJet.Pt() - jet_pt_nom) * sin(recoJet.Phi());
  }
  else if(systname.find("_Scale_") != std::string::npos){ // in this case jec_param = + jes_delta
  	met_px_shifted = met_px_shifted - ( (jet_pt_nom * (1 + jec_param)) - jet_pt_nom) * cos(recoJet.Phi());
  	met_py_shifted = met_py_shifted - ( (jet_pt_nom * (1 + jec_param)) - jet_pt_nom) * sin(recoJet.Phi());
  }

  // Add this to the systematics vector
  systVec.at(syst)->SetPxPyPzE(met_px_shifted, met_py_shifted, systVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_shifted,2) + pow(met_py_shifted,2)));

}


void Met::propagateJER(TLorentzVector recoJet, double const& jer_sf_nom, double const& jer_sf_shift, int syst){
  
  if(systVec.at(syst) == nullptr) return;

  // update shifted values
  met_px_jerShifted = systVec.at(syst)->Px();
  met_py_jerShifted = systVec.at(syst)->Py();

  // std::cout << "syst = " << syst << ", systVec.at(" << syst << ").Px() = " << systVec.at(syst)->Px() << ", py = " << systVec.at(syst)->Py() << ", pt = " << systVec.at(syst)->Pt() << std::endl;
  // std::cout << "met_px_jerShifted = " << met_px_jerShifted << ". met_px = " << met_px << ", met_py_jerShifted = " << met_py_jerShifted << ", met_py = " << met_py << std::endl;
  double jet_pt_nom = jer_sf_nom * recoJet.Pt();
  // make sure this is positive
  if(jet_pt_nom < 0.0) jet_pt_nom *= -1.0;

  // Check that the jet pt is greater than the thrshold for unclustered energy in order to do the propagation
  if(jet_pt_nom <= jet_unclEnThreshold) return;

  // met_px_nom = met_px_nom - (jet_pt_nom - recoJet.Pt()) * cos(recoJet.Phi());
  // met_px_nom = met_py_nom - (jet_pt_nom - recoJet.Pt()) * sin(recoJet.Phi());

  // update those values we want for the systematic calculated, this also includes the nominal value.
  if(syst == 0){
    met_px_jerShifted = met_px_jerShifted - (jet_pt_nom - recoJet.Pt()) * cos(recoJet.Phi());
    met_py_jerShifted = met_py_jerShifted - (jet_pt_nom - recoJet.Pt()) * sin(recoJet.Phi());
  }
  else{
    met_px_jerShifted = met_px_jerShifted - (jer_sf_shift * recoJet.Pt() - jet_pt_nom) * cos(recoJet.Phi());
    met_py_jerShifted = met_py_jerShifted - (jer_sf_shift * recoJet.Pt() - jet_pt_nom) * sin(recoJet.Phi());
  }
  // Add this to the systematics vector
  systVec.at(syst)->SetPxPyPzE(met_px_jerShifted, met_py_jerShifted, systVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_jerShifted,2) + pow(met_py_jerShifted,2)));
  

  // std::cout << "met_px_jerShifted = " << met_px_jerShifted - (jer_sf_shift * recoJet.Pt() - jet_pt_nom) * cos(recoJet.Phi()) << ". met_px = " << met_px << ", met_py_jerShifted = " << met_py_jerShifted - (jer_sf_shift * recoJet.Pt() - jet_pt_nom) * sin(recoJet.Phi()) << ", met_py = " << met_py << std::endl;
  // std::cout << "syst = " << syst << ", systVec.at(" << syst << ").Px() = " << systVec.at(syst)->Px() << ", py = " << systVec.at(syst)->Py() << ", pt = " << systVec.at(syst)->Pt() << std::endl;

}



void Met::propagateJES(TLorentzVector recoJet, double const& jer_sf_nom, double const& jes_delta, double const& jes_sigma, int syst){
  
  if(systVec.at(syst) == nullptr) return;

  met_px_jesShifted = systVec.at(syst)->Px();
  met_py_jesShifted = systVec.at(syst)->Py();

  double jet_pt_nom = jer_sf_nom * recoJet.Pt();
  if(jet_pt_nom < 0.0) jet_pt_nom *= -1.0;

  // Check that the jet pt is greater than the thrshold for unclustered energy in order to do the propagation
  if(jet_pt_nom <= jet_unclEnThreshold) return;

  // update those values we want for the systematic calculated, this also includes the nominal value.
  met_px_jesShifted = met_px_jesShifted - ( (jet_pt_nom * (1 + jes_sigma*jes_delta)) - jet_pt_nom) * cos(recoJet.Phi());
  met_py_jesShifted = met_py_jesShifted - ( (jet_pt_nom * (1 + jes_sigma*jes_delta)) - jet_pt_nom) * sin(recoJet.Phi());

  // Add this to the systematics vector
  systVec.at(syst)->SetPxPyPzE(met_px_jesShifted, met_py_jesShifted, systVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_jesShifted,2) + pow(met_py_jesShifted,2)));

}



void Met::update(PartStats& stats, Jet& jet, int syst=0){
  ///Calculates met from values from each file plus smearing and treating muons as neutrinos
  if(systVec.at(syst) == nullptr) return;

  //Final update of the MET: systVec.at(0) contains the met_nom. If syst=0 then don't update any further this met.
  if(syst ==0){
    finalmet_px_shifted = systVec.at(0)->Px();
    finalmet_py_shifted = systVec.at(0)->Py();
  }
  else{
    finalmet_px_shifted = systVec.at(syst)->Px() + (systVec.at(0)->Px() - met_px);
    finalmet_py_shifted = systVec.at(syst)->Py() + (systVec.at(0)->Py() - met_py);
  }
  // std::cout << "met_px (update 1) = " << met_px << ". met_py (update 2) = " << met_py << std::endl;
  // std::cout << "systVec.at(" << syst << ")->Px() = " << systVec.at(syst)->Px() << ", systVec.at(" << syst << ")->Py() = " << systVec.at(syst)->Py() << std::endl;
  // std::cout << "systVec.at(0)->Px() = " << systVec.at(0)->Px() << ", systVec.at(0)->Py() = " << systVec.at(0)->Py() << std::endl;
  // std::cout << "finalmet_px_shifted = " << finalmet_px_shifted << ", finalmet_py_shifted = " << finalmet_py_shifted << std::endl;

  // Set the vector of the current systematics with the final values
  systVec.at(syst)->SetPxPyPzE(finalmet_px_shifted, finalmet_py_shifted, systVec.at(syst)->Pz(), TMath::Sqrt(pow(met_px_jesShifted,2) + pow(met_py_jesShifted,2)));

  // std::cout << "2: ===== systVec.at(" << syst << ")->Px() = " << systVec.at(syst)->Px() << ", systVec.at(" << syst << ")->Py() = " << systVec.at(syst)->Py() << std::endl;

/*
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

*/
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
  if( systVec.at(syst) == nullptr) {
    cur_P = systVec.at(0);
    activeSystematic = 0;
  } else {
    cur_P = systVec[syst];
    activeSystematic=syst;
  }
}


void Met::unBranch() {
  BOOM->SetBranchStatus((GenName+"*").c_str(), 0);
}

