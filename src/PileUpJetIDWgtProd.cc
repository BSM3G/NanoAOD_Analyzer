#include "PileUpJetIDWgtProd.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>


PileUpJetIDWgtProd::PileUpJetIDWgtProd(const std::string& datapath, const std::string& year, const int& WP){
	
  static std::map<int, std::string> WP_map = {
    {0, "T"},
    {1, "M"},
    {2, "L"}
};

  std::string WP_value = WP_map[WP];

  //Set the root file to pull from
  effcyMap_rootFile = Form("%sPUJetIDWeightEffcy.root", datapath.data());
  SFMap_rootFile = Form("%sPUJetIDWeightScaleFactors.root", datapath.data());
  
  //Set the histrogram to pull from within the root file
  effcyMap = Form("h2_eff_mc%s_%s", year.data(), WP_value.data());
  SFMap = Form("h2_eff_sf%s_%s", year.data(), WP_value.data());
  mistag_effcyMap = Form("h2_mistag_mc%s_%s", year.data(), WP_value.data());
  mistag_SFMap = Form("h2_mistag_sf%s_%s", year.data(), WP_value.data());

  //Create TFile from the efficiency map
  TFile* file_effcyMap = new TFile(effcyMap_rootFile);
  if(!file_effcyMap or file_effcyMap->IsZombie()){
    std::cerr << std::endl << "ERROR! Failed to open input file '" << file_effcyMap << "'!" << std::endl;
	}

  //Get the histogram to pull efficiency and mistag values from
  h_PileUpJetIDEffcymap = dynamic_cast<TH2F*>(file_effcyMap->Get(effcyMap));
  h_PileUpJetIDEffcymap->SetDirectory(nullptr);
  h_PileUpJetIDEffcyMistag_map = dynamic_cast<TH2F*>(file_effcyMap->Get(mistag_effcyMap));
  h_PileUpJetIDEffcyMistag_map->SetDirectory(nullptr);
  file_effcyMap->Close();
  delete file_effcyMap;

  //Create TFile from the scale factor map
  TFile* file_SFMap = new TFile(SFMap_rootFile);
  if(!file_SFMap or file_SFMap->IsZombie()){
    std::cerr << std::endl << "ERROR! Failed to open input file '" << file_SFMap << "!" << std::endl;
  }

  //Get the histogram to pull scale factor values from
  h_PileUpJetIDSFmap = dynamic_cast<TH2F*>(file_SFMap->Get(SFMap));
  h_PileUpJetIDSFmap->SetDirectory(nullptr);
  h_PileUpJetIDSFMistag_map = dynamic_cast<TH2F*>(file_SFMap->Get(mistag_SFMap));
  h_PileUpJetIDSFMistag_map->SetDirectory(nullptr);
  file_SFMap->Close();
  delete file_SFMap;
}


  //This function will fill the Passing/Failing_Jets_Data/MC vectors which hold the effcy and SF values needed to calculate the weights. The actual calculation happens in producePUJetIDWeights.  
float PileUpJetIDWgtProd::getPUJetIDWeights(Jet& jets, std::vector<int> passing_jets, std::vector<int> failing_jets){
  // First we get values for jets which pass PU ID.
  std::cout << "jet collection size is " << jets.size() << std::endl; 
  std::cout << "passing jet collection size is " << passing_jets.size() << std::endl;
  std::cout << "failing jet collection size is " << failing_jets.size() << std::endl;
  
  for(size_t i = 0; i < passing_jets.size(); i++){
  
    // Get the Lorentz vector for the corresponding jet
    TLorentzVector jetP4 = jets.p4(passing_jets[i]);
    float eta_jet = jetP4.Eta();
    float pt_jet = jetP4.Pt();
    int matchedGenJetIndex = jets.genJetIdx[passing_jets[i]];
    std::cout << "Passing Jet " << i << " has pt " << pt_jet << std::endl;
    std::cout << "Passing Jet "<< i <<" has eta " << eta_jet << std::endl;
    std::cout << "Passing matchedGenJetIndex is: " << matchedGenJetIndex << std::endl;

    // Check that it is in the affected regions
    if(pt_jet < 15.0 || pt_jet > 50.0) continue;
    if(abs(jetP4.Eta()) < -5.0 || abs(jetP4.Eta()) > 5.0) continue;	

    if (matchedGenJetIndex < 0){ //These are real PU Jets which passed PU Jet ID & were therefore "mistagged" (the naming convention here is confusing), therefore they pull from the mistag histograms.
    JetEffcyValue = getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDEffcyMistag_map);
    JetSFValue =  getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDSFMistag_map); 
    }
    
    else{//These are real hard scatter jets which were correctly identified.
    JetEffcyValue = getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDEffcymap);
    JetSFValue =  getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDSFmap);
	}

    std::cout << "Passing Jet " << i << " has Effcy value " << JetEffcyValue << std::endl;
    std::cout << "Passing Jet "<< i <<" has SF value " << JetSFValue << std::endl;

    //Push values into appropriate vectors.
    MC_PU_values.push_back(JetEffcyValue);
    Data_PU_values.push_back(JetSFValue*JetEffcyValue);
    	}

// Now we do the same but for jets which fail PU ID.
  for(size_t i = 0; i < failing_jets.size(); i++){
  
    // Get the Lorentz vector for the corresponding jet
    TLorentzVector jetP4 = jets.p4(failing_jets[i]);
    float eta_jet = jetP4.Eta();
    float pt_jet = jetP4.Pt();
    int matchedGenJetIndex = jets.genJetIdx[failing_jets[i]];
    std::cout << "Failing Jet " << i << " has pt " << pt_jet << std::endl;
    std::cout << "Failing Jet "<< i <<" has eta " << eta_jet << std::endl;
    std::cout << "Failing matchedGenJetIndex is: " << matchedGenJetIndex << std::endl;

    // Check that it is in the affected regions
    if(pt_jet < 15.0 || pt_jet > 50.0) continue;
    if(abs(jetP4.Eta()) < -5.0 || abs(jetP4.Eta()) > 5.0) continue;	
    
    if (matchedGenJetIndex < 0){//Genuine PU Jets that fail PU Jet ID, they get entered as 1-mistag_efficiency.
    JetEffcyValue = getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDEffcyMistag_map);
    JetSFValue =  getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDSFMistag_map);
    }

    else{//Hard scatter jets which failed PU ID, these get entered as 1-efficiency.
    JetEffcyValue = getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDEffcymap);
    JetSFValue =  getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDSFmap);
    }

    std::cout << "Failing Jet " << i << " has Effcy value " << JetEffcyValue << std::endl;
    std::cout << "Failing Jet "<< i <<" has SF value " << JetSFValue << std::endl;

    //Push values into appropriate vectors.
    MC_PU_values.push_back(1.0-JetEffcyValue);
    Data_PU_values.push_back(1.0-(JetSFValue*JetEffcyValue));
	}
  
  float final_weight = producePUJetIDWeights(Data_PU_values, MC_PU_values); 

  MC_PU_values.clear();
  Data_PU_values.clear();

  return final_weight;

}


float PileUpJetIDWgtProd::getPileUpEffcyOrSF(float eta, float pt, TH2F* h_PUWeightMap){
  //This function will open h_PUWeightmap and pull the appropriate efficiency.

  if(h_PUWeightMap == nullptr){
    std::cout << "Prefiring map not found, setting prefiring rate to 0" << std::endl;
    return 0.0;
  }

  int thebin = h_PUWeightMap->FindBin(pt,eta);
  float bin_value = h_PUWeightMap->GetBinContent(thebin);
  
  return bin_value;

}

//This function will combine the effcy/SF values held in Data_PU_values & MC_PU_values and produce the PUJetIDweight
float PileUpJetIDWgtProd::producePUJetIDWeights(std::vector<float> Data_PU_values, std::vector<float> MC_PU_values){

   if(!(Data_PU_values.size() == MC_PU_values.size())){
    std::cerr << std::endl << "ERROR! Pile-Up Jet vectors not same length!" << std::endl;
    }

   float Prob_Data = 1.0;
   float Prob_MC = 1.0;

   for(size_t i = 0; i < Data_PU_values.size(); i++){
     
     Prob_Data *= Data_PU_values[i];
     Prob_MC *= MC_PU_values[i];
   }
   
   float PUJetIDWeight = Prob_Data/Prob_MC;

   return PUJetIDWeight;
 }
