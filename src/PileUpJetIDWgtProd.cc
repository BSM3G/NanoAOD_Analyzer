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

  //Create TFile from the efficiency map
  TFile* file_effcyMap = new TFile(effcyMap_rootFile);
  if(!file_effcyMap or file_effcyMap->IsZombie()){
    std::cerr << std::endl << "ERROR! Failed to open input file '" << file_effcyMap << "'!" << std::endl;
	}

  //Get the histogram to pull efficiency values from
  h_PileUpJetIDEffcymap = dynamic_cast<TH2F*>(file_effcyMap->Get(effcyMap));
  h_PileUpJetIDEffcymap->SetDirectory(nullptr);
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
  file_SFMap->Close();
  delete file_SFMap;
}


  //This function will fill the Passing/Failing_Jets_Data/MC vectors which hold the effcy and SF values needed to calculate the weights. The actual calculation happens in producePUJetIDWeights.  
float PileUpJetIDWgtProd::getPUJetIDWeights(Jet& jets, std::vector<int> passing_jets, std::vector<int> failing_jets){
  // First we get values for jets which pass PU ID.
  for(size_t i = 0; i < passing_jets.size(); i++){
  
    // Get the Lorentz vector for the corresponding jet
    TLorentzVector jetP4 = jets.p4(passing_jets[i]);
    float eta_jet = jetP4.Eta();
    float pt_jet = jetP4.Pt();

    // Check that it is in the affected regions
    if(pt_jet < 15.0 || pt_jet > 50.0) continue;
    if(abs(jetP4.Eta()) < -5.0 || abs(jetP4.Eta()) > 5.0) continue;	
    
    //Get actual effcy/SF value for pt & eta value from histogram 
    float JetEffcyValue = getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDEffcymap);
    float JetSFValue =  getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDSFmap);

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

    // Check that it is in the affected regions
    if(pt_jet < 15.0 || pt_jet > 50.0) continue;
    if(abs(jetP4.Eta()) < -5.0 || abs(jetP4.Eta()) > 5.0) continue;	
    
    //Get actual effcy/SF value for pt & eta value from histogram 
    float JetEffcyValue = getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDEffcymap);
    float JetSFValue =  getPileUpEffcyOrSF(eta_jet, pt_jet, h_PileUpJetIDSFmap);

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
