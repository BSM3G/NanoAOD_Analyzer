#include "PileUpJetIDWgtProd.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>


PileUpJetIDWgtProd::PileUpJetIDWgtProd(const std::string& datapath, const std::string& year){
	
	static std::map<std::string, std::string> yearmap = {
		{"2016", "2016BtoH"},
		{"2017", "2017BtoF"},
		{"2018", "2018AtoD"}
	};

	PileUpJetID_dataera_ = yearmap[year];

	//prefiringRateSystUnc_ = 0.2; // Recommended: https://lathomas.web.cern.ch/lathomas/TSGStuff/L1Prefiring/PrefiringMaps_2016and2017/

	PileUpJetID_jetmapname = Form("PileUpJetID_%s", PileUpJetID_dataera_.data());

	filename_PileUpJetID_jetmap = Form("%s/%s.root", datapath.data(), PileUpJetID_jetmapname.Data());

	TFile* file_PileUpJetID_map = new TFile(filename_PileUpJetID_jetmap);
	if(!file_PileUpJetID_map or file_PileUpJetID_map->IsZombie()){
		std::cerr << std::endl << "ERROR! Failed to open input file '" << file_PileUpJetID_map << "'!" << std::endl;
	}


	h_PileUpJetIDmap_jet = dynamic_cast<TH2F*>(file_PileUpJetID_map->Get(PileUpJetID_jetmapname));
	h_PileUpJetIDmap_jet->SetDirectory(nullptr);
	file_PileUpJetID_map->Close();
	delete file_PileUpJetID_map;

}

void PileUpJetIDWgtProd::produceWeights(Jet& jets){

	// Probability for the event NOT to prefire, computed with the prefiring maps per object.
	// Up and down values correspond to the resulting value when shifting up/down all prefiring rates in prefiring maps

	for(const auto var: {variations::central, variations::up, variations::down}) {

		// Now apply the prefiring maps to jets in the affected regions.
		for(size_t i = 0; i < jets.size(); i++){
			// Get the Lorentz vector for the corresponding jet
			TLorentzVector jetP4 = jets.RecoP4(i);

			float pt_jet = jetP4.Pt();

			// Check that it is in the affected regions
			if(pt_jet < 15.0 || pt_jet > 50.0) continue;
			if(abs(jetP4.Eta()) < -5.0 || abs(jetP4.Eta()) > 5.0) continue;

			//float nonprefiringprobfromoverlappingjet = 1.0 - getPileUpJetIDRate(jetP4.Eta(), pt_jet, h_prefmap_jet, var);
//
//			if(nonprefiringprobfromoverlappingphotons == 1.0){
//				// Case 1: there are no overlapping photons
//				nonPrefiringProbability[var] *= nonprefiringprobfromoverlappingjet;
//			}
//			else if(nonprefiringprobfromoverlappingphotons > nonprefiringprobfromoverlappingjet){
//				// Case 2: if overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
//				if(nonprefiringprobfromoverlappingphotons != 0.){
//					nonPrefiringProbability[var] *= nonprefiringprobfromoverlappingjet / nonprefiringprobfromoverlappingphotons;
//				} else {
//					nonPrefiringProbability[var] = 0.0;
//				}
//			}
//			else{
//				// Case 3: if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight and do nothing.
			}
		}
	}
}

float PileUpJetIDWgtProd::getPileUpJetIDRate(double eta, double pt, TH2F* h_PileUpJetIDmap, variations var){

	if(h_PileUpJetIDmap == nullptr){
		std::cout << "PileUp Jet ID map not found, setting prefiring rate to 0" << std::endl;
		return 0.0;
	}

	// Check pt is not above map overflow
	int nbinsy = h_PileUpJetIDmap->GetNbinsY();
	float maxy = h_PileUpJetIDmap->GetYaxis()->GetBinLowEdge(nbinsy + 1);
	if(pt >= maxy)
		pt = maxy - 0.01;
	int thebin = h_PileUpJetIDmap->FindBin(eta,pt);

	float PileUpJetIDrate = h_PileUpJetIDmap->GetBinContent(thebin);

	float statuncty = h_PileUpJetIDmap->GetBinError(thebin);
	float systuncty = PileUpJetIDRateSystUnc_ * PileUpJetIDrate;

	if(var == up)
		PileUpJetIDrate = std::min(1., PileUpJetIDrate + sqrt(pow(statuncty,2) + pow(systuncty, 2)));
	if(var == down)
		PileUpJetIDrate = std::max(0., PileUpJetIDrate - sqrt(pow(statuncty,2) + pow(systuncty, 2)));

	return PileUpJetIDrate;

}

float PileUpJetIDWgtProd::getPileUpJetIDWeight(const std::string& syst_name){

	if(syst_name == "")
		return nonPrefiringProbability[0];
	else if(syst_name == "Up")
		return nonPrefiringProbability[1];
	else if(syst_name == "Down")
		return nonPrefiringProbability[2];

	return 1.0;
}

void PileUpJetIDWgtProd::resetWeights(){

	for(const auto var: {variations::central, variations::up, variations::down}){
		nonPrefiringProbability[var] = 1.0;
	}

}
