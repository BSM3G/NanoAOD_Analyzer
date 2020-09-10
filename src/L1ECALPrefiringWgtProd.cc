#include "L1ECALPrefiringWgtProd.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>


L1ECALPrefiringWgtProd::L1ECALPrefiringWgtProd(const std::string& datapath, const std::string& year, const bool useEMpt){
	
	static std::map<std::string, std::string> yearmap = {
		{"2016", "2016BtoH"},
		{"2017", "2017BtoF"}
	};

	dataera_ = yearmap[year];

	useEMpt_ = useEMpt;
	prefiringRateSystUnc_ = 0.2; // Recommended: https://lathomas.web.cern.ch/lathomas/TSGStuff/L1Prefiring/PrefiringMaps_2016and2017/

	photonmapname = Form("L1prefiring_photonpt_%s", dataera_.data());
	jetmapname = (useEMpt_) ? Form("L1prefiring_jetpt_%s", dataera_.data()) : Form("L1prefiring_jetempt_%s", dataera_.data());

	filename_photonmap = Form("%s/%s.root", datapath.data(), photonmapname.Data());
	filename_jetmap = Form("%s/%s.root", datapath.data(), jetmapname.Data());

	TFile* file_photonprefiringmap = new TFile(filename_photonmap);

	if(!file_photonprefiringmap or file_photonprefiringmap->IsZombie()){
		std::cerr << std::endl << "ERROR! Failed to open input file '" << file_photonprefiringmap << "'!" << std::endl;
	}

	TFile* file_jetprefiringmap = new TFile(filename_jetmap);
	if(!file_jetprefiringmap or file_jetprefiringmap->IsZombie()){
		std::cerr << std::endl << "ERROR! Failed to open input file '" << file_photonprefiringmap << "'!" << std::endl;
	}

	h_prefmap_photon = dynamic_cast<TH2F*>(file_photonprefiringmap->Get(photonmapname));
	h_prefmap_photon->SetDirectory(nullptr);
	file_photonprefiringmap->Close();
	delete file_photonprefiringmap;

	h_prefmap_jet = dynamic_cast<TH2F*>(file_jetprefiringmap->Get(jetmapname));
	h_prefmap_jet->SetDirectory(nullptr);
	file_jetprefiringmap->Close();
	delete file_jetprefiringmap;

}

void L1ECALPrefiringWgtProd::produceWeights(Photon& photons, Jet& jets){

	// Probability for the event NOT to prefire, computed with the prefiring maps per object.
	// Up and down values correspond to the resulting value when shifting up/down all prefiring rates in prefiring maps

	for(const auto var: {variations::central, variations::up, variations::down}) {

		// Start looping over photons
		for(size_t i = 0; i < photons.size(); i++){

			// Get Lorentz vector of the corresponding photon
			TLorentzVector photonP4 = photons.RecoP4(i);

			// Check that it is in the affected region of pt and eta
			if(photonP4.Pt() < 20.0) continue;
			if(abs(photonP4.Eta()) < 2.0 || abs(photonP4.Eta()) > 3.0) continue;

			// get the prefiring rate
			float prefiringprob_gamma = getPrefiringRate(photonP4.Eta(), photonP4.Pt(), h_prefmap_photon, var);
			nonPrefiringProbability[var] *= (1.0 - prefiringprob_gamma);
		}

		// Now apply the prefiring maps to jets in the affected regions.
		for(size_t i = 0; i < jets.size(); i++){
			// Get the Lorentz vector for the corresponding jet
			TLorentzVector jetP4 = jets.RecoP4(i);

			float pt_jet = jetP4.Pt();

			// Check that it is in the affected regions
			if(pt_jet < 20.0) continue;
			if(abs(jetP4.Eta()) < 2.0 || abs(jetP4.Eta()) > 3.0) continue;

			// Loop over photons to remove overlap
			float nonprefiringprobfromoverlappingphotons = 1.0;

			for(size_t j = 0; j < photons.size(); j++){

				// Get Lorentz vector of the corresponding photon
				TLorentzVector photon_forjetP4 = photons.RecoP4(j);

				// Check that it is in the affected region of pt and eta
				if(photon_forjetP4.Pt() < 20.0) continue;
				if(abs(photon_forjetP4.Eta()) < 2.0 || abs(photon_forjetP4.Eta()) > 3.0) continue;

				// Check if there is an overlap, if not, skip to the next photon.
				float photonjet_deltaR = jetP4.DeltaR(photon_forjetP4);
				if(photonjet_deltaR > 0.4) continue;

				float prefiringprob_gamma_overlap = getPrefiringRate(photon_forjetP4.Eta(), photon_forjetP4.Pt(), h_prefmap_photon, var);
				nonprefiringprobfromoverlappingphotons *= (1.0 - prefiringprob_gamma_overlap);
			}

			// useEMpt = true if one wants to use maps parametrized vs Jet EM pt instead of pt
			if(useEMpt_) 
				pt_jet *= (jets.neutralEmEnergyFraction[i] + jets.chargedEmEnergyFraction[i]);

			float nonprefiringprobfromoverlappingjet = 1.0 - getPrefiringRate(jetP4.Eta(), pt_jet, h_prefmap_jet, var);

			if(nonprefiringprobfromoverlappingphotons == 1.0){
				// Case 1: there are no overlapping photons
				nonPrefiringProbability[var] *= nonprefiringprobfromoverlappingjet;
			}
			else if(nonprefiringprobfromoverlappingphotons > nonprefiringprobfromoverlappingjet){
				// Case 2: if overlapping photons have a non prefiring rate larger than the jet, then replace these weights by the jet one
				if(nonprefiringprobfromoverlappingphotons != 0.){
					nonPrefiringProbability[var] *= nonprefiringprobfromoverlappingjet / nonprefiringprobfromoverlappingphotons;
				} else {
					nonPrefiringProbability[var] = 0.0;
				}
			}
			else{
				// Case 3: if overlapping photons have a non prefiring rate smaller than the jet, don't consider the jet in the event weight and do nothing.
			}
		}
	}
}

float L1ECALPrefiringWgtProd::getPrefiringRate(double eta, double pt, TH2F* h_prefmap, variations var){

	if(h_prefmap == nullptr){
		std::cout << "Prefiring map not found, setting prefiring rate to 0" << std::endl;
		return 0.0;
	}

	// Check pt is not above map overflow
	int nbinsy = h_prefmap->GetNbinsY();
	float maxy = h_prefmap->GetYaxis()->GetBinLowEdge(nbinsy + 1);
	if(pt >= maxy)
		pt = maxy - 0.01;
	int thebin = h_prefmap->FindBin(eta,pt);

	float prefrate = h_prefmap->GetBinContent(thebin);

	float statuncty = h_prefmap->GetBinError(thebin);
	float systuncty = prefiringRateSystUnc_ * prefrate;

	if(var == up)
		prefrate = std::min(1., prefrate + sqrt(pow(statuncty,2) + pow(systuncty, 2)));
	if(var == down)
		prefrate = std::max(0., prefrate - sqrt(pow(statuncty,2) + pow(systuncty, 2)));

	return prefrate;

}

float L1ECALPrefiringWgtProd::getPrefiringWeight(const std::string& syst_name){

	if(syst_name == "")
		return nonPrefiringProbability[0];
	else if(syst_name == "Up")
		return nonPrefiringProbability[1];
	else if(syst_name == "Down")
		return nonPrefiringProbability[2];

	return 1.0;
}

void L1ECALPrefiringWgtProd::resetWeights(){

	for(const auto var: {variations::central, variations::up, variations::down}){
		nonPrefiringProbability[var] = 1.0;
	}

}