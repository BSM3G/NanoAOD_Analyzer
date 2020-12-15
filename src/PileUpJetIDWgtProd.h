#ifndef PileUpJetIDWgtProd_h
#define PileUpJetIDWgtProd_h

/* -*- C++ -*-
 *
 * Class: PileUpJetIDWgtProd
 * Based on recipe by Laurent Thomas @ CMS
 * Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
 * Source: 
 * Goal: calculate event weights for Pile Up Jet ID. 
 * Adapted by: Dale Julson.
 * Date: Dec 1, 2020
*/

#include <iostream>
#include <memory>
#include <TFile.h>   // TFile
#include <TString.h> // Form
#include <string>    // std::string
#include <vector>    // std::vector
#include <map>       // std::map
#include <TH2.h>
#include "Particle.h"

class PileUpJetIDWgtProd{
	public:
		PileUpJetIDWgtProd() { };
		PileUpJetIDWgtProd(const std::string& datapath, const std::string& year, const std::string& WP);
		~PileUpJetIDWgtProd() { };

		void getPUJetIDWeights(Jet& passing_jets, Jet& failing_jets);
		//void resetPileUpJetIDWeights();
		float producePUJetIDWeights(std::vector<float> Data_PU_values, std::vector<float> MC_PU_values);
		
	private:

		float getPileUpEffcyOrSF(float eta, float pt, TH2F* h_PUWeightMap);
		TString effcyMap_rootFile;
		TString SFMap_rootFile;
		TString effcyMap;
		TString SFMap;
		TH2F* h_PileUpJetIDEffcymap;
		TH2F* h_PileUpJetIDSFmap;
		std::vector<float> Data_PU_values;
		std::vector<float> MC_PU_values;
		std::string WP_value;
		float PUJetIDWeight;

		//double PileUpJetIDRateSystUnc_;
		
};

#endif
