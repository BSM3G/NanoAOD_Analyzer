#ifndef L1ECALPrefiringWgtProd_h
#define L1ECALPrefiringWgtProd_h

/* -*- C++ -*-
 *
 * Class: L1ECALPrefiringWgtProd
 * Based on recipe by Laurent Thomas @ CMS
 * Twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe 
 * Source: https://github.com/cms-sw/cmssw/blob/8706dbe8a09e7e1314f2127288cfc39051851eea/PhysicsTools/PatUtils/plugins/L1ECALPrefiringWeightProducer.cc
 * Goal: calculate event weights for 2016+2017 to mitigate the L1 ECAL prefiring issue. 
 * Adapted by: Brenda Fabela.
 * Date: Sept 10, 2020
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

enum variations {central = 0, up, down};

class L1ECALPrefiringWgtProd{
	public:
		L1ECALPrefiringWgtProd() { };
		L1ECALPrefiringWgtProd(const std::string& datapath, const std::string& year, const bool useEMpt=false);
		~L1ECALPrefiringWgtProd() { };

		void produceWeights(Photon& photons, Jet& jets);
		void resetWeights();
		float getPrefiringWeight(const std::string& syst_name="");

	private:

		float getPrefiringRate(double eta, double pt, TH2F* h_prefmap, variations var);
		float nonPrefiringProbability[3] = {1.0, 1.0, 1.0}; // 0: central, 1: up, 2: down

		TString photonmapname;
		TString jetmapname;
		TString filename_photonmap;
		TString filename_jetmap;
		TH2F* h_prefmap_photon;
		TH2F* h_prefmap_jet;
		bool useEMpt_;
		std::string dataera_;
		double prefiringRateSystUnc_;

};

#endif
