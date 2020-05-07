#ifndef JETRECALIBRATOR
#define JETRECALIBRATOR

#include <unordered_map>
#include <TLorentzVector.h>
#include <string>
#include <TH1.h>
#include <TGraph.h>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>

// Include files from CMSSW
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/Utilities.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"

class JetRecalibrator {

	public:
		JetRecalibrator();
		JetRecalibrator(const std::string& path, const std::string& globalTag, const std::string& jetFlavor, bool& doResidualJECs, const std::string& jecPathFile, int& upToLevel=3, 
			bool& calculateSeparateCorrections=false, bool& calculateTypeIMETCorr=false, double& jetPtThreshold=15.0, double& skipEMfractionThreshold=0.9, 
			bool& skipMuons=true);

		float getCorrection(TLorentzVector& jet4vec, float& jet_area, float& jet_rawFactor, float& rho, float& delta = 0.0);
		TLorentzVector correctedP4(TLorentzVector& jet4vec, float& corr_factor, float& jet_rawFactor);

		FactorizedJetCorrector JetCorrector;
		JetCorrectionUncertainty JetUncertainty;
		bool nouncertainty = false;

	private:
		JetCorrectorParameters L1JetParameters;
		JetCorrectorParameters L2JetParameters;
		JetCorrectorParameters L3JetParameters;
		JetCorrectorParameters ResidualJetParameters;

		std::vector <JetCorrectorParameters> vPar;
		std::vector <JetCorrectorParameters> vParL1;
		std::vector <JetCorrectorParameters> vParL2;
		std::vector <JetCorrectorParameters> vParL3;
		std::vector <JetCorrectorParameters> vParL3Residuals;

		std::map<std::string, FactorizedJetCorrector> separateJetCorrectors;

};

#endif