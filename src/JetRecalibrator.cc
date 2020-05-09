// https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
// https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/JetReCalibrator.py
// https://github.com/yhaddad/nanoAOD-tools/blob/4d7f94fe7f4d319a760a79c9548a05d017885fb4/python/postprocessing/modules/jme/JetReCalibrator.py
#include "JetRecalibrator.h"
#include <sstream>
#include <cmath>
#include <algorithm>

JetRecalibrator::JetRecalibrator(){};

JetRecalibrator::JetRecalibrator(const std::string path, const std::string globalTag, const std::string jetFlavor, bool doResidualJECs, int upToLevel, 
			bool calculateSeparateCorrections, bool calculateTypeIMETCorr){

	path_str = path;
	globalTag_str = globalTag;
	jetFlavor_str = jetFlavor;
	doResJECs = doResidualJECs;
	corrLevel = upToLevel;
	separateCorr = calculateSeparateCorrections;
	calcT1MetCorr = calculateTypeIMETCorr;


};

float JetRecalibrator::getCorrection(TLorentzVector jet4vec, float jet_area, float jet_rawFactor, float rho, float delta){

	// Create a corrector object that applies the L1,L2,L3 and possibly the residual corrections to the jets.
	// If configured to do so, it will also compute the type1 MET corrections

	FactorizedJetCorrector jetcorrector;
	JetCorrectionUncertainty jetuncertainty;

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

	// Make base corrections
	L1JetParameters = JetCorrectorParameters((path_str+globalTag_str+"_L1FastJet_"+jetFlavor_str+".txt").c_str(), "");
	L2JetParameters = JetCorrectorParameters((path_str+globalTag_str+"_L2Relative_"+jetFlavor_str+".txt").c_str(), "");
	L3JetParameters = JetCorrectorParameters((path_str+globalTag_str+"_L3Absolute_"+jetFlavor_str+".txt").c_str(), "");

	vPar.push_back(L1JetParameters);

	if(corrLevel >= 2) vPar.push_back(L2JetParameters);
	if(corrLevel >= 3) vPar.push_back(L3JetParameters);

	// Add residuals if needed
	if(doResJECs){
		ResidualJetParameters = JetCorrectorParameters((path_str+globalTag_str+"_L2L3Residual_"+jetFlavor_str+".txt").c_str());
		vPar.push_back(ResidualJetParameters);
	}

	// Construct a FactorizedJetCorrector object
	FactorizedJetCorrector* jetcorrector = new FactorizedJetCorrector(vPar);

	std::ifstream uncty(path_str+globalTag_str+"_Uncertainty_"+jetFlavor_str+".txt");
	std::ifstream unctyfake(path_str+"Uncertainty_FAKE.txt");
	if(uncty.good()) {
		jetuncertainty = JetCorrectionUncertainty(path_str+globalTag_str+"_Uncertainty_"+jetFlavor_str+".txt");
	}
	else if (unctyfake.good()){
		jetuncertainty = JetCorrectionUncertainty(path_str+"/Uncertainty_FAKE.txt");
	}
	else{
		std::cout << "Missing JEC uncertainty file " << (path_str+globalTag_str).c_str() << "_Uncertainty_" << jetFlavor_str << ".txt, so jet energy uncertainties will not be available." << std::endl;
		nouncertainty = true;
	}

	if(separateCorr || calcT1MetCorr){

		vParL1.push_back(L1JetParameters);
		separateJetCorrectors["L1"] = FactorizedJetCorrector(vParL1);

		if(corrLevel >= 2 && separateCorr){
			vParL2.push_back(L1JetParameters);
			vParL2.push_back(L2JetParameters);

			separateJetCorrectors["L1L2"] = FactorizedJetCorrector(vParL2);
		}
		if(corrLevel >= 3 && separateCorr){
			vParL3.push_back(L1JetParameters);
			vParL3.push_back(L2JetParameters);
			vParL3.push_back(L3JetParameters);

			separateJetCorrectors["L1L2L3"] = FactorizedJetCorrector(vParL3);
		}
		if(doResJECs && separateCorr){
			vParL3Residuals.push_back(L1JetParameters);
			vParL3Residuals.push_back(L2JetParameters);
			vParL3Residuals.push_back(L3JetParameters);
			vParL3Residuals.push_back(ResidualJetParameters);

			separateJetCorrectors["L1L2L3Res"] = FactorizedJetCorrector(vParL3Residuals);
		}
	}

	jetcorrector.setJetEta(jet4vec.Eta());
	jetcorrector.setJetPt(jet4vec.Pt() * (1 - jet_rawFactor));
	jetcorrector.setJetA(jet_area);
	jetcorrector.setRho(rho);

	float corr = jetcorrector.getCorrection();
	float jetEnergyCorrUncertainty = 1.0;

	if(delta != 0){
		if(nouncertainty){
			std::cout << "Jet energy scale uncertainty shifts requested, but not available." << std::endl;
		}
		else{
			jetuncertainty.setJetEta(jet4vec.Eta());
			jetuncertainty.setJetPt(corr * jet4vec.Pt() * jet_rawFactor);

			try{
				jetEnergyCorrUncertainty = jetuncertainty.getUncertainty(true);
			}
			catch(std::runtime_error& err){
				std::cout << "Caught " << err.what() << " when getting uncertainty for jet of pt = %.1f" << corr * jet4vec.Pt() * jet_rawFactor << ", eta = %.2f" << jet4vec.Eta() << std::endl;
				jetEnergyCorrUncertainty = 0.5;
			}
			corr *= std::max(0, 1.0+delta+jetEnergyCorrUncertainty);
		}
	}

	return corr;
}

TLorentzVector JetRecalibrator::correctedP4(TLorentzVector jet4vec, float corr_factor, float jet_rawFactor){

	double raw = 1.0 - jet_rawFactor;
	if(corr_factor <= 0.0) return jet4vec;

	double newpt = jet4vec.Pt() * raw * corr_factor;
	double newmass = jet4vec.M() * raw * corr_factor;

	TLorentzVector correctedJet4Vec;
	correctedJet4Vec.SetPtEtaPhiM(newpt, jet4vec.Eta(), jet4vec.Phi(), newmass);

	return correctedJet4Vec;

}