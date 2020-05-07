// https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py
// https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/JetReCalibrator.py
// https://github.com/yhaddad/nanoAOD-tools/blob/4d7f94fe7f4d319a760a79c9548a05d017885fb4/python/postprocessing/modules/jme/JetReCalibrator.py
#include "JetRecalibrator.h"
#include <sstream>
#include <cmath>
#include <algorithm>

JetRecalibrator::JetRecalibrator(){};

JetRecalibrator::JetRecalibrator(const std::string& globalTag, const std::string& jetFlavor, bool& doResidualJECs, const std::string& jecPathFile, const std::string& jecPathFile, int& upToLevel, 
			bool& calculateSeparateCorrections, bool& calculateTypeIMETCorr, double& jetPtThreshold, double& skipEMfractionThreshold, 
			bool& skipMuons){

	// Create a corrector object that applies the L1,L2,L3 and possibly the residual corrections to the jets.
	// If configured to do so, it will also compute the type1 MET corrections

	// Make base corrections
	L1JetParameters = JetCorrectorParameters((path+"/"+globalTag+"_L1FastJet_"+jetFlavor+".txt").c_str(), "");
	L2JetParameters = JetCorrectorParameters((path+"/"+globalTag+"_L2Relative_"+jetFlavor+".txt").c_str(), "");
	L3JetParameters = JetCorrectorParameters((path+"/"+globalTag+"_L3Absolute_"+jetFlavor+".txt").c_str(), "");

	vPar.push_back(L1JetParameters);

	if(upToLevel >= 2) vPar.push_back(L2JetParameters);
	if(upToLevel >= 3) vPar.push_back(L3JetParameters);

	// Add residuals if needed
	if(doResidualJECs){
		ResidualJetParameters = JetCorrectorParameters((path+"/"+globalTag+"_L2L3Residual_"+jetFlavor+".txt").c_str());
		vPar.push_back(ResidualJetParameters);
	}

	// Construct a FactorizedJetCorrector object
	JetCorrector = FactorizedJetCorrector(vPar);

	std::ifstream uncty(path+"/"+globalTag+"_Uncertainty_"+jetFlavor+".txt");
	std::ifstream unctyfake(path+"/Uncertainty_FAKE.txt");
	if(uncty.good()) {
		JetUncertainty = JetCorrectionUncertainty(path+"/"+globalTag+"_Uncertainty_"+jetFlavor+".txt");
	}
	else if (unctyfake.good()){
		JetUncertainty = JetCorrectionUncertainty(path+"/Uncertainty_FAKE.txt");
	}
	else{
		std::cout << "Missing JEC uncertainty file " << path << "/" << globalTag << "_Uncertainty_" << jetFlavor << ".txt, so jet energy uncertainties will not be available." << std::endl;
		nouncertainty = true;
	}

	if(calculateSeparateCorrections || calculateTypeIMETCorr){

		vParL1.push_back(L1JetParameters);
		separateJetCorrectors["L1"] = FactorizedJetCorrector(vParL1);

		if(upToLevel >= 2 && calculateSeparateCorrections){
			vParL2.push_back(L1JetParameters);
			vParL2.push_back(L2JetParameters);

			separateJetCorrectors["L1L2"] = FactorizedJetCorrector(vParL2);
		}
		if(upToLevel >= 3 && calculateSeparateCorrections){
			vParL3.push_back(L1JetParameters);
			vParL3.push_back(L2JetParameters);
			vParL3.push_back(L3JetParameters);

			separateJetCorrectors["L1L2L3"] = FactorizedJetCorrector(vParL3);
		}
		if(doResidualJECs && calculateSeparateCorrections){
			vParL3Residuals.push_back(L1JetParameters);
			vParL3Residuals.push_back(L2JetParameters);
			vParL3Residuals.push_back(L3JetParameters);
			vParL3Residuals.push_back(ResidualJetParameters);

			separateJetCorrectors["L1L2L3Res"] = FactorizedJetCorrector(vParL3Residuals);
		}
	}

};

float JetRecalibrator::getCorrection(TLorentzVector& jet4vec, float& jet_area, float& jet_rawFactor, float& rho, float& delta){

	JetCorrector.setJetEta(jet4vec.Eta());
	JetCorrector.setJetPt(jet4vec.Pt() * (1 - jet_rawFactor));
	JetCorrector.setJetA(jet_area);
	JetCorrector.setRho(rho);

	float corr = JetCorrector.getCorrection();
	float jetEnergyCorrUncertainty = 1.0;

	if(delta != 0){
		if(nouncertainty){
			std::cout << "Jet energy scale uncertainty shifts requested, but not available." << std::endl;
		}
		else{
			JetUncertainty.setJetEta(jet4vec.Eta());
			JetUncertainty.setJetPt(corr * jet4vec.Pt() * jet_rawFactor);

			try{
				jetEnergyCorrUncertainty = JetUncertainty.getUncertainty(true);
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

TLorentzVector JetRecalibrator::correctedP4(TLorentzVector& jet4vec, float& corr_factor, float& jet_rawFactor){

	double raw = 1.0 - jet_rawFactor;
	if(corr_factor <= 0.0) return jet4vec;

	double newpt = jet4vec.Pt() * raw * corr_factor;
	double newmass = jet4vec.M() * raw * corr_factor;

	TLorentzVector correctedJet4Vec;
	correctedJet4Vec.SetPtEtaPhiM(newpt, jet4vec.Eta(), jet4vec.Phi(), newmass);

	return correctedJet4Vec;

}