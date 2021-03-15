#include "JetScaleResolution.h"
#include <sstream>
#include <cmath>



JetScaleResolution::JetScaleResolution(){
}

JetScaleResolution::JetScaleResolution(const std::string& scalefilename, const std::string& parttype, const std::string& resolutionfile, const std::string& sfresolutionfile){
    // InitScale(scalefilename, parttype);
    // InitResolution(resolutionfile, sfresolutionfile);
    /*
    resfilename = resolutionfile;
    resfunctyfilename = sfresolutionfile;
    scaleunctyfilename = scalefilename;
    type_name = parttype;
    */
    jetcorrectorparams = JetCorrectorParameters(scalefilename.c_str(),parttype.c_str());
    jesUncertainty = new JetCorrectionUncertainty(jetcorrectorparams);

    jer = JME::JetResolution(resolutionfile.c_str());
    jerSF_and_Uncertainty = JME::JetResolutionScaleFactor(sfresolutionfile.c_str());
}

double JetScaleResolution::GetScaleDelta(const double& recojetPt, const double& recojetEta){
    if(!(recojetPt > 0.0)){
        return recojetPt;
    }

    //Call jesuncertainty with Pt being that obtained after applying nominal jet energy resolution scale factor
    // jet_pt_nom = jer_sf_nom * recojet.Pt();
    jesUncertainty->setJetPt(recojetPt);
    jesUncertainty->setJetEta(recojetEta);

    double delta = jesUncertainty->getUncertainty(true);

    return delta;
}

// double JetScaleResolution::GetSmearValsPtSF(const TLorentzVector& recojet, const TLorentzVector& genjet, double rho, double sigmares){
std::vector<float> JetScaleResolution::GetSmearValsPtSF(const TLorentzVector& recojet, const TLorentzVector& genjet, double rho){


    if(!(recojet.Pt() > 0.0)){
        //return recojet.Pt();
        std::vector<float> nosmearvalues(3, 1.0);
        return nosmearvalues;
    }

    JME::JetParameters params_sf_and_uncertainty;
    JME::JetParameters params_resolution;

    // CV: define enums to access JER scale factors and uncertainties (cf. CondFormats/JetMETObjects/interface/JetResolutionObject.h)
    // int index_nominal = 0, index_shift_down = 1, index_shift_up = 2;

    //float jet_pt_sf_and_uncertainty[3] = { };
    std::vector<float> jet_pt_sf_and_uncertainty;

    for(int index = 0; index < 3; index++){
        Variation var = (Variation) index;
        params_sf_and_uncertainty.setJetEta(recojet.Eta());
        params_sf_and_uncertainty.setJetPt(recojet.Pt()); // new, june 11, 2020.
        jet_pt_sf_and_uncertainty.push_back(jerSF_and_Uncertainty.getScaleFactor(params_sf_and_uncertainty, var));
    }

    //float smear_values[3] = { };
    std::vector<float> smear_values;

    if(genjet != TLorentzVector(0.,0.,0.,0.)){
        for(int index = 0; index < 3; index++){
            double smearFactor = 0.0;

            // Case 1: we have a "good" generator level jet matched to the reconstructed jet

            double deltaPt = recojet.Pt() - genjet.Pt();
            smearFactor = 1.0 + (jet_pt_sf_and_uncertainty.at(index) - 1.0) * (deltaPt / recojet.Pt());

            // check that the smeared jet energy remains positive
            // as the direction of the jet would change ("flip") otherwise - and this is not what we want

            if(smearFactor * recojet.Pt() < 1.0e-2){
                smearFactor = 1.0e-2;
            }

            //smear_values[index] = smearFactor;
            smear_values.push_back(smearFactor);
        }
    }
    else{
        // Call params_resolution:
        params_resolution.setJetPt(recojet.Pt());
        params_resolution.setJetEta(recojet.Eta());
        params_resolution.setRho(rho);

        float jet_pt_resolution = jer.getResolution(params_resolution);
        double rand = rnd->Gaus(0.0, jet_pt_resolution);

        for(int index = 0; index < 3; index++){

            double smearFactor = 0.0;

            if(jet_pt_sf_and_uncertainty[index] > 1.0){
                // Case 2: we don't have a generator level jet. Smear jet pT using a random Gaussian variation
                smearFactor =  1.0 + rand * sqrt( pow(jet_pt_sf_and_uncertainty.at(index), 2) - 1.0);
            }
            else{
                // Case 3: we cannot smear this jet, as we don't have a generator level jet and the resolution in data is better than the resolution in the simulation,
                // so we would need to randomly "unsmear" the jet, which is impossible
                smearFactor = 1.0;
            }

            // Check that the smeared jet energy remains positive, as the direction of the jet would change ("flip") otherwise, and this is not what we want
            if(smearFactor * recojet.Pt() < 1.0e-2){
                smearFactor = 1.0e-2;
            }

            smear_values.push_back(smearFactor);

        }
    }

    return smear_values;
}

bool JetScaleResolution::resolution_matching(const TLorentzVector& recojet, const TLorentzVector& genjet, double rho){

  JME::JetParameters params_resolution;

  params_resolution.setJetPt(recojet.Pt());
  params_resolution.setJetEta(recojet.Eta());
  params_resolution.setRho(rho);

  double resolution = jer.getResolution(params_resolution);

  return abs(recojet.Pt() - genjet.Pt()) < 3.0 * resolution * recojet.Pt();
}

std::vector<std::string> string_split(const std::string& in, const std::vector<std::string> splits)
{
    std::vector<std::pair<size_t, size_t> > positions;
    positions.push_back(std::pair<size_t, size_t>(0, 0));
    for(size_t s = 0 ; s < splits.size() ; ++s)
    {
        size_t lastpos = 0;
        while(lastpos < in.size())
        {
            lastpos = in.find(splits[s], lastpos);
            if(lastpos == std::string::npos)
            {
                break;
            }
            else
            {
                positions.push_back(std::pair<size_t, size_t>(lastpos, splits[s].size()));
                //lastpos += splits[s].size()+1;
                lastpos += splits[s].size();
            }
        }

    }
    positions.push_back(std::pair<size_t, size_t>(in.size(), 0));
    sort(positions.begin(), positions.end(), [](const std::pair<size_t, size_t>& A, const std::pair<size_t, size_t>& B){return A.first < B.first;});
    std::vector<std::string> result;
    for(size_t p = 0 ; p < positions.size()-1 ; ++p)
    {
        size_t begin = positions[p].first + positions[p].second;
        size_t end = positions[p+1].first;
        if(end != begin)result.push_back(in.substr(begin, end-begin));
    }
    return result;
}
