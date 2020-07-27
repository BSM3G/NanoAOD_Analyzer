#ifndef TauFESTool_h
#define TauFESTool_h

/*
 * @class TauIDSFTool
 *
 * Class to retrieve tau ID SFs.
 *  - pT-dependent SFs for MVAoldDM2017v2
 *  - DM-dependent SFs for MVAoldDM2017v2
 *  - eta-dependent SFs for anti-lepton discriminators
 * Source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
 * Inspiration from TauTriggerSFs/src/TauTriggerSFs2017.cc
 *
 * @author Izaak Neutelings
 * @date July 2019
 *
 */

#include <TFile.h>   // TFile
#include <TH1.h>     // TH1
#include <TF1.h>     // TF1
#include <TString.h> // Form
#include <string>    // std::string
#include <vector>    // std::vector
#include <map>       // std::map
#include <stdlib.h>  // getenv
#include <functional>

class TauFESTool {

  protected:

    std::map<const std::string,const TF1*> func;
    TH1* hist;
    [[noreturn]] void disabled() const;

  public:

    std::string ID;
    float pt_low;
    float pt_high;
    std::vector<int> DMs;
    std::vector<int> genmatches;
    std::map< std::string, std::map<int, std::vector<float> > > FESs;



    TauFESTool(const std::string& datapath, const std::string& year, const std::string& id="MVAoldDM2017v2");
    ~TauFESTool() { }

    float getFES(double eta,  int dm, int genmatch, const std::string& unc="");
    float getFES(double eta,  int dm,               const std::string& unc="");

};

#endif