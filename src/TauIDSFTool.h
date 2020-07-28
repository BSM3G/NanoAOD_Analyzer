#ifndef TauIDSFTool_h
#define TauIDSFTool_h

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

class TauIDSFTool {

  protected:

    std::map<const std::string,const TF1*> func;
    TH1* hist;
    [[noreturn]] void disabled() const;

  public:

    std::string ID;
    std::string WP;
    std::vector<int> DMs;
    std::vector<int> genmatches;
    bool isVsPT  = false;
    bool isVsDM  = false;
    bool isVsEta = false;

    TauIDSFTool() {};
    TauIDSFTool(const std::string& datapath, const std::string& year, const std::string& id="MVAoldDM2017v2", const std::string& wp="Tight", const bool dm=false, const bool embedding=false);
    ~TauIDSFTool() { }

    float getSFvsPT( double pt,          int genmatch, const std::string& unc="");
    float getSFvsPT( double pt,                        const std::string& unc="");
    float getSFvsDM( double pt,  int dm, int genmatch, const std::string& unc="") const;
    float getSFvsDM( double pt,  int dm,               const std::string& unc="") const;
    float getSFvsEta(double eta,         int genmatch, const std::string& unc="") const;

};

#endif // TauIDSFTool_h


#ifndef TauESTool_h
#define TauESTool_h

/*
 * @class TauESTool
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

class TauESTool {

  protected:

    std::map<const std::string,const TF1*> func;
    TH1* hist_lowpt;
    TH1* hist_highpt;

  public:

    std::string ID;
    float pt_low;
    float pt_high;
    std::vector<int> DMs;

    TauESTool() { };
    TauESTool(const std::string& datapath, const std::string& year, const std::string& id="MVAoldDM2017v2");
    ~TauESTool() { }

    float getTES(       double pt,  int dm, int genmatch, const std::string& unc="");
    float getTES(       double pt,  int dm,               const std::string& unc="");
    float getTES_highpt(            int dm, int genmatch, const std::string& unc="");
    float getTES_highpt(            int dm,               const std::string& unc="");

};

#endif // TauESTool_h


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
#include <TGraphAsymmErrors.h> // TGraphAsymmErrors
#include <string>    // std::string
#include <vector>    // std::vector
#include <map>       // std::map
#include <stdlib.h>  // getenv
#include <functional>

class TauFESTool {

  protected:

    std::map<const std::string,const TF1*> func;
    TH1* hist;

  public:

    std::string ID;
    float pt_low;
    float pt_high;
    std::vector<int> DMs;
    std::vector<int> genmatches;
    std::map< std::string, std::map<int, std::vector<float> > > FESs;


    TauFESTool() { };
    TauFESTool(const std::string& datapath, const std::string& year, const std::string& id="DeepTau2017v2p1VSjet");
    ~TauFESTool() { }

    float getFES(double eta,  int dm, int genmatch, const std::string& unc="");
    float getFES(double eta,  int dm,               const std::string& unc="");

};

#endif // TauFESTool_h