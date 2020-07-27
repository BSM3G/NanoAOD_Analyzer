#include "TauFESTool.h"
#include <iostream> // std::cerr, std::endl
#include <iomanip>
#include <assert.h> // assert



TFile* ensureTFile(const TString filename, bool verbose=false){
  if(verbose)
    std::cout << "Opening " << filename << std::endl;
  TFile* file = new TFile(filename);
  if(!file or file->IsZombie()) {
    std::cerr << std::endl << "ERROR! Failed to open input file = '" << filename << "'!" << std::endl;
    assert(0);
  }
  return file;
}

TH1* extractTH1(const TFile* file, const std::string& histname){
  TH1* hist = dynamic_cast<TH1*>((const_cast<TFile*>(file))->Get(histname.data()));
  if(!hist){
    std::cerr << std::endl << "ERROR! Failed to load histogram = '" << histname << "' from input file!" << std::endl;
    assert(0);
  }
  return hist;
}

const TF1* extractTF1(const TFile* file, const std::string& funcname){
  const TF1* function = dynamic_cast<TF1*>((const_cast<TFile*>(file))->Get(funcname.data()));
  if(!function){
    std::cerr << std::endl << "ERROR! Failed to load function = '" << funcname << "' from input file!" << std::endl;
    assert(0);
  }
  return function;
}



void TauFESTool::disabled() const {
  std::cerr << std::endl << "ERROR! Method has been disabled! isVsPT = "<<isVsPT<<", isVsDM = "
            << isVsDM<<", isVsEta = "<<isVsEta<< std::endl;
  assert(0);
}


TauFESTool::TauFESTool(const std::string& datapath, const std::string& year, const std::string& id): ID(id){

  bool verbose = false;
  //std::string datapath                = Form("%s/src/TauPOG/TauIDSFs/data",getenv("CMSSW_BASE"));
  std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco"};
  std::vector<std::string> antiJetIDs = {"MVAoldDM2017v2","DeepTau2017v2p1VSjet"};
  std::vector<std::string> antiEleIDs = {"antiEleMVA6",   "DeepTau2017v2p1VSe"};
  std::vector<std::string> antiMuIDs  = {"antiMu3",       "DeepTau2017v2p1VSmu"};

  if(std::find(years.begin(),years.end(),year)==years.end()){
    std::cerr << std::endl << "ERROR! '"<<year<<"' is not a valid year! Please choose from ";
     std::vector<std::string>::iterator it = years.begin();
    for(it=years.begin(); it!=years.end(); it++){
      if(it!=years.begin()) std::cerr << ", ";
      std::cerr << *it;
    }
    std::cerr << std::endl;
    assert(0);
  }

  if(std::find(antiJetIDs.begin(),antiJetIDs.end(),ID)!=antiJetIDs.end()){

    filename = Form("%s/TauFES_eta-dm_%s_%s.root",datapath.data(),ID.data(),year.data());

    TFile* file = ensureTFile(filename,verbose);
    TGraphAsymmErrors *graph = dynamic_cast<TGraphAsymmErrors*>((const_cast<TFile*>(file))->Get("fes"));

    DMs = {0, 1};

    int i = 0;
    // std::map <int, std::vector<float> > innermap;
    std::vector<std::string> regions = {"barrel", "endcap"};

    for(size_t i = 0; i < regions.size(); i++){
      std::string region = regions.at(i);

      for(size_t j = 0; j < DMs.size(); j++){

        int dm = DMs.at(j);
        float y = graph->GetY()[j];
        float yup = graph->GetErrorYhigh(j);
        float ylow = graph->GetErrorYlow(j);

        FESs[region][dm] = {y-ylow, y, y + yup}; 
      }
    }

    genmatches = {1,3};
    file->Close();
  }
}


float TauFESTool::getFES(double eta,  int dm, int genmatch, const std::string& unc){
  if( (std::find(DMs.begin(),DMs.end(),dm)!=DMs.end()) && (std::find(genmatches.begin(),genmatches.end(),genmatch) != genmatches.end()) ){

    std::string region;
    if(abs(eta) < 1.5) 
      region = "barrel";
    else
      region = "endcap";

    std::vector<float> fesVector = FESs[region][dm];
    float fes = 1.0;

    if(unc == "Up"){
      fes = fesVector[2];
    }
    else if(unc == "Down"){
      fes = fesVector[0];
    }
    else{
      fes = fesVector[1];
    }

    return fes;
  }

  return 1.0;

}

float TauFESTool::getFES(double eta, int dm, const std::string& unc=""){
  return getFES(eta, dm, 1, unc);
}
