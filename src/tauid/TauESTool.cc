#include "TauESTool.h"
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



void TauESTool::disabled() const {
  std::cerr << std::endl << "ERROR! Method has been disabled! isVsPT = "<<isVsPT<<", isVsDM = "
            << isVsDM<<", isVsEta = "<<isVsEta<< std::endl;
  assert(0);
}


TauESTool::TauESTool(const std::string& datapath, const std::string& year, const std::string& id): ID(id){

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

    filename_lowpt = Form("%s/TauES_dm_%s_%s.root",datapath.data(),ID.data(),year.data());
    filename_highpt = Form("%s/TauES_dm_%s_%s_ptgt100.root",datapath.data(),ID.data(),year.data());

    TFile* file_lowpt = ensureTFile(filename_lowpt,verbose);
    TFile* file_highpt = ensureTFile(filename_highpt,verbose);

    hist_lowpt = extractTH1(file_lowpt,"tes");
    hist_highpt = extractTH1(file_highpt,"tes");

    hist_lowpt->SetDirectory(nullptr);
    hist_highpt->SetDirectory(nullptr);

    pt_low = 34.0;    // average pT in Z-> tautau measurement (incl. in DM)
    pt_high = 170.0;  // average pT in W* -> taunu measurement (incl. in DM)
    
    DMs    = {0,1,10};
    if (ID.find("oldDM") == std::string::npos)
    {
        DMs.push_back(11);
    }

    file_lowpt->Close();
    file_highpt->Close();
  }
}


float TauESTool::getTES(double pt,  int dm, int genmatch, const std::string& unc){
  if(std::find(DMs.begin(),DMs.end(),dm)!=DMs.end()){
    if(genmatch==5){
      Int_t bin = hist_lowpt->GetXaxis()->FindBin(dm);
      float tes  = static_cast<float>(hist_lowpt->GetBinContent(bin));
      float err;

      if(!unc.empty()){
        if(pt >= pt_high){ // high pt
          float bin_high = hist_highpt->GetXaxis()->FindBin(dm);
          err = hist_highpt->GetXaxis()->FindBin(bin_high);
        }
        else if(pt > pt_low){ // linearly interpolate between low and high pT
          float bin_high = hist_highpt->GetXaxis()->FindBin(dm);
          float err_high = hist_highpt->GetBinError(bin_high);
          float err_low = hist_lowpt->GetBinError(bin);
          err = err_low + (err_high-err_low)/(pt_high-pt_low)*(pt-pt_low);
        }
        else{ // low pT
          err = hist_lowpt->GetBinError(bin);
        }
      }

      if(unc=="Up"){
        tes += err;
      }
      else if(unc=="Down")
        tes -= err;
      return tes;
    }
    return 1.0;
  }
}

float TauESTool::getTES(double pt, int dm, const std::string& unc=""){
  return getTES(pt, dm, 5, unc);
}

float TauESTool::getTES_highpt(int dm, int genmatch, const std::string& unc){
  if(std::find(DMs.begin(),DMs.end(),dm)!=DMs.end()){
    if(genmatch==5){
      Int_t bin = hist_highpt->GetXaxis()->FindBin(dm);
      float tes  = static_cast<float>(hist_lowpt->GetBinContent(bin));
      if(unc=="Up"){
        tes += hist_highpt->GetBinError(bin);
      }
      else if(unc=="Down")
        tes -= hist_highpt->GetBinError(bin);
      return tes;
    }
    return 1.0;
  }
}

float TauESTool::getTES_highpt(int dm, const std::string& unc){
  return getTES_highpt(dm,5,unc);
}
