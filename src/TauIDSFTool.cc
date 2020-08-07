#include "TauIDSFTool.h"
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



void TauIDSFTool::disabled() const {
  std::cerr << std::endl << "ERROR! Method has been disabled! isVsPT = "<<isVsPT<<", isVsDM = "
            << isVsDM<<", isVsEta = "<<isVsEta<< std::endl;
  assert(0);
}

// ------ TauIDSFTool ------ //
TauIDSFTool::TauIDSFTool(const std::string& datapath, const std::string& year, const std::string& id, const std::string& wp, const bool dm, const bool embedding): ID(id), WP(wp){

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
    if(dm){
      TString filename;
      if (embedding) {
          if (ID.find("oldDM") != std::string::npos)
          {
             std::cerr << "Scale factors for embedded samples are not provided for the MVA IDs." << std::endl;
             assert(0);
          }
          filename = Form("%s/TauID_SF_dm_%s_%s_EMB.root",datapath.data(),ID.data(),year.data());
      }
      else {
          filename = Form("%s/TauID_SF_dm_%s_%s.root",datapath.data(),ID.data(),year.data());
      }
      TFile* file = ensureTFile(filename,verbose);
      hist = extractTH1(file,WP);
      hist->SetDirectory(nullptr);
      file->Close();
      delete file;
      DMs    = {0,1,10};
      if (ID.find("oldDM") == std::string::npos)
      {
          DMs.push_back(11);
      }
      isVsDM = true;
    }else{
      TString filename;
      if (embedding) {
          if (ID.find("oldDM") != std::string::npos)
          {
             std::cerr << "Scale factors for embedded samples are not provided for the MVA IDs." << std::endl;
             assert(0);
          }
          filename = Form("%s/TauID_SF_pt_%s_%s_EMB.root",datapath.data(),ID.data(),year.data());
      }
      else {
          filename = Form("%s/TauID_SF_pt_%s_%s.root",datapath.data(),ID.data(),year.data());
      }
      TFile* file = ensureTFile(filename,verbose);
      func[""]     = extractTF1(file,Form("%s_cent",WP.data()));
      func["Up"]   = extractTF1(file,Form("%s_up",  WP.data()));
      func["Down"] = extractTF1(file,Form("%s_down",WP.data()));
      file->Close();
      delete file;
      isVsPT = true;
    }
  }else if(std::find(antiEleIDs.begin(),antiEleIDs.end(),ID)!=antiEleIDs.end()){
      if (embedding){
          std::cerr << "SF for ID " << ID << " not available for the embedded samples!" << std::endl;
          assert(0);
      }
      TString filename = Form("%s/TauID_SF_eta_%s_%s.root",datapath.data(),ID.data(),year.data());
      TFile* file = ensureTFile(filename,verbose);
      hist = extractTH1(file,WP);
      hist->SetDirectory(nullptr);
      file->Close();
      delete file;
      genmatches = {1,3};
      isVsEta    = true;
  }else if(std::find(antiMuIDs.begin(),antiMuIDs.end(),ID)!=antiMuIDs.end()){
      if (embedding){
          std::cerr << "SF for ID " << ID << " not available for the embedded samples!" << std::endl;
          assert(0);
      }
      TString filename = Form("%s/TauID_SF_eta_%s_%s.root",datapath.data(),ID.data(),year.data());
      TFile* file = ensureTFile(filename,verbose);
      hist = extractTH1(file,WP);
      hist->SetDirectory(nullptr);
      file->Close();
      delete file;
      genmatches = {2,4};
      isVsEta    = true;
  }else{
      std::cerr << "Did not recognize tau ID '" << ID << "'!" << std::endl;
      assert(0);
  }
}



float TauIDSFTool::getSFvsPT(double pt, int genmatch, const std::string& unc){
  if(!isVsPT) disabled();
  if(genmatch==5){
    float SF = static_cast<float>(func[unc]->Eval(pt));
    return SF;
  }
  return 1.0;
}

float TauIDSFTool::getSFvsPT(double pt, const std::string& unc){
  return getSFvsPT(pt,5,unc);
}



float TauIDSFTool::getSFvsDM(double pt, int dm, int genmatch, const std::string& unc) const{
  if(!isVsDM) disabled();
  if(std::find(DMs.begin(),DMs.end(),dm)!=DMs.end() or pt<=40){
    if(genmatch==5){
      Int_t bin = hist->GetXaxis()->FindBin(dm);
      float SF  = static_cast<float>(hist->GetBinContent(bin));
      if(unc=="Up")
        SF += hist->GetBinError(bin);
      else if(unc=="Down")
        SF -= hist->GetBinError(bin);
      return SF;
    }
    return 1.0;
  }
  return 0.0;
}

float TauIDSFTool::getSFvsDM(double pt, int dm, const std::string& unc) const{
  return getSFvsDM(pt,dm,5,unc);
}



float TauIDSFTool::getSFvsEta(double eta, int genmatch, const std::string& unc) const{
  if(!isVsEta) disabled();
  if(std::find(genmatches.begin(),genmatches.end(),genmatch)!=genmatches.end()){
    Int_t bin = hist->GetXaxis()->FindBin(eta);
    float SF  = static_cast<float>(hist->GetBinContent(bin));
    if(unc=="Up")
      SF += hist->GetBinError(bin);
    else if(unc=="Down")
      SF -= hist->GetBinError(bin);
    return SF;
  }
  return 1.0;
}


// ------- TauESTool ------- //

TauESTool::TauESTool(const std::string& datapath, const std::string& year, const std::string& id): ID(id){

  bool verbose = false;
  std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco"};
  std::vector<std::string> antiJetIDs = {"MVAoldDM2017v2","DeepTau2017v2p1VSjet"};

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

    TString filename_lowpt = Form("%s/TauES_dm_%s_%s.root",datapath.data(),ID.data(),year.data());
    TString filename_highpt = Form("%s/TauES_dm_%s_%s_ptgt100.root",datapath.data(),ID.data(),year.data());

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

      if(!unc.empty()){
        float err;

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

        if(unc=="Up"){
          tes += err;
        }
        else if(unc=="Down"){
          tes -= err;
        }
      }
      
      return tes;
    }
    
    return 1.0;
  }

  return 1.0;
}

float TauESTool::getTES(double pt, int dm, const std::string& unc){
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
  return 1.0;
}

float TauESTool::getTES_highpt(int dm, const std::string& unc){
  return getTES_highpt(dm,5,unc);
}

// ------- TauFESTool ------- //

TauFESTool::TauFESTool(const std::string& datapath, const std::string& year, const std::string& id): ID(id){

  bool verbose = false;
  std::vector<std::string> years      = {"2016Legacy","2017ReReco","2018ReReco"};
  std::vector<std::string> antiEleIDs = {"DeepTau2017v2p1VSe"};

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

  if(ID.find("DeepTau2017v2p1VSe") != std::string::npos){

    TString filename = Form("%s/TauFES_eta-dm_%s_%s.root",datapath.data(),ID.data(),year.data());

    TFile* file = ensureTFile(filename,verbose);
    TGraphAsymmErrors *graph = dynamic_cast<TGraphAsymmErrors* >((const_cast<TFile*>(file))->Get("fes"));

    DMs = {0, 1};
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

float TauFESTool::getFES(double eta, int dm, const std::string& unc){
  return getFES(eta, dm, 1, unc);
}

