
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"
#include "BTagUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "THStack.h"
#include "TPad.h"
#include "TEfficiency.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TBranch.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

double        Jet1pt;
double        Jet2pt;
double        Jet1eta;
double        Jet2eta;
double        Jet0daughter0_pt;
double        Jet0daughter1_pt;
double        Jet1daughter0_pt;
double        Jet1daughter1_pt;
double        Jet0daughter0_eta;
double        Jet0daughter1_eta;
double        Jet1daughter0_eta;
double        Jet1daughter1_eta;
double        Jet0_bjetness;
double        Jet1_bjetness;
double        Jet0_daughter0_bjetness;
double        Jet0_daughter1_bjetness;
double        Jet1_daughter0_bjetness;
double        Jet1_daughter1_bjetness;
double        Jet0_deltaR;
double        Jet1_deltaR;

void ChooseJets(std::string* JetChoice, float jetPt[4], float jetEta[4], float jetcsvBtag[4], float jetFlavour[4], int& njets);

int main(int argc, char** argv)
{

  // Set style options
  setTDRStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(1110); 
  //gStyle->SetOptStat(0000); 
  gStyle->SetOptFit(1); 

  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> Analyses::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");
  std::string inputSFReader = gConfigParser -> readStringOption("Input::inputSFReader");
  std::string inputEffReader = gConfigParser -> readStringOption("Input::inputEffReader");
 
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");
  
  std::map<int,TChain*> ntu;
  std::map<int, TFile*> Files;

  char trees[500];
  FILE *f_trees;
  
  f_trees = fopen(inputList.c_str(),"r");
  int pos = 0;
  int pos_total = 0;

  while(fscanf(f_trees,"%s \n", trees) !=EOF ){
    std::string TREES = std::string(trees);
 
    if(TREES.find("#") != std::string::npos) continue;
    std::cout << "\nReading input: " << trees << std::endl;
    ntu[pos] = new TChain(inputTree.c_str());
    ntu[pos] -> Add(trees);
  
    pos++;
    
  }

  pos_total = pos;
  
 
  TFile* outFile = new TFile(outputName.c_str(),"RECREATE");
  outFile->cd();

  TTree* outTree = (TTree*)ntu[0]->CloneTree();
  //outTree->SetDirectory(0);
  
  for(int ii = 0; ii < pos_total; ii++){
      /*ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("Jet1pt",1);
      ntu[ii]->SetBranchStatus("Jet2pt",1);
      ntu[ii]->SetBranchStatus("Jet1eta",1);
      ntu[ii]->SetBranchStatus("Jet2eta",1);
      ntu[ii]->SetBranchStatus("Jet0_deltaR",1);
      ntu[ii]->SetBranchStatus("Jet1_deltaR",1);
      ntu[ii]->SetBranchStatus("Jet0daughter0_pt",1);
      ntu[ii]->SetBranchStatus("Jet0daughter1_pt",1);
      ntu[ii]->SetBranchStatus("Jet1daughter0_pt",1);
      ntu[ii]->SetBranchStatus("Jet1daughter1_pt",1);
      ntu[ii]->SetBranchStatus("Jet0daughter0_eta",1);
      ntu[ii]->SetBranchStatus("Jet0daughter1_eta",1);
      ntu[ii]->SetBranchStatus("Jet1daughter0_eta",1);
      ntu[ii]->SetBranchStatus("Jet1daughter1_eta",1);
      ntu[ii]->SetBranchStatus("Jet0_bjetness",1);
      ntu[ii]->SetBranchStatus("Jet1_bjetness",1);
      ntu[ii]->SetBranchStatus("Jet0_daughter0_bjetness",1);
      ntu[ii]->SetBranchStatus("Jet0_daughter1_bjetness",1);
      ntu[ii]->SetBranchStatus("Jet1_daughter0_bjetness",1);
      ntu[ii]->SetBranchStatus("Jet1_daughter1_bjetness",1);*/
      ntu[ii]->SetBranchAddress("Jet1pt",&Jet1pt);
      ntu[ii]->SetBranchAddress("Jet2pt",&Jet2pt);
      ntu[ii]->SetBranchAddress("Jet1eta",&Jet1eta);
      ntu[ii]->SetBranchAddress("Jet2eta",&Jet2eta);
      ntu[ii]->SetBranchAddress("Jet0_deltaR",&Jet0_deltaR);
      ntu[ii]->SetBranchAddress("Jet1_deltaR",&Jet1_deltaR);
      ntu[ii]->SetBranchAddress("Jet0daughter0_pt",&Jet0daughter0_pt);
      ntu[ii]->SetBranchAddress("Jet0daughter1_pt",&Jet0daughter1_pt);
      ntu[ii]->SetBranchAddress("Jet1daughter0_pt",&Jet1daughter0_pt);
      ntu[ii]->SetBranchAddress("Jet1daughter1_pt",&Jet1daughter1_pt);
      ntu[ii]->SetBranchAddress("Jet0daughter0_eta",&Jet0daughter0_eta);
      ntu[ii]->SetBranchAddress("Jet0daughter1_eta",&Jet0daughter1_eta);
      ntu[ii]->SetBranchAddress("Jet1daughter0_eta",&Jet1daughter0_eta);
      ntu[ii]->SetBranchAddress("Jet1daughter1_eta",&Jet1daughter1_eta);
      ntu[ii]->SetBranchAddress("Jet0_bjetness",&Jet0_bjetness);
      ntu[ii]->SetBranchAddress("Jet1_bjetness",&Jet1_bjetness);
      ntu[ii]->SetBranchAddress("Jet0_daughter0_bjetness",&Jet0_daughter0_bjetness);
      ntu[ii]->SetBranchAddress("Jet0_daughter1_bjetness",&Jet0_daughter1_bjetness);
      ntu[ii]->SetBranchAddress("Jet1_daughter0_bjetness",&Jet1_daughter0_bjetness);
      ntu[ii]->SetBranchAddress("Jet1_daughter1_bjetness",&Jet1_daughter1_bjetness);
  }
 
  int eventID; 
  double evweight_btagSF_sj1_sj2;
  double evweight_btagSF_errorup_sj1_sj2;
  double evweight_btagSF_errordown_sj1_sj2;  
  double evweight_btagSF_errorup_sj1d_sj2;
  double evweight_btagSF_errordown_sj1d_sj2;  
  double evweight_btagSF_errorup_sj1_sj2d;
  double evweight_btagSF_errordown_sj1_sj2d;  
  double evweight_btagSF_errorup_sj1d_sj2d;
  double evweight_btagSF_errordown_sj1d_sj2d;  
  double evweight_btagSF_fj1_sj2;
  double evweight_btagSF_errorup_fj1_sj2;
  double evweight_btagSF_errordown_fj1_sj2;  
  double evweight_btagSF_errorup_fj1_sj2d;
  double evweight_btagSF_errordown_fj1_sj2d;  
  double evweight_btagSF_sj1_fj2;
  double evweight_btagSF_errorup_sj1_fj2;
  double evweight_btagSF_errordown_sj1_fj2;  
  double evweight_btagSF_errorup_sj1d_fj2;
  double evweight_btagSF_errordown_sj1d_fj2;  
  double evweight_btagSF_fj1_fj2;
  double evweight_btagSF_errorup_fj1_fj2;
  double evweight_btagSF_errordown_fj1_fj2; 
  
  TBranch* b_eventID = outTree -> Branch("eventID",&eventID,"eventID/I"           );
  TBranch* b_evweight_btagSF_sj1_sj2 = outTree -> Branch("evweight_btagSF_sj1_sj2",&evweight_btagSF_sj1_sj2,"evweight_btagSF_sj1_sj2/D"           );
  TBranch* b_evweight_btagSF_errorup_sj1_sj2 = outTree -> Branch("evweight_btagSF_errorup_sj1_sj2",&evweight_btagSF_errorup_sj1_sj2,"evweight_btagSF_errorup_sj1_sj2/D"  );
  TBranch* b_evweight_btagSF_errordown_sj1_sj2 = outTree -> Branch("evweight_btagSF_errordown_sj1_sj2",&evweight_btagSF_errordown_sj1_sj2,"evweight_btagSF_errordown_sj1_sj2/D");
  TBranch* b_evweight_btagSF_errorup_sj1d_sj2 = outTree -> Branch("evweight_btagSF_errorup_sj1d_sj2",&evweight_btagSF_errorup_sj1d_sj2,"evweight_btagSF_errorup_sj1d_sj2/D"  );
  TBranch* b_evweight_btagSF_errordown_sj1d_sj2 = outTree -> Branch("evweight_btagSF_errordown_sj1d_sj2",&evweight_btagSF_errordown_sj1d_sj2,"evweight_btagSF_errordown_sj1d_sj2/D");
  TBranch* b_evweight_btagSF_errorup_sj1_sj2d = outTree -> Branch("evweight_btagSF_errorup_sj1_sj2d",&evweight_btagSF_errorup_sj1_sj2d,"evweight_btagSF_errorup_sj1_sj2d/D"  );
  TBranch* b_evweight_btagSF_errordown_sj1_sj2d = outTree -> Branch("evweight_btagSF_errordown_sj1_sj2d",&evweight_btagSF_errordown_sj1_sj2d,"evweight_btagSF_errordown_sj1_sj2d/D");
  TBranch* b_evweight_btagSF_errorup_sj1d_sj2d = outTree -> Branch("evweight_btagSF_errorup_sj1d_sj2d",&evweight_btagSF_errorup_sj1d_sj2d,"evweight_btagSF_errorup_sj1d_sj2d/D"  );
  TBranch* b_evweight_btagSF_errordown_sj1d_sj2d = outTree -> Branch("evweight_btagSF_errordown_sj1d_sj2d",&evweight_btagSF_errordown_sj1d_sj2d,"evweight_btagSF_errordown_sj1d_sj2d/D");
  TBranch* b_evweight_btagSF_fj1_sj2 = outTree -> Branch("evweight_btagSF_fj1_sj2",&evweight_btagSF_fj1_sj2,"evweight_btagSF_fj1_sj2/D"           );
  TBranch* b_evweight_btagSF_errorup_fj1_sj2 = outTree -> Branch("evweight_btagSF_errorup_fj1_sj2",&evweight_btagSF_errorup_fj1_sj2,"evweight_btagSF_errorup_fj1_sj2/D"  );
  TBranch* b_evweight_btagSF_errordown_fj1_sj2 = outTree -> Branch("evweight_btagSF_errordown_fj1_sj2",&evweight_btagSF_errordown_fj1_sj2,"evweight_btagSF_errordown_fj1_sj2/D");
  TBranch* b_evweight_btagSF_errorup_fj1_sj2d = outTree -> Branch("evweight_btagSF_errorup_fj1_sj2d",&evweight_btagSF_errorup_fj1_sj2d,"evweight_btagSF_errorup_fj1_sj2d/D"  );
  TBranch* b_evweight_btagSF_errordown_fj1_sj2d = outTree -> Branch("evweight_btagSF_errordown_fj1_sj2d",&evweight_btagSF_errordown_fj1_sj2d,"evweight_btagSF_errordown_fj1_sj2d/D");
  TBranch* b_evweight_btagSF_sj1_fj2 = outTree -> Branch("evweight_btagSF_sj1_fj2",&evweight_btagSF_sj1_fj2,"evweight_btagSF_sj1_fj2/D"           );
  TBranch* b_evweight_btagSF_errorup_sj1_fj2 = outTree -> Branch("evweight_btagSF_errorup_sj1_fj2",&evweight_btagSF_errorup_sj1_fj2,"evweight_btagSF_errorup_sj1_fj2/D"  );
  TBranch* b_evweight_btagSF_errordown_sj1_fj2 = outTree -> Branch("evweight_btagSF_errordown_sj1_fj2",&evweight_btagSF_errordown_sj1_fj2,"evweight_btagSF_errordown_sj1_fj2/D");
  TBranch* b_evweight_btagSF_errorup_sj1d_fj2 = outTree -> Branch("evweight_btagSF_errorup_sj1d_fj2",&evweight_btagSF_errorup_sj1d_fj2,"evweight_btagSF_errorup_sj1d_fj2/D"  );
  TBranch* b_evweight_btagSF_errordown_sj1d_fj2 = outTree -> Branch("evweight_btagSF_errordown_sj1d_fj2",&evweight_btagSF_errordown_sj1d_fj2,"evweight_btagSF_errordown_sj1d_fj2/D");
  TBranch* b_evweight_btagSF_fj1_fj2 = outTree -> Branch("evweight_btagSF_fj1_fj2",&evweight_btagSF_fj1_fj2,"evweight_btagSF_fj1_fj2/D"           );
  TBranch* b_evweight_btagSF_errorup_fj1_fj2 = outTree -> Branch("evweight_btagSF_errorup_fj1_fj2",&evweight_btagSF_errorup_fj1_fj2,"evweight_btagSF_errorup_fj1_fj2/D"  );
  TBranch* b_evweight_btagSF_errordown_fj1_fj2 = outTree -> Branch("evweight_btagSF_errordown_fj1_fj2",&evweight_btagSF_errordown_fj1_fj2,"evweight_btagSF_errordown_fj1_fj2/D");

  float jetPt[4];
  float jetEta[4];
  float jetcsvBtag[4];
  float jetFlavour[4];
  
  int njets = 4;

  float jetBtagSF[4] = {0.,0.,0.,0.};
  float jetBtagSFErrorUP[4] = {1.,1.,1.,1.};
  float jetBtagSFErrorDown[4] = {1.,1.,1.,1.};
  float jetBtagEff[4] = {1.,1.,1.,1.};
  float jetBtagEffError[4] = {1.,1.,1.,1.};

  TLorentzVector* p4 = new TLorentzVector();

  BtagSFReader* SFReader = new BtagSFReader(inputSFReader);
  
  std::map<int,BtagEfficiencyReader*> EffReader;

  std::string WP = "loose";
  std::string JetChoice = "";
  
  for(int jj = 0; jj < 9; jj++){
         
    if(jj == 0) JetChoice = "sj1_sj2";
    else if(jj == 1) JetChoice  = "sj1d_sj2";
    else if(jj == 2) JetChoice  = "sj1_sj2d";
    else if(jj == 3) JetChoice  = "sj1d_sj2d";
    else if(jj == 4) JetChoice  = "fj1_sj2";
    else if(jj == 5) JetChoice  = "fj1_sj2d";
    else if(jj == 6) JetChoice  = "sj1_fj2";
    else if(jj == 7) JetChoice  = "sj1d_fj2";
    else if(jj == 8) JetChoice  = "fj1_fj2";


    EffReader[jj] = new BtagEfficiencyReader();
    EffReader[jj]->Init(std::string(inputEffReader+"_"+JetChoice+".root"));
    //std::cout << "BtagEff = " << std::string(inputEffReader+"_"+JetChoice+".root") << std::endl;
  }

  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){

          if(ientry%1000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry << std::endl;
          ntu[nn]->GetEntry(ientry);

          eventID = ientry;
   
          evweight_btagSF_sj1_sj2 = 0.;
          evweight_btagSF_errorup_sj1_sj2 = 0.;
          evweight_btagSF_errordown_sj1_sj2 = 0.; 
          evweight_btagSF_errorup_sj1d_sj2 = 0.;
          evweight_btagSF_errordown_sj1d_sj2 = 0.;
          evweight_btagSF_errorup_sj1_sj2d = 0.;
          evweight_btagSF_errordown_sj1_sj2d = 0.; 
          evweight_btagSF_errorup_sj1d_sj2d = 0.;
          evweight_btagSF_errordown_sj1d_sj2d = 0.; 
          evweight_btagSF_fj1_sj2 = 0.;
          evweight_btagSF_errorup_fj1_sj2 = 0.;
          evweight_btagSF_errordown_fj1_sj2 = 0.;
          evweight_btagSF_errorup_fj1_sj2d = 0.;
          evweight_btagSF_errordown_fj1_sj2d = 0.;
          evweight_btagSF_sj1_fj2 = 0.;
          evweight_btagSF_errorup_sj1_fj2 = 0.;
          evweight_btagSF_errordown_sj1_fj2 = 0.;
          evweight_btagSF_errorup_sj1d_fj2 = 0.;
          evweight_btagSF_errordown_sj1d_fj2 = 0.;
          evweight_btagSF_fj1_fj2 = 0.;
          evweight_btagSF_errorup_fj1_fj2 = 0.;
          evweight_btagSF_errordown_fj1_fj2 = 0.;
           
          float weight = 0.;
          float weight_errorup = 0.;
          float weight_errorup_naive = 0.;
          float weight_errordown = 0.;
          float weight_errordown_naive = 0.;
 
          for(int jj = 0; jj < 9; jj++){
         
           if(jj == 0) JetChoice = "sj1_sj2";
           else if(jj == 1) JetChoice  = "sj1d_sj2";
           else if(jj == 2) JetChoice  = "sj1_sj2d";
           else if(jj == 3) JetChoice  = "sj1d_sj2d";
           else if(jj == 4) JetChoice  = "fj1_sj2";
           else if(jj == 5) JetChoice  = "fj1_sj2d";
           else if(jj == 6) JetChoice  = "sj1_fj2";
           else if(jj == 7) JetChoice  = "sj1d_fj2";
           else if(jj == 8) JetChoice  = "fj1_fj2";

           ChooseJets(&JetChoice,jetPt,jetEta,jetcsvBtag,jetFlavour,njets);

           bool isBad = false;
           for(int ii = 0; ii < njets; ii++)
               if(jetPt[ii] == 0.) isBad = true;
           if(isBad == true) continue;
              
           for(int ii = 0; ii < njets; ii++){
              p4->SetPtEtaPhiE(jetPt[ii],jetEta[ii],0.,0.);  
              jetBtagSF[ii] = SFReader->getSF(p4,jetFlavour[ii],WP);
              jetBtagSFErrorUP[ii] = SFReader->getSFErrorUp(p4,jetFlavour[ii],WP);
              jetBtagSFErrorDown[ii] = SFReader->getSFErrorDown(p4,jetFlavour[ii],WP);
              jetBtagEff[ii] = EffReader[jj]->getBtagEfficiency(p4,WP,jetFlavour[ii]);
              jetBtagEffError[ii] = EffReader[jj]->getBtagEfficiencyError(p4,WP,jetFlavour[ii]);  
           }

           if(JetChoice.find("sj1d") != std::string::npos){
              jetBtagSFErrorUP[0] = 2*jetBtagSFErrorUP[0]; jetBtagSFErrorUP[1] = 2*jetBtagSFErrorUP[1];
              jetBtagSFErrorDown[0] = 2*jetBtagSFErrorDown[0]; jetBtagSFErrorDown[1] = 2*jetBtagSFErrorDown[1]; 
           }  
           if(JetChoice.find("sj2d") != std::string::npos){
              if(njets == 4){
                 jetBtagSFErrorUP[2] = 2*jetBtagSFErrorUP[2]; jetBtagSFErrorUP[3] = 2*jetBtagSFErrorUP[3];
                 jetBtagSFErrorDown[2] = 2*jetBtagSFErrorDown[2]; jetBtagSFErrorDown[3] = 2*jetBtagSFErrorDown[3]; 
              } 
              else if(njets == 3){
                 jetBtagSFErrorUP[1] = 2*jetBtagSFErrorUP[1]; jetBtagSFErrorUP[2] = 2*jetBtagSFErrorUP[2];
                 jetBtagSFErrorDown[1] = 2*jetBtagSFErrorDown[1]; jetBtagSFErrorDown[2] = 2*jetBtagSFErrorDown[2]; 
              }
           }  

           if(njets == 4){
             weight = eventWeight_4jets(WP,jetBtagSF[0],jetBtagSF[1],jetBtagSF[2],jetBtagSF[3],jetBtagEff[0],jetBtagEff[1],jetBtagEff[2],jetBtagEff[3],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2],jetcsvBtag[3]);
             weight_errorup = eventWeight_error_4jets(WP,jetBtagSF[0],jetBtagSFErrorUP[0],jetBtagSF[1],jetBtagSFErrorUP[1],jetBtagSF[2],jetBtagSFErrorUP[2],jetBtagSF[3],jetBtagSFErrorUP[3],jetBtagEff[0],jetBtagEffError[0],jetBtagEff[1],jetBtagEffError[1],jetBtagEff[2],jetBtagEffError[2],jetBtagEff[3],jetBtagEffError[3],jetFlavour[0],jetFlavour[1],jetFlavour[2],jetFlavour[3],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2],jetcsvBtag[3]);
              weight_errorup_naive = eventWeight_NAIVE_4jets_UP(WP,jetBtagSF[0],jetBtagSFErrorUP[0],jetBtagSF[1],jetBtagSFErrorUP[1],jetBtagSF[2],jetBtagSFErrorUP[2],jetBtagSF[3],jetBtagSFErrorUP[3],jetBtagEff[0],jetBtagEff[1],jetBtagEff[2],jetBtagEff[3],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2],jetcsvBtag[3]);
              weight_errordown = eventWeight_error_4jets(WP,jetBtagSF[0],jetBtagSFErrorDown[0],jetBtagSF[1],jetBtagSFErrorDown[1],jetBtagSF[2],jetBtagSFErrorDown[2],jetBtagSF[3],jetBtagSFErrorDown[3],jetBtagEff[0],jetBtagEffError[0],jetBtagEff[1],jetBtagEffError[1],jetBtagEff[2],jetBtagEffError[2],jetBtagEff[3],jetBtagEffError[3],jetFlavour[0],jetFlavour[1],jetFlavour[2],jetFlavour[3],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2],jetcsvBtag[3]);
              weight_errordown_naive = eventWeight_NAIVE_4jets_DOWN(WP,jetBtagSF[0],jetBtagSFErrorDown[0],jetBtagSF[1],jetBtagSFErrorDown[1],jetBtagSF[2],jetBtagSFErrorDown[2],jetBtagSF[3],jetBtagSFErrorDown[3],jetBtagEff[0],jetBtagEff[1],jetBtagEff[2],jetBtagEff[3],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2],jetcsvBtag[3]);
            }
            else if(njets == 3){
              weight = eventWeight_3jets(WP,jetBtagSF[0],jetBtagSF[1],jetBtagSF[2],jetBtagEff[0],jetBtagEff[1],jetBtagEff[2],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2]);
             weight_errorup = eventWeight_error_3jets(WP,jetBtagSF[0],jetBtagSFErrorUP[0],jetBtagSF[1],jetBtagSFErrorUP[1],jetBtagSF[2],jetBtagSFErrorUP[2],jetBtagEff[0],jetBtagEffError[0],jetBtagEff[1],jetBtagEffError[1],jetBtagEff[2],jetBtagEffError[2],jetFlavour[0],jetFlavour[1],jetFlavour[2],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2]);
             weight_errorup_naive = eventWeight_NAIVE_3jets_UP(WP,jetBtagSF[0],jetBtagSFErrorUP[0],jetBtagSF[1],jetBtagSFErrorUP[1],jetBtagSF[2],jetBtagSFErrorUP[2],jetBtagEff[0],jetBtagEff[1],jetBtagEff[2],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2]);
             weight_errordown = eventWeight_error_3jets(WP,jetBtagSF[0],jetBtagSFErrorDown[0],jetBtagSF[1],jetBtagSFErrorDown[1],jetBtagSF[2],jetBtagSFErrorDown[2],jetBtagEff[0],jetBtagEffError[0],jetBtagEff[1],jetBtagEffError[1],jetBtagEff[2],jetBtagEffError[2],jetFlavour[0],jetFlavour[1],jetFlavour[2],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2]);
             weight_errordown_naive = eventWeight_NAIVE_3jets_DOWN(WP,jetBtagSF[0],jetBtagSFErrorDown[0],jetBtagSF[1],jetBtagSFErrorDown[1],jetBtagSF[2],jetBtagSFErrorDown[2],jetBtagEff[0],jetBtagEff[1],jetBtagEff[2],jetcsvBtag[0],jetcsvBtag[1],jetcsvBtag[2]);
            }
            else if(njets == 2){
             weight = eventWeight_2jets(WP,jetBtagSF[0],jetBtagSF[1],jetBtagEff[0],jetBtagEff[1],jetcsvBtag[0],jetcsvBtag[1]);
             weight_errorup = eventWeight_error_2jets(WP,jetBtagSF[0],jetBtagSFErrorUP[0],jetBtagSF[1],jetBtagSFErrorUP[1],jetBtagEff[0],jetBtagEffError[0],jetBtagEff[1],jetBtagEffError[1],jetFlavour[0],jetFlavour[1],jetcsvBtag[0],jetcsvBtag[1]);
             weight_errorup_naive = eventWeight_NAIVE_2jets_UP(WP,jetBtagSF[0],jetBtagSFErrorUP[0],jetBtagSF[1],jetBtagSFErrorUP[1],jetBtagEff[0],jetBtagEff[1],jetcsvBtag[0],jetcsvBtag[1]);
             weight_errordown = eventWeight_error_2jets(WP,jetBtagSF[0],jetBtagSFErrorDown[0],jetBtagSF[1],jetBtagSFErrorDown[1],jetBtagEff[0],jetBtagEffError[0],jetBtagEff[1],jetBtagEffError[1],jetFlavour[0],jetFlavour[1],jetcsvBtag[0],jetcsvBtag[1]);
             weight_errordown_naive = eventWeight_NAIVE_2jets_DOWN(WP,jetBtagSF[0],jetBtagSFErrorDown[0],jetBtagSF[1],jetBtagSFErrorDown[1],jetBtagEff[0],jetBtagEff[1],jetcsvBtag[0],jetcsvBtag[1]);
            }
            
            float signUP = (weight_errorup_naive-weight)/fabs(weight_errorup_naive-weight); 
            float signDOWN = (weight_errordown_naive-weight)/fabs(weight_errordown_naive-weight);  
            weight_errorup = weight_errorup*signUP; 
            weight_errordown = weight_errordown*signDOWN;

            if(jj == 0){
               evweight_btagSF_sj1_sj2 = (double)weight;
               evweight_btagSF_errorup_sj1_sj2 = (double)weight_errorup;
               evweight_btagSF_errordown_sj1_sj2 = (double)weight_errordown;
            }
            else if(jj == 1){
               evweight_btagSF_errorup_sj1d_sj2 = (double)weight_errorup;
               evweight_btagSF_errordown_sj1d_sj2 = (double)weight_errordown;
            }
            else if(jj == 2){
               evweight_btagSF_errorup_sj1_sj2d = (double)weight_errorup;
               evweight_btagSF_errordown_sj1_sj2d = (double)weight_errordown;
            }
            else if(jj == 3){
               evweight_btagSF_errorup_sj1d_sj2d = (double)weight_errorup;
               evweight_btagSF_errordown_sj1d_sj2d = (double)weight_errordown;
            }
            else if(jj == 4){
               evweight_btagSF_fj1_sj2 = (double)weight;
               evweight_btagSF_errorup_fj1_sj2 = (double)weight_errorup;
               evweight_btagSF_errordown_fj1_sj2 = (double)weight_errordown;
            }
            else if(jj == 5){
               evweight_btagSF_errorup_fj1_sj2d = (double)weight_errorup;
               evweight_btagSF_errordown_fj1_sj2d = (double)weight_errordown;
            }
            else if(jj == 6){
               evweight_btagSF_sj1_fj2 = (double)weight;
               evweight_btagSF_errorup_sj1_fj2 = (double)weight_errorup;
               evweight_btagSF_errordown_sj1_fj2 = (double)weight_errordown;
            }
            else if(jj == 7){
               evweight_btagSF_errorup_sj1d_fj2 = (double)weight_errorup;
               evweight_btagSF_errordown_sj1d_fj2 = (double)weight_errordown;
            }
            else if(jj == 8){
               evweight_btagSF_fj1_fj2 = (double)weight;
               evweight_btagSF_errorup_fj1_fj2 = (double)weight_errorup;
               evweight_btagSF_errordown_fj1_fj2 = (double)weight_errordown;
            }

           } 

           b_eventID->Fill();

           b_evweight_btagSF_sj1_sj2->Fill();    
           b_evweight_btagSF_errorup_sj1_sj2->Fill();    
           b_evweight_btagSF_errordown_sj1_sj2->Fill();    
       
           b_evweight_btagSF_errorup_sj1d_sj2->Fill(); 
           b_evweight_btagSF_errordown_sj1d_sj2->Fill();

           b_evweight_btagSF_errorup_sj1_sj2d->Fill();    
           b_evweight_btagSF_errordown_sj1_sj2d->Fill();

           b_evweight_btagSF_errorup_sj1d_sj2d->Fill();    
           b_evweight_btagSF_errordown_sj1d_sj2d->Fill(); 

           b_evweight_btagSF_fj1_sj2->Fill();    
           b_evweight_btagSF_errorup_fj1_sj2->Fill();    
           b_evweight_btagSF_errordown_fj1_sj2->Fill();

           b_evweight_btagSF_errorup_fj1_sj2d->Fill();     
           b_evweight_btagSF_errordown_fj1_sj2d->Fill();

           b_evweight_btagSF_sj1_fj2->Fill();    
           b_evweight_btagSF_errorup_sj1_fj2->Fill();   
           b_evweight_btagSF_errordown_sj1_fj2->Fill();
           
           b_evweight_btagSF_errorup_sj1d_fj2->Fill();    
           b_evweight_btagSF_errordown_sj1d_fj2->Fill();
            
           b_evweight_btagSF_fj1_fj2->Fill();   
           b_evweight_btagSF_errorup_fj1_fj2->Fill();  
           b_evweight_btagSF_errordown_fj1_fj2->Fill();       
      } 
  }

  outFile->cd();
  outTree->Write();
  outFile->Close();
 
}

void ChooseJets(std::string* JetChoice, float jetPt[4], float jetEta[4], float jetcsvBtag[4], float jetFlavour[4], int& njets){

     if(*JetChoice == "sj1_sj2"){
          
          njets = 4;
          
          jetPt[0] = (float)Jet0daughter0_pt;
          jetPt[1] = (float)Jet0daughter1_pt;
          jetPt[2] = (float)Jet1daughter0_pt;
          jetPt[3] = (float)Jet1daughter1_pt;

          jetEta[0] = (float)Jet0daughter0_eta;
          jetEta[1] = (float)Jet0daughter1_eta;
          jetEta[2] = (float)Jet1daughter0_eta;
          jetEta[3] = (float)Jet1daughter1_eta;
 
          jetcsvBtag[0] = (float)Jet0_daughter0_bjetness;
          jetcsvBtag[1] = (float)Jet0_daughter1_bjetness;
          jetcsvBtag[2] = (float)Jet1_daughter0_bjetness;
          jetcsvBtag[3] = (float)Jet1_daughter1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;
          jetFlavour[3] = 5;  

     }
     else if(*JetChoice == "sj1d_sj2"){
          
          njets = 4;
          
          jetPt[0] = (float)Jet0daughter0_pt;
          jetPt[1] = (float)Jet0daughter1_pt;
          jetPt[2] = (float)Jet1daughter0_pt;
          jetPt[3] = (float)Jet1daughter1_pt;

          jetEta[0] = (float)Jet0daughter0_eta;
          jetEta[1] = (float)Jet0daughter1_eta;
          jetEta[2] = (float)Jet1daughter0_eta;
          jetEta[3] = (float)Jet1daughter1_eta;
 
          jetcsvBtag[0] = (float)Jet0_daughter0_bjetness;
          jetcsvBtag[1] = (float)Jet0_daughter1_bjetness;
          jetcsvBtag[2] = (float)Jet1_daughter0_bjetness;
          jetcsvBtag[3] = (float)Jet1_daughter1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;
          jetFlavour[3] = 5;  

     }
     else if(*JetChoice == "sj1_sj2d"){
          
          njets = 4;
          
          jetPt[0] = (float)Jet0daughter0_pt;
          jetPt[1] = (float)Jet0daughter1_pt;
          jetPt[2] = (float)Jet1daughter0_pt;
          jetPt[3] = (float)Jet1daughter1_pt;

          jetEta[0] = (float)Jet0daughter0_eta;
          jetEta[1] = (float)Jet0daughter1_eta;
          jetEta[2] = (float)Jet1daughter0_eta;
          jetEta[3] = (float)Jet1daughter1_eta;
 
          jetcsvBtag[0] = (float)Jet0_daughter0_bjetness;
          jetcsvBtag[1] = (float)Jet0_daughter1_bjetness;
          jetcsvBtag[2] = (float)Jet1_daughter0_bjetness;
          jetcsvBtag[3] = (float)Jet1_daughter1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;
          jetFlavour[3] = 5;  

     }
     if(*JetChoice == "sj1d_sj2d"){
          
          njets = 4;
          
          jetPt[0] = (float)Jet0daughter0_pt;
          jetPt[1] = (float)Jet0daughter1_pt;
          jetPt[2] = (float)Jet1daughter0_pt;
          jetPt[3] = (float)Jet1daughter1_pt;

          jetEta[0] = (float)Jet0daughter0_eta;
          jetEta[1] = (float)Jet0daughter1_eta;
          jetEta[2] = (float)Jet1daughter0_eta;
          jetEta[3] = (float)Jet1daughter1_eta;
 
          jetcsvBtag[0] = (float)Jet0_daughter0_bjetness;
          jetcsvBtag[1] = (float)Jet0_daughter1_bjetness;
          jetcsvBtag[2] = (float)Jet1_daughter0_bjetness;
          jetcsvBtag[3] = (float)Jet1_daughter1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;
          jetFlavour[3] = 5;  

     }
     else if(*JetChoice == "fj1_sj2"){
          
          njets = 3;
          
          jetPt[0] = (float)Jet1pt;
          jetPt[1] = (float)Jet1daughter0_pt;
          jetPt[2] = (float)Jet1daughter1_pt;
          
          jetEta[0] = (float)Jet1eta;
          jetEta[1] = (float)Jet1daughter0_eta;
          jetEta[2] = (float)Jet1daughter1_eta;
          
          jetcsvBtag[0] = (float)Jet0_bjetness;
          jetcsvBtag[1] = (float)Jet1_daughter0_bjetness;
          jetcsvBtag[2] = (float)Jet1_daughter1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;

     }
     else if(*JetChoice == "fj1_sj2d"){
          
          njets = 3;
          
          jetPt[0] = (float)Jet1pt;
          jetPt[1] = (float)Jet1daughter0_pt;
          jetPt[2] = (float)Jet1daughter1_pt;
          
          jetEta[0] = (float)Jet1eta;
          jetEta[1] = (float)Jet1daughter0_eta;
          jetEta[2] = (float)Jet1daughter1_eta;
          
          jetcsvBtag[0] = (float)Jet0_bjetness;
          jetcsvBtag[1] = (float)Jet1_daughter0_bjetness;
          jetcsvBtag[2] = (float)Jet1_daughter1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;

     }
     else if(*JetChoice == "sj1_fj2"){
          
          njets = 3;
          
          jetPt[0] = (float)Jet0daughter0_pt;
          jetPt[1] = (float)Jet0daughter1_pt;
          jetPt[2] = (float)Jet2pt;
          
          jetEta[0] = (float)Jet0daughter0_eta;
          jetEta[1] = (float)Jet0daughter1_eta;
          jetEta[2] = (float)Jet2eta;
          
          jetcsvBtag[0] = (float)Jet0_daughter0_bjetness;
          jetcsvBtag[1] = (float)Jet0_daughter1_bjetness;
          jetcsvBtag[2] = (float)Jet1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;

     }
     else if(*JetChoice == "sj1d_fj2"){
          
          njets = 3;
          
          jetPt[0] = (float)Jet0daughter0_pt;
          jetPt[1] = (float)Jet0daughter1_pt;
          jetPt[2] = (float)Jet2pt;
          
          jetEta[0] = (float)Jet0daughter0_eta;
          jetEta[1] = (float)Jet0daughter1_eta;
          jetEta[2] = (float)Jet2eta;
          
          jetcsvBtag[0] = (float)Jet0_daughter0_bjetness;
          jetcsvBtag[1] = (float)Jet0_daughter1_bjetness;
          jetcsvBtag[2] = (float)Jet1_bjetness;

          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          jetFlavour[2] = 5;

     }
     else if(*JetChoice == "fj1_fj2"){
          
          njets = 2;
          
          jetPt[0] = (float)Jet1pt;
          jetPt[1] = (float)Jet2pt;
         
          jetEta[0] = (float)Jet1eta;
          jetEta[1] = (float)Jet2eta;
          
          jetcsvBtag[0] = (float)Jet0_bjetness;
          jetcsvBtag[1] = (float)Jet1_bjetness;
          
          jetFlavour[0] = 5;
          jetFlavour[1] = 5;
          
     }

}
