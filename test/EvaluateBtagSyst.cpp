
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

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

void drawPlotUP(TGraph* corr, std::string outputPlot, std::string cat_x, std::string cat_y);
void drawPlotDOWN(TGraph* corr, std::string outputPlot, std::string cat_x, std::string cat_y);

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
  
  std::string inputFile = gConfigParser -> readStringOption("Input::inputFile");
  float xMin = gConfigParser -> readFloatOption("Input::inputxMin");
  float xMax = gConfigParser -> readFloatOption("Input::inputxMax");
  
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");
  std::string outputPlotName = gConfigParser -> readStringOption("Output::outputPlotName");
  
  std::map<int,TChain*> ntu;
  
  
  ntu[0] = new TChain("TCVARS");
  ntu[0] -> Add(inputFile.c_str());
  
  int pos_total = 1;

  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file" << ii << " is empty" << std::endl; 
         return -1;
      }
  }

  float mtot;
  float evWeight;
  float weightBtagSF,weightBtagSFerrUp,weightBtagSFerrDown; 
  int cut_based_ct;
  int jet1_flavour;
  float jet1_csvBtag, jet1_btagSF_M, jet1_btagSFErrorUp_M, jet1_btagSFErrorDown_M, jet1_btagEff_M, jet1_btagEffError_M;
  float jet1_pt,jet1_eta,jet1_phi,jet1_e;
  int jet2_flavour;
  float jet2_csvBtag, jet2_btagSF_M, jet2_btagSFErrorUp_M, jet2_btagSFErrorDown_M, jet2_btagEff_M, jet2_btagEffError_M;
  float jet2_pt,jet2_eta,jet2_phi,jet2_e;
  
  int n_1btag = 0;
  int n_2btag = 0;
  int n_nob_1btag = 0;
  int n_nob_2btag = 0;

  std::vector<int> njets;
  
  for(int ii = 0; ii < pos_total; ii++){
      ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("evWeight",1);
      ntu[ii]->SetBranchStatus("mtot",1);
      ntu[ii]->SetBranchStatus("jet1_flavour",1);
      ntu[ii]->SetBranchStatus("jet1_csvBtag",1);
      ntu[ii]->SetBranchStatus("jet1_btagSF_M",1);
      ntu[ii]->SetBranchStatus("jet1_btagSFErrorUp_M",1);
      ntu[ii]->SetBranchStatus("jet1_btagSFErrorDown_M",1);
      ntu[ii]->SetBranchStatus("jet1_btagEff_M",1);
      ntu[ii]->SetBranchStatus("jet1_btagEffError_M",1);
      ntu[ii]->SetBranchStatus("jet1_pt",1);
      ntu[ii]->SetBranchStatus("jet1_eta",1);
      ntu[ii]->SetBranchStatus("jet1_phi",1);
      ntu[ii]->SetBranchStatus("jet1_e",1);
      ntu[ii]->SetBranchStatus("jet2_flavour",1);
      ntu[ii]->SetBranchStatus("jet2_csvBtag",1);
      ntu[ii]->SetBranchStatus("jet2_btagSF_M",1);
      ntu[ii]->SetBranchStatus("jet2_btagSFErrorUp_M",1);
      ntu[ii]->SetBranchStatus("jet2_btagSFErrorDown_M",1);
      ntu[ii]->SetBranchStatus("jet2_btagEff_M",1);
      ntu[ii]->SetBranchStatus("jet2_btagEffError_M",1);
      ntu[ii]->SetBranchStatus("jet2_pt",1);
      ntu[ii]->SetBranchStatus("jet2_eta",1);
      ntu[ii]->SetBranchStatus("jet2_phi",1);
      ntu[ii]->SetBranchStatus("jet2_e",1);
      ntu[ii]->SetBranchStatus("cut_based_ct",1);
      ntu[ii]->SetBranchStatus("weightBtagSF",1);
      ntu[ii]->SetBranchStatus("weightBtagSFerrUp",1);
      ntu[ii]->SetBranchStatus("weightBtagSFerrDown",1);
      ntu[ii]->SetBranchAddress("evWeight",&evWeight);
      ntu[ii]->SetBranchAddress("mtot",&mtot);
      ntu[ii]->SetBranchAddress("jet1_flavour",&jet1_flavour);
      ntu[ii]->SetBranchAddress("jet1_csvBtag",&jet1_csvBtag);
      ntu[ii]->SetBranchAddress("jet1_btagSF_M",&jet1_btagSF_M);
      ntu[ii]->SetBranchAddress("jet1_btagSFErrorUp_M",&jet1_btagSFErrorUp_M);
      ntu[ii]->SetBranchAddress("jet1_btagSFErrorDown_M",&jet1_btagSFErrorDown_M);
      ntu[ii]->SetBranchAddress("jet1_btagEff_M",&jet1_btagEff_M);
      ntu[ii]->SetBranchAddress("jet1_btagEffError_M",&jet1_btagEffError_M);
      ntu[ii]->SetBranchAddress("jet1_pt",&jet1_pt);
      ntu[ii]->SetBranchAddress("jet1_eta",&jet1_eta);
      ntu[ii]->SetBranchAddress("jet1_phi",&jet1_phi);
      ntu[ii]->SetBranchAddress("jet1_e",&jet1_e);
      ntu[ii]->SetBranchAddress("jet2_flavour",&jet2_flavour);
      ntu[ii]->SetBranchAddress("jet2_csvBtag",&jet2_csvBtag);
      ntu[ii]->SetBranchAddress("jet2_btagSF_M",&jet2_btagSF_M);
      ntu[ii]->SetBranchAddress("jet2_btagSFErrorUp_M",&jet2_btagSFErrorUp_M);
      ntu[ii]->SetBranchAddress("jet2_btagSFErrorDown_M",&jet2_btagSFErrorDown_M);
      ntu[ii]->SetBranchAddress("jet2_btagEff_M",&jet2_btagEff_M);
      ntu[ii]->SetBranchAddress("jet2_btagEffError_M",&jet2_btagEffError_M);
      ntu[ii]->SetBranchAddress("jet2_pt",&jet2_pt);
      ntu[ii]->SetBranchAddress("jet2_eta",&jet2_eta);
      ntu[ii]->SetBranchAddress("jet2_phi",&jet2_phi);
      ntu[ii]->SetBranchAddress("jet2_e",&jet2_e);
      ntu[ii]->SetBranchAddress("cut_based_ct",&cut_based_ct);
      ntu[ii]->SetBranchAddress("weightBtagSF",&weightBtagSF);
      ntu[ii]->SetBranchAddress("weightBtagSFerrUp",&weightBtagSFerrUp);
      ntu[ii]->SetBranchAddress("weightBtagSFerrDown",&weightBtagSFerrDown);
  }
  
  TH1F* h_4body_0btag = new TH1F("h_4body_0btag","h_4body_0btag",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag = new TH1F("h_4body_1btag","h_4body_1btag",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag = new TH1F("h_4body_2btag","h_4body_2btag",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_up = new TH1F("h_4body_0btag_up","h_4body_0btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_up = new TH1F("h_4body_1btag_up","h_4body_1btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_up = new TH1F("h_4body_2btag_up","h_4body_2btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_down = new TH1F("h_4body_0btag_down","h_4body_0btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_down = new TH1F("h_4body_1btag_down","h_4body_1btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_down = new TH1F("h_4body_2btag_down","h_4body_2btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_05up = new TH1F("h_4body_0btag_05up","h_4body_0btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_05up = new TH1F("h_4body_1btag_05up","h_4body_1btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_05up = new TH1F("h_4body_2btag_05up","h_4body_2btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_05down = new TH1F("h_4body_0btag_05down","h_4body_0btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_05down = new TH1F("h_4body_1btag_05down","h_4body_1btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_05down = new TH1F("h_4body_2btag_05down","h_4body_2btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_15up = new TH1F("h_4body_0btag_15up","h_4body_0btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_15up = new TH1F("h_4body_1btag_15up","h_4body_1btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_15up = new TH1F("h_4body_2btag_15up","h_4body_2btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_15down = new TH1F("h_4body_0btag_15down","h_4body_0btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_15down = new TH1F("h_4body_1btag_15down","h_4body_1btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_15down = new TH1F("h_4body_2btag_15down","h_4body_2btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_20up = new TH1F("h_4body_0btag_20up","h_4body_0btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_20up = new TH1F("h_4body_1btag_20up","h_4body_1btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_20up = new TH1F("h_4body_2btag_20up","h_4body_2btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_20down = new TH1F("h_4body_0btag_20down","h_4body_0btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_20down = new TH1F("h_4body_1btag_20down","h_4body_1btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_20down = new TH1F("h_4body_2btag_20down","h_4body_2btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_25up = new TH1F("h_4body_0btag_25up","h_4body_0btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_25up = new TH1F("h_4body_1btag_25up","h_4body_1btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_25up = new TH1F("h_4body_2btag_25up","h_4body_2btag_up",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_0btag_25down = new TH1F("h_4body_0btag_25down","h_4body_0btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_1btag_25down = new TH1F("h_4body_1btag_25down","h_4body_1btag_down",int((xMax-xMin)/4.),xMin,xMax);
  TH1F* h_4body_2btag_25down = new TH1F("h_4body_2btag_25down","h_4body_2btag_down",int((xMax-xMin)/4.),xMin,xMax);

  TH1F* h_Mjj_nob_1btag = new TH1F("h_Mjj_nob_1btag","h_Mjj_nob_1btag",50,75.,175.);
  TH1F* h_Mjj_1b_1btag = new TH1F("h_Mjj_1b_1btag","h_Mjj_1b_1btag",50,75.,175.); 
  TH1F* h_Mjj_2b_1btag = new TH1F("h_Mjj_2b_1btag","h_Mjj_2b_1btag",50,75.,175.); 
  TH1F* h_Mjj_nob_2btag = new TH1F("h_Mjj_nob_2btag","h_Mjj_nob_2btag",50,75.,175.); 
  TH1F* h_Mjj_1b_2btag = new TH1F("h_Mjj_1b_2btag","h_Mjj_1b_2btag",50,75.,175.); 
  TH1F* h_Mjj_2b_2btag = new TH1F("h_Mjj_2b_2btag","h_Mjj_2b_2btag",50,75.,175.); 

  TH1F* h_weight_1btag = new TH1F("h_weight_1btag","h_weight_1btag",400,-2.,2.); 
  TH1F* h_weight_2btag = new TH1F("h_weight_2btag","h_weight_2btag",400,-2.,2.); 
  TH1F* h_weightUP_1btag = new TH1F("h_weightUP_1btag","h_weightUP_1btag",400,-2.,2.); 
  TH1F* h_weightUP_2btag = new TH1F("h_weightUP_2btag","h_weightUP_2btag",400,-2.,2.); 
  TH1F* h_weightDOWN_1btag = new TH1F("h_weightDOWN_1btag","h_weightDOWN_1btag",400,-2.,2.); 
  TH1F* h_weightDOWN_2btag = new TH1F("h_weightDOWN_2btag","h_weightDOWN_2btag",400,-2.,2.); 
  
  TH1F* h_jet1pt_1btag = new TH1F("h_jet1pt_1btag","h_jet1pt_1btag",500,0.,500.); 
  TH1F* h_jet2pt_1btag = new TH1F("h_jet2pt_1btag","h_jet2pt_1btag",500,0.,500.); 
  TH1F* h_jet1pt_2btag = new TH1F("h_jet1pt_2btag","h_jet1pt_2btag",500,0.,500.); 
  TH1F* h_jet2pt_2btag = new TH1F("h_jet2pt_2btag","h_jet2pt_2btag",500,0.,500.); 
  
   
  TLorentzVector* jet1_p4 = new TLorentzVector();
  TLorentzVector* jet2_p4 = new TLorentzVector();
  TLorentzVector* sumP4 = new TLorentzVector(); 

  TGraph* corrUp_1b2b = new TGraph();
  TGraph* corrDown_1b2b = new TGraph();
  TGraph* corrUp_0b1b = new TGraph();
  TGraph* corrDown_0b1b = new TGraph();
  TGraph* corrUp_0b2b = new TGraph();
  TGraph* corrDown_0b2b = new TGraph();

  int btag1_n = 0;
  int btag1_n_up_plus = 0;
  int btag1_n_up_minus = 0;
  int btag1_n_down_plus = 0;
  int btag1_n_down_minus = 0;

  float sumWeight_1btagUP = 0.;
  float sumWeight_1btagDOWN = 0.;

  int btag2_n = 0;
  int btag2_n_up_plus = 0;
  int btag2_n_up_minus = 0;
  int btag2_n_down_plus = 0;
  int btag2_n_down_minus = 0;
  
  float sumWeight_2btagUP = 0.;
  float sumWeight_2btagDOWN = 0.;
  
  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%10000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

          float weight_tmp = eventWeight_2jets("medium", jet1_btagSF_M, jet2_btagSF_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          float weight_errorUp_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, jet1_btagSFErrorUp_M, jet2_btagSF_M, jet2_btagSFErrorUp_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_errorDown_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, jet1_btagSFErrorDown_M, jet2_btagSF_M, jet2_btagSFErrorDown_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);

          float weight_NAIVE_UP = eventWeight_NAIVE_2jets_UP("medium", jet1_btagSF_M, jet1_btagSFErrorUp_M, jet2_btagSF_M, jet2_btagSFErrorUp_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          float weight_NAIVE_DOWN = eventWeight_NAIVE_2jets_DOWN("medium", jet1_btagSF_M, jet1_btagSFErrorDown_M, jet2_btagSF_M, jet2_btagSFErrorDown_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);


          float weight_NAIVE_UP_tagged = eventWeight_NAIVE_2jets_UP_tagged("medium", jet1_btagSF_M, jet1_btagSFErrorUp_M, jet2_btagSF_M, jet2_btagSFErrorUp_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          float weight_NAIVE_UP_untagged = eventWeight_NAIVE_2jets_UP_untagged("medium", jet1_btagSF_M, jet1_btagSFErrorUp_M, jet2_btagSF_M, jet2_btagSFErrorUp_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          float weight_NAIVE_DOWN_tagged = eventWeight_NAIVE_2jets_DOWN_tagged("medium", jet1_btagSF_M, jet1_btagSFErrorDown_M, jet2_btagSF_M, jet2_btagSFErrorDown_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);
          float weight_NAIVE_DOWN_untagged = eventWeight_NAIVE_2jets_DOWN_untagged("medium", jet1_btagSF_M, jet1_btagSFErrorDown_M, jet2_btagSF_M, jet2_btagSFErrorDown_M, jet1_btagEff_M, jet2_btagEff_M, jet1_csvBtag, jet2_csvBtag);

          float signUP = (weight_NAIVE_UP-weight_tmp)/fabs(weight_NAIVE_UP-weight_tmp); 
          float signDOWN = (weight_NAIVE_DOWN-weight_tmp)/fabs(weight_NAIVE_DOWN-weight_tmp); 
          
          float signUP_old = signUP;
          float signDOWN_old = signDOWN;

          if(fabs(jet1_flavour) == 5 && fabs(jet2_flavour) != 5 && fabs(jet2_flavour) != 4){
             signUP = (weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp)/fabs(weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp);
             signDOWN = (weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp)/fabs(weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp);
          }
          if(fabs(jet1_flavour) == 4 && fabs(jet2_flavour) != 5 && fabs(jet2_flavour) != 4){
             signUP = (weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp)/fabs(weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp);
             signDOWN = (weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp)/fabs(weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp);
          } 
          if(fabs(jet2_flavour) == 5 && fabs(jet1_flavour) != 5 && fabs(jet1_flavour) != 4){
             signUP = (weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp)/fabs(weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp);
             signDOWN = (weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp)/fabs(weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp);
          }
          if(fabs(jet2_flavour) == 4 && fabs(jet1_flavour) != 5 && fabs(jet1_flavour) != 4){
             signUP = (weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp)/fabs(weight_NAIVE_UP_tagged+weight_NAIVE_UP_untagged-2*weight_tmp);
             signDOWN = (weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp)/fabs(weight_NAIVE_DOWN_tagged+weight_NAIVE_DOWN_untagged-2*weight_tmp);
          }
 
          float weight_05errorUp_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 0.5*jet1_btagSFErrorUp_M, jet2_btagSF_M, 0.5*jet2_btagSFErrorUp_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_15errorUp_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 1.5*jet1_btagSFErrorUp_M, jet2_btagSF_M, 1.5*jet2_btagSFErrorUp_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_20errorUp_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 2.0*jet1_btagSFErrorUp_M, jet2_btagSF_M, 2.0*jet2_btagSFErrorUp_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_25errorUp_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 2.5*jet1_btagSFErrorUp_M, jet2_btagSF_M, 2.5*jet2_btagSFErrorUp_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          
          float weight_05errorDown_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 0.5*jet1_btagSFErrorDown_M, jet2_btagSF_M, 0.5*jet2_btagSFErrorDown_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_15errorDown_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 1.5*jet1_btagSFErrorDown_M, jet2_btagSF_M, 1.5*jet2_btagSFErrorDown_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_20errorDown_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 2.0*jet1_btagSFErrorDown_M, jet2_btagSF_M, 2.0*jet2_btagSFErrorDown_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);
          float weight_25errorDown_tmp = eventWeight_error_2jets("medium",jet1_btagSF_M, 2.5*jet1_btagSFErrorDown_M, jet2_btagSF_M, 2.5*jet2_btagSFErrorDown_M, jet1_btagEff_M, jet1_btagEffError_M, jet2_btagEff_M, jet2_btagEffError_M, jet1_flavour, jet2_flavour, jet1_csvBtag, jet2_csvBtag);

          jet1_p4->SetPtEtaPhiE(jet1_pt,jet1_eta,jet1_phi,jet1_e);
          jet2_p4->SetPtEtaPhiE(jet2_pt,jet2_eta,jet2_phi,jet2_e);
          *sumP4 = *jet1_p4 + *jet2_p4;
          
          if(cut_based_ct == 1){

             btag1_n++;
             if(signUP > 0) btag1_n_up_plus++;
             if(signUP < 0) btag1_n_up_minus++;
             if(signDOWN > 0) btag1_n_down_plus++;
             if(signDOWN < 0) btag1_n_down_minus++;

             h_jet1pt_1btag->Fill(jet1_pt);
             h_jet2pt_1btag->Fill(jet2_pt);
           
             h_4body_1btag_up->Fill(mtot,evWeight/(1*weightBtagSF)*(weightBtagSF+signUP*weightBtagSFerrUp));
             h_4body_1btag->Fill(mtot,evWeight/1.);
             h_4body_1btag_down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weightBtagSFerrDown));
             h_4body_1btag_05up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_05errorUp_tmp));
             h_4body_1btag_05down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_05errorDown_tmp));
             h_4body_1btag_15up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_15errorUp_tmp));
             h_4body_1btag_15down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_15errorDown_tmp));
             h_4body_1btag_20up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_20errorUp_tmp));
             h_4body_1btag_20down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_20errorDown_tmp));
             h_4body_1btag_25up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_25errorUp_tmp));
             h_4body_1btag_25down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_25errorDown_tmp));

             h_weightUP_1btag->Fill(weightBtagSF+signUP*weightBtagSFerrUp);
             h_weight_1btag->Fill(weightBtagSF);
             h_weightDOWN_1btag->Fill(weightBtagSF+signDOWN*weightBtagSFerrDown);

             if(fabs(jet1_flavour) != 5 ) n_nob_1btag++;
             if(fabs(jet2_flavour) != 5 ) n_nob_1btag++;
             n_1btag++;
             n_1btag++;

             if((fabs(jet1_flavour) != 5 && fabs(jet2_flavour) != 5)) h_Mjj_nob_1btag->Fill(sumP4->M());
             if((fabs(jet1_flavour) != 5 && fabs(jet2_flavour) == 5) || (fabs(jet1_flavour) == 5 && fabs(jet2_flavour) != 5)) h_Mjj_1b_1btag->Fill(sumP4->M());
             if((fabs(jet1_flavour) == 5 && fabs(jet2_flavour) == 5)) h_Mjj_2b_1btag->Fill(sumP4->M());
          }

          if(cut_based_ct == 0){

             
             btag2_n++;
             if(signUP > 0) btag2_n_up_plus++;
             if(signUP < 0) btag2_n_up_minus++;
             if(signDOWN > 0) btag2_n_down_plus++;
             if(signDOWN < 0) btag2_n_down_minus++;

             h_jet1pt_2btag->Fill(jet1_pt);
             h_jet2pt_2btag->Fill(jet2_pt);

             h_4body_2btag_up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weightBtagSFerrUp));
             h_4body_2btag->Fill(mtot,evWeight/1.);
             h_4body_2btag_down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weightBtagSFerrDown));
             h_4body_2btag_05up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_05errorUp_tmp));
             h_4body_2btag_05down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_05errorDown_tmp));
             h_4body_2btag_15up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_15errorUp_tmp));
             h_4body_2btag_15down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_15errorDown_tmp));
             h_4body_2btag_20up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_20errorUp_tmp));
             h_4body_2btag_20down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_20errorDown_tmp));
             h_4body_2btag_25up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_25errorUp_tmp));
             h_4body_2btag_25down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_25errorDown_tmp));

             h_weightUP_2btag->Fill(weightBtagSF+signUP*weightBtagSFerrUp);
             h_weight_2btag->Fill(weightBtagSF);
             h_weightDOWN_2btag->Fill(weightBtagSF+signDOWN*weightBtagSFerrDown);

             if(fabs(jet1_flavour) != 5 ) n_nob_2btag++;
             if(fabs(jet2_flavour) != 5 ) n_nob_2btag++;
             n_2btag++;
             n_2btag++;

             if((fabs(jet1_flavour) != 5 && fabs(jet2_flavour) != 5)) h_Mjj_nob_2btag->Fill(sumP4->M()); 
             if((fabs(jet1_flavour) != 5 && fabs(jet2_flavour) == 5) || (fabs(jet1_flavour) == 5 && fabs(jet2_flavour) != 5)) h_Mjj_1b_2btag->Fill(sumP4->M());
             if((fabs(jet1_flavour) == 5 && fabs(jet2_flavour) == 5)) h_Mjj_2b_2btag->Fill(sumP4->M());
          }

          if(cut_based_ct == 2){

             h_4body_0btag_up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weightBtagSFerrUp));
             h_4body_0btag->Fill(mtot,evWeight/1.);
             h_4body_0btag_down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weightBtagSFerrDown));
             h_4body_0btag_05up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_05errorUp_tmp));
             h_4body_0btag_05down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_05errorDown_tmp));
             h_4body_0btag_15up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_15errorUp_tmp));
             h_4body_0btag_15down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_15errorDown_tmp));
             h_4body_0btag_20up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_20errorUp_tmp));
             h_4body_0btag_20down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_20errorDown_tmp));
             h_4body_0btag_25up->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signUP*weight_25errorUp_tmp));
             h_4body_0btag_25down->Fill(mtot,evWeight/(1.*weightBtagSF)*(weightBtagSF+signDOWN*weight_25errorDown_tmp));
          }
      } 
  }

  std::cout << "N non-b jet 1btag: " << n_nob_1btag << std::endl;
  std::cout << "N jet 1btag: " << n_1btag << std::endl;
  std::cout << "N non-b jet 2btag: " << n_nob_2btag << std::endl;
  std::cout << "N jet 2btag: " << n_2btag << std::endl;

  std::cout << "Syst 1btag: " << (h_4body_1btag_up->Integral()-h_4body_1btag_down->Integral())/(2*h_4body_1btag->Integral())*100 << std::endl;
  std::cout << "Syst 2btag: " << (h_4body_2btag_up->Integral()-h_4body_2btag_down->Integral())/(2*h_4body_2btag->Integral())*100 << std::endl;

  std::cout << "signUP > 0: " << float(btag1_n_up_plus)/float(btag1_n)*100 << " - " << float(btag2_n_up_plus)/float(btag2_n)*100 << std::endl;
  std::cout << "signUP < 0: " << float(btag1_n_up_minus)/float(btag1_n)*100 << " - " << float(btag2_n_up_minus)/float(btag2_n)*100 << std::endl;
  std::cout << "signDOWN > 0: " << float(btag1_n_down_plus)/float(btag1_n)*100 << " - " << float(btag2_n_down_plus)/float(btag2_n)*100 << std::endl;
  std::cout << "signDOWN < 0: " << float(btag1_n_down_minus)/float(btag1_n)*100 << " - " << float(btag2_n_down_minus)/float(btag2_n)*100 << std::endl;

  float Nevents_mean_0btag = h_4body_0btag->Integral();
  float Nevents_mean_1btag = h_4body_1btag->Integral();
  float Nevents_mean_2btag = h_4body_2btag->Integral();

  float Nevents_05up_0btag = h_4body_0btag_05up->Integral();
  float Nevents_05up_1btag = h_4body_1btag_05up->Integral();
  float Nevents_05up_2btag = h_4body_2btag_05up->Integral();
  float Nevents_05down_0btag = h_4body_0btag_05down->Integral();
  float Nevents_05down_1btag = h_4body_1btag_05down->Integral();
  float Nevents_05down_2btag = h_4body_2btag_05down->Integral();

  float Nevents_up_0btag = h_4body_0btag_up->Integral();
  float Nevents_up_1btag = h_4body_1btag_up->Integral();
  float Nevents_up_2btag = h_4body_2btag_up->Integral();
  float Nevents_down_0btag = h_4body_0btag_down->Integral();
  float Nevents_down_1btag = h_4body_1btag_down->Integral();
  float Nevents_down_2btag = h_4body_2btag_down->Integral();

  float Nevents_15up_0btag = h_4body_0btag_15up->Integral();
  float Nevents_15up_1btag = h_4body_1btag_15up->Integral();
  float Nevents_15up_2btag = h_4body_2btag_15up->Integral();
  float Nevents_15down_0btag = h_4body_0btag_15down->Integral();
  float Nevents_15down_1btag = h_4body_1btag_15down->Integral();
  float Nevents_15down_2btag = h_4body_2btag_15down->Integral();

  float Nevents_20up_0btag = h_4body_0btag_20up->Integral();
  float Nevents_20up_1btag = h_4body_1btag_20up->Integral();
  float Nevents_20up_2btag = h_4body_2btag_20up->Integral();
  float Nevents_20down_0btag = h_4body_0btag_20down->Integral();
  float Nevents_20down_1btag = h_4body_1btag_20down->Integral();
  float Nevents_20down_2btag = h_4body_2btag_20down->Integral();

  float Nevents_25up_0btag = h_4body_0btag_25up->Integral();
  float Nevents_25up_1btag = h_4body_1btag_25up->Integral();
  float Nevents_25up_2btag = h_4body_2btag_25up->Integral();
  float Nevents_25down_0btag = h_4body_0btag_25down->Integral();
  float Nevents_25down_1btag = h_4body_1btag_25down->Integral();
  float Nevents_25down_2btag = h_4body_2btag_25down->Integral();

  corrDown_1b2b->SetPoint(0,Nevents_25down_1btag-Nevents_mean_1btag,Nevents_25down_2btag-Nevents_mean_2btag);
  corrDown_1b2b->SetPoint(1,Nevents_20down_1btag-Nevents_mean_1btag,Nevents_20down_2btag-Nevents_mean_2btag);
  corrDown_1b2b->SetPoint(2,Nevents_15down_1btag-Nevents_mean_1btag,Nevents_15down_2btag-Nevents_mean_2btag);
  corrDown_1b2b->SetPoint(3,Nevents_down_1btag-Nevents_mean_1btag,Nevents_down_2btag-Nevents_mean_2btag);
  corrDown_1b2b->SetPoint(4,Nevents_05down_1btag-Nevents_mean_1btag,Nevents_05down_2btag-Nevents_mean_2btag);
  corrDown_1b2b->SetPoint(5,0.,0.);

  corrDown_0b1b->SetPoint(0,Nevents_25down_0btag-Nevents_mean_0btag,Nevents_25down_1btag-Nevents_mean_1btag);
  corrDown_0b1b->SetPoint(1,Nevents_20down_0btag-Nevents_mean_0btag,Nevents_20down_1btag-Nevents_mean_1btag);
  corrDown_0b1b->SetPoint(2,Nevents_15down_0btag-Nevents_mean_0btag,Nevents_15down_1btag-Nevents_mean_1btag);
  corrDown_0b1b->SetPoint(3,Nevents_down_0btag-Nevents_mean_0btag,Nevents_down_1btag-Nevents_mean_1btag);
  corrDown_0b1b->SetPoint(4,Nevents_05down_0btag-Nevents_mean_0btag,Nevents_05down_1btag-Nevents_mean_1btag);
  corrDown_0b1b->SetPoint(5,0.,0.);
  
  corrDown_0b2b->SetPoint(0,Nevents_25down_0btag-Nevents_mean_0btag,Nevents_25down_2btag-Nevents_mean_2btag);
  corrDown_0b2b->SetPoint(1,Nevents_20down_0btag-Nevents_mean_0btag,Nevents_20down_2btag-Nevents_mean_2btag);
  corrDown_0b2b->SetPoint(2,Nevents_15down_0btag-Nevents_mean_0btag,Nevents_15down_2btag-Nevents_mean_2btag);
  corrDown_0b2b->SetPoint(3,Nevents_down_0btag-Nevents_mean_0btag,Nevents_down_2btag-Nevents_mean_2btag);
  corrDown_0b2b->SetPoint(4,Nevents_05down_0btag-Nevents_mean_0btag,Nevents_05down_2btag-Nevents_mean_2btag);
  corrDown_0b2b->SetPoint(5,0.,0.);
 
  corrUp_1b2b->SetPoint(0,0.,0.);
  corrUp_1b2b->SetPoint(1,Nevents_05up_1btag-Nevents_mean_1btag,Nevents_05up_2btag-Nevents_mean_2btag);
  corrUp_1b2b->SetPoint(2,Nevents_up_1btag-Nevents_mean_1btag,Nevents_up_2btag-Nevents_mean_2btag);
  corrUp_1b2b->SetPoint(3,Nevents_15up_1btag-Nevents_mean_1btag,Nevents_15up_2btag-Nevents_mean_2btag);
  corrUp_1b2b->SetPoint(4,Nevents_20up_1btag-Nevents_mean_1btag,Nevents_20up_2btag-Nevents_mean_2btag);
  corrUp_1b2b->SetPoint(5,Nevents_25up_1btag-Nevents_mean_1btag,Nevents_25up_2btag-Nevents_mean_2btag);
  
  corrUp_0b1b->SetPoint(0,0.,0.);
  corrUp_0b1b->SetPoint(1,Nevents_05up_0btag-Nevents_mean_0btag,Nevents_05up_1btag-Nevents_mean_1btag);
  corrUp_0b1b->SetPoint(2,Nevents_up_0btag-Nevents_mean_0btag,Nevents_up_1btag-Nevents_mean_1btag);
  corrUp_0b1b->SetPoint(3,Nevents_15up_0btag-Nevents_mean_0btag,Nevents_15up_1btag-Nevents_mean_1btag);
  corrUp_0b1b->SetPoint(4,Nevents_20up_0btag-Nevents_mean_0btag,Nevents_20up_1btag-Nevents_mean_1btag);
  corrUp_0b1b->SetPoint(5,Nevents_25up_0btag-Nevents_mean_0btag,Nevents_25up_1btag-Nevents_mean_1btag);
   
  corrUp_0b2b->SetPoint(0,0.,0.);
  corrUp_0b2b->SetPoint(1,Nevents_05up_0btag-Nevents_mean_0btag,Nevents_05up_2btag-Nevents_mean_2btag);
  corrUp_0b2b->SetPoint(2,Nevents_up_0btag-Nevents_mean_0btag,Nevents_up_2btag-Nevents_mean_2btag);
  corrUp_0b2b->SetPoint(3,Nevents_15up_0btag-Nevents_mean_0btag,Nevents_15up_2btag-Nevents_mean_2btag);
  corrUp_0b2b->SetPoint(4,Nevents_20up_0btag-Nevents_mean_0btag,Nevents_20up_2btag-Nevents_mean_2btag);
  corrUp_0b2b->SetPoint(5,Nevents_25up_0btag-Nevents_mean_0btag,Nevents_25up_2btag-Nevents_mean_2btag);
   
  
  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  
  output->cd();

  h_4body_1btag_up->Write();
  h_4body_1btag->Write();
  h_4body_1btag_down->Write();
  
  h_4body_2btag_up->Write();
  h_4body_2btag->Write();
  h_4body_2btag_down->Write();

  h_Mjj_nob_1btag->Write();
  h_Mjj_1b_1btag->Write();
  h_Mjj_2b_1btag->Write();
  h_Mjj_nob_2btag->Write();
  h_Mjj_1b_2btag->Write();
  h_Mjj_2b_2btag->Write();

  h_weightUP_1btag->Write();
  h_weight_1btag->Write();
  h_weightDOWN_1btag->Write();

  h_weightUP_2btag->Write();
  h_weight_2btag->Write();
  h_weightDOWN_2btag->Write(); 

  h_jet1pt_1btag->Write();
  h_jet2pt_1btag->Write();
  h_jet1pt_2btag->Write();
  h_jet2pt_2btag->Write();
  
  output->Close();

  drawPlotUP(corrUp_1b2b,outputPlotName+std::string("_1b2b_up"), std::string("1btag"), std::string("2btag"));
  drawPlotDOWN(corrDown_1b2b,outputPlotName+std::string("_1b2b_down"), std::string("1btag"), std::string("2btag"));

  drawPlotUP(corrUp_0b1b,outputPlotName+std::string("_0b1b_up"), std::string("0btag"), std::string("1btag"));
  drawPlotDOWN(corrDown_0b1b,outputPlotName+std::string("_0b1b_down"), std::string("0btag"), std::string("1btag"));
   
  drawPlotUP(corrUp_0b2b,outputPlotName+std::string("_0b2b_up"), std::string("0btag"), std::string("2btag"));
  drawPlotDOWN(corrDown_0b2b,outputPlotName+std::string("_0b2b_down"), std::string("0btag"), std::string("2btag"));

}

void drawPlotUP(TGraph* corr, std::string outputPlot, std::string cat_x, std::string cat_y){
    
  TCanvas* cplot = new TCanvas("gplot", "gplot",100,100,1100,600);
  cplot->cd();
  
  corr ->GetXaxis()-> SetTitle(("#DeltaEvents_"+cat_x).c_str());
  corr ->GetYaxis()-> SetTitle(("#DeltaEvents_"+cat_y).c_str());
  
  TLatex *latex0 = new TLatex(corr->GetX()[1], corr->GetY()[1],"+0.5#sigma");
  latex0->SetTextSize(0.03);
  latex0->SetTextColor(kBlue);
  corr->GetListOfFunctions()->Add(latex0);
  TLatex *latex1 = new TLatex(corr->GetX()[2], corr->GetY()[2],"+1.0#sigma");
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kBlue);
  corr->GetListOfFunctions()->Add(latex1);
  TLatex *latex2 = new TLatex(corr->GetX()[3], corr->GetY()[3],"+1.5#sigma");
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kBlue);
  corr->GetListOfFunctions()->Add(latex2);
  TLatex *latex3 = new TLatex(corr->GetX()[4], corr->GetY()[4],"+2.0#sigma");
  latex3->SetTextSize(0.03);
  latex3->SetTextColor(kBlue);
  corr->GetListOfFunctions()->Add(latex3);
  TLatex *latex4 = new TLatex(corr->GetX()[5], corr->GetY()[5],"+2.5#sigma");
  latex4->SetTextSize(0.03);
  latex4->SetTextColor(kBlue);
  corr->GetListOfFunctions()->Add(latex4);

  corr->Draw("ap");

  cplot -> Print((outputPlot+".png").c_str(),"png");
  cplot -> Print((outputPlot+".pdf").c_str(),"pdf");

  delete cplot;
}


void drawPlotDOWN(TGraph* corr, std::string outputPlot, std::string cat_x, std::string cat_y){
    
  TCanvas* cplot = new TCanvas("gplot", "gplot",100,100,1100,600);
  cplot->cd();
  
  corr ->GetXaxis()-> SetTitle(("#DeltaEvents_"+cat_x).c_str());
  corr ->GetYaxis()-> SetTitle(("#DeltaEvents_"+cat_y).c_str());
  
  TLatex *latex0 = new TLatex(corr->GetX()[0], corr->GetY()[0],"-2.5#sigma");
  latex0->SetTextSize(0.03);
  latex0->SetTextColor(kRed);
  corr->GetListOfFunctions()->Add(latex0);
  TLatex *latex1 = new TLatex(corr->GetX()[1], corr->GetY()[1],"-2.0#sigma");
  latex1->SetTextSize(0.03);
  latex1->SetTextColor(kRed);
  corr->GetListOfFunctions()->Add(latex1);
  TLatex *latex2 = new TLatex(corr->GetX()[2], corr->GetY()[2],"-1.5#sigma");
  latex2->SetTextSize(0.03);
  latex2->SetTextColor(kRed);
  corr->GetListOfFunctions()->Add(latex2);
  TLatex *latex3 = new TLatex(corr->GetX()[3], corr->GetY()[3],"-1.0#sigma");
  latex3->SetTextSize(0.03);
  latex3->SetTextColor(kRed);
  corr->GetListOfFunctions()->Add(latex3);
  TLatex *latex4 = new TLatex(corr->GetX()[4], corr->GetY()[4],"-0.5#sigma");
  latex4->SetTextSize(0.03);
  latex4->SetTextColor(kRed);
  corr->GetListOfFunctions()->Add(latex4);

  corr->Draw("ap");

  cplot -> Print((outputPlot+".png").c_str(),"png");
  cplot -> Print((outputPlot+".pdf").c_str(),"pdf");

  delete cplot;
}

