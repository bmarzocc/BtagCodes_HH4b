
#include "ConfigParser.h"
#include "ParserUtils.h"
#include "setTDRStyle.h"

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
  if(argc != 3)
  {
    std::cerr << ">>>>> Analyses::usage: " << argv[0] << " configFileName" << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputList = gConfigParser -> readStringOption("Input::inputList");
  std::string inputTree = gConfigParser -> readStringOption("Input::inputTree");
 
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");

  char* inputJetChoice = argv[2];
  std::string(JetChoice) = inputJetChoice;

  std::cout << "inputList      = " << inputList << std::endl;
  std::cout << "inputTree      = " << inputTree << std::endl;
  std::cout << "outputName     = " << outputName << std::endl; 
  std::cout << "inputJetChoice = " << inputJetChoice << std::endl;
  
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

  for(int ii = 0; ii < pos_total; ii++){
      ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("Jet1pt",1);
      ntu[ii]->SetBranchStatus("Jet2pt",1);
      ntu[ii]->SetBranchStatus("Jet1eta",1);
      ntu[ii]->SetBranchStatus("Jet2eta",1);
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
      ntu[ii]->SetBranchStatus("Jet1_daughter1_bjetness",1);
      ntu[ii]->SetBranchAddress("Jet1pt",&Jet1pt);
      ntu[ii]->SetBranchAddress("Jet2pt",&Jet2pt);
      ntu[ii]->SetBranchAddress("Jet1eta",&Jet1eta);
      ntu[ii]->SetBranchAddress("Jet2eta",&Jet2eta);
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

  float ptmin[21] = {20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 600., 800.,900.,1000.,1200.,1500.};
  float etamin[8] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5};

  TH2F* h2_BTaggingEff_Denom_b_L = new TH2F("h2_BTaggingEff_Denom_b_L","h2_BTaggingEff_Denom_b_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_b_L = new TH2F("h2_BTaggingEff_Num_b_L","h2_BTaggingEff_Num_b_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_b_L = new TH2F("h2_BTaggingEff_b_L","h2_BTaggingEff_b_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Denom_b_M = new TH2F("h2_BTaggingEff_Denom_b_M","h2_BTaggingEff_Denom_b_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_b_M = new TH2F("h2_BTaggingEff_Num_b_M","h2_BTaggingEff_Num_b_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_b_M = new TH2F("h2_BTaggingEff_b_M","h2_BTaggingEff_b_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Denom_b_T = new TH2F("h2_BTaggingEff_Denom_b_T","h2_BTaggingEff_Denom_b_T",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_b_T = new TH2F("h2_BTaggingEff_Num_b_T","h2_BTaggingEff_Num_b_T",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_b_T = new TH2F("h2_BTaggingEff_b_T","h2_BTaggingEff_b_T",20, ptmin, 7,etamin);

  TH2F* h2_BTaggingEff_Denom_c_L = new TH2F("h2_BTaggingEff_Denom_c_L","h2_BTaggingEff_Denom_c_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_c_L = new TH2F("h2_BTaggingEff_Num_c_L","h2_BTaggingEff_Num_c_L",20, ptmin,7,etamin);
  TH2F* h2_BTaggingEff_c_L = new TH2F("h2_BTaggingEff_c_L","h2_BTaggingEff_c_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Denom_c_M = new TH2F("h2_BTaggingEff_Denom_c_M","h2_BTaggingEff_Denom_c_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_c_M = new TH2F("h2_BTaggingEff_Num_c_M","h2_BTaggingEff_Num_c_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_c_M = new TH2F("h2_BTaggingEff_c_M","h2_BTaggingEff_c_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Denom_c_T = new TH2F("h2_BTaggingEff_Denom_c_T","h2_BTaggingEff_Denom_c_T",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_c_T = new TH2F("h2_BTaggingEff_Num_c_T","h2_BTaggingEff_Num_c_T",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_c_T = new TH2F("h2_BTaggingEff_c_T","h2_BTaggingEff_c_T",20, ptmin, 7,etamin);

  TH2F* h2_BTaggingEff_Denom_udsg_L = new TH2F("h2_BTaggingEff_Denom_udsg_L","h2_BTaggingEff_Denom_udsg_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_udsg_L = new TH2F("h2_BTaggingEff_Num_udsg_L","h2_BTaggingEff_Num_udsg_L",20, ptmin,7,etamin);
  TH2F* h2_BTaggingEff_udsg_L = new TH2F("h2_BTaggingEff_udsg_L","h2_BTaggingEff_udsg_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Denom_udsg_M = new TH2F("h2_BTaggingEff_Denom_udsg_M","h2_BTaggingEff_Denom_udsg_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_udsg_M = new TH2F("h2_BTaggingEff_Num_udsg_M","h2_BTaggingEff_Num_udsg_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_udsg_M = new TH2F("h2_BTaggingEff_udsg_M","h2_BTaggingEff_udsg_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Denom_udsg_T = new TH2F("h2_BTaggingEff_Denom_udsg_T","h2_BTaggingEff_Denom_c_T",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_Num_udsg_T = new TH2F("h2_BTaggingEff_Num_udsg_T","h2_BTaggingEff_Num_udsg_T",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_udsg_T = new TH2F("h2_BTaggingEff_udsg_T","h2_BTaggingEff_udsg_T",20, ptmin, 7,etamin);

  float jetPt[4];
  float jetEta[4];
  float jetcsvBtag[4];
  float jetFlavour[4]; 
  int njets = 4;

  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%1000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry); 

          ChooseJets(&JetChoice,jetPt,jetEta,jetcsvBtag,jetFlavour,njets);
        
          float puRe = 1.;

          for(int ii = 0; ii < njets; ii++){

              if(abs(jetFlavour[ii]) == 5){
                 
                 h2_BTaggingEff_Denom_b_L->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.244) h2_BTaggingEff_Num_b_L->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                   
                 h2_BTaggingEff_Denom_b_M->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.679) h2_BTaggingEff_Num_b_M->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 
                 h2_BTaggingEff_Denom_b_T->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.898) h2_BTaggingEff_Num_b_T->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
              }
 
              if(abs(jetFlavour[ii]) == 4){
                 
                 h2_BTaggingEff_Denom_c_L->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.244) h2_BTaggingEff_Num_c_L->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                   
                 h2_BTaggingEff_Denom_c_M->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.679) h2_BTaggingEff_Num_c_M->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 
                 h2_BTaggingEff_Denom_c_T->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.898) h2_BTaggingEff_Num_c_T->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
              }

              if(abs(jetFlavour[ii]) == 1 || abs(jetFlavour[ii]) == 2 || abs(jetFlavour[ii]) == 3 || abs(jetFlavour[ii]) == 21){
                 
                 h2_BTaggingEff_Denom_udsg_L->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.244) h2_BTaggingEff_Num_udsg_L->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                   
                 h2_BTaggingEff_Denom_udsg_M->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.679) h2_BTaggingEff_Num_udsg_M->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 
                 h2_BTaggingEff_Denom_udsg_T->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
                 if(jetcsvBtag[ii] > 0.898) h2_BTaggingEff_Num_udsg_T->Fill(jetPt[ii],fabs(jetEta[ii]),puRe);
              }
 
          }

      } 
  }
  
  h2_BTaggingEff_Num_b_L->Sumw2();
  h2_BTaggingEff_Denom_b_L->Sumw2();
  h2_BTaggingEff_b_L->Sumw2();
  h2_BTaggingEff_Num_c_L->Sumw2();
  h2_BTaggingEff_Denom_c_L->Sumw2();
  h2_BTaggingEff_c_L->Sumw2();
  h2_BTaggingEff_Num_udsg_L->Sumw2();
  h2_BTaggingEff_Denom_udsg_L->Sumw2();
  h2_BTaggingEff_udsg_L->Sumw2();
  h2_BTaggingEff_Num_b_M->Sumw2();
  h2_BTaggingEff_Denom_b_M->Sumw2();
  h2_BTaggingEff_b_M->Sumw2();
  h2_BTaggingEff_Num_c_M->Sumw2();
  h2_BTaggingEff_Denom_c_M->Sumw2();
  h2_BTaggingEff_c_M->Sumw2();
  h2_BTaggingEff_Num_udsg_M->Sumw2();
  h2_BTaggingEff_Denom_udsg_M->Sumw2();
  h2_BTaggingEff_udsg_M->Sumw2();
  h2_BTaggingEff_Num_b_T->Sumw2();
  h2_BTaggingEff_Denom_b_T->Sumw2();
  h2_BTaggingEff_b_T->Sumw2();
  h2_BTaggingEff_Num_c_T->Sumw2();
  h2_BTaggingEff_Denom_c_T->Sumw2();
  h2_BTaggingEff_c_T->Sumw2();
  h2_BTaggingEff_Num_udsg_T->Sumw2();
  h2_BTaggingEff_Denom_udsg_T->Sumw2();
  h2_BTaggingEff_udsg_T->Sumw2();

  h2_BTaggingEff_b_L->Divide(h2_BTaggingEff_Num_b_L,h2_BTaggingEff_Denom_b_L,1,1,"B");
  h2_BTaggingEff_b_M->Divide(h2_BTaggingEff_Num_b_M,h2_BTaggingEff_Denom_b_M,1,1,"B");
  h2_BTaggingEff_b_T->Divide(h2_BTaggingEff_Num_b_T,h2_BTaggingEff_Denom_b_T,1,1,"B");

  h2_BTaggingEff_c_L->Divide(h2_BTaggingEff_Num_c_L,h2_BTaggingEff_Denom_c_L,1,1,"B");
  h2_BTaggingEff_c_M->Divide(h2_BTaggingEff_Num_c_M,h2_BTaggingEff_Denom_c_M,1,1,"B");
  h2_BTaggingEff_c_T->Divide(h2_BTaggingEff_Num_c_T,h2_BTaggingEff_Denom_c_T,1,1,"B");
  
  h2_BTaggingEff_udsg_L->Divide(h2_BTaggingEff_Num_udsg_L,h2_BTaggingEff_Denom_udsg_L,1,1,"B");
  h2_BTaggingEff_udsg_M->Divide(h2_BTaggingEff_Num_udsg_M,h2_BTaggingEff_Denom_udsg_M,1,1,"B");
  h2_BTaggingEff_udsg_T->Divide(h2_BTaggingEff_Num_udsg_T,h2_BTaggingEff_Denom_udsg_T,1,1,"B");
  
  TFile* output = new TFile((outputName+"_"+JetChoice+".root").c_str(),"RECREATE");
  
  output->cd();

  h2_BTaggingEff_Denom_b_L->Write();
  h2_BTaggingEff_Num_b_L->Write();
  h2_BTaggingEff_b_L->Write();
  h2_BTaggingEff_Denom_b_M->Write();
  h2_BTaggingEff_Num_b_M->Write();
  h2_BTaggingEff_b_M->Write();
  h2_BTaggingEff_Denom_b_T->Write();
  h2_BTaggingEff_Num_b_T->Write();
  h2_BTaggingEff_b_T->Write();
  h2_BTaggingEff_Denom_c_L->Write();
  h2_BTaggingEff_Num_c_L->Write();
  h2_BTaggingEff_c_L->Write();
  h2_BTaggingEff_Denom_c_M->Write();
  h2_BTaggingEff_Num_c_M->Write();
  h2_BTaggingEff_c_M->Write();
  h2_BTaggingEff_Denom_c_T->Write();
  h2_BTaggingEff_Num_c_T->Write();
  h2_BTaggingEff_c_T->Write();
  h2_BTaggingEff_Denom_udsg_L->Write();
  h2_BTaggingEff_Num_udsg_L->Write();
  h2_BTaggingEff_udsg_L->Write();
  h2_BTaggingEff_Denom_udsg_M->Write();
  h2_BTaggingEff_Num_udsg_M->Write();
  h2_BTaggingEff_udsg_M->Write();
  h2_BTaggingEff_Denom_udsg_T->Write();
  h2_BTaggingEff_Num_udsg_T->Write();
  h2_BTaggingEff_udsg_T->Write();

  output->Close();

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


