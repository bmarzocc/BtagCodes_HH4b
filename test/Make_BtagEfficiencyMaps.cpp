
#include "ConfigParser.h"
#include "ParserUtils.h"

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

bool passCutBasedJetId(const TLorentzVector* jetP4,const TLorentzVector* phoP4,const bool& pfloose, const float& betastarclassic, const float& dR2Mean, const int& nvtx);

int main(int argc, char** argv)
{

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

  TLorentzVector p4;

  for(int ii = 0; ii < pos_total; ii++){
      if(ntu[ii]->GetEntries() == 0 )
      {
         std::cout << "Error: input file" << ii << " is empty" << std::endl; 
         return -1;
      }
  }

  int            vtx_std_n;
  int            jet_algoPF1_n;
  int            jet_algoPF1_flavour[200];
  TClonesArray  *jet_algoPF1_p4;
  float          jet_algoPF1_csvBtag[200];  
  float          jet_algoPF1_dR2Mean[200]; 
  float          jet_algoPF1_betaStarClassic[200];
  bool           jet_algoPF1_pfloose[200];
  TClonesArray  *pho_p4;
  int            dipho_n;
  int            dipho_leadind[20];   
  int            dipho_subleadind[20]; 
  int            dipho_vtxind[20]; 
  int            vtx_std_sel; 

  jet_algoPF1_p4 = 0;
  pho_p4 = 0;
  
  for(int ii = 0; ii < pos_total; ii++){
      ntu[ii]->SetBranchStatus("*",0);
      ntu[ii]->SetBranchStatus("vtx_std_n",1); 
      ntu[ii]->SetBranchStatus("jet_algoPF1_n",1);
      ntu[ii]->SetBranchStatus("jet_algoPF1_flavour",1);
      ntu[ii]->SetBranchStatus("jet_algoPF1_p4",1);
      ntu[ii]->SetBranchStatus("jet_algoPF1_csvBtag",1); 
      ntu[ii]->SetBranchStatus("jet_algoPF1_dR2Mean",1);
      ntu[ii]->SetBranchStatus("jet_algoPF1_betaStarClassic",1); 
      ntu[ii]->SetBranchStatus("jet_algoPF1_pfloose",1); 
      ntu[ii]->SetBranchStatus("pho_p4",1); 
      ntu[ii]->SetBranchStatus("dipho_n",1); 
      ntu[ii]->SetBranchStatus("dipho_leadind",1); 
      ntu[ii]->SetBranchStatus("dipho_subleadind",1); 
      ntu[ii]->SetBranchStatus("dipho_vtxind",1); 
      ntu[ii]->SetBranchStatus("vtx_std_sel",1); 
      ntu[ii]->SetBranchAddress("vtx_std_n",&vtx_std_n); 
      ntu[ii]->SetBranchAddress("jet_algoPF1_n",&jet_algoPF1_n);
      ntu[ii]->SetBranchAddress("jet_algoPF1_flavour",&jet_algoPF1_flavour);
      ntu[ii]->SetBranchAddress("jet_algoPF1_p4",&jet_algoPF1_p4);
      ntu[ii]->SetBranchAddress("jet_algoPF1_csvBtag",&jet_algoPF1_csvBtag); 
      ntu[ii]->SetBranchAddress("jet_algoPF1_dR2Mean",&jet_algoPF1_dR2Mean);
      ntu[ii]->SetBranchAddress("jet_algoPF1_betaStarClassic",&jet_algoPF1_betaStarClassic); 
      ntu[ii]->SetBranchAddress("jet_algoPF1_pfloose",&jet_algoPF1_pfloose); 
      ntu[ii]->SetBranchAddress("pho_p4",&pho_p4);
      ntu[ii]->SetBranchAddress("dipho_n",&dipho_n);
      ntu[ii]->SetBranchAddress("dipho_leadind",&dipho_leadind);
      ntu[ii]->SetBranchAddress("dipho_subleadind",&dipho_subleadind);
      ntu[ii]->SetBranchAddress("dipho_vtxind",&dipho_vtxind);
      ntu[ii]->SetBranchAddress("vtx_std_sel",&vtx_std_sel);
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

  std::map<int,TLorentzVector*> jetP4;
  TLorentzVector* phoP4;

  for(int nn = 0; nn < pos_total; nn++){
      for(int ientry = 0; ientry < ntu[nn]->GetEntries(); ientry++){
          if(ientry%1000==0) std::cout<<"--- Reading file_" << nn << " entry = "<< ientry <<std::endl;
          ntu[nn]->GetEntry(ientry);

          for(int jj = 0; jj < dipho_n; jj++){
              if(dipho_vtxind[jj] != vtx_std_sel) continue;
              phoP4 = (TLorentzVector*)pho_p4->ConstructedAt(dipho_leadind[jj]);
          }

          for(int ii = 0; ii < jet_algoPF1_n; ii++){

              jetP4[ii] = (TLorentzVector*)jet_algoPF1_p4->ConstructedAt(ii);

              if(passCutBasedJetId(jetP4[ii],phoP4,jet_algoPF1_pfloose[ii],jet_algoPF1_betaStarClassic[ii],jet_algoPF1_dR2Mean[ii],vtx_std_n) == false) continue;
               
              if(abs(jet_algoPF1_flavour[ii]) == 5){
                 
                 h2_BTaggingEff_Denom_b_L->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.244) h2_BTaggingEff_Num_b_L->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                   
                 h2_BTaggingEff_Denom_b_M->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.679) h2_BTaggingEff_Num_b_M->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 
                 h2_BTaggingEff_Denom_b_T->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.898) h2_BTaggingEff_Num_b_T->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
              }
 
              if(abs(jet_algoPF1_flavour[ii]) == 4){
                 
                 h2_BTaggingEff_Denom_c_L->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.244) h2_BTaggingEff_Num_c_L->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                   
                 h2_BTaggingEff_Denom_c_M->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.679) h2_BTaggingEff_Num_c_M->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 
                 h2_BTaggingEff_Denom_c_T->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.898) h2_BTaggingEff_Num_c_T->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
              }

              if(abs(jet_algoPF1_flavour[ii]) == 1 || abs(jet_algoPF1_flavour[ii]) == 2 || abs(jet_algoPF1_flavour[ii]) == 3 || abs(jet_algoPF1_flavour[ii]) == 21){
                 
                 h2_BTaggingEff_Denom_udsg_L->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.244) h2_BTaggingEff_Num_udsg_L->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                   
                 h2_BTaggingEff_Denom_udsg_M->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.679) h2_BTaggingEff_Num_udsg_M->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 
                 h2_BTaggingEff_Denom_udsg_T->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
                 if(jet_algoPF1_csvBtag[ii] > 0.898) h2_BTaggingEff_Num_udsg_T->Fill(jetP4[ii]->Pt(),fabs(jetP4[ii]->Eta()));
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

  TFile* output = new TFile(outputName.c_str(),"RECREATE");
  
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

bool passCutBasedJetId(const TLorentzVector* jetP4,const TLorentzVector* phoP4,const bool& pfloose, const float& betastarclassic, const float& dR2Mean, const int& nvtx) {

  bool isGood = true;
  float eta = jetP4->Eta();
  
  //if(jetP4->DeltaR(*phoP4) < 0.5) isGood = false;

  if(pfloose == false) isGood = false;

  if ( fabs(eta) < 2.5 ) {
    if ( betastarclassic > 0.2 * log( nvtx - 0.64) )  isGood = false;
    if (dR2Mean > 0.06)                               isGood = false;
  } 
  else if (fabs(eta) < 3.){
    if ( dR2Mean > 0.05)  isGood =false;
  } 
  else {
    if ( dR2Mean > 0.055) isGood =false;
  }
  
  return isGood;
}
















