
#include "ConfigParser.h"
#include "ParserUtils.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
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
#include "THStack.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>

int main(int argc, char** argv)
{
  // Input parameters
  
  std::cout << "\n*******************************************************************************************************************" << std::endl;
  std::cout << "arcg: " << argc << std::endl;
  
  //Check if all nedeed arguments to parse are there
  if(argc != 2)
  {
    std::cerr << ">>>>> Analyses::usage: " << argv[0] << " configFileName " << std::endl ;
    return 1;
  }
 
  /// Parse the config file
  parseConfigFile (argv[1]) ;
  
  std::string inputDir = gConfigParser -> readStringOption("Input::inputDir");
  
  std::string outputName = gConfigParser -> readStringOption("Output::outputName");

  system(("ls "+inputDir+" | grep root | awk '{print}' > inputFileAdd.dat").c_str() );
  
  std::map<int, TFile*> Files;
  std::map<int, std::map<int, TH2F*> > Histos;
  std::map<int, TH2F* > FinalHistos;
  
  std::map<int,std::string> FileName;
  std::map<int,std::string> HistoName;
  
  char file[1000];
  FILE *f_file;
  f_file = fopen("inputFileAdd.dat","r");
  
  int file_pos = 0;
  int file_tot = 0;
  int histo_tot = 0;

  
  while(fscanf(f_file,"%s \n", file) !=EOF ){
    std::string FILE = std::string(file);

    if(FILE.find("#") != std::string::npos) continue;
    std::cout << "\nReading File: " << file << std::endl;
    Files[file_pos] = TFile::Open((inputDir+"/"+file).c_str());
    file_pos++;
  }
  
  file_tot = file_pos;

 
  HistoName[0] = std::string("h2_BTaggingEff_Denom_b_L");
  HistoName[1] = std::string("h2_BTaggingEff_Num_b_L");
  HistoName[2] = std::string("h2_BTaggingEff_Denom_b_M");
  HistoName[3] = std::string("h2_BTaggingEff_Num_b_M");
  HistoName[4] = std::string("h2_BTaggingEff_Denom_b_T");
  HistoName[5] = std::string("h2_BTaggingEff_Num_b_T");
  HistoName[6] = std::string("h2_BTaggingEff_Denom_c_L");
  HistoName[7] = std::string("h2_BTaggingEff_Num_c_L");
  HistoName[8] = std::string("h2_BTaggingEff_Denom_c_M");
  HistoName[9] = std::string("h2_BTaggingEff_Num_c_M");
  HistoName[10] = std::string("h2_BTaggingEff_Denom_c_T");
  HistoName[11] = std::string("h2_BTaggingEff_Num_c_T");
  HistoName[12] = std::string("h2_BTaggingEff_Denom_udsg_L");
  HistoName[13] = std::string("h2_BTaggingEff_Num_udsg_L");
  HistoName[14] = std::string("h2_BTaggingEff_Denom_udsg_M");
  HistoName[15] = std::string("h2_BTaggingEff_Num_udsg_M");
  HistoName[16] = std::string("h2_BTaggingEff_Denom_udsg_T");
  HistoName[17] = std::string("h2_BTaggingEff_Num_udsg_T");
  
  histo_tot = 18;

  float ptmin[21] = {20., 30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500., 600., 800.,900.,1000.,1200.,1500.};
  float etamin[8] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5};

  TH2F* h2_BTaggingEff_b_L = new TH2F("h2_BTaggingEff_b_L","h2_BTaggingEff_b_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_b_M = new TH2F("h2_BTaggingEff_b_M","h2_BTaggingEff_b_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_b_T = new TH2F("h2_BTaggingEff_b_T","h2_BTaggingEff_b_T",20, ptmin, 7,etamin);

  TH2F* h2_BTaggingEff_c_L = new TH2F("h2_BTaggingEff_c_L","h2_BTaggingEff_c_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_c_M = new TH2F("h2_BTaggingEff_c_M","h2_BTaggingEff_c_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_c_T = new TH2F("h2_BTaggingEff_c_T","h2_BTaggingEff_c_T",20, ptmin, 7,etamin);

  TH2F* h2_BTaggingEff_udsg_L = new TH2F("h2_BTaggingEff_udsg_L","h2_BTaggingEff_udsg_L",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_udsg_M = new TH2F("h2_BTaggingEff_udsg_M","h2_BTaggingEff_udsg_M",20, ptmin, 7,etamin);
  TH2F* h2_BTaggingEff_udsg_T = new TH2F("h2_BTaggingEff_udsg_T","h2_BTaggingEff_udsg_T",20, ptmin, 7,etamin);


  if (file_tot == 0 )
  {
    std::cout << "Error: give at least one file!" << std::endl; 
    return -1;
  }
  
  for(int ii = 0; ii < file_tot; ii++)
      for(int jj = 0; jj < histo_tot; jj++)
          Histos[ii][jj] = (TH2F*)Files[ii]->Get(HistoName[jj].c_str());
  
  for(int jj = 0; jj < histo_tot; jj++){
      FinalHistos[jj] = (TH2F*)Histos[0][jj]->Clone((HistoName[jj]).c_str());
      //FinalHistos[jj]->Sumw2();
      for(int ii = 1; ii < file_tot; ii++)
          FinalHistos[jj]->Add(Histos[ii][jj]);
  }

  h2_BTaggingEff_b_L->Divide(FinalHistos[1],FinalHistos[0],1,1,"B");
  h2_BTaggingEff_b_M->Divide(FinalHistos[3],FinalHistos[2],1,1,"B");
  h2_BTaggingEff_b_T->Divide(FinalHistos[5],FinalHistos[4],1,1,"B");

  h2_BTaggingEff_c_L->Divide(FinalHistos[7],FinalHistos[6],1,1,"B");
  h2_BTaggingEff_c_M->Divide(FinalHistos[9],FinalHistos[8],1,1,"B");
  h2_BTaggingEff_c_T->Divide(FinalHistos[11],FinalHistos[10],1,1,"B");
  
  h2_BTaggingEff_udsg_L->Divide(FinalHistos[13],FinalHistos[12],1,1,"B");
  h2_BTaggingEff_udsg_M->Divide(FinalHistos[15],FinalHistos[14],1,1,"B");
  h2_BTaggingEff_udsg_T->Divide(FinalHistos[17],FinalHistos[16],1,1,"B");

  TFile* outputFile = new TFile((outputName).c_str(), "RECREATE");
  outputFile->cd();

  for(int ii = 0; ii < histo_tot; ii++)
      FinalHistos[ii]->Write();

      h2_BTaggingEff_b_L->Write();
      h2_BTaggingEff_b_M->Write();
      h2_BTaggingEff_b_T->Write();
      h2_BTaggingEff_c_L->Write();
      h2_BTaggingEff_c_M->Write();
      h2_BTaggingEff_c_T->Write();
      h2_BTaggingEff_udsg_L->Write();
      h2_BTaggingEff_udsg_M->Write();
      h2_BTaggingEff_udsg_T->Write(); 

  outputFile->Close();

  system("rm inputFileAdd.dat");
          
}


















