#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<bitset>

#include <TSystem.h>
#include <string>
#include <sstream> //string stream
#include <sys/stat.h> //get the status of files. "st_"
#include <sys/types.h>
#include <unistd.h> //UNIx STanDard Header file
#include <climits> //char limits

#include"TROOT.h"
#include"TDirectory.h"
#include"TFile.h"
#include"TTree.h"
#include"TCut.h"
#include"TCutG.h"
#include"TMath.h"
#include"TString.h"
#include"TEnv.h"
#include"TBits.h"

using namespace std;
using namespace TMath;

int main(int argc, char *argv[]){

  time_t start, stop;
  time(&start);

  Int_t FileNumber = TString(argv[1]).Atoi();
  printf("\n%s %d %s\n\n","=== Execute ana_beam for RUN",FileNumber," ===" );
  //===== Load input file =================================================

  TString FileName = Form("/home/koiwai/analysis/rootfiles/unpacked/run%04d.root",FileNumber);
  TFile *infile = TFile::Open(FileName);
  
  TTree *caltr;
  infile->GetObject("caltr",caltr);

  printf("%-20s %s \n","Input data file:",FileName.Data());


  //=== Load Env file ===
  TEnv *env_set = new TEnv("/home/koiwai/analysis/conversion_settings.prm");
  printf("%-20s %s \n","Setting file:","/home/koiwai/analysis/conversion_settings.prm");

  
  //===== in tree variables ===============================================
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;

  //=== Plastic ===
  Double_t plaF3_Q, plaF5_Q, plaF7_Q;
  Double_t plaF3_T, plaF5_T, plaF7_T;
  Double_t plaF3_dT, plaF5_dT, plaF7_dT;

  Double_t plaF3_QL, plaF3_QR, plaF5_QL, plaF5_QR, plaF7_QL, plaF7_QR; //newly added.
  Double_t plaF3_TL, plaF3_TR, plaF5_TL, plaF5_TR, plaF7_TL, plaF7_TR; //newly added.

  //=== PPAC ===
  Double_t ppacF31A_X, ppacF31A_Y, ppacF31B_X, ppacF31B_Y, ppacF32A_X, ppacF32A_Y, ppacF32B_X, ppacF32B_Y;
  Double_t ppacF51A_X, ppacF51A_Y, ppacF51B_X, ppacF51B_Y, ppacF52A_X, ppacF52A_Y, ppacF52B_X, ppacF52B_Y;
  Double_t ppacF71A_X, ppacF71A_Y, ppacF71B_X, ppacF71B_Y, ppacF72A_X, ppacF72A_Y, ppacF72B_X, ppacF72B_Y;
  
  Int_t ppacF31A_X_T1, ppacF31A_X_T2, ppacF31A_Y_T1, ppacF31A_Y_T2;
  Int_t ppacF31B_X_T1, ppacF31B_X_T2, ppacF31B_Y_T1, ppacF31B_Y_T2;
  Int_t ppacF32A_X_T1, ppacF32A_X_T2, ppacF32A_Y_T1, ppacF32A_Y_T2;
  Int_t ppacF32B_X_T1, ppacF32B_X_T2, ppacF32B_Y_T1, ppacF32B_Y_T2;
  Int_t ppacF51A_X_T1, ppacF51A_X_T2, ppacF51A_Y_T1, ppacF51A_Y_T2;
  Int_t ppacF51B_X_T1, ppacF51B_X_T2, ppacF51B_Y_T1, ppacF51B_Y_T2;
  Int_t ppacF52A_X_T1, ppacF52A_X_T2, ppacF52A_Y_T1, ppacF52A_Y_T2;
  Int_t ppacF52B_X_T1, ppacF52B_X_T2, ppacF52B_Y_T1, ppacF52B_Y_T2;
  Int_t ppacF71A_X_T1, ppacF71A_X_T2, ppacF71A_Y_T1, ppacF71A_Y_T2;
  Int_t ppacF71B_X_T1, ppacF71B_X_T2, ppacF71B_Y_T1, ppacF71B_Y_T2;
  Int_t ppacF72A_X_T1, ppacF72A_X_T2, ppacF72A_Y_T1, ppacF72A_Y_T2;
  Int_t ppacF72B_X_T1, ppacF72B_X_T2, ppacF72B_Y_T1, ppacF72B_Y_T2;

  //=== F7IC === 
  Double_t icF7_E;
  
  Int_t icF7_raw[6];

  //=== SBT ===
  Double_t sbt1_Q, sbt2_Q, sbt1_T, sbt2_T, sbt1_dT, sbt2_dT;
  Double_t sbt1_QL, sbt1_QR, sbt1_TL, sbt1_TR;
  Double_t sbt2_QL, sbt2_QR, sbt2_TL, sbt2_TR;

  Double_t sbt1_Tslew,  sbt2_Tslew;
  Double_t sbt1_dTslew, sbt2_dTslew;
  
  //=== Set value ===
  caltr->SetBranchAddress("RunNumber",&RunNumber);
  caltr->SetBranchAddress("EventNumber",&EventNumber);
  
  caltr->SetBranchAddress("plaF3_Q",&plaF3_Q);
  caltr->SetBranchAddress("plaF5_Q",&plaF5_Q);
  caltr->SetBranchAddress("plaF7_Q",&plaF7_Q);
  caltr->SetBranchAddress("plaF3_T",&plaF3_T);
  caltr->SetBranchAddress("plaF5_T",&plaF5_T);
  caltr->SetBranchAddress("plaF7_T",&plaF7_T);
  caltr->SetBranchAddress("plaF3_dT",&plaF3_dT);
  caltr->SetBranchAddress("plaF5_dT",&plaF5_dT);
  caltr->SetBranchAddress("plaF7_dT",&plaF7_dT);
  
  caltr->SetBranchAddress("plaF3_QL",&plaF3_QL);
  caltr->SetBranchAddress("plaF3_QR",&plaF3_QR);
  caltr->SetBranchAddress("plaF5_QL",&plaF5_QL);
  caltr->SetBranchAddress("plaF5_QR",&plaF5_QR);
  caltr->SetBranchAddress("plaF7_QL",&plaF7_QL);
  caltr->SetBranchAddress("plaF7_QR",&plaF7_QR);
  
  caltr->SetBranchAddress("plaF3_TL",&plaF3_TL);
  caltr->SetBranchAddress("plaF3_TR",&plaF3_TR);
  caltr->SetBranchAddress("plaF5_TL",&plaF5_TL);
  caltr->SetBranchAddress("plaF5_TR",&plaF5_TR);
  caltr->SetBranchAddress("plaF7_TL",&plaF7_TL);
  caltr->SetBranchAddress("plaF7_TR",&plaF7_TR);
  
  caltr->SetBranchAddress("ppacF31A_X",&ppacF31A_X);
  caltr->SetBranchAddress("ppacF31A_Y",&ppacF31A_Y);
  caltr->SetBranchAddress("ppacF31B_X",&ppacF31B_X);
  caltr->SetBranchAddress("ppacF31B_Y",&ppacF31B_Y);
  caltr->SetBranchAddress("ppacF32A_X",&ppacF32A_X);
  caltr->SetBranchAddress("ppacF32A_Y",&ppacF32A_Y);
  caltr->SetBranchAddress("ppacF32B_X",&ppacF32B_X);
  caltr->SetBranchAddress("ppacF32B_Y",&ppacF32B_Y);
  
  caltr->SetBranchAddress("ppacF51A_X",&ppacF51A_X);
  caltr->SetBranchAddress("ppacF51A_Y",&ppacF51A_Y);
  caltr->SetBranchAddress("ppacF51B_X",&ppacF51B_X);
  caltr->SetBranchAddress("ppacF51B_Y",&ppacF51B_Y);
  caltr->SetBranchAddress("ppacF52A_X",&ppacF52A_X);
  caltr->SetBranchAddress("ppacF52A_Y",&ppacF52A_Y);
  caltr->SetBranchAddress("ppacF52B_X",&ppacF52B_X);
  caltr->SetBranchAddress("ppacF52B_Y",&ppacF52B_Y);
  
  caltr->SetBranchAddress("ppacF71A_X",&ppacF71A_X);
  caltr->SetBranchAddress("ppacF71A_Y",&ppacF71A_Y);
  caltr->SetBranchAddress("ppacF71B_X",&ppacF71B_X);
  caltr->SetBranchAddress("ppacF71B_Y",&ppacF71B_Y);
  caltr->SetBranchAddress("ppacF72A_X",&ppacF72A_X);
  caltr->SetBranchAddress("ppacF72A_Y",&ppacF72A_Y);
  caltr->SetBranchAddress("ppacF72B_X",&ppacF72B_X);
  caltr->SetBranchAddress("ppacF72B_Y",&ppacF72B_Y);
  
  caltr->SetBranchAddress("ppacF31A_X_T1",&ppacF31A_X_T1);
  caltr->SetBranchAddress("ppacF31A_X_T2",&ppacF31A_X_T2);
  caltr->SetBranchAddress("ppacF31A_Y_T1",&ppacF31A_Y_T1);
  caltr->SetBranchAddress("ppacF31A_Y_T2",&ppacF31A_Y_T2);
  caltr->SetBranchAddress("ppacF31B_X_T1",&ppacF31B_X_T1);
  caltr->SetBranchAddress("ppacF31B_X_T2",&ppacF31B_X_T2);
  caltr->SetBranchAddress("ppacF31B_Y_T1",&ppacF31B_Y_T1);
  caltr->SetBranchAddress("ppacF31B_Y_T2",&ppacF31B_Y_T2);
  caltr->SetBranchAddress("ppacF32A_X_T1",&ppacF32A_X_T1);
  caltr->SetBranchAddress("ppacF32A_X_T2",&ppacF32A_X_T2);
  caltr->SetBranchAddress("ppacF32A_Y_T1",&ppacF32A_Y_T1);
  caltr->SetBranchAddress("ppacF32A_Y_T2",&ppacF32A_Y_T2);
  caltr->SetBranchAddress("ppacF32B_X_T1",&ppacF32B_X_T1);
  caltr->SetBranchAddress("ppacF32B_X_T2",&ppacF32B_X_T2);
  caltr->SetBranchAddress("ppacF32B_Y_T1",&ppacF32B_Y_T1);
  caltr->SetBranchAddress("ppacF32B_Y_T2",&ppacF32B_Y_T2);
  caltr->SetBranchAddress("ppacF51A_X_T1",&ppacF51A_X_T1);
  caltr->SetBranchAddress("ppacF51A_X_T2",&ppacF51A_X_T2);
  caltr->SetBranchAddress("ppacF51A_Y_T1",&ppacF51A_Y_T1);
  caltr->SetBranchAddress("ppacF51A_Y_T2",&ppacF51A_Y_T2);
  caltr->SetBranchAddress("ppacF51B_X_T1",&ppacF51B_X_T1);
  caltr->SetBranchAddress("ppacF51B_X_T2",&ppacF51B_X_T2);
  caltr->SetBranchAddress("ppacF51B_Y_T1",&ppacF51B_Y_T1);
  caltr->SetBranchAddress("ppacF51B_Y_T2",&ppacF51B_Y_T2);
  caltr->SetBranchAddress("ppacF52A_X_T1",&ppacF52A_X_T1);
  caltr->SetBranchAddress("ppacF52A_X_T2",&ppacF52A_X_T2);
  caltr->SetBranchAddress("ppacF52A_Y_T1",&ppacF52A_Y_T1);
  caltr->SetBranchAddress("ppacF52A_Y_T2",&ppacF52A_Y_T2);
  caltr->SetBranchAddress("ppacF52B_X_T1",&ppacF52B_X_T1);
  caltr->SetBranchAddress("ppacF52B_X_T2",&ppacF52B_X_T2);
  caltr->SetBranchAddress("ppacF52B_Y_T1",&ppacF52B_Y_T1);
  caltr->SetBranchAddress("ppacF52B_Y_T2",&ppacF52B_Y_T2);
  caltr->SetBranchAddress("ppacF71A_X_T1",&ppacF71A_X_T1);
  caltr->SetBranchAddress("ppacF71A_X_T2",&ppacF71A_X_T2);
  caltr->SetBranchAddress("ppacF71A_Y_T1",&ppacF71A_Y_T1);
  caltr->SetBranchAddress("ppacF71A_Y_T2",&ppacF71A_Y_T2);
  caltr->SetBranchAddress("ppacF71B_X_T1",&ppacF71B_X_T1);
  caltr->SetBranchAddress("ppacF71B_X_T2",&ppacF71B_X_T2);
  caltr->SetBranchAddress("ppacF71B_Y_T1",&ppacF71B_Y_T1);
  caltr->SetBranchAddress("ppacF71B_Y_T2",&ppacF71B_Y_T2);
  caltr->SetBranchAddress("ppacF72A_X_T1",&ppacF72A_X_T1);
  caltr->SetBranchAddress("ppacF72A_X_T2",&ppacF72A_X_T2);
  caltr->SetBranchAddress("ppacF72A_Y_T1",&ppacF72A_Y_T1);
  caltr->SetBranchAddress("ppacF72A_Y_T2",&ppacF72A_Y_T2);
  caltr->SetBranchAddress("ppacF72B_X_T1",&ppacF72B_X_T1);
  caltr->SetBranchAddress("ppacF72B_X_T2",&ppacF72B_X_T2);
  caltr->SetBranchAddress("ppacF72B_Y_T1",&ppacF72B_Y_T1);
  caltr->SetBranchAddress("ppacF72B_Y_T2",&ppacF72B_Y_T2);
  
  caltr->SetBranchAddress("icF7_E",&icF7_E);
  caltr->SetBranchAddress("icF7_raw",icF7_raw);

  caltr->SetBranchAddress("sbt1_Q",&sbt1_Q);
  caltr->SetBranchAddress("sbt1_T",&sbt1_T);
  caltr->SetBranchAddress("sbt1_dT",&sbt1_dT);
  caltr->SetBranchAddress("sbt2_Q",&sbt2_Q);
  caltr->SetBranchAddress("sbt2_T",&sbt2_T);
  caltr->SetBranchAddress("sbt2_dT",&sbt2_dT);
  caltr->SetBranchAddress("sbt1_QL",&sbt1_QL);
  caltr->SetBranchAddress("sbt1_QR",&sbt1_QR);
  caltr->SetBranchAddress("sbt1_TL",&sbt1_TL);
  caltr->SetBranchAddress("sbt1_TR",&sbt1_TR);
  caltr->SetBranchAddress("sbt2_QL",&sbt2_QL);
  caltr->SetBranchAddress("sbt2_QR",&sbt2_QR);
  caltr->SetBranchAddress("sbt2_TL",&sbt2_TL);
  caltr->SetBranchAddress("sbt2_TR",&sbt2_TR);

  caltr->SetBranchAddress("sbt1_Tslew",&sbt1_Tslew);
  caltr->SetBranchAddress("sbt2_Tslew",&sbt2_Tslew);    
    
  //===== Load CUT files ==================================================
  //=== Plastic (graphical cut)===
  TString placutfilename = env_set->GetValue("placut","");
  TFile *fcutpla = TFile::Open(placutfilename);
  printf("%-20s %-25s \n","Cut file for Pla:",placutfilename.Data());
  TCutG *cplaF3 = (TCutG*)fcutpla->Get("CPLAF3");
  TCutG *cplaF5 = (TCutG*)fcutpla->Get("CPLAF5");
  TCutG *cplaF7 = (TCutG*)fcutpla->Get("CPLAF7");

  /*
  //=== Chrage change @ F5 (graphical cut) ===
  TFile *cutfileF5Qchange = new TFile("/home/koiwai/analysis/cutfiles/cut_bigripsZvsF5Qchange.root");
  TCutG *cF5Qchange = (TCutG*)cutfileF5Qchange->Get("CUTG");
  */
  
  //=== PPAC Tsum gate ===
  TString ppaccutfilename = env_set->GetValue("ppaccut","");
  ifstream finppac;
  finppac.open(ppaccutfilename);
  printf("%-20s %-25s \n","Cut file for PPAC:",ppaccutfilename.Data());
  if(finppac.fail()){
    //cout << "Error: file is not found." << endl;
    printf("Error: file %s is not found.\n",ppaccutfilename.Data());
    return 1;
  }
  string dummyppac[24];
  Double_t cppac_low[24], cppac_up[24];
  
  for(Int_t cPPAC_index = 0;cPPAC_index<24;++cPPAC_index)
    finppac >> dummyppac[cPPAC_index] >> cppac_low[cPPAC_index] >> cppac_up[cPPAC_index];

  //=== IC gate ===
  TString iccutfilename = env_set->GetValue("iccut","");
  ifstream finic;
  finic.open(iccutfilename);
  printf("%-20s %-25s \n","Cut file for IC:",iccutfilename.Data());
  if(finic.fail()){
    //cout << "Error: file is not found." << endl;
    printf("Error: file %s is not found.\n",iccutfilename.Data());
    return 1;
  }

  Double_t cIC_low[5], cIC_up[5];
  
  for(Int_t cIC_index = 0;cIC_index<5;++cIC_index)
    finic >> cIC_low[cIC_index] >> cIC_up[cIC_index];

  
  //===== Load .dat files =====
  TEnv *env = new TEnv(env_set->GetValue("geometrydata",""));
  printf("%-20s %s \n","Geometry file:",env_set->GetValue("geometrydata",""));
  TEnv *env_pla  = new TEnv("/home/koiwai/analysis/db/pla2pos_Z.dat");
  TEnv *env_q = new TEnv("/home/koiwai/analysis/db/plaQ2f7icE.dat");
  
  //===== Create output file/tree =========================================
  TString AnaFileName = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNumber);
  TFile *anafile = new TFile(AnaFileName,"recreate");

  TTree *anatrB = new TTree("anatrB","anatrB");

  printf("\n%-20s %s \n\n","Output file:",AnaFileName.Data());

  //===== Declear const.s =================================================
  //double pla3pos[2];
  //double pla7pos[2];
  //double pla7z[3];
  //double pla3q2e[2];
  //double pla7q2e[2];
  //double pla13_1q2e[2];
  //double pla13_2q2e[2];
  //for(Int_t i=0;i<2;++i){
    //pla3pos[i] = env_pla->GetValue(Form("pla3pos%d",i),0.0);
    //pla7pos[i] = env_pla->GetValue(Form("pla7pos%d",i),0.0);
    //pla7z[i] = env_pla->GetValue(Form("pla7Z[%d]",i),0.0);
    //pla3q2e[i]    = env_q->GetValue(Form("pla3q2e%d",i),0.0);
    //pla7q2e[i]    = env_q->GetValue(Form("pla7q2e%d",i),0.0);
    //pla13_1q2e[i] = env_q->GetValue(Form("pla13_1q2e%d",i),0.0);
    //pla13_2q2e[i] = env_q->GetValue(Form("pla13_2q2e%d",i),0.0);
  //}

  double DistF3F7 = 46568.; //[mm]
  double DistF3F5 = 23284.;
  double DistF5F7 = 23284.;
  double DistF7F13 = env->GetValue("Dist_F7SBT",0.0);
  double DistF3F13 = env->GetValue("Dist_F3F13",0.0);
  double OffsetF3F7 = env->GetValue("offsetF3F7",292.379); //[nsec]
  double OffsetF3F5 = env->GetValue("offsetF3F5",159.572);
  double OffsetF5F7 = env->GetValue("offsetF5F7",132.807);
  double OffsetF7F13 = env->GetValue("offsetF7SBT",588.109);
  double OffsetF3F13 = env->GetValue("offsetF3F13",0.0);
  double Ionpair = 4.866; //[keV]
  double m_e = 511.; //[keV]
  double m_u = 931.49432; //[MeV]
  double clight = 299.792458; //[mm/nsec]
  double zetBR_c0 = env->GetValue("zetBR_c0",0.0);
  double zetBR_c1 = env->GetValue("zetBR_c1",0.0);
  double DistF3PPAC = 890.; //[mm]
  double DistF5PPAC = 650.;
  double DistF7PPAC = 945.;
  //=== Optical matrix ===
  //from matrix/mat1.mat 
  double XXF3F5 = 0.917467; //[]
  double XAF3F5 = -0.00520039; //[mm/mrad]
  double XDF3F5 = 31.6051; //[mm/%]
  double AXF3F5 = -0.018721; //[mrad/mm]
  double AAF3F5 = 1.09006; //[]
  double ADF3F5 = -0.0129966; //[mrad/%]
  double YYF3F5 = 1.07852;
  double YBF3F5 = 0.0325887;
  double BYF3F5 = 0.331125;
  double BBF3F5 = 0.937205;
  //from matrix/mat2.mat
  double XXF5F7 = 1.09101;  //[]
  double XAF5F7 = 0.020415; //[mm/mrad]
  //double XDF5F7 = -34.4457; //[mm/%]
  double XDF5F7 = -33.9457; //[mm/%]
  double AXF5F7 = -0.0172247; //[mrad/mm]
  double AAF5F7 = 0.916262; //[]
  double ADF5F7 = 0.590371; //[mrad/%]
  double YYF5F7 = 0.939059;
  double YBF5F7 = 0.0255208;
  double BYF5F7 = 0.33606;
  double BBF5F7 = 1.07403;
  //=== Central Brho from run sheet ===
  double Brho0F3F5 = 7.0621; //Tm
  double Brho0F5F7 = 6.854; //Tm
  
  //===== Define ana variables ============================================
  Int_t EventNum, RunNum;
  
  //=== for Z ===
  Double_t tofF3F7, tofF3F5, tofF5F7, tofF7F13, tofF3F13;
  Double_t vF3F7, vF3F5, vF5F7, vF7F13, vF3F13;
  Double_t betaF3F7, betaF3F5, betaF5F7, betaF7F13, betaF3F13;
  Double_t gammaF3F7, gammaF3F5, gammaF5F7, gammaF7F13, gammaF3F13;
  Double_t zetBRraw;
  Double_t zetBR, zetBR313, zetBR37;

  Int_t tsum_f31ax, tsum_f31bx, tsum_f32ax, tsum_f32bx;
  Int_t tsum_f31ay, tsum_f31by, tsum_f32ay, tsum_f32by;
  Int_t tsum_f51ax, tsum_f51bx, tsum_f52ax, tsum_f52bx;
  Int_t tsum_f51ay, tsum_f51by, tsum_f52ay, tsum_f52by;
  Int_t tsum_f71ax, tsum_f71bx, tsum_f72ax, tsum_f72bx;
  Int_t tsum_f71ay, tsum_f71by, tsum_f72ay, tsum_f72by;
  
  

  //=== for A/Q ===
  Double_t F3X, F3Y, F3A, F3B; //X, Y:[mm]
  Double_t F5X, F5Y, F5A, F5B; //A, B:[mrad]
  Double_t F7X, F7Y, F7A, F7B;

  Double_t F31_X, F31_Y, F32_X, F32_Y;
  Double_t F51_X, F51_Y, F52_X, F52_Y;
  Double_t F71_X, F71_Y, F72_X, F72_Y;

  Double_t deltaF3F5, deltaF5F7;
  Double_t deltaF3F5_X, deltaF3F5_A;
  Double_t brhoF3F5, brhoF5F7;
  Double_t aoqF3F13, aoqF5F7, aoqF3F5;
  Double_t recoF3A, recodeltaF3F5_X, recodeltaF3F5_A;
  Double_t aoqBR;

  Double_t reco1f7x, reco2f7x;
  Double_t reco1deltaF5F7, reco2deltaF5F7, reco1brhoF5F7, reco2brhoF5F7;
  Double_t reco1aoq57, reco2aoq57;
  
  Bool_t BG_flag; //flag for background
  Bool_t plaflag[3], ppacflag[3];
  Bool_t icflag;
  
  Int_t f71flag, f72flag;



  
  //======
  anatrB->Branch("EventNumber",&EventNum);
  anatrB->Branch("RunNumber",&RunNum);
  
  anatrB->Branch("tofF3F7",&tofF3F7);
  anatrB->Branch("tofF3F5",&tofF3F5);
  anatrB->Branch("tofF5F7",&tofF5F7);
  anatrB->Branch("tofF7F13",&tofF7F13);
  anatrB->Branch("tofF3F13",&tofF3F13);
  anatrB->Branch("vF3F7",&vF3F7);
  anatrB->Branch("vF3F5",&vF3F5);
  anatrB->Branch("vF5F7",&vF5F7);
  anatrB->Branch("vF7F13",&vF7F13);
  anatrB->Branch("vF3F13",&vF3F13);
  anatrB->Branch("betaF3F7",&betaF3F7);
  anatrB->Branch("betaF3F5",&betaF3F5);
  anatrB->Branch("betaF5F7",&betaF5F7);
  anatrB->Branch("betaF7F13",&betaF7F13);
  anatrB->Branch("betaF3F13",&betaF3F13);
  anatrB->Branch("gammaF3F7",&gammaF3F7);
  anatrB->Branch("gammaF3F5",&gammaF3F5);
  anatrB->Branch("gammaF5F7",&gammaF5F7);
  anatrB->Branch("gammaF7F13",&gammaF7F13);
  anatrB->Branch("gammaF3F13",&gammaF3F13);
  anatrB->Branch("zetBRraw",&zetBRraw);
  anatrB->Branch("zetBR",&zetBR);
  anatrB->Branch("zetBR313",&zetBR313);
  anatrB->Branch("zetBR37",&zetBR37);
 
  anatrB->Branch("F3X",&F3X);
  anatrB->Branch("F3Y",&F3Y);
  anatrB->Branch("F3A",&F3A);
  anatrB->Branch("F3B",&F3B);
  anatrB->Branch("F5X",&F5X);
  anatrB->Branch("F5Y",&F5Y);
  anatrB->Branch("F5A",&F5A);
  anatrB->Branch("F5B",&F5B);
  anatrB->Branch("F7X",&F7X);
  anatrB->Branch("F7Y",&F7Y);
  anatrB->Branch("F7A",&F7A);
  anatrB->Branch("F7B",&F7B);
  anatrB->Branch("deltaF3F5",&deltaF3F5);
  anatrB->Branch("deltaF5F7",&deltaF5F7);
  anatrB->Branch("brhoF3F5",&brhoF3F5);
  anatrB->Branch("brhoF5F7",&brhoF5F7);
  anatrB->Branch("aoqF3F13",&aoqF3F13);
  anatrB->Branch("aoqF5F7",&aoqF5F7);
  anatrB->Branch("aoqF3F5",&aoqF3F5);
  anatrB->Branch("recoF3A",&recoF3A);
  anatrB->Branch("recodeltaF3F5_A",&recodeltaF3F5_A);
  anatrB->Branch("recodeltaF3F5_X",&recodeltaF3F5_X);
  anatrB->Branch("deltaF3F5_A",&deltaF3F5_A);
  anatrB->Branch("deltaF3F5_X",&deltaF3F5_X);
  anatrB->Branch("aoqBR",&aoqBR);

  anatrB->Branch("reco1f7x",&reco1f7x);
  anatrB->Branch("reco2f7x",&reco2f7x);

  anatrB->Branch("reco1aoq57",&reco1aoq57);
  anatrB->Branch("reco2aoq57",&reco2aoq57);
  
  anatrB->Branch("BG_flag",&BG_flag,"BG_flag/O");
  anatrB->Branch("plaflag",plaflag,"plaflag[3]/O");
  anatrB->Branch("ppacflag",ppacflag,"ppacflag[3]/O");
  anatrB->Branch("icflag",&icflag,"icflag/O");
  
  anatrB->Branch("f71flag",&f71flag);
  anatrB->Branch("f72flag",&f72flag);

  infile->cd();

 //===== Begin LOOP ======================================================

  //cout << "Start conversion." << endl;
  printf("Conversion START!\n");
  
  int nEntry = caltr->GetEntries();
  //for(int iEntry=0;iEntry<nEntry;++iEntry){
    for(int iEntry=0;iEntry<10;++iEntry){
    caltr->GetEntry(iEntry);

    if(iEntry%100 == 0){
      std::clog << iEntry/1000 << "k / " << nEntry/1000 << "k events treated..." << "\r";
    }

    EventNum = EventNumber;
    RunNum = RunNumber;
    
    //=== Initialization ===
    BG_flag = kTRUE;
    for(int flag_index=0;flag_index<3;flag_index++){
      plaflag[flag_index] = kTRUE;
      ppacflag[flag_index] = kTRUE;
    }
    icflag = kTRUE;
    
    f71flag = 0;
    f72flag = 0;
    tofF3F7 = TMath::Sqrt(-1);
    tofF3F5 = TMath::Sqrt(-1);
    tofF5F7 = TMath::Sqrt(-1);
    tofF7F13 = TMath::Sqrt(-1);
    tofF3F13 = TMath::Sqrt(-1);
    vF3F7 = TMath::Sqrt(-1);
    vF3F5 = TMath::Sqrt(-1);
    vF5F7 = TMath::Sqrt(-1);
    vF7F13 = TMath::Sqrt(-1);
    vF3F13 = TMath::Sqrt(-1);
    betaF3F7 = TMath::Sqrt(-1);
    betaF3F5 = TMath::Sqrt(-1);
    betaF5F7 = TMath::Sqrt(-1);
    betaF7F13 = TMath::Sqrt(-1);
    betaF3F13 = TMath::Sqrt(-1);
    gammaF3F7 = TMath::Sqrt(-1);
    gammaF3F5 = TMath::Sqrt(-1);
    gammaF5F7 = TMath::Sqrt(-1);
    gammaF7F13 = TMath::Sqrt(-1);
    gammaF3F13 = TMath::Sqrt(-1);

    
    zetBRraw= TMath::Sqrt(-1);
    zetBR = TMath::Sqrt(-1);
    zetBR313 = TMath::Sqrt(-1);

    tsum_f31ax = TMath::Sqrt(-1);
    tsum_f31bx = TMath::Sqrt(-1);
    tsum_f32ax = TMath::Sqrt(-1);
    tsum_f32bx = TMath::Sqrt(-1);
    tsum_f31ay = TMath::Sqrt(-1);
    tsum_f31by = TMath::Sqrt(-1);
    tsum_f32ay = TMath::Sqrt(-1);
    tsum_f32by = TMath::Sqrt(-1);
    tsum_f51ax = TMath::Sqrt(-1);
    tsum_f51bx = TMath::Sqrt(-1);
    tsum_f52ax = TMath::Sqrt(-1);
    tsum_f52bx = TMath::Sqrt(-1);
    tsum_f51ay = TMath::Sqrt(-1);
    tsum_f51by = TMath::Sqrt(-1);
    tsum_f52ay = TMath::Sqrt(-1);
    tsum_f52by = TMath::Sqrt(-1);
    tsum_f71ax = TMath::Sqrt(-1);
    tsum_f71bx = TMath::Sqrt(-1);
    tsum_f72ax = TMath::Sqrt(-1);
    tsum_f72bx = TMath::Sqrt(-1);
    tsum_f71ay = TMath::Sqrt(-1);
    tsum_f71by = TMath::Sqrt(-1);
    tsum_f72ay = TMath::Sqrt(-1);
    tsum_f72by = TMath::Sqrt(-1);
    
    F3X = TMath::Sqrt(-1);
    F3Y = TMath::Sqrt(-1);
    F3A = TMath::Sqrt(-1);
    F3B = TMath::Sqrt(-1);
    F5X = TMath::Sqrt(-1);
    F5Y = TMath::Sqrt(-1);
    F5A = TMath::Sqrt(-1);
    F5B = TMath::Sqrt(-1);
    F7X = TMath::Sqrt(-1);
    F7Y = TMath::Sqrt(-1);
    F7A = TMath::Sqrt(-1);
    F7B = TMath::Sqrt(-1);

    F31_X = TMath::Sqrt(-1);
    F31_Y = TMath::Sqrt(-1);
    F32_X = TMath::Sqrt(-1);
    F32_Y = TMath::Sqrt(-1);
    F51_X = TMath::Sqrt(-1);
    F51_Y = TMath::Sqrt(-1);
    F52_X = TMath::Sqrt(-1);
    F52_Y = TMath::Sqrt(-1);
    F71_X = TMath::Sqrt(-1);
    F71_Y = TMath::Sqrt(-1);
    F72_X = TMath::Sqrt(-1);
    F72_Y = TMath::Sqrt(-1);

    deltaF3F5 = TMath::Sqrt(-1);
    deltaF5F7 = TMath::Sqrt(-1);
    brhoF3F5 = TMath::Sqrt(-1);
    brhoF5F7 = TMath::Sqrt(-1);
    aoqF3F13 = TMath::Sqrt(-1);
    aoqF5F7 = TMath::Sqrt(-1);
    aoqF3F5 = TMath::Sqrt(-1);
    recoF3A = TMath::Sqrt(-1);
    recodeltaF3F5_A = TMath::Sqrt(-1);
    recodeltaF3F5_X = TMath::Sqrt(-1);
    deltaF3F5_A = TMath::Sqrt(-1);
    deltaF3F5_X = TMath::Sqrt(-1);
    aoqBR = TMath::Sqrt(-1);

    reco1f7x = TMath::Sqrt(-1);
    reco2f7x = TMath::Sqrt(-1);

    reco1aoq57 = TMath::Sqrt(-1);
    reco2aoq57 = TMath::Sqrt(-1);

    reco1deltaF5F7 = TMath::Sqrt(-1);
    reco2deltaF5F7 = TMath::Sqrt(-1);
    reco1brhoF5F7 = TMath::Sqrt(-1);
    reco2brhoF5F7 = TMath::Sqrt(-1);

    

    //=== Calculation ===
    tofF3F7    = plaF7_T - plaF3_T + OffsetF3F7;    
    tofF3F5    = plaF5_T - plaF3_T + OffsetF3F5;    
    tofF5F7    = plaF7_T - plaF5_T + OffsetF5F7;    
    tofF7F13   = sbt1_T  - plaF7_T + OffsetF7F13;
    tofF3F13   = sbt1_T  - plaF3_T + OffsetF3F13;
    vF3F7      = DistF3F7/tofF3F7;
    vF3F5      = DistF3F5/tofF3F5;
    vF5F7      = DistF5F7/tofF5F7;
    vF7F13     = DistF7F13/tofF7F13;
    vF3F13     = DistF3F13/tofF3F13;
    betaF3F7   = vF3F7/clight;
    betaF3F5   = vF3F5/clight;
    betaF5F7   = vF5F7/clight;
    betaF7F13  = vF7F13/clight;
    betaF3F13  = vF3F13/clight;
    gammaF3F7  = 1/TMath::Sqrt(1.-betaF3F7*betaF3F7);
    gammaF3F5  = 1/TMath::Sqrt(1.-betaF3F5*betaF3F5);
    gammaF5F7  = 1/TMath::Sqrt(1.-betaF5F7*betaF5F7);
    gammaF7F13 = 1/TMath::Sqrt(1.-betaF7F13*betaF7F13);
    gammaF3F13 = 1/TMath::Sqrt(1.-betaF3F13*betaF3F13);
    
    double dev37  = TMath::Log(2*m_e*vF3F7*vF3F7/Ionpair)-TMath::Log(1-betaF3F7*betaF3F7)-betaF3F7*betaF3F7;
    double dev313 = TMath::Log(2*m_e*vF3F13*vF3F13/Ionpair)-TMath::Log(1-betaF3F13*betaF3F13)-betaF3F13*betaF3F13;

    zetBR37  = vF3F7  * TMath::Sqrt(icF7_E/dev37);
    zetBR313 = vF3F13 * TMath::Sqrt(icF7_E/dev313);

    zetBRraw = zetBR313;
    
    zetBR = zetBR_c1 * zetBRraw + zetBR_c0;

    /*
    //@@@ tried to deduce Z with pla charge @@@
    Double_t de;
    //square average
    de = F7IC_E*(pla3q2e[0]+pla3q2e[1]*F3_Charge)*(pla7q2e[0]+pla7q2e[1]*F7_Charge)*(pla13_1q2e[0]+pla13_1q2e[1]*SBT1_Charge)*(pla13_2q2e[0]+pla13_2q2e[1]*SBT2_Charge);
    de = pow(de,1./5.); 

    //sum average
    de = (F7IC_E + (pla3q2e[0]+pla3q2e[1]*F3_Charge) + (pla7q2e[0]+pla7q2e[1]*F7_Charge) + (pla13_1q2e[0]+pla13_1q2e[1]*SBT1_Charge) + (pla13_2q2e[0]+pla13_2q2e[1]*SBT2_Charge))/5.;
    
    zetplaic = vF3F7 * TMath::Sqrt(de/(TMath::Log(2*m_e*vF3F7*vF3F7/Ionpair)-TMath::Log(1-betaF3F7*betaF3F7)-betaF3F7*betaF3F7));
    */

    bitset<4> f3x(0), f3y(0), f5x(0), f5y(0), f7x(0), f7y(0); // bit[2B 2A 1B 1A];

    // ===== TSum gate =====
    tsum_f31ax = ppacF31A_X_T1 + ppacF31A_X_T2;  
    tsum_f31bx = ppacF31B_X_T1 + ppacF31B_X_T2;
    tsum_f32ax = ppacF32A_X_T1 + ppacF32A_X_T2;
    tsum_f32bx = ppacF32B_X_T1 + ppacF32B_X_T2;
    tsum_f31ay = ppacF31A_Y_T1 + ppacF31A_Y_T2;
    tsum_f31by = ppacF31B_Y_T1 + ppacF31B_Y_T2;
    tsum_f32ay = ppacF32A_Y_T1 + ppacF32A_Y_T2;
    tsum_f32by = ppacF32B_Y_T1 + ppacF32B_Y_T2;
    tsum_f51ax = ppacF51A_X_T1 + ppacF51A_X_T2;
    tsum_f51bx = ppacF51B_X_T1 + ppacF51B_X_T2;
    tsum_f52ax = ppacF52A_X_T1 + ppacF52A_X_T2;
    tsum_f52bx = ppacF52B_X_T1 + ppacF52B_X_T2;
    tsum_f51ay = ppacF51A_Y_T1 + ppacF51A_Y_T2;
    tsum_f51by = ppacF51B_Y_T1 + ppacF51B_Y_T2;
    tsum_f52ay = ppacF52A_Y_T1 + ppacF52A_Y_T2;
    tsum_f52by = ppacF52B_Y_T1 + ppacF52B_Y_T2;
    tsum_f71ax = ppacF71A_X_T1 + ppacF71A_X_T2;
    tsum_f71bx = ppacF71B_X_T1 + ppacF71B_X_T2;
    tsum_f72ax = ppacF72A_X_T1 + ppacF72A_X_T2;
    tsum_f72bx = ppacF72B_X_T1 + ppacF72B_X_T2;
    tsum_f71ay = ppacF71A_Y_T1 + ppacF71A_Y_T2;
    tsum_f71by = ppacF71B_Y_T1 + ppacF71B_Y_T2;
    tsum_f72ay = ppacF72A_Y_T1 + ppacF72A_Y_T2;
    tsum_f72by = ppacF72B_Y_T1 + ppacF72B_Y_T2;
    
    if(cppac_low[0]  < tsum_f31ax && tsum_f31ax < cppac_up[0])  f3x.set(0);
    if(cppac_low[1]  < tsum_f31bx && tsum_f31bx < cppac_up[1])  f3x.set(1);
    if(cppac_low[2]  < tsum_f32ax && tsum_f32ax < cppac_up[2])  f3x.set(2);
    if(cppac_low[3]  < tsum_f32bx && tsum_f32bx < cppac_up[3])  f3x.set(3);
    if(cppac_low[4]  < tsum_f31ay && tsum_f31ay < cppac_up[4])  f3y.set(0);
    if(cppac_low[5]  < tsum_f31by && tsum_f31by < cppac_up[5])  f3y.set(1);
    if(cppac_low[6]  < tsum_f32ay && tsum_f32ay < cppac_up[6])  f3y.set(2);
    if(cppac_low[7]  < tsum_f32by && tsum_f32by < cppac_up[7])  f3y.set(3);
    if(cppac_low[8]  < tsum_f51ax && tsum_f51ax < cppac_up[8])  f5x.set(0);
    if(cppac_low[9]  < tsum_f51bx && tsum_f51bx < cppac_up[9])  f5x.set(1);
    if(cppac_low[10] < tsum_f52ax && tsum_f52ax < cppac_up[10]) f5x.set(2);
    if(cppac_low[11] < tsum_f52bx && tsum_f52bx < cppac_up[11]) f5x.set(3);
    if(cppac_low[12] < tsum_f51ay && tsum_f51ay < cppac_up[12]) f5y.set(0);
    if(cppac_low[13] < tsum_f51by && tsum_f51by < cppac_up[13]) f5y.set(1);
    if(cppac_low[14] < tsum_f52ay && tsum_f52ay < cppac_up[14]) f5y.set(2);
    if(cppac_low[15] < tsum_f52by && tsum_f52by < cppac_up[15]) f5y.set(3);
    if(cppac_low[16] < tsum_f71ax && tsum_f71ax < cppac_up[16]) f7x.set(0);
    if(cppac_low[17] < tsum_f71bx && tsum_f71bx < cppac_up[17]) f7x.set(1);
    if(cppac_low[18] < tsum_f72ax && tsum_f72ax < cppac_up[18]) f7x.set(2);
    if(cppac_low[19] < tsum_f72bx && tsum_f72bx < cppac_up[19]) f7x.set(3);
    if(cppac_low[20] < tsum_f71ay && tsum_f71ay < cppac_up[20]) f7y.set(0);
    if(cppac_low[21] < tsum_f71by && tsum_f71by < cppac_up[21]) f7y.set(1);
    if(cppac_low[22] < tsum_f72ay && tsum_f72ay < cppac_up[22]) f7y.set(2);
    if(cppac_low[23] < tsum_f72by && tsum_f72by < cppac_up[23]) f7y.set(3);

    cout << f3x << endl;

    //=== F31 X ===
    if(f3x[0]&f3x[1]) F31_X = (ppacF31A_X + ppacF31B_X)/2.;
    else if(f3x[0]&!f3x[1])   F31_X = ppacF31A_X;
    else if(f3x[1]&!f3x[0])   F31_X = ppacF31B_X;
    else{             //F31_X = pla3pos[1]*plaF3_dT + pla3pos[0];
                      ppacflag[0] = kFALSE;
    }

    //=== F32 X ===
    if(f3x[2]&f3x[3]) {
      F32_X = (ppacF32A_X + ppacF32B_X)/2.;
      cout << "tsum f32AX " << tsum_f32ax << endl;
      cout << "tsum f32BX " << tsum_f32bx << endl; 
      cout << "ppacF32A_X " << ppacF32A_X << endl;
      cout << "ppacF32B_X " << ppacF32B_X << endl;
      cout << "F32_X      " << F32_X      << endl;    
    }
    else if(f3x[2]&!f3x[3])   {
      F32_X = ppacF32A_X;
      //cout << F32_X << endl;
    }
    else if(f3x[3]&!f3x[2])   {
      F32_X = ppacF32B_X;
      //cout << F32_X << endl;
    }
    else{             //F32_X = pla3pos[1]*plaF3_dT + pla3pos[0];
                      ppacflag[0] = kFALSE;
    }

    //=== F3 X ===
    F3X = (F31_X + F32_X)/2.;

    //=== F31 Y ===
    if(f3y[0]&f3y[1]) F31_Y = (ppacF31A_Y + ppacF31B_Y)/2.;
    else if(f3y[0])   F31_Y = ppacF31A_Y;
    else if(f3y[1])   F31_Y = ppacF31B_Y;
    
    //=== F32 Y ===
    if(f3y[2]&f3y[3]) F32_Y = (ppacF32A_Y + ppacF32B_Y)/2.;
    else if(f3y[2])   F32_Y = ppacF32A_Y;
    else if(f3y[3])   F32_Y = ppacF32B_Y;

    //=== F3 Y ===
    F3Y = (F31_Y + F32_Y)/2.;

    //=== F51 X ===
    if(f5x[0]&f5x[1]) F51_X = (ppacF51A_X + ppacF51B_X)/2.;
    else if(f5x[0])   F51_X = ppacF51A_X;
    else if(f5x[1])   F51_X = ppacF51B_X;
    else              ppacflag[1] = kFALSE;

    //=== F52 X ===
    if(f5x[2]&f5x[3]) F52_X = (ppacF52A_X + ppacF52B_X)/2.;
    else if(f5x[2])   F52_X = ppacF52A_X;
    else if(f5x[3])   F52_X = ppacF52B_X;
    else              ppacflag[1] = kFALSE;

    //=== F5 X ===
    F5X = (F51_X + F52_X)/2.;

    //=== F51 Y ===
    if(f5y[0]&f5y[1]) F51_Y = (ppacF51A_Y + ppacF51B_Y)/2.;
    else if(f5y[0])   F51_Y = ppacF51A_Y;
    else if(f5y[1])   F51_Y = ppacF51B_Y;

    //=== F52 Y ===
    if(f5y[2]&f5y[3]) F52_Y = (ppacF52A_Y + ppacF52B_Y)/2.;
    else if(f5y[2])   F52_Y = ppacF52A_Y;
    else if(f5y[3])   F52_Y = ppacF52B_Y;
    
    //=== F5 Y ===
    F5Y = (F51_Y + F52_Y)/2.;

    //=== F71 X ===
    if(f7x[0]&f7x[1]) F71_X = (ppacF71A_X + ppacF71B_X)/2.;
    else if(f7x[0])   F71_X = ppacF71A_X;
    else if(f7x[1])   F71_X = ppacF71B_X;
    else{             //F71_X = pla7pos[1]*plaF7_dT + pla7pos[0];
      ppacflag[2] = kFALSE;
    }

    //=== F72 X ===
    if(f7x[2]&f7x[3]) F72_X = (ppacF72A_X + ppacF72B_X)/2.;
    else if(f7x[2])   F72_X = ppacF72A_X;
    else if(f7x[3])   F72_X = ppacF72B_X;
    else{             //F71_X = pla7pos[1]*plaF7_dT + pla7pos[0];
      ppacflag[2] = kFALSE;
    }

    //=== F7 X ===
    F7X = (F71_X + F72_X)/2.;

    //=== F71 Y ===
    if(f7y[0]&f7y[1]) F71_Y = (ppacF71A_Y + ppacF71B_Y)/2.;
    else if(f7y[0])   F71_Y = ppacF71A_Y;
    else if(f7y[1])   F71_Y = ppacF71B_Y;

    //=== F72 Y ===
    if(f7y[2]&f7y[3]) F72_Y = (ppacF72A_Y + ppacF72B_Y)/2.;
    else if(f7y[2])   F72_Y = ppacF72A_Y;
    else if(f7y[3])   F72_Y = ppacF72B_Y;
 
    //=== F7 Y ===
    F7Y = (F71_Y + F72_Y)/2.;
    

    if((f3x[0]|f3x[1])&(f3x[2]|f3x[3])) F3A = 1000.*TMath::ATan((F31_X-F32_X)/DistF3PPAC);
    F3B = 1000.*TMath::ATan((F31_Y-F32_Y)/DistF3PPAC);
   
    F5A = 1000.*TMath::ATan((F51_X-F52_X)/DistF5PPAC);
    F5B = 1000.*TMath::ATan((F51_Y-F52_Y)/DistF5PPAC);
    
    F7A = 1000.*TMath::ATan((F71_X-F72_X)/DistF7PPAC);
    F7B = 1000.*TMath::ATan((F71_Y-F72_Y)/DistF7PPAC);

    //reco1f7x = F7_TimeDiff*22.1925 + 22.1925*2.60488;
    //reco2f7x = F7_TimeDiff*10.05 + 10.05*2.60488;

    recoF3A = (ADF3F5*F5X - XDF3F5*F5A - (ADF3F5*XXF3F5 - XDF3F5*AXF3F5)*F3X)/(ADF3F5*XAF3F5 - XDF3F5*AAF3F5);
    
    deltaF3F5 = (F5X - XXF3F5*F3X - XAF3F5*recoF3A)/XDF3F5; //[%] 
    deltaF5F7 = (F7X - XXF5F7*F5X - XAF5F7*F5A)/XDF5F7;

    reco1deltaF5F7 = (reco1f7x - XXF5F7*F5X - XAF5F7*F5A)/XDF5F7;
    reco2deltaF5F7 = (reco2f7x - XXF5F7*F5X - XAF5F7*F5A)/XDF5F7;

    brhoF3F5 = Brho0F3F5*(1 + deltaF3F5*0.01);
    brhoF5F7 = Brho0F5F7*(1 + deltaF5F7*0.01);

    reco1brhoF5F7 = Brho0F5F7*(1 + reco1deltaF5F7*0.01);
    reco2brhoF5F7 = Brho0F5F7*(1 + reco2deltaF5F7*0.01);
    
    //aoqF3F5 = brhoF3F5*clight/m_u/betaF3F5/gammaF3F5 + 0.0017*F3X + 0.00005*F5X + 0.000002*F5X*F5X;
    //aoqF5F7 = brhoF5F7*clight/m_u/betaF5F7/gammaF5F7 + 0.0006*F5A + 0.00015*F7A;

    aoqF3F5 = brhoF3F5*clight/m_u/betaF3F13/gammaF3F13 + 0.0006*F3X;
    aoqF3F13 = brhoF5F7*clight/m_u/betaF3F13/gammaF3F13;
    aoqF5F7 = brhoF5F7*clight/m_u/betaF7F13/gammaF7F13;

    reco1aoq57 = reco1brhoF5F7*clight/m_u/betaF3F13/gammaF3F13;
    reco2aoq57 = reco2brhoF5F7*clight/m_u/betaF3F13/gammaF3F13;
    
    //aoqBR = aoqF5F7 + 0.0002*F7X + 0.00025*F5A + 0.00065*F7A + 0.00009*F5X - 0.000002*F5X*F5X + 0.000003*F7Y*F7Y + 0.000007*F7B*F7B;
    //aoqBR = aoqF5F7 + 0.0005*F7X - 0.00001*F7X*F7X + 0.00005*F7A - 0.00003*F5X + 0.0000002*F5X*F5X;
    //aoqBR = aoqF3F13 + 0.0006*F7X + 0.00007*F7A + 0.000000*1*F5X*F5X - 0.00015*F3X;
    //aoqBR = reco1aoq57 + 0.007*(F7_TimeDiff+2.6) - 0.0001*F7A + 0.0000001*F5X*F5X -0.00001*F5X;

    aoqBR = reco2aoq57 + 0.0000001*F5X*F5X -0.00001*F5X;

    //===== Cut by graphical cut ==========================================================
    //if(!cF3pla->IsInside(F3_TR-F3_TL,log(F3_QL/F3_QR))) BG_flag = 1;
    //if(!cF5pla->IsInside(F5_TR-F5_TL,log(F5_QL/F5_QR))) BG_flag = 1;
    //if(!cF7pla->IsInside(F7_TR-F7_TL,log(F7_QL/F7_QR))) BG_flag = 1;
    //if(!cF5Qchange->IsInside(aoqF3F5/aoqF5F7,zetBR)) BG_flag = 4; after finalizing the PID
    //if(cBR56Ca->IsInside(aoqBR,zetBR)) BR56Ca = 1;
    //if(cBR53Ca->IsInside(aoqBR,zetBR)) BR53Ca = 1;
    //if(cBR51K->IsInside(aoqBR,zetBR))  BR51K  = 1;
    //if(cBR56Sc->IsInside(aoqBR,zetBR)) BR56Sc = 1;
    
    if(!cplaF3->IsInside(plaF3_TR-plaF3_TL,log(plaF3_QL/plaF3_QR))){
      BG_flag = kFALSE;
      plaflag[0] = kFALSE;
    }
    if(!cplaF5->IsInside(plaF5_TR-plaF5_TL,log(plaF5_QL/plaF5_QR))){
      BG_flag = kFALSE;
      plaflag[1] = kFALSE;
    }
    if(!cplaF7->IsInside(plaF7_TR-plaF7_TL,log(plaF7_QL/plaF7_QR))){
      BG_flag = kFALSE;
      plaflag[2] = kFALSE;
    }

    if(cIC_low[0] > icF7_raw[0]-icF7_raw[1] || icF7_raw[0]-icF7_raw[1] > cIC_up[0]) icflag = kFALSE;
    if(cIC_low[1] > icF7_raw[1]-icF7_raw[2] || icF7_raw[1]-icF7_raw[2] > cIC_up[1]) icflag = kFALSE;
    if(cIC_low[2] > icF7_raw[2]-icF7_raw[3] || icF7_raw[2]-icF7_raw[3] > cIC_up[2]) icflag = kFALSE;
    if(cIC_low[3] > icF7_raw[3]-icF7_raw[4] || icF7_raw[3]-icF7_raw[4] > cIC_up[3]) icflag = kFALSE;
    if(cIC_low[4] > icF7_raw[4]-icF7_raw[5] || icF7_raw[4]-icF7_raw[5] > cIC_up[4]) icflag = kFALSE;

    if(!icflag) BG_flag = kFALSE;


    
  
    anatrB->Fill();
  }
  anafile->cd();
  anatrB->Write();
  anafile->Close();

  time(&stop);
  printf("Elapsed time: %.1f seconds\n",difftime(stop,start));

  //cout << nEntry/1000 << "k events have been treated." << endl;
  //cout << "Conversion Finished!" << endl;
  //cout << endl;
  printf("%d k events have been treated.\n",nEntry/1000);
  printf("Conversion Finished!\n\n");
  
  
  return 0;
}
