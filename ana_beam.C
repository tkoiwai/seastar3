#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>

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

using namespace std;
using namespace TMath;

int main(int argc, char *argv[]){

  Int_t FileNumber = TString(argv[1]).Atoi();

  //===== Load input file =================================================

  TString FileName = Form("/home/koiwai/analysis/rootfiles/all/run%04d_ALL.root",FileNumber);
  TFile *infile = TFile::Open(FileName);
  
  TTree *caltr;
  infile->GetObject("caltr",caltr);

  
  //===== in tree variables ===============================================
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;

  //=== Plastic ===
  Double_t F3_Charge, F5_Charge, F7_Charge;
  Double_t F3_Time, F5_Time, F7_Time;
  Double_t F3_TimeDiff, F5_TimeDiff, F7_TimeDiff;

  Double_t F3_QL, F3_QR, F5_QL, F5_QR, F7_QL, F7_QR;
  Double_t F3_TL, F3_TR, F5_TL, F5_TR, F7_TL, F7_TR;

  //=== PPAC ===
  Double_t F31A_X, F31A_Y, F31B_X, F31B_Y, F32A_X, F32A_Y, F32B_X, F32B_Y;
  Double_t F51A_X, F51A_Y, F51B_X, F51B_Y, F52A_X, F52A_Y, F52B_X, F52B_Y;
  Double_t F71A_X, F71A_Y, F71B_X, F71B_Y, F72A_X, F72A_Y, F72B_X, F72B_Y;
  
  Int_t F31A_X_T1, F31A_X_T2, F31A_Y_T1, F31A_Y_T2;
  Int_t F31B_X_T1, F31B_X_T2, F31B_Y_T1, F31B_Y_T2;
  Int_t F32A_X_T1, F32A_X_T2, F32A_Y_T1, F32A_Y_T2;
  Int_t F32B_X_T1, F32B_X_T2, F32B_Y_T1, F32B_Y_T2;
  Int_t F51A_X_T1, F51A_X_T2, F51A_Y_T1, F51A_Y_T2;
  Int_t F51B_X_T1, F51B_X_T2, F51B_Y_T1, F51B_Y_T2;
  Int_t F52A_X_T1, F52A_X_T2, F52A_Y_T1, F52A_Y_T2;
  Int_t F52B_X_T1, F52B_X_T2, F52B_Y_T1, F52B_Y_T2;
  Int_t F71A_X_T1, F71A_X_T2, F71A_Y_T1, F71A_Y_T2;
  Int_t F71B_X_T1, F71B_X_T2, F71B_Y_T1, F71B_Y_T2;
  Int_t F72A_X_T1, F72A_X_T2, F72A_Y_T1, F72A_Y_T2;
  Int_t F72B_X_T1, F72B_X_T2, F72B_Y_T1, F72B_Y_T2;

  //=== F7IC === 
  Double_t F7IC_E;
  
  Int_t F7IC_raw[6];

  //=== SBT ===
  Double_t SBT1_Charge, SBT2_Charge;
  Double_t SBT1_Time, SBT2_Time;
  Double_t SBT1_TimeDiff, SBT2_TimeDiff;
  
  //=== Set value ===
  caltr->SetBranchAddress("RunNumber",&RunNumber);
  caltr->SetBranchAddress("EventNumber",&EventNumber);
  
  caltr->SetBranchAddress("F3_Charge",&F3_Charge);
  caltr->SetBranchAddress("F5_Charge",&F5_Charge);
  caltr->SetBranchAddress("F7_Charge",&F7_Charge);
  caltr->SetBranchAddress("F3_Time",&F3_Time);
  caltr->SetBranchAddress("F5_Time",&F5_Time);
  caltr->SetBranchAddress("F7_Time",&F7_Time);
  caltr->SetBranchAddress("F3_TimeDiff",&F3_TimeDiff);
  caltr->SetBranchAddress("F5_TimeDiff",&F5_TimeDiff);
  caltr->SetBranchAddress("F7_TimeDiff",&F7_TimeDiff);
  
  caltr->SetBranchAddress("F3_QL",&F3_QL);
  caltr->SetBranchAddress("F3_QR",&F3_QR);
  caltr->SetBranchAddress("F5_QL",&F5_QL);
  caltr->SetBranchAddress("F5_QR",&F5_QR);
  caltr->SetBranchAddress("F7_QL",&F7_QL);
  caltr->SetBranchAddress("F7_QR",&F7_QR);
  
  caltr->SetBranchAddress("F3_TL",&F3_TL);
  caltr->SetBranchAddress("F3_TR",&F3_TR);
  caltr->SetBranchAddress("F5_TL",&F5_TL);
  caltr->SetBranchAddress("F5_TR",&F5_TR);
  caltr->SetBranchAddress("F7_TL",&F7_TL);
  caltr->SetBranchAddress("F7_TR",&F7_TR);
  
  caltr->SetBranchAddress("F31A_X",&F31A_X);
  caltr->SetBranchAddress("F31A_Y",&F31A_Y);
  caltr->SetBranchAddress("F31B_X",&F31B_X);
  caltr->SetBranchAddress("F31B_Y",&F31B_Y);
  caltr->SetBranchAddress("F32A_X",&F32A_X);
  caltr->SetBranchAddress("F32A_Y",&F32A_Y);
  caltr->SetBranchAddress("F32B_X",&F32B_X);
  caltr->SetBranchAddress("F32B_Y",&F32B_Y);
  
  caltr->SetBranchAddress("F51A_X",&F51A_X);
  caltr->SetBranchAddress("F51A_Y",&F51A_Y);
  caltr->SetBranchAddress("F51B_X",&F51B_X);
  caltr->SetBranchAddress("F51B_Y",&F51B_Y);
  caltr->SetBranchAddress("F52A_X",&F52A_X);
  caltr->SetBranchAddress("F52A_Y",&F52A_Y);
  caltr->SetBranchAddress("F52B_X",&F52B_X);
  caltr->SetBranchAddress("F52B_Y",&F52B_Y);
  
  caltr->SetBranchAddress("F71A_X",&F71A_X);
  caltr->SetBranchAddress("F71A_Y",&F71A_Y);
  caltr->SetBranchAddress("F71B_X",&F71B_X);
  caltr->SetBranchAddress("F71B_Y",&F71B_Y);
  caltr->SetBranchAddress("F72A_X",&F72A_X);
  caltr->SetBranchAddress("F72A_Y",&F72A_Y);
  caltr->SetBranchAddress("F72B_X",&F72B_X);
  caltr->SetBranchAddress("F72B_Y",&F72B_Y);
  
  caltr->SetBranchAddress("F31A_X_T1",&F31A_X_T1);
  caltr->SetBranchAddress("F31A_X_T2",&F31A_X_T2);
  caltr->SetBranchAddress("F31A_Y_T1",&F31A_Y_T1);
  caltr->SetBranchAddress("F31A_Y_T2",&F31A_Y_T2);
  caltr->SetBranchAddress("F31B_X_T1",&F31B_X_T1);
  caltr->SetBranchAddress("F31B_X_T2",&F31B_X_T2);
  caltr->SetBranchAddress("F31B_Y_T1",&F31B_Y_T1);
  caltr->SetBranchAddress("F31B_Y_T2",&F31B_Y_T2);
  caltr->SetBranchAddress("F32A_X_T1",&F32A_X_T1);
  caltr->SetBranchAddress("F32A_X_T2",&F32A_X_T2);
  caltr->SetBranchAddress("F32A_Y_T1",&F32A_Y_T1);
  caltr->SetBranchAddress("F32A_Y_T2",&F32A_Y_T2);
  caltr->SetBranchAddress("F32B_X_T1",&F32B_X_T1);
  caltr->SetBranchAddress("F32B_X_T2",&F32B_X_T2);
  caltr->SetBranchAddress("F32B_Y_T1",&F32B_Y_T1);
  caltr->SetBranchAddress("F32B_Y_T2",&F32B_Y_T2);
  caltr->SetBranchAddress("F51A_X_T1",&F51A_X_T1);
  caltr->SetBranchAddress("F51A_X_T2",&F51A_X_T2);
  caltr->SetBranchAddress("F51A_Y_T1",&F51A_Y_T1);
  caltr->SetBranchAddress("F51A_Y_T2",&F51A_Y_T2);
  caltr->SetBranchAddress("F51B_X_T1",&F51B_X_T1);
  caltr->SetBranchAddress("F51B_X_T2",&F51B_X_T2);
  caltr->SetBranchAddress("F51B_Y_T1",&F51B_Y_T1);
  caltr->SetBranchAddress("F51B_Y_T2",&F51B_Y_T2);
  caltr->SetBranchAddress("F52A_X_T1",&F52A_X_T1);
  caltr->SetBranchAddress("F52A_X_T2",&F52A_X_T2);
  caltr->SetBranchAddress("F52A_Y_T1",&F52A_Y_T1);
  caltr->SetBranchAddress("F52A_Y_T2",&F52A_Y_T2);
  caltr->SetBranchAddress("F52B_X_T1",&F52B_X_T1);
  caltr->SetBranchAddress("F52B_X_T2",&F52B_X_T2);
  caltr->SetBranchAddress("F52B_Y_T1",&F52B_Y_T1);
  caltr->SetBranchAddress("F52B_Y_T2",&F52B_Y_T2);
  caltr->SetBranchAddress("F71A_X_T1",&F71A_X_T1);
  caltr->SetBranchAddress("F71A_X_T2",&F71A_X_T2);
  caltr->SetBranchAddress("F71A_Y_T1",&F71A_Y_T1);
  caltr->SetBranchAddress("F71A_Y_T2",&F71A_Y_T2);
  caltr->SetBranchAddress("F71B_X_T1",&F71B_X_T1);
  caltr->SetBranchAddress("F71B_X_T2",&F71B_X_T2);
  caltr->SetBranchAddress("F71B_Y_T1",&F71B_Y_T1);
  caltr->SetBranchAddress("F71B_Y_T2",&F71B_Y_T2);
  caltr->SetBranchAddress("F72A_X_T1",&F72A_X_T1);
  caltr->SetBranchAddress("F72A_X_T2",&F72A_X_T2);
  caltr->SetBranchAddress("F72A_Y_T1",&F72A_Y_T1);
  caltr->SetBranchAddress("F72A_Y_T2",&F72A_Y_T2);
  caltr->SetBranchAddress("F72B_X_T1",&F72B_X_T1);
  caltr->SetBranchAddress("F72B_X_T2",&F72B_X_T2);
  caltr->SetBranchAddress("F72B_Y_T1",&F72B_Y_T1);
  caltr->SetBranchAddress("F72B_Y_T2",&F72B_Y_T2);
  
  caltr->SetBranchAddress("F7IC_E",&F7IC_E);
  caltr->SetBranchAddress("F7IC_raw",F7IC_raw);

  caltr->SetBranchAddress("SBT1_Charge",&SBT1_Charge);
  caltr->SetBranchAddress("SBT2_Charge",&SBT2_Charge);
  caltr->SetBranchAddress("SBT1_Time",&SBT1_Time);
  caltr->SetBranchAddress("SBT2_Time",&SBT2_Time);
  caltr->SetBranchAddress("SBT1_TimeDiff",&SBT1_TimeDiff);
  caltr->SetBranchAddress("SBT2_TimeDiff",&SBT2_TimeDiff);

  
  //===== Load CUT files ==================================================
  //=== Plastic (graphical cut)===
  TFile *cutfileF3pla = new TFile("/home/koiwai/analysis/cutfiles/cutf3pla.root");
  TCutG *cF3pla = (TCutG*)cutfileF3pla->Get("f3pla");
  TFile *cutfileF5pla = new TFile("/home/koiwai/analysis/cutfiles/cutf5pla.root");
  TCutG *cF5pla = (TCutG*)cutfileF5pla->Get("f5pla");
  TFile *cutfileF7pla = new TFile("/home/koiwai/analysis/cutfiles/cutf7pla.root");
  TCutG *cF7pla = (TCutG*)cutfileF7pla->Get("f7pla");

  /*
  //=== Chrage change @ F5 (graphical cut) ===
  TFile *cutfileF5Qchange = new TFile("/home/koiwai/analysis/cutfiles/cut_bigripsZvsF5Qchange.root");
  TCutG *cF5Qchange = (TCutG*)cutfileF5Qchange->Get("CUTG");
  */
  
  //=== PPAC Tsum gate ===
  ifstream fin;
  fin.open("/home/koiwai/analysis/cutfiles/cut_PPAC_Tsum.dat");
  if(fin.fail()){
    cout << "Error: file is not found." << endl;
    return 1;
  }
  
  Int_t cPPAC_Tsum_low[24], cPPAC_Tsum_up[24];
  
  for(Int_t cPPAC_index = 0;cPPAC_index<24;++cPPAC_index){
    fin >> cPPAC_Tsum_low[cPPAC_index] >> cPPAC_Tsum_up[cPPAC_index];
  }

  //=== 56Ca gate ===
  TFile *BRpid = new TFile("/home/koiwai/analysis/cutfiles/BRpid.root");
  TCutG *cBR56Ca = (TCutG*)BRpid->Get("BR56Ca");
  TCutG *cBR53Ca = (TCutG*)BRpid->Get("BR53Ca");
  TCutG *cBR51K  = (TCutG*)BRpid->Get("BR51K");
  
  //===== Load .dat files =====
  
  TEnv *env = new TEnv("/home/koiwai/analysis/db/geometry_psp17.dat");
  TEnv *env_pla  = new TEnv("/home/koiwai/analysis/db/pla2pos_Z.dat");
  TEnv *env_q = new TEnv("/home/koiwai/analysis/db/plaQ2f7icE.dat");
  
  //===== Create output file/tree =========================================
  TString AnaFileName = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNumber);
  TFile *anafile = new TFile(AnaFileName,"recreate");

  TTree *anatrB = new TTree("anatrB","anatrB");

  //===== Declear const.s =================================================
  double pla3pos[2];
  double pla7pos[2];
  //double pla7z[3];
  double pla3q2e[2];
  double pla7q2e[2];
  double pla13_1q2e[2];
  double pla13_2q2e[2];
  for(Int_t i=0;i<2;++i){
    pla3pos[i] = env_pla->GetValue(Form("pla3pos%d",i),0.0);
    pla7pos[i] = env_pla->GetValue(Form("pla7pos%d",i),0.0);
    //pla7z[i] = env_pla->GetValue(Form("pla7Z[%d]",i),0.0);
    pla3q2e[i]    = env_q->GetValue(Form("pla3q2e%d",i),0.0);
    pla7q2e[i]    = env_q->GetValue(Form("pla7q2e%d",i),0.0);
    pla13_1q2e[i] = env_q->GetValue(Form("pla13_1q2e%d",i),0.0);
    pla13_2q2e[i] = env_q->GetValue(Form("pla13_2q2e%d",i),0.0);
  }

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
  Double_t zetBR, zetplaic;

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

  //=== Imported from intree ===
  Double_t anaF3_Time, anaF5_Time, anaF7_Time;
  Double_t anaF3_TimeDiff, anaF5_TimeDiff, anaF7_TimeDiff;
  
  Int_t BG_flag; //flag for background
  Int_t BR56Ca,BR53Ca,BR51K;


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
  anatrB->Branch("zetplaic",&zetplaic);
 
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
  
  anatrB->Branch("anaF3_Time",&anaF3_Time); 
  anatrB->Branch("anaF5_Time",&anaF5_Time);
  anatrB->Branch("anaF7_Time",&anaF7_Time);
  anatrB->Branch("anaF3_TimeDiff",&anaF3_TimeDiff);
  anatrB->Branch("anaF5_TimeDiff",&anaF5_TimeDiff);
  anatrB->Branch("anaF7_TimeDiff",&anaF7_TimeDiff);
  
  
  anatrB->Branch("BG_flag",&BG_flag);
  anatrB->Branch("BR56Ca",&BR56Ca);
  anatrB->Branch("BR53Ca",&BR53Ca);
  anatrB->Branch("BR51K",&BR51K);
  
  infile->cd();


 //===== Begin LOOP ======================================================
  int nEntry = caltr->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    //for(int iEntry=0;iEntry<10;++iEntry){
    caltr->GetEntry(iEntry);

    if(iEntry%100 == 0){
      std::clog << iEntry/1000 << "k events treated..." << "\r";
    }

    EventNum = EventNumber;
    RunNum = RunNumber;
    
    //=== Initialization ===
    BG_flag = 0;
    BR56Ca = 0;
    BR53Ca = 0;
    BR51K  = 0;
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
    zetplaic = TMath::Sqrt(-1);
    
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

    anaF3_Time = TMath::Sqrt(-1);
    anaF5_Time = TMath::Sqrt(-1);
    anaF7_Time = TMath::Sqrt(-1);

    //=== Calculation ===
    tofF3F7 = F7_Time - F3_Time + OffsetF3F7;
    tofF3F5 = F5_Time - F3_Time + OffsetF3F5;
    tofF5F7 = F7_Time - F5_Time + OffsetF5F7;
    tofF7F13 = SBT1_Time - F7_Time + OffsetF7F13;
    tofF3F13 = SBT1_Time - F3_Time + OffsetF3F13;
    vF3F7 = DistF3F7/tofF3F7;
    vF3F5 = DistF3F5/tofF3F5;
    vF5F7 = DistF5F7/tofF5F7;
    vF7F13 = DistF7F13/tofF7F13;
    vF3F13 = DistF3F13/tofF3F13;

    betaF3F7 = vF3F7/clight;
    betaF3F5 = vF3F5/clight;
    betaF5F7 = vF5F7/clight;
    betaF7F13 = vF7F13/clight;
    betaF3F13 = vF3F13/clight;
    gammaF3F7 = 1/TMath::Sqrt(1.-betaF3F7*betaF3F7);
    gammaF3F5 = 1/TMath::Sqrt(1.-betaF3F5*betaF3F5);
    gammaF5F7 = 1/TMath::Sqrt(1.-betaF5F7*betaF5F7);
    gammaF7F13 = 1/TMath::Sqrt(1.-betaF7F13*betaF7F13);
    gammaF3F13 = 1/TMath::Sqrt(1.-betaF3F13*betaF3F13);

    zetBRraw = vF7F13 * TMath::Sqrt(F7IC_E/(TMath::Log(2*m_e*vF7F13*vF7F13/Ionpair)-TMath::Log(1-betaF7F13*betaF7F13)-betaF7F13*betaF7F13));
    
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
    
    //=== F3 Tsum gate ===
    if((cPPAC_Tsum_low[0]<F31A_X_T1+F31A_X_T2&&F31A_X_T1+F31A_X_T2<cPPAC_Tsum_up[0])&&(cPPAC_Tsum_low[1]<F31B_X_T1+F31B_X_T2&&F31B_X_T1+F31B_X_T2<cPPAC_Tsum_up[1])){
      F31_X = (F31A_X + F31B_X)/2.;
    }else if((cPPAC_Tsum_low[0]<F31A_X_T1+F31A_X_T2&&F31A_X_T1+F31A_X_T2<cPPAC_Tsum_up[0])&&!(cPPAC_Tsum_low[1]<F31B_X_T1+F31B_X_T2&&F31B_X_T1+F31B_X_T2<cPPAC_Tsum_up[1])){
      F31_X = F31A_X;
    }else if(!(cPPAC_Tsum_low[0]<F31A_X_T1+F31A_X_T2&&F31A_X_T1+F31A_X_T2<cPPAC_Tsum_up[0])&&(cPPAC_Tsum_low[1]<F31B_X_T1+F31B_X_T2&&F31B_X_T1+F31B_X_T2<cPPAC_Tsum_up[1])){
      F31_X = F31B_X;
    }else{
      F31_X = pla3pos[1]*F3_TimeDiff + pla3pos[0];
      //BG_flag = 2;
    }
    if((cPPAC_Tsum_low[2]<F32A_X_T1+F32A_X_T2&&F32A_X_T1+F32A_X_T2<cPPAC_Tsum_up[2])&&(cPPAC_Tsum_low[3]<F32B_X_T1+F32B_X_T2&&F32B_X_T1+F32B_X_T2<cPPAC_Tsum_up[3])){
      F32_X = (F32A_X + F32B_X)/2.;
    }else if((cPPAC_Tsum_low[2]<F32A_X_T1+F32A_X_T2&&F32A_X_T1+F32A_X_T2<cPPAC_Tsum_up[2])&&!(cPPAC_Tsum_low[3]<F32B_X_T1+F32B_X_T2&&F32B_X_T1+F32B_X_T2<cPPAC_Tsum_up[3])){
      F32_X = F32A_X;
    }else if(!(cPPAC_Tsum_low[2]<F32A_X_T1+F32A_X_T2&&F32A_X_T1+F32A_X_T2<cPPAC_Tsum_up[2])&&(cPPAC_Tsum_low[3]<F32B_X_T1+F32B_X_T2&&F32B_X_T1+F32B_X_T2<cPPAC_Tsum_up[3])){
      F32_X = F32B_X;
    }else{
      F32_X = pla3pos[1]*F3_TimeDiff + pla3pos[0];
      //BG_flag = 2;
    }
     if((cPPAC_Tsum_low[4]<F31A_Y_T1+F31A_Y_T2&&F31A_Y_T1+F31A_Y_T2<cPPAC_Tsum_up[4])&&(cPPAC_Tsum_low[5]<F31B_Y_T1+F31B_Y_T2&&F31B_Y_T1+F31B_Y_T2<cPPAC_Tsum_up[5])){
      F31_Y = (F31A_Y + F31B_Y)/2.;
    }else if((cPPAC_Tsum_low[4]<F31A_Y_T1+F31A_Y_T2&&F31A_Y_T1+F31A_Y_T2<cPPAC_Tsum_up[4])&&!(cPPAC_Tsum_low[5]<F31B_Y_T1+F31B_Y_T2&&F31B_Y_T1+F31B_Y_T2<cPPAC_Tsum_up[5])){
      F31_Y = F31A_Y;
    }else if(!(cPPAC_Tsum_low[4]<F31A_Y_T1+F31A_Y_T2&&F31A_Y_T1+F31A_Y_T2<cPPAC_Tsum_up[4])&&(cPPAC_Tsum_low[5]<F31B_Y_T1+F31B_Y_T2&&F31B_Y_T1+F31B_Y_T2<cPPAC_Tsum_up[5])){
      F31_Y = F31B_Y;
    }else{
       //BG_flag = 2;
    }
    if((cPPAC_Tsum_low[6]<F32A_Y_T1+F32A_Y_T2&&F32A_Y_T1+F32A_Y_T2<cPPAC_Tsum_up[6])&&(cPPAC_Tsum_low[7]<F32B_Y_T1+F32B_Y_T2&&F32B_Y_T1+F32B_Y_T2<cPPAC_Tsum_up[7])){
      F32_Y = (F32A_Y + F32B_Y)/2.;
    }else if((cPPAC_Tsum_low[6]<F32A_Y_T1+F32A_Y_T2&&F32A_Y_T1+F32A_Y_T2<cPPAC_Tsum_up[6])&&!(cPPAC_Tsum_low[7]<F32B_Y_T1+F32B_Y_T2&&F32B_Y_T1+F32B_Y_T2<cPPAC_Tsum_up[7])){
      F32_Y = F32A_Y;
    }else if(!(cPPAC_Tsum_low[6]<F32A_Y_T1+F32A_Y_T2&&F32A_Y_T1+F32A_Y_T2<cPPAC_Tsum_up[6])&&(cPPAC_Tsum_low[7]<F32B_Y_T1+F32B_Y_T2&&F32B_Y_T1+F32B_Y_T2<cPPAC_Tsum_up[7])){
      F32_Y = F32B_Y;
    }else{
      //BG_flag = 2;
    }
    //=== F5 ===
     if((cPPAC_Tsum_low[8]<F51A_X_T1+F51A_X_T2&&F51A_X_T1+F51A_X_T2<cPPAC_Tsum_up[8])&&(cPPAC_Tsum_low[9]<F51B_X_T1+F51B_X_T2&&F51B_X_T1+F51B_X_T2<cPPAC_Tsum_up[9])){
      F51_X = (F51A_X + F51B_X)/2.;
    }else if((cPPAC_Tsum_low[8]<F51A_X_T1+F51A_X_T2&&F51A_X_T1+F51A_X_T2<cPPAC_Tsum_up[8])&&!(cPPAC_Tsum_low[9]<F51B_X_T1+F51B_X_T2&&F51B_X_T1+F51B_X_T2<cPPAC_Tsum_up[9])){
      F51_X = F51A_X;
    }else if(!(cPPAC_Tsum_low[8]<F51A_X_T1+F51A_X_T2&&F51A_X_T1+F51A_X_T2<cPPAC_Tsum_up[8])&&(cPPAC_Tsum_low[9]<F51B_X_T1+F51B_X_T2&&F51B_X_T1+F51B_X_T2<cPPAC_Tsum_up[9])){
      F51_X = F51B_X;
    }else{
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[10]<F52A_X_T1+F52A_X_T2&&F52A_X_T1+F52A_X_T2<cPPAC_Tsum_up[10])&&(cPPAC_Tsum_low[11]<F52B_X_T1+F52B_X_T2&&F52B_X_T1+F52B_X_T2<cPPAC_Tsum_up[11])){
      F52_X = (F52A_X + F52B_X)/2.;
    }else if((cPPAC_Tsum_low[10]<F52A_X_T1+F52A_X_T2&&F52A_X_T1+F52A_X_T2<cPPAC_Tsum_up[10])&&!(cPPAC_Tsum_low[11]<F52B_X_T1+F52B_X_T2&&F52B_X_T1+F52B_X_T2<cPPAC_Tsum_up[11])){
      F52_X = F52A_X;
    }else if(!(cPPAC_Tsum_low[10]<F52A_X_T1+F52A_X_T2&&F52A_X_T1+F52A_X_T2<cPPAC_Tsum_up[10])&&(cPPAC_Tsum_low[11]<F52B_X_T1+F52B_X_T2&&F52B_X_T1+F52B_X_T2<cPPAC_Tsum_up[11])){
      F52_X = F52B_X;
    }else{
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[12]<F51A_Y_T1+F51A_Y_T2&&F51A_Y_T1+F51A_Y_T2<cPPAC_Tsum_up[12])&&(cPPAC_Tsum_low[13]<F51B_Y_T1+F51B_Y_T2&&F51B_Y_T1+F51B_Y_T2<cPPAC_Tsum_up[13])){
      F51_Y = (F51A_Y + F51B_Y)/2.;
    }else if((cPPAC_Tsum_low[12]<F51A_Y_T1+F51A_Y_T2&&F51A_Y_T1+F51A_Y_T2<cPPAC_Tsum_up[12])&&!(cPPAC_Tsum_low[13]<F51B_Y_T1+F51B_Y_T2&&F51B_Y_T1+F51B_Y_T2<cPPAC_Tsum_up[13])){
      F51_Y = F51A_Y;
    }else if(!(cPPAC_Tsum_low[12]<F51A_Y_T1+F51A_Y_T2&&F51A_Y_T1+F51A_Y_T2<cPPAC_Tsum_up[12])&&(cPPAC_Tsum_low[13]<F51B_Y_T1+F51B_Y_T2&&F51B_Y_T1+F51B_Y_T2<cPPAC_Tsum_up[13])){
      F51_Y = F51B_Y;
    }else{
      //BG_flag = 2;
    }
    if((cPPAC_Tsum_low[14]<F52A_Y_T1+F52A_Y_T2&&F52A_Y_T1+F52A_Y_T2<cPPAC_Tsum_up[14])&&(cPPAC_Tsum_low[15]<F52B_Y_T1+F52B_Y_T2&&F52B_Y_T1+F52B_Y_T2<cPPAC_Tsum_up[15])){
      F52_Y = (F52A_Y + F52B_Y)/2.;
    }else if((cPPAC_Tsum_low[14]<F52A_Y_T1+F52A_Y_T2&&F52A_Y_T1+F52A_Y_T2<cPPAC_Tsum_up[14])&&!(cPPAC_Tsum_low[15]<F52B_Y_T1+F52B_Y_T2&&F52B_Y_T1+F52B_Y_T2<cPPAC_Tsum_up[15])){
      F52_Y = F52A_Y;
    }else if(!(cPPAC_Tsum_low[14]<F52A_Y_T1+F52A_Y_T2&&F52A_Y_T1+F52A_Y_T2<cPPAC_Tsum_up[14])&&(cPPAC_Tsum_low[15]<F52B_Y_T1+F52B_Y_T2&&F52B_Y_T1+F52B_Y_T2<cPPAC_Tsum_up[15])){
      F52_Y = F52B_Y;
    }else{
      //BG_flag = 2;
    }
    //=== F7 ===
     if((cPPAC_Tsum_low[16]<F71A_X_T1+F71A_X_T2&&F71A_X_T1+F71A_X_T2<cPPAC_Tsum_up[16])&&(cPPAC_Tsum_low[17]<F71B_X_T1+F71B_X_T2&&F71B_X_T1+F71B_X_T2<cPPAC_Tsum_up[17])){
      F71_X = (F71A_X + F71B_X)/2.;
    }else if((cPPAC_Tsum_low[16]<F71A_X_T1+F71A_X_T2&&F71A_X_T1+F71A_X_T2<cPPAC_Tsum_up[16])&&!(cPPAC_Tsum_low[17]<F71B_X_T1+F71B_X_T2&&F71B_X_T1+F71B_X_T2<cPPAC_Tsum_up[17])){
      F71_X = F71A_X;
    }else if(!(cPPAC_Tsum_low[16]<F71A_X_T1+F71A_X_T2&&F71A_X_T1+F71A_X_T2<cPPAC_Tsum_up[16])&&(cPPAC_Tsum_low[17]<F71B_X_T1+F71B_X_T2&&F71B_X_T1+F71B_X_T2<cPPAC_Tsum_up[17])){
      F71_X = F71B_X;
    }else{
       F71_X = pla7pos[1] * F7_TimeDiff + pla7pos[0];
       //BG_flag = 2;
    }
    if((cPPAC_Tsum_low[18]<F72A_X_T1+F72A_X_T2&&F72A_X_T1+F72A_X_T2<cPPAC_Tsum_up[18])&&(cPPAC_Tsum_low[19]<F72B_X_T1+F72B_X_T2&&F72B_X_T1+F72B_X_T2<cPPAC_Tsum_up[19])){
      F72_X = (F72A_X + F72B_X)/2.;
    }else if((cPPAC_Tsum_low[18]<F72A_X_T1+F72A_X_T2&&F72A_X_T1+F72A_X_T2<cPPAC_Tsum_up[18])&&!(cPPAC_Tsum_low[19]<F72B_X_T1+F72B_X_T2&&F72B_X_T1+F72B_X_T2<cPPAC_Tsum_up[19])){
      F72_X = F72A_X;
    }else if(!(cPPAC_Tsum_low[18]<F72A_X_T1+F72A_X_T2&&F72A_X_T1+F72A_X_T2<cPPAC_Tsum_up[18])&&(cPPAC_Tsum_low[19]<F72B_X_T1+F72B_X_T2&&F72B_X_T1+F72B_X_T2<cPPAC_Tsum_up[19])){
      F72_X = F72B_X;
    }else{
      F72_X = pla7pos[1] * F7_TimeDiff + pla7pos[0];
      //BG_flag = 2;
    }
     if((cPPAC_Tsum_low[20]<F71A_Y_T1+F71A_Y_T2&&F71A_Y_T1+F71A_Y_T2<cPPAC_Tsum_up[20])&&(cPPAC_Tsum_low[21]<F71B_Y_T1+F71B_Y_T2&&F71B_Y_T1+F71B_Y_T2<cPPAC_Tsum_up[21])){
      F71_Y = (F71A_Y + F71B_Y)/2.;
    }else if((cPPAC_Tsum_low[20]<F71A_Y_T1+F71A_Y_T2&&F71A_Y_T1+F71A_Y_T2<cPPAC_Tsum_up[20])&&!(cPPAC_Tsum_low[21]<F71B_Y_T1+F71B_Y_T2&&F71B_Y_T1+F71B_Y_T2<cPPAC_Tsum_up[21])){
      F71_Y = F71A_Y;
    }else if(!(cPPAC_Tsum_low[20]<F71A_Y_T1+F71A_Y_T2&&F71A_Y_T1+F71A_Y_T2<cPPAC_Tsum_up[20])&&(cPPAC_Tsum_low[21]<F71B_Y_T1+F71B_Y_T2&&F71B_Y_T1+F71B_Y_T2<cPPAC_Tsum_up[21])){
      F71_Y = F71B_Y;
    }else{
       //BG_flag = 2;
    }
    if((cPPAC_Tsum_low[22]<F72A_Y_T1+F72A_Y_T2&&F72A_Y_T1+F72A_Y_T2<cPPAC_Tsum_up[22])&&(cPPAC_Tsum_low[23]<F72B_Y_T1+F72B_Y_T2&&F72B_Y_T1+F72B_Y_T2<cPPAC_Tsum_up[23])){
      F72_Y = (F72A_Y + F72B_Y)/2.;
    }else if((cPPAC_Tsum_low[22]<F72A_Y_T1+F72A_Y_T2&&F72A_Y_T1+F72A_Y_T2<cPPAC_Tsum_up[22])&&!(cPPAC_Tsum_low[23]<F72B_Y_T1+F72B_Y_T2&&F72B_Y_T1+F72B_Y_T2<cPPAC_Tsum_up[23])){
      F72_Y = F72A_Y;
    }else if(!(cPPAC_Tsum_low[22]<F72A_Y_T1+F72A_Y_T2&&F72A_Y_T1+F72A_Y_T2<cPPAC_Tsum_up[22])&&(cPPAC_Tsum_low[23]<F72B_Y_T1+F72B_Y_T2&&F72B_Y_T1+F72B_Y_T2<cPPAC_Tsum_up[23])){
      F72_Y = F72B_Y;
    }else{
      //BG_flag = 2;
    }
    
    
    F3X = (F31_X + F32_X)/2.;
    F3Y = (F31_Y + F32_Y)/2.;
    if(F31_X!=F32_X) F3A = 1000.*TMath::ATan((F31_X-F32_X)/DistF3PPAC);
    F3B = 1000.*TMath::ATan((F31_Y-F32_Y)/DistF3PPAC);
   
    F5X = (F51_X + F52_X)/2.;
    F5Y = (F51_Y + F52_Y)/2.;
    F5A = 1000.*TMath::ATan((F51_X-F52_X)/DistF5PPAC);
    F5B = 1000.*TMath::ATan((F51_Y-F52_Y)/DistF5PPAC);
    
    F7X = (F71_X + F72_X)/2.;
    F7Y = (F71_Y + F72_Y)/2.;
    F7A = 1000.*TMath::ATan((F71_X-F72_X)/DistF7PPAC);
    F7B = 1000.*TMath::ATan((F71_Y-F72_Y)/DistF7PPAC);

    recoF3A = (ADF3F5*F5X - XDF3F5*F5A - (ADF3F5*XXF3F5 - XDF3F5*AXF3F5)*F3X)/(ADF3F5*XAF3F5 - XDF3F5*AAF3F5);
    
    deltaF3F5 = (F5X - XXF3F5*F3X - XAF3F5*recoF3A)/XDF3F5; //[%] 
    deltaF5F7 = (F7X - XXF5F7*F5X - XAF5F7*F5A)/XDF5F7;

    brhoF3F5 = Brho0F3F5*(1 + deltaF3F5*0.01);
    brhoF5F7 = Brho0F5F7*(1 + deltaF5F7*0.01);

    //aoqF3F5 = brhoF3F5*clight/m_u/betaF3F5/gammaF3F5 + 0.0017*F3X + 0.00005*F5X + 0.000002*F5X*F5X;
    //aoqF5F7 = brhoF5F7*clight/m_u/betaF5F7/gammaF5F7 + 0.0006*F5A + 0.00015*F7A;

    aoqF3F5 = brhoF3F5*clight/m_u/betaF3F13/gammaF3F13;
    aoqF3F13 = brhoF5F7*clight/m_u/betaF3F13/gammaF3F13;
    aoqF5F7 = brhoF5F7*clight/m_u/betaF7F13/gammaF7F13;
    
    //aoqBR = aoqF5F7 + 0.0002*F7X + 0.00025*F5A + 0.00065*F7A + 0.00009*F5X - 0.000002*F5X*F5X + 0.000003*F7Y*F7Y + 0.000007*F7B*F7B;
    //aoqBR = aoqF5F7 + 0.0005*F7X - 0.00001*F7X*F7X + 0.00005*F7A - 0.00003*F5X + 0.0000002*F5X*F5X;

    aoqBR = aoqF3F13 + 0.0006*F7X + 0.0000*F7A + 0.000000*1*F5X*F5X - 0.00015*F3X;
    
    anaF3_Time = F3_Time;
    anaF5_Time = F5_Time;
    anaF7_Time = F7_Time;
    anaF3_TimeDiff = F3_TimeDiff;
    anaF5_TimeDiff = F5_TimeDiff;
    anaF7_TimeDiff = F7_TimeDiff;

    //===== Cut by graphical cut ==========================================================
    if(!cF3pla->IsInside(F3_TR-F3_TL,log(F3_QL/F3_QR))) BG_flag = 1;
    if(!cF5pla->IsInside(F5_TR-F5_TL,log(F5_QL/F5_QR))) BG_flag = 1;
    if(!cF7pla->IsInside(F7_TR-F7_TL,log(F7_QL/F7_QR))) BG_flag = 1;
    //if(!cF5Qchange->IsInside(aoqF3F5/aoqF5F7,zetBR)) BG_flag = 4; after finalizing the PID
    if(cBR56Ca->IsInside(aoqBR,zetBR)) BR56Ca = 1;
    if(cBR53Ca->IsInside(aoqBR,zetBR)) BR53Ca = 1;
    if(cBR51K->IsInside(aoqBR,zetBR)) BR51K  = 1;
  
    anatrB->Fill();
  }
  anafile->cd();
  anatrB->Write();
  anafile->Close();

  return 0;
}
