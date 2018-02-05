#include<stdio.h>
#include<iostream>
#include<math.h>

#include"TROOT.h"
#include"TDirectory.h"
#include"TFile.h"
#include"TTree.h"
#include"TCut.h"
#include"TCutG.h"
#include"TMath.h"

using namespace std;

void ana_beam(){

  //=====Load input file=================================================
  TFile *infile = TFile::Open("rootfiles/run0056/run0056_BEAM.root");
  TTree *caltreeB;
  infile->GetObject("caltreeB",caltreeB);

  //=====in tree variables===============================================
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;

  //Plastic
  Double_t F3_Charge, F5_Charge, F7_Charge;
  Double_t F3_Time, F5_Time, F7_Time;
  Double_t F3_TimeDiff, F5_TimeDiff, F7_TimeDiff;

  Double_t F3_QL, F3_QR, F5_QL, F5_QR, F7_QL, F7_QR;
  Double_t F3_TL, F3_TR, F5_TL, F5_TR, F7_TL, F7_TR;

  //PPAC
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
  
  //IC
  Double_t F7IC_E;
  
  Int_t F7IC_raw[6];
  
  
  caltreeB->SetBranchAddress("RunNumber",&RunNumber);
  caltreeB->SetBranchAddress("EventNumber",&EventNumber);
  
  caltreeB->SetBranchAddress("F3_Charge",&F3_Charge);
  caltreeB->SetBranchAddress("F5_Charge",&F5_Charge);
  caltreeB->SetBranchAddress("F7_Charge",&F7_Charge);
  caltreeB->SetBranchAddress("F3_Time",&F3_Time);
  caltreeB->SetBranchAddress("F5_Time",&F5_Time);
  caltreeB->SetBranchAddress("F7_Time",&F7_Time);
  caltreeB->SetBranchAddress("F3_TimeDiff",&F3_TimeDiff);
  caltreeB->SetBranchAddress("F5_TimeDiff",&F5_TimeDiff);
  caltreeB->SetBranchAddress("F7_TimeDiff",&F7_TimeDiff);
  
  caltreeB->SetBranchAddress("F3_QL",&F3_QL);
  caltreeB->SetBranchAddress("F3_QR",&F3_QR);
  caltreeB->SetBranchAddress("F5_QL",&F5_QL);
  caltreeB->SetBranchAddress("F5_QR",&F5_QR);
  caltreeB->SetBranchAddress("F7_QL",&F7_QL);
  caltreeB->SetBranchAddress("F7_QR",&F7_QR);
  
  caltreeB->SetBranchAddress("F3_TL",&F3_TL);
  caltreeB->SetBranchAddress("F3_TR",&F3_TR);
  caltreeB->SetBranchAddress("F5_TL",&F5_TL);
  caltreeB->SetBranchAddress("F5_TR",&F5_TR);
  caltreeB->SetBranchAddress("F7_TL",&F7_TL);
  caltreeB->SetBranchAddress("F7_TR",&F7_TR);
  
  caltreeB->SetBranchAddress("F31A_X",&F31A_X);
  caltreeB->SetBranchAddress("F31A_Y",&F31A_Y);
  caltreeB->SetBranchAddress("F31B_X",&F31B_X);
  caltreeB->SetBranchAddress("F31B_Y",&F31B_Y);
  caltreeB->SetBranchAddress("F32A_X",&F32A_X);
  caltreeB->SetBranchAddress("F32A_Y",&F32A_Y);
  caltreeB->SetBranchAddress("F32B_X",&F32B_X);
  caltreeB->SetBranchAddress("F32B_Y",&F32B_Y);
  
  caltreeB->SetBranchAddress("F51A_X",&F51A_X);
  caltreeB->SetBranchAddress("F51A_Y",&F51A_Y);
  caltreeB->SetBranchAddress("F51B_X",&F51B_X);
  caltreeB->SetBranchAddress("F51B_Y",&F51B_Y);
  caltreeB->SetBranchAddress("F52A_X",&F52A_X);
  caltreeB->SetBranchAddress("F52A_Y",&F52A_Y);
  caltreeB->SetBranchAddress("F52B_X",&F52B_X);
  caltreeB->SetBranchAddress("F52B_Y",&F52B_Y);
  
  caltreeB->SetBranchAddress("F71A_X",&F71A_X);
  caltreeB->SetBranchAddress("F71A_Y",&F71A_Y);
  caltreeB->SetBranchAddress("F71B_X",&F71B_X);
  caltreeB->SetBranchAddress("F71B_Y",&F71B_Y);
  caltreeB->SetBranchAddress("F72A_X",&F72A_X);
  caltreeB->SetBranchAddress("F72A_Y",&F72A_Y);
  caltreeB->SetBranchAddress("F72B_X",&F72B_X);
  caltreeB->SetBranchAddress("F72B_Y",&F72B_Y);
  
  caltreeB->SetBranchAddress("F31A_X_T1",&F31A_X_T1);
  caltreeB->SetBranchAddress("F31A_X_T2",&F31A_X_T2);
  caltreeB->SetBranchAddress("F31A_Y_T1",&F31A_Y_T1);
  caltreeB->SetBranchAddress("F31A_Y_T2",&F31A_Y_T2);
  caltreeB->SetBranchAddress("F31B_X_T1",&F31B_X_T1);
  caltreeB->SetBranchAddress("F31B_X_T2",&F31B_X_T2);
  caltreeB->SetBranchAddress("F31B_Y_T1",&F31B_Y_T1);
  caltreeB->SetBranchAddress("F31B_Y_T2",&F31B_Y_T2);
  caltreeB->SetBranchAddress("F32A_X_T1",&F32A_X_T1);
  caltreeB->SetBranchAddress("F32A_X_T2",&F32A_X_T2);
  caltreeB->SetBranchAddress("F32A_Y_T1",&F32A_Y_T1);
  caltreeB->SetBranchAddress("F32A_Y_T2",&F32A_Y_T2);
  caltreeB->SetBranchAddress("F32B_X_T1",&F32B_X_T1);
  caltreeB->SetBranchAddress("F32B_X_T2",&F32B_X_T2);
  caltreeB->SetBranchAddress("F32B_Y_T1",&F32B_Y_T1);
  caltreeB->SetBranchAddress("F32B_Y_T2",&F32B_Y_T2);
  caltreeB->SetBranchAddress("F51A_X_T1",&F51A_X_T1);
  caltreeB->SetBranchAddress("F51A_X_T2",&F51A_X_T2);
  caltreeB->SetBranchAddress("F51A_Y_T1",&F51A_Y_T1);
  caltreeB->SetBranchAddress("F51A_Y_T2",&F51A_Y_T2);
  caltreeB->SetBranchAddress("F51B_X_T1",&F51B_X_T1);
  caltreeB->SetBranchAddress("F51B_X_T2",&F51B_X_T2);
  caltreeB->SetBranchAddress("F51B_Y_T1",&F51B_Y_T1);
  caltreeB->SetBranchAddress("F51B_Y_T2",&F51B_Y_T2);
  caltreeB->SetBranchAddress("F52A_X_T1",&F52A_X_T1);
  caltreeB->SetBranchAddress("F52A_X_T2",&F52A_X_T2);
  caltreeB->SetBranchAddress("F52A_Y_T1",&F52A_Y_T1);
  caltreeB->SetBranchAddress("F52A_Y_T2",&F52A_Y_T2);
  caltreeB->SetBranchAddress("F52B_X_T1",&F52B_X_T1);
  caltreeB->SetBranchAddress("F52B_X_T2",&F52B_X_T2);
  caltreeB->SetBranchAddress("F52B_Y_T1",&F52B_Y_T1);
  caltreeB->SetBranchAddress("F52B_Y_T2",&F52B_Y_T2);
  caltreeB->SetBranchAddress("F71A_X_T1",&F71A_X_T1);
  caltreeB->SetBranchAddress("F71A_X_T2",&F71A_X_T2);
  caltreeB->SetBranchAddress("F71A_Y_T1",&F71A_Y_T1);
  caltreeB->SetBranchAddress("F71A_Y_T2",&F71A_Y_T2);
  caltreeB->SetBranchAddress("F71B_X_T1",&F71B_X_T1);
  caltreeB->SetBranchAddress("F71B_X_T2",&F71B_X_T2);
  caltreeB->SetBranchAddress("F71B_Y_T1",&F71B_Y_T1);
  caltreeB->SetBranchAddress("F71B_Y_T2",&F71B_Y_T2);
  caltreeB->SetBranchAddress("F72A_X_T1",&F72A_X_T1);
  caltreeB->SetBranchAddress("F72A_X_T2",&F72A_X_T2);
  caltreeB->SetBranchAddress("F72A_Y_T1",&F72A_Y_T1);
  caltreeB->SetBranchAddress("F72A_Y_T2",&F72A_Y_T2);
  caltreeB->SetBranchAddress("F72B_X_T1",&F72B_X_T1);
  caltreeB->SetBranchAddress("F72B_X_T2",&F72B_X_T2);
  caltreeB->SetBranchAddress("F72B_Y_T1",&F72B_Y_T1);
  caltreeB->SetBranchAddress("F72B_Y_T2",&F72B_Y_T2);
  
  caltreeB->SetBranchAddress("F7IC_E",&F7IC_E);
  caltreeB->SetBranchAddress("F7IC_raw",F7IC_raw);


  
  //=====Load CUT files==================================================
  TFile *cutfileF3pla = new TFile("cutfiles/cut_F3pla.root");
  TCutG *cF3pla = (TCutG*)cutfileF3pla->Get("F3pla");
  TFile *cutfileF7pla = new TFile("cutfiles/cut_F7pla.root");
  TCutG *cF7pla = (TCutG*)cutfileF7pla->Get("F7pla");

  TFile *cutfileF31A_X_Tsum = new TFile("cutfiles/cut_F31A_X_Tsum.root");
  TCutG *cF31A_X_Tsum = (TCutG*)cutfileF31A_X_Tsum->Get("F31A_X_Tsum");
  TFile *cutfileF31A_Y_Tsum = new TFile("cutfiles/cut_F31A_Y_Tsum.root");
  TCutG *cF31A_Y_Tsum = (TCutG*)cutfileF31A_Y_Tsum->Get("F31A_Y_Tsum");
  TFile *cutfileF31B_X_Tsum = new TFile("cutfiles/cut_F31B_X_Tsum.root");
  TCutG *cF31B_X_Tsum = (TCutG*)cutfileF31B_X_Tsum->Get("F31B_X_Tsum");
  TFile *cutfileF31B_Y_Tsum = new TFile("cutfiles/cut_F31B_Y_Tsum.root");
  TCutG *cF31B_Y_Tsum = (TCutG*)cutfileF31B_Y_Tsum->Get("F31B_Y_Tsum");
  TFile *cutfileF32A_X_Tsum = new TFile("cutfiles/cut_F32A_X_Tsum.root");
  TCutG *cF32A_X_Tsum = (TCutG*)cutfileF32A_X_Tsum->Get("F32A_X_Tsum");
  TFile *cutfileF32A_Y_Tsum = new TFile("cutfiles/cut_F32A_Y_Tsum.root");
  TCutG *cF32A_Y_Tsum = (TCutG*)cutfileF32A_Y_Tsum->Get("F32A_Y_Tsum");
  TFile *cutfileF32B_X_Tsum = new TFile("cutfiles/cut_F32B_X_Tsum.root");
  TCutG *cF32B_X_Tsum = (TCutG*)cutfileF32B_X_Tsum->Get("F32B_X_Tsum");
  TFile *cutfileF32B_Y_Tsum = new TFile("cutfiles/cut_F32B_Y_Tsum.root");
  TCutG *cF32B_Y_Tsum = (TCutG*)cutfileF32B_Y_Tsum->Get("F32B_Y_Tsum");
  TFile *cutfileF51A_X_Tsum = new TFile("cutfiles/cut_F51A_X_Tsum.root");
  TCutG *cF51A_X_Tsum = (TCutG*)cutfileF51A_X_Tsum->Get("F51A_X_Tsum");
  TFile *cutfileF51A_Y_Tsum = new TFile("cutfiles/cut_F51A_Y_Tsum.root");
  TCutG *cF51A_Y_Tsum = (TCutG*)cutfileF51A_Y_Tsum->Get("F51A_Y_Tsum");
  TFile *cutfileF51B_X_Tsum = new TFile("cutfiles/cut_F51B_X_Tsum.root");
  TCutG *cF51B_X_Tsum = (TCutG*)cutfileF51B_X_Tsum->Get("F51B_X_Tsum");
  TFile *cutfileF51B_Y_Tsum = new TFile("cutfiles/cut_F51B_Y_Tsum.root");
  TCutG *cF51B_Y_Tsum = (TCutG*)cutfileF51B_Y_Tsum->Get("F51B_Y_Tsum");
  TFile *cutfileF52A_X_Tsum = new TFile("cutfiles/cut_F52A_X_Tsum.root");
  TCutG *cF52A_X_Tsum = (TCutG*)cutfileF52A_X_Tsum->Get("F52A_X_Tsum");
  TFile *cutfileF52A_Y_Tsum = new TFile("cutfiles/cut_F52A_Y_Tsum.root");
  TCutG *cF52A_Y_Tsum = (TCutG*)cutfileF52A_Y_Tsum->Get("F52A_Y_Tsum");
  TFile *cutfileF52B_X_Tsum = new TFile("cutfiles/cut_F52B_X_Tsum.root");
  TCutG *cF52B_X_Tsum = (TCutG*)cutfileF52B_X_Tsum->Get("F52B_X_Tsum");
  TFile *cutfileF52B_Y_Tsum = new TFile("cutfiles/cut_F52B_Y_Tsum.root");
  TCutG *cF52B_Y_Tsum = (TCutG*)cutfileF52B_Y_Tsum->Get("F52B_Y_Tsum");
  TFile *cutfileF71A_X_Tsum = new TFile("cutfiles/cut_F71A_X_Tsum.root");
  TCutG *cF71A_X_Tsum = (TCutG*)cutfileF71A_X_Tsum->Get("F71A_X_Tsum");
  TFile *cutfileF71A_Y_Tsum = new TFile("cutfiles/cut_F71A_Y_Tsum.root");
  TCutG *cF71A_Y_Tsum = (TCutG*)cutfileF71A_Y_Tsum->Get("F71A_Y_Tsum");
  TFile *cutfileF71B_X_Tsum = new TFile("cutfiles/cut_F71B_X_Tsum.root");
  TCutG *cF71B_X_Tsum = (TCutG*)cutfileF71B_X_Tsum->Get("F71B_X_Tsum");
  TFile *cutfileF71B_Y_Tsum = new TFile("cutfiles/cut_F71B_Y_Tsum.root");
  TCutG *cF71B_Y_Tsum = (TCutG*)cutfileF71B_Y_Tsum->Get("F71B_Y_Tsum");
  TFile *cutfileF72A_X_Tsum = new TFile("cutfiles/cut_F72A_X_Tsum.root");
  TCutG *cF72A_X_Tsum = (TCutG*)cutfileF72A_X_Tsum->Get("F72A_X_Tsum");
  TFile *cutfileF72A_Y_Tsum = new TFile("cutfiles/cut_F72A_Y_Tsum.root");
  TCutG *cF72A_Y_Tsum = (TCutG*)cutfileF72A_Y_Tsum->Get("F72A_Y_Tsum");
  TFile *cutfileF72B_X_Tsum = new TFile("cutfiles/cut_F72B_X_Tsum.root");
  TCutG *cF72B_X_Tsum = (TCutG*)cutfileF72B_X_Tsum->Get("F72B_X_Tsum");
  TFile *cutfileF72B_Y_Tsum = new TFile("cutfiles/cut_F72B_Y_Tsum.root");
  TCutG *cF72B_Y_Tsum = (TCutG*)cutfileF72B_Y_Tsum->Get("F72B_Y_Tsum");

  //=====Create output file/tree=========================================
  TFile *anafile = new TFile("rootfiles/ana_beam.root","RECREATE");
  TTree *anatreeB = new TTree("anatreeB","anatreeB");

  //=====Declear const.s=================================================
  double DistF3F7 = 46568.; //mm
  double DistF3F5 = 23284.;
  double DistF5F7 = 23284.;
  double OffsetF3F7 = 292.379; //nsec
  double OffsetF3F5 = 159.572;
  double OffsetF5F7 = 132.807;
  double Ionpair = 4.866; //keV
  double m_e = 511.; //keV
  double m_u = 931.494; //MeV
  double clight = 299.792458; //mm/nsec
  double zet_c1 = 0.0374179; //slope
  double zet_c2 = -4.2354; //const.
  double DistF3PPAC = 890.; //mm
  double DistF5PPAC = 650.;
  double DistF7PPAC = 945.;
  //===optical matrix===
  double XXF3F5 = 0.917467; //from matrix/mat1.mat //[]
  double XAF3F5 = -0.018721; //[mm/mrad]
  double XDF3F5 = -0.031663; //[mm/%]
  double AXF3F5 = -0.00520039; //[mrad/mm]
  double AAF3F5 = 1.09006; //[]
  double ADF3F5 = 1.88156; //[mrad/%]
  double XXF5F7 = 1.09101; //from matrix/mat2.mat //[]
  double XAF5F7 = -0.0172247; //[mm/mrad]
  double XDF5F7 = -0.00277197; //[mm/%]
  double AXF5F7 = 0.020415; //[mrad/mm]
  double AAF5F7 = 0.916262; //[]
  double ADF5F7 = -1.72349; //[mrad/%]
  //===central Brho===
  double Brho0F3F5 = 7.0621; //Tm
  double Brho0F5F7 = 6.854; //Tm
  
  //=====Define ana variables============================================
  //===for Z===
  Double_t vF3F7, vF3F5, vF5F7;
  Double_t betaF3F7, betaF3F5, betaF5F7;
  Double_t gammaF3F7, gammaF3F5, gammaF5F7;
  Double_t raw_zet;
  Double_t zet;

  //===for A/Q===
  Double_t F3X, F3Y, F3A, F3B; //X, Y:[mm]
  Double_t F5X, F5Y, F5A, F5B; //A, B:[mrad]
  Double_t F7X, F7Y, F7A, F7B;

  Double_t F31_X, F31_Y, F32_X, F32_Y;
  Double_t F51_X, F51_Y, F52_X, F52_Y;
  Double_t F71_X, F71_Y, F72_X, F72_Y;

  Double_t deltaF3F5, deltaF5F7;
  Double_t brhoF3F5, brhoF5F7;
  Double_t aoqF3F5, aoqF5F7;
  Double_t recoF3A;

  //===imported from intree===
  Double_t anaF3_Time, anaF5_Time, anaF7_Time;
  
  Int_t BG_flag; //flag for background

  anatreeB->Branch("vF3F7",&vF3F7);
  anatreeB->Branch("vF3F5",&vF3F5);
  anatreeB->Branch("vF5F7",&vF5F7);
  anatreeB->Branch("betaF3F7",&betaF3F7);
  anatreeB->Branch("betaF3F5",&betaF3F5);
  anatreeB->Branch("betaF5F7",&betaF5F7);
  anatreeB->Branch("gammaF3F7",&gammaF3F7);
  anatreeB->Branch("gammaF3F5",&gammaF3F5);
  anatreeB->Branch("gammaF5F7",&gammaF5F7);
  anatreeB->Branch("zet",&zet);
 
  anatreeB->Branch("F3X",&F3X);
  anatreeB->Branch("F3Y",&F3Y);
  anatreeB->Branch("F3A",&F3A);
  anatreeB->Branch("F3B",&F3B);
  anatreeB->Branch("F5X",&F5X);
  anatreeB->Branch("F5Y",&F5Y);
  anatreeB->Branch("F5A",&F5A);
  anatreeB->Branch("F5B",&F5B);
  anatreeB->Branch("F7X",&F7X);
  anatreeB->Branch("F7Y",&F7Y);
  anatreeB->Branch("F7A",&F7A);
  anatreeB->Branch("F7B",&F7B);
  anatreeB->Branch("deltaF3F5",&deltaF3F5);
  anatreeB->Branch("deltaF5F7",&deltaF5F7);
  anatreeB->Branch("brhoF3F5",&brhoF3F5);
  anatreeB->Branch("brhoF5F7",&brhoF5F7);
  anatreeB->Branch("aoqF3F5",&aoqF3F5);
  anatreeB->Branch("aoqF5F7",&aoqF5F7);
  anatreeB->Branch("recoF3A",&recoF3A);

  anatreeB->Branch("anaF3_Time",&anaF3_Time); 
  anatreeB->Branch("anaF5_Time",&anaF5_Time);
  anatreeB->Branch("anaF7_Time",&anaF7_Time);  
  
  anatreeB->Branch("BG_flag",&BG_flag);

  //=====Define BIT FLAGs for PPACs===============================================
  Int_t flag_1AX = 0x0001;
  Int_t flag_1BX = 0x0002;
  Int_t flag_2AX = 0x0004;
  Int_t flag_2BX = 0x0008;
  Int_t flag_1AY = 0x0010;
  Int_t flag_1BY = 0x0020;
  Int_t flag_2AY = 0x0040;
  Int_t flag_2BY = 0x0080;

  Int_t FLAG_F3 = 0x0000;
  Int_t FLAG_F5 = 0x0000;
  Int_t FLAG_F7 = 0x0000;
  
  infile->cd();
  //=====Begin LOOP======================================================
  int nEntry = caltreeB->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    caltreeB->GetEntry(iEntry);

    //===Initialization===
    BG_flag = 0;
    vF3F7 = TMath::Sqrt(-1);
    vF3F5 = TMath::Sqrt(-1);
    vF5F7 = TMath::Sqrt(-1);
    betaF3F7 = TMath::Sqrt(-1);
    betaF3F5 = TMath::Sqrt(-1);
    betaF5F7 = TMath::Sqrt(-1);
    gammaF3F7 = TMath::Sqrt(-1);
    gammaF3F5 = TMath::Sqrt(-1);
    gammaF5F7 = TMath::Sqrt(-1);

    raw_zet= TMath::Sqrt(-1);
    zet = TMath::Sqrt(-1);
    
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
    aoqF3F5 = TMath::Sqrt(-1);
    aoqF5F7 = TMath::Sqrt(-1);
    recoF3A = TMath::Sqrt(-1);

    anaF3_Time = TMath::Sqrt(-1);
    anaF5_Time = TMath::Sqrt(-1);
    anaF7_Time = TMath::Sqrt(-1);

    //===Calculation===
    vF3F7 = DistF3F7/(F7_Time - F3_Time + OffsetF3F7);
    vF3F5 = DistF3F5/(F5_Time - F3_Time + OffsetF3F5);
    vF5F7 = DistF5F7/(F7_Time - F5_Time + OffsetF5F7);    

    betaF3F7 = vF3F7/clight;
    betaF3F5 = vF3F5/clight;
    betaF5F7 = vF5F7/clight;
    gammaF3F7 = 1/TMath::Sqrt(1-betaF3F7*betaF3F7);
    gammaF3F5 = 1/TMath::Sqrt(1-betaF3F5*betaF3F5);
    gammaF5F7 = 1/TMath::Sqrt(1-betaF5F7*betaF5F7);

    raw_zet = vF3F7 * TMath::Sqrt(F7IC_E/(TMath::Log(2*m_e*vF3F7*vF3F7/Ionpair)-TMath::Log(1-betaF3F7*betaF3F7)-betaF3F7*betaF3F7));

    //zet = raw_zet;
    zet = zet_c1 * raw_zet + zet_c2;

    /*//F3
    if(cF31A_X_Tsum->IsInside(F31A_X_T1+F31A_X_T2,F31A_X)
       &&cF31B_X_Tsum->IsInside(F31B_X_T1+F31B_X_T2,F31B_X)){
      F31_X = (F31A_X + F31B_X)/2.;
    }else if(cF31A_X_Tsum->IsInside(F31A_X_T1+F31A_X_T2,F31A_X)
	     &&!cF31B_X_Tsum->IsInside(F31B_X_T1+F31B_X_T2,F31B_X)){
      F31_X = F31A_X;
    }else if(!cF31A_X_Tsum->IsInside(F31A_X_T1+F31A_X_T2,F31A_X)
	     &&cF31B_X_Tsum->IsInside(F31B_X_T1+F31B_X_T2,F31B_X)){
      F31_X = F31B_X;
    }
    if(cF31A_Y_Tsum->IsInside(F31A_Y_T1+F31A_Y_T2,F31A_Y)
       &&cF31B_Y_Tsum->IsInside(F31B_Y_T1+F31B_Y_T2,F31B_Y)){
      F31_Y = (F31A_Y + F31B_Y)/2.;
    }else if(cF31A_Y_Tsum->IsInside(F31A_Y_T1+F31A_Y_T2,F31A_Y)
	     &&!cF31B_Y_Tsum->IsInside(F31B_Y_T1+F31B_Y_T2,F31B_Y)){
      F31_Y = F31A_Y;
    }else if(!cF31A_Y_Tsum->IsInside(F31A_Y_T1+F31A_Y_T2,F31A_Y)
	     &&cF31B_Y_Tsum->IsInside(F31B_Y_T1+F31B_Y_T2,F31B_Y)){
      F31_Y = F31B_Y;
    }
    if(cF32A_X_Tsum->IsInside(F32A_X_T1+F32A_X_T2,F32A_X)
       &&cF32B_X_Tsum->IsInside(F32B_X_T1+F32B_X_T2,F32B_X)){
      F32_X = (F32A_X + F32B_X)/2.;
    }else if(cF32A_X_Tsum->IsInside(F32A_X_T1+F32A_X_T2,F32A_X)
	     &&!cF32B_X_Tsum->IsInside(F32B_X_T1+F32B_X_T2,F32B_X)){
      F32_X = F32A_X;
    }else if(!cF32A_X_Tsum->IsInside(F32A_X_T1+F32A_X_T2,F32A_X)
	     &&cF32B_X_Tsum->IsInside(F32B_X_T1+F32B_X_T2,F32B_X)){
      F32_X = F32B_X;
    }
    if(cF32A_Y_Tsum->IsInside(F32A_Y_T1+F32A_Y_T2,F32A_Y)
       &&cF32B_Y_Tsum->IsInside(F32B_Y_T1+F32B_Y_T2,F32B_Y)){
      F32_Y = (F32A_Y + F32B_Y)/2.;
    }else if(cF32A_Y_Tsum->IsInside(F32A_Y_T1+F32A_Y_T2,F32A_Y)
	     &&!cF32B_Y_Tsum->IsInside(F32B_Y_T1+F32B_Y_T2,F32B_Y)){
      F32_Y = F32A_Y;
    }else if(!cF32A_Y_Tsum->IsInside(F32A_Y_T1+F32A_Y_T2,F32A_Y)
	     &&cF32B_Y_Tsum->IsInside(F32B_Y_T1+F32B_Y_T2,F32B_Y)){
      F32_Y = F32B_Y;
      }*/

    //=====TEST of BIT FLAG for F3=======-
    if(cF31A_X_Tsum->IsInside(F31A_X_T1+F31A_X_T2,F31A_X)){
      FLAG_F3 = (FLAG_F3 | flag_1AX);
    }
    if(cF31B_X_Tsum->IsInside(F31B_X_T1+F31B_X_T2,F31B_X)){
      FLAG_F3 = (FLAG_F3 | flag_1BX);
    }
    if(cF32A_X_Tsum->IsInside(F32A_X_T1+F32A_X_T2,F32A_X)){
      FLAG_F3 = (FLAG_F3 | flag_2AX);
    }
    if(cF32B_X_Tsum->IsInside(F32B_X_T1+F32B_X_T2,F32B_X)){
      FLAG_F3 = (FLAG_F3 | flag_2BX);
    }
    if(cF31A_Y_Tsum->IsInside(F31A_Y_T1+F31A_Y_T2,F31A_Y)){
      FLAG_F3 = (FLAG_F3 | flag_1AY);
    }
    if(cF31B_Y_Tsum->IsInside(F31B_Y_T1+F31B_Y_T2,F31B_Y)){
      FLAG_F3 = (FLAG_F3 | flag_1BY);
    }
    if(cF32A_Y_Tsum->IsInside(F32A_Y_T1+F32A_Y_T2,F32A_Y)){
      FLAG_F3 = (FLAG_F3 | flag_2AY);
    }
    if(cF32B_Y_Tsum->IsInside(F32B_Y_T1+F32B_Y_T2,F32B_Y)){
      FLAG_F3 = (FLAG_F3 | flag_2BY);
    }
    
    

    
    //F5
    if(cF51A_X_Tsum->IsInside(F51A_X_T1+F51A_X_T2,F51A_X)
       &&cF51B_X_Tsum->IsInside(F51B_X_T1+F51B_X_T2,F51B_X)){
      F51_X = (F51A_X + F51B_X)/2.;
    }else if(cF51A_X_Tsum->IsInside(F51A_X_T1+F51A_X_T2,F51A_X)
	     &&!cF51B_X_Tsum->IsInside(F51B_X_T1+F51B_X_T2,F51B_X)){
      F51_X = F51A_X;
    }else if(!cF51A_X_Tsum->IsInside(F51A_X_T1+F51A_X_T2,F51A_X)
	     &&cF51B_X_Tsum->IsInside(F51B_X_T1+F51B_X_T2,F51B_X)){
      F51_X = F51B_X;
    }
    if(cF51A_Y_Tsum->IsInside(F51A_Y_T1+F51A_Y_T2,F51A_Y)
       &&cF51B_Y_Tsum->IsInside(F51B_Y_T1+F51B_Y_T2,F51B_Y)){
      F51_Y = (F51A_Y + F51B_Y)/2.;
    }else if(cF51A_Y_Tsum->IsInside(F51A_Y_T1+F51A_Y_T2,F51A_Y)
	     &&!cF51B_Y_Tsum->IsInside(F51B_Y_T1+F51B_Y_T2,F51B_Y)){
      F51_Y = F51A_Y;
    }else if(!cF51A_Y_Tsum->IsInside(F51A_Y_T1+F51A_Y_T2,F51A_Y)
	     &&cF51B_Y_Tsum->IsInside(F51B_Y_T1+F51B_Y_T2,F51B_Y)){
      F51_Y = F51B_Y;
    }
    if(cF52A_X_Tsum->IsInside(F52A_X_T1+F52A_X_T2,F52A_X)
       &&cF52B_X_Tsum->IsInside(F52B_X_T1+F52B_X_T2,F52B_X)){
      F52_X = (F52A_X + F52B_X)/2.;
    }else if(cF52A_X_Tsum->IsInside(F52A_X_T1+F52A_X_T2,F52A_X)
	     &&!cF52B_X_Tsum->IsInside(F52B_X_T1+F52B_X_T2,F52B_X)){
      F52_X = F52A_X;
    }else if(!cF52A_X_Tsum->IsInside(F52A_X_T1+F52A_X_T2,F52A_X)
	     &&cF52B_X_Tsum->IsInside(F52B_X_T1+F52B_X_T2,F52B_X)){
      F52_X = F52B_X;
    }
    if(cF52A_Y_Tsum->IsInside(F52A_Y_T1+F52A_Y_T2,F52A_Y)
       &&cF52B_Y_Tsum->IsInside(F52B_Y_T1+F52B_Y_T2,F52B_Y)){
      F52_Y = (F52A_Y + F52B_Y)/2.;
    }else if(cF52A_Y_Tsum->IsInside(F52A_Y_T1+F52A_Y_T2,F52A_Y)
	     &&!cF52B_Y_Tsum->IsInside(F52B_Y_T1+F52B_Y_T2,F52B_Y)){
      F52_Y = F52A_Y;
    }else if(!cF52A_Y_Tsum->IsInside(F52A_Y_T1+F52A_Y_T2,F52A_Y)
	     &&cF52B_Y_Tsum->IsInside(F52B_Y_T1+F52B_Y_T2,F52B_Y)){
      F52_Y = F52B_Y;
    }

    //F7
    if(cF71A_X_Tsum->IsInside(F71A_X_T1+F71A_X_T2,F71A_X)
       &&cF71B_X_Tsum->IsInside(F71B_X_T1+F71B_X_T2,F71B_X)){
      F71_X = (F71A_X + F71B_X)/2.;
    }else if(cF71A_X_Tsum->IsInside(F71A_X_T1+F71A_X_T2,F71A_X)
	     &&!cF71B_X_Tsum->IsInside(F71B_X_T1+F71B_X_T2,F71B_X)){
      F71_X = F71A_X;
    }else if(!cF71A_X_Tsum->IsInside(F71A_X_T1+F71A_X_T2,F71A_X)
	     &&cF71B_X_Tsum->IsInside(F71B_X_T1+F71B_X_T2,F71B_X)){
      F71_X = F71B_X;
    }
    if(cF71A_Y_Tsum->IsInside(F71A_Y_T1+F71A_Y_T2,F71A_Y)
       &&cF71B_Y_Tsum->IsInside(F71B_Y_T1+F71B_Y_T2,F71B_Y)){
      F71_Y = (F71A_Y + F71B_Y)/2.;
    }else if(cF71A_Y_Tsum->IsInside(F71A_Y_T1+F71A_Y_T2,F71A_Y)
	     &&!cF71B_Y_Tsum->IsInside(F71B_Y_T1+F71B_Y_T2,F71B_Y)){
      F71_Y = F71A_Y;
    }else if(!cF71A_Y_Tsum->IsInside(F71A_Y_T1+F71A_Y_T2,F71A_Y)
	     &&cF71B_Y_Tsum->IsInside(F71B_Y_T1+F71B_Y_T2,F71B_Y)){
      F71_Y = F71B_Y;
    }
    if(cF72A_X_Tsum->IsInside(F72A_X_T1+F72A_X_T2,F72A_X)
       &&cF72B_X_Tsum->IsInside(F72B_X_T1+F72B_X_T2,F72B_X)){
      F72_X = (F72A_X + F72B_X)/2.;
    }else if(cF72A_X_Tsum->IsInside(F72A_X_T1+F72A_X_T2,F72A_X)
	     &&!cF72B_X_Tsum->IsInside(F72B_X_T1+F72B_X_T2,F72B_X)){
      F72_X = F72A_X;
    }else if(!cF72A_X_Tsum->IsInside(F72A_X_T1+F72A_X_T2,F72A_X)
	     &&cF72B_X_Tsum->IsInside(F72B_X_T1+F72B_X_T2,F72B_X)){
      F72_X = F72B_X;
    }
    if(cF72A_Y_Tsum->IsInside(F72A_Y_T1+F72A_Y_T2,F72A_Y)
       &&cF72B_Y_Tsum->IsInside(F72B_Y_T1+F72B_Y_T2,F72B_Y)){
      F72_Y = (F72A_Y + F72B_Y)/2.;
    }else if(cF72A_Y_Tsum->IsInside(F72A_Y_T1+F72A_Y_T2,F72A_Y)
	     &&!cF72B_Y_Tsum->IsInside(F72B_Y_T1+F72B_Y_T2,F72B_Y)){
      F72_Y = F72A_Y;
    }else if(!cF72A_Y_Tsum->IsInside(F72A_Y_T1+F72A_Y_T2,F72A_Y)
	     &&cF72B_Y_Tsum->IsInside(F72B_Y_T1+F72B_Y_T2,F72B_Y)){
      F72_Y = F72B_Y;
    }

   
    if((FLAG_F3==(flag_1AX|flag_1BX|flag_2AX|flag_2BX))!=0){
      F3X = (F31A_X + F31B_X + F32A_X + F32B_X)/4.;
    }
    


    
    
    //F3X = (F31_X + F32_X)/2.;
    F3Y = (F31_Y + F32_Y)/2.;
    F3A = 1000.*(F31_X - F32_X)/TMath::Abs(F31_X - F32_X)*TMath::ACos(DistF3PPAC/TMath::Sqrt(DistF3PPAC*DistF3PPAC + (F31_X - F32_X)*(F31_X - F32_X)));
    F3B = 1000.*(F31_Y - F32_Y)/TMath::Abs(F31_Y - F32_Y)*TMath::ACos(DistF3PPAC/TMath::Sqrt(DistF3PPAC*DistF3PPAC + (F31_Y - F32_Y)*(F31_Y - F32_Y)));
    
    F5X = (F51_X + F52_X)/2.;
    F5Y = (F51_Y + F52_Y)/2.;
    F5A = 1000.*(F51_X - F52_X)/TMath::Abs(F51_X - F52_X)*TMath::ACos(DistF5PPAC/TMath::Sqrt(DistF5PPAC*DistF5PPAC + (F51_X - F52_X)*(F51_X - F52_X)));
    F5B = 1000.*(F51_Y - F52_Y)/TMath::Abs(F51_Y - F52_Y)*TMath::ACos(DistF5PPAC/TMath::Sqrt(DistF5PPAC*DistF5PPAC + (F51_Y - F52_Y)*(F51_Y - F52_Y)));
    
    F7X = (F71_X + F72_X)/2.;
    F7Y = (F71_Y + F72_Y)/2.;
    F7A = 1000.*(F71_X - F72_X)/TMath::Abs(F71_X - F72_X)*TMath::ACos(DistF7PPAC/TMath::Sqrt(DistF7PPAC*DistF7PPAC + (F71_X - F72_X)*(F71_X - F72_X)));
    F7B = 1000.*(F71_Y - F72_Y)/TMath::Abs(F71_Y - F72_Y)*TMath::ACos(DistF7PPAC/TMath::Sqrt(DistF7PPAC*DistF7PPAC + (F71_Y - F72_Y)*(F71_Y - F72_Y)));

    deltaF3F5 = (F5A - AXF3F5*F3X - AAF3F5*F3A)/ADF3F5;
    deltaF5F7 = (F5A - AXF5F7*F5X - AAF5F7*F5A)/ADF5F7;

    brhoF3F5 = Brho0F3F5*(1 + deltaF3F5);
    brhoF5F7 = Brho0F5F7*(1 + deltaF5F7);

    aoqF3F5 = brhoF3F5/betaF3F5/gammaF3F5*clight/m_u;
    aoqF5F7 = brhoF5F7/betaF5F7/gammaF5F7*clight/m_u;

    recoF3A = (ADF3F5*F5X - XDF3F5*F5A - (ADF3F5*XXF3F5 - XDF3F5*AXF3F5)*F3X)/(ADF3F5*XAF3F5 - XDF3F5*AAF3F5);
    
    anaF3_Time = F3_Time;
    anaF5_Time = F5_Time;
    anaF7_Time = F7_Time;
    
 
    if(!cF3pla->IsInside(F3_TR-F3_TL,log(F3_QL/F3_QR))
       ||!cF7pla->IsInside(F7_TR-F7_TL,log(F7_QL/F7_QR))
       /*&&(cF31A_X_Tsum->IsInside(F31A_X_T1+F31A_X_T2,F31A_X)
	||cF31B_X_Tsum->IsInside(F31B_X_T1+F31B_X_T2,F31B_X))
       &&(cF31A_Y_Tsum->IsInside(F31A_Y_T1+F31A_Y_T2,F31A_Y)
        ||cF31B_Y_Tsum->IsInside(F31B_Y_T1+F31B_Y_T2,F31B_Y))
       &&(cF32A_X_Tsum->IsInside(F32A_X_T1+F32A_X_T2,F32A_X)
        ||cF32B_X_Tsum->IsInside(F32B_X_T1+F32B_X_T2,F32B_X))
       &&(cF32A_Y_Tsum->IsInside(F32A_Y_T1+F32A_Y_T2,F32A_Y)
        ||cF32B_Y_Tsum->IsInside(F32B_Y_T1+F32B_Y_T2,F32B_Y))
       &&(cF71A_X_Tsum->IsInside(F71A_X_T1+F71A_X_T2,F71A_X)
        ||cF71B_X_Tsum->IsInside(F71B_X_T1+F71B_X_T2,F71B_X))
       &&(cF71A_Y_Tsum->IsInside(F71A_Y_T1+F71A_Y_T2,F71A_Y)
        ||cF71B_Y_Tsum->IsInside(F71B_Y_T1+F71B_Y_T2,F71B_Y))
       &&(cF72A_X_Tsum->IsInside(F72A_X_T1+F72A_X_T2,F72A_X)
        ||cF72B_X_Tsum->IsInside(F72B_X_T1+F72B_X_T2,F72B_X))
       &&(cF72A_Y_Tsum->IsInside(F72A_Y_T1+F72A_Y_T2,F72A_Y)
       ||cF72B_Y_Tsum->IsInside(F72B_Y_T1+F72B_Y_T2,F72B_Y))*/
       ){
      BG_flag = 1;
    }
    anatreeB->Fill();
  }
  anafile->cd();
  anatreeB->Write();
  anafile->Close();
}
