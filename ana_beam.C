#include<stdio.h>
#include<iostream>
#include<fstream>
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

  //===== Load input file =================================================
  TFile *infile = TFile::Open("rootfiles/run0056/run0056_BEAM.root");
  TTree *caltreeB;
  infile->GetObject("caltreeB",caltreeB);

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
  
  //=== Set value ===
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


  
  //===== Load CUT files ==================================================
  //=== Plastic (graphical cut)===
  TFile *cutfileF3pla = new TFile("cutfiles/cut_F3pla.root");
  TCutG *cF3pla = (TCutG*)cutfileF3pla->Get("F3pla");
  TFile *cutfileF7pla = new TFile("cutfiles/cut_F7pla.root");
  TCutG *cF7pla = (TCutG*)cutfileF7pla->Get("F7pla");

  //=== Chrage change @ F5 (graphical cut) ===
  TFile *cutfileF5Qchange = new TFile("cutfiles/cut_bigripsZvsF5Qchange.root");
  TCutG *cF5Qchange = (TCutG*)cutfileF5Qchange->Get("CUTG");

  //=== PPAC Tsum gate ===
  ifstream fin;
  fin.open("cutfiles/cut_PPAC_Tsum.dat");
  if(fin.fail()){
    cout << "Error: file is not found." << endl;
    return;
  }
  
  Int_t cPPAC_Tsum_low[24], cPPAC_Tsum_up[24];
  
  for(Int_t cPPAC_index = 0;cPPAC_index<24;++cPPAC_index){
    fin >> cPPAC_Tsum_low[cPPAC_index] >> cPPAC_Tsum_up[cPPAC_index];
  }
  
  //===== Create output file/tree =========================================
  TFile *anafile = new TFile("rootfiles/ana_beam.root","RECREATE");
  TTree *anatreeB = new TTree("anatreeB","anatreeB");

  //===== Declear const.s =================================================
  double DistF3F7 = 46568.; //[mm]
  double DistF3F5 = 23284.;
  double DistF5F7 = 23284.;
  double OffsetF3F7 = 292.379; //[nsec]
  double OffsetF3F5 = 159.572;
  //double OffsetF3F5 = 160.572;
  double OffsetF5F7 = 132.807;
  double Ionpair = 4.866; //[keV]
  double m_e = 511.; //[keV]
  double m_u = 931.49432; //[MeV]
  double clight = 299.792458; //[mm/nsec]
  double zetBR_c1 = 0.0374179; //slope
  double zetBR_c2 = -4.2354; //const.
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
  double XDF5F7 = -34.4457; //[mm/%]
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
  //=== for Z ===
  Double_t tofF3F7;
  Double_t vF3F7, vF3F5, vF5F7;
  Double_t betaF3F7, betaF3F5, betaF5F7;
  Double_t gammaF3F7, gammaF3F5, gammaF5F7;
  Double_t raw_zetBR;
  Double_t zetBR;

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
  Double_t aoqF3F5, aoqF5F7;
  Double_t recoF3A, recodeltaF3F5_X, recodeltaF3F5_A;
  Double_t aoqBR;

  //=== Imported from intree ===
  Double_t anaF3_Time, anaF5_Time, anaF7_Time;
  
  Int_t BG_flag; //flag for background

  //======
  anatreeB->Branch("tofF3F7",&tofF3F7);
  anatreeB->Branch("vF3F7",&vF3F7);
  anatreeB->Branch("vF3F5",&vF3F5);
  anatreeB->Branch("vF5F7",&vF5F7);
  anatreeB->Branch("betaF3F7",&betaF3F7);
  anatreeB->Branch("betaF3F5",&betaF3F5);
  anatreeB->Branch("betaF5F7",&betaF5F7);
  anatreeB->Branch("gammaF3F7",&gammaF3F7);
  anatreeB->Branch("gammaF3F5",&gammaF3F5);
  anatreeB->Branch("gammaF5F7",&gammaF5F7);
  anatreeB->Branch("zetBR",&zetBR);
 
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
  anatreeB->Branch("recodeltaF3F5_A",&recodeltaF3F5_A);
  anatreeB->Branch("recodeltaF3F5_X",&recodeltaF3F5_X);
  anatreeB->Branch("deltaF3F5_A",&deltaF3F5_A);
  anatreeB->Branch("deltaF3F5_X",&deltaF3F5_X);
  anatreeB->Branch("aoqBR",&aoqBR);
  
  anatreeB->Branch("anaF3_Time",&anaF3_Time); 
  anatreeB->Branch("anaF5_Time",&anaF5_Time);
  anatreeB->Branch("anaF7_Time",&anaF7_Time);

  anatreeB->Branch("BG_flag",&BG_flag);
  
  infile->cd();

  //===== Begin LOOP ======================================================
  int nEntry = caltreeB->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    caltreeB->GetEntry(iEntry);

    //=== Initialization ===
    BG_flag = 0;
    tofF3F7 = TMath::Sqrt(-1);
    vF3F7 = TMath::Sqrt(-1);
    vF3F5 = TMath::Sqrt(-1);
    vF5F7 = TMath::Sqrt(-1);
    betaF3F7 = TMath::Sqrt(-1);
    betaF3F5 = TMath::Sqrt(-1);
    betaF5F7 = TMath::Sqrt(-1);
    gammaF3F7 = TMath::Sqrt(-1);
    gammaF3F5 = TMath::Sqrt(-1);
    gammaF5F7 = TMath::Sqrt(-1);

    raw_zetBR= TMath::Sqrt(-1);
    zetBR = TMath::Sqrt(-1);
    
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
    vF3F7 = DistF3F7/(F7_Time - F3_Time + OffsetF3F7);
    vF3F5 = DistF3F5/(F5_Time - F3_Time + OffsetF3F5);
    vF5F7 = DistF5F7/(F7_Time - F5_Time + OffsetF5F7);    

    betaF3F7 = vF3F7/clight;
    betaF3F5 = vF3F5/clight;
    betaF5F7 = vF5F7/clight;
    gammaF3F7 = 1/TMath::Sqrt(1-betaF3F7*betaF3F7);
    gammaF3F5 = 1/TMath::Sqrt(1-betaF3F5*betaF3F5);
    gammaF5F7 = 1/TMath::Sqrt(1-betaF5F7*betaF5F7);

    raw_zetBR = vF3F7 * TMath::Sqrt(F7IC_E/(TMath::Log(2*m_e*vF3F7*vF3F7/Ionpair)-TMath::Log(1-betaF3F7*betaF3F7)-betaF3F7*betaF3F7));

    //zetBR = raw_zetBR;
    zetBR = zetBR_c1 * raw_zetBR + zetBR_c2;
  
  
   
    //=== F3 Tsum gate ===
    if((cPPAC_Tsum_low[0]<F31A_X_T1+F31A_X_T2&&F31A_X_T1+F31A_X_T2<cPPAC_Tsum_up[0])&&(cPPAC_Tsum_low[1]<F31B_X_T1+F31B_X_T2&&F31B_X_T1+F31B_X_T2<cPPAC_Tsum_up[1])){
      F31_X = (F31A_X + F31B_X)/2.;
    }else if((cPPAC_Tsum_low[0]<F31A_X_T1+F31A_X_T2&&F31A_X_T1+F31A_X_T2<cPPAC_Tsum_up[0])&&!(cPPAC_Tsum_low[1]<F31B_X_T1+F31B_X_T2&&F31B_X_T1+F31B_X_T2<cPPAC_Tsum_up[1])){
      F31_X = F31A_X;
    }else if(!(cPPAC_Tsum_low[0]<F31A_X_T1+F31A_X_T2&&F31A_X_T1+F31A_X_T2<cPPAC_Tsum_up[0])&&(cPPAC_Tsum_low[1]<F31B_X_T1+F31B_X_T2&&F31B_X_T1+F31B_X_T2<cPPAC_Tsum_up[1])){
      F31_X = F31B_X;
    }else{
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[2]<F32A_X_T1+F32A_X_T2&&F32A_X_T1+F32A_X_T2<cPPAC_Tsum_up[2])&&(cPPAC_Tsum_low[3]<F32B_X_T1+F32B_X_T2&&F32B_X_T1+F32B_X_T2<cPPAC_Tsum_up[3])){
      F32_X = (F32A_X + F32B_X)/2.;
    }else if((cPPAC_Tsum_low[2]<F32A_X_T1+F32A_X_T2&&F32A_X_T1+F32A_X_T2<cPPAC_Tsum_up[2])&&!(cPPAC_Tsum_low[3]<F32B_X_T1+F32B_X_T2&&F32B_X_T1+F32B_X_T2<cPPAC_Tsum_up[3])){
      F32_X = F32A_X;
    }else if(!(cPPAC_Tsum_low[2]<F32A_X_T1+F32A_X_T2&&F32A_X_T1+F32A_X_T2<cPPAC_Tsum_up[2])&&(cPPAC_Tsum_low[3]<F32B_X_T1+F32B_X_T2&&F32B_X_T1+F32B_X_T2<cPPAC_Tsum_up[3])){
      F32_X = F32B_X;
    }else{
      BG_flag = 2;
    }
     if((cPPAC_Tsum_low[4]<F31A_Y_T1+F31A_Y_T2&&F31A_Y_T1+F31A_Y_T2<cPPAC_Tsum_up[4])&&(cPPAC_Tsum_low[5]<F31B_Y_T1+F31B_Y_T2&&F31B_Y_T1+F31B_Y_T2<cPPAC_Tsum_up[5])){
      F31_Y = (F31A_Y + F31B_Y)/2.;
    }else if((cPPAC_Tsum_low[4]<F31A_Y_T1+F31A_Y_T2&&F31A_Y_T1+F31A_Y_T2<cPPAC_Tsum_up[4])&&!(cPPAC_Tsum_low[5]<F31B_Y_T1+F31B_Y_T2&&F31B_Y_T1+F31B_Y_T2<cPPAC_Tsum_up[5])){
      F31_Y = F31A_Y;
    }else if(!(cPPAC_Tsum_low[4]<F31A_Y_T1+F31A_Y_T2&&F31A_Y_T1+F31A_Y_T2<cPPAC_Tsum_up[4])&&(cPPAC_Tsum_low[5]<F31B_Y_T1+F31B_Y_T2&&F31B_Y_T1+F31B_Y_T2<cPPAC_Tsum_up[5])){
      F31_Y = F31B_Y;
    }else{
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[6]<F32A_Y_T1+F32A_Y_T2&&F32A_Y_T1+F32A_Y_T2<cPPAC_Tsum_up[6])&&(cPPAC_Tsum_low[7]<F32B_Y_T1+F32B_Y_T2&&F32B_Y_T1+F32B_Y_T2<cPPAC_Tsum_up[7])){
      F32_Y = (F32A_Y + F32B_Y)/2.;
    }else if((cPPAC_Tsum_low[6]<F32A_Y_T1+F32A_Y_T2&&F32A_Y_T1+F32A_Y_T2<cPPAC_Tsum_up[6])&&!(cPPAC_Tsum_low[7]<F32B_Y_T1+F32B_Y_T2&&F32B_Y_T1+F32B_Y_T2<cPPAC_Tsum_up[7])){
      F32_Y = F32A_Y;
    }else if(!(cPPAC_Tsum_low[6]<F32A_Y_T1+F32A_Y_T2&&F32A_Y_T1+F32A_Y_T2<cPPAC_Tsum_up[6])&&(cPPAC_Tsum_low[7]<F32B_Y_T1+F32B_Y_T2&&F32B_Y_T1+F32B_Y_T2<cPPAC_Tsum_up[7])){
      F32_Y = F32B_Y;
    }else{
      BG_flag = 2;
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
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[14]<F52A_Y_T1+F52A_Y_T2&&F52A_Y_T1+F52A_Y_T2<cPPAC_Tsum_up[14])&&(cPPAC_Tsum_low[15]<F52B_Y_T1+F52B_Y_T2&&F52B_Y_T1+F52B_Y_T2<cPPAC_Tsum_up[15])){
      F52_Y = (F52A_Y + F52B_Y)/2.;
    }else if((cPPAC_Tsum_low[14]<F52A_Y_T1+F52A_Y_T2&&F52A_Y_T1+F52A_Y_T2<cPPAC_Tsum_up[14])&&!(cPPAC_Tsum_low[15]<F52B_Y_T1+F52B_Y_T2&&F52B_Y_T1+F52B_Y_T2<cPPAC_Tsum_up[15])){
      F52_Y = F52A_Y;
    }else if(!(cPPAC_Tsum_low[14]<F52A_Y_T1+F52A_Y_T2&&F52A_Y_T1+F52A_Y_T2<cPPAC_Tsum_up[14])&&(cPPAC_Tsum_low[15]<F52B_Y_T1+F52B_Y_T2&&F52B_Y_T1+F52B_Y_T2<cPPAC_Tsum_up[15])){
      F52_Y = F52B_Y;
    }else{
      BG_flag = 2;
    }
    //=== F7 ===
     if((cPPAC_Tsum_low[16]<F71A_X_T1+F71A_X_T2&&F71A_X_T1+F71A_X_T2<cPPAC_Tsum_up[16])&&(cPPAC_Tsum_low[17]<F71B_X_T1+F71B_X_T2&&F71B_X_T1+F71B_X_T2<cPPAC_Tsum_up[17])){
      F71_X = (F71A_X + F71B_X)/2.;
    }else if((cPPAC_Tsum_low[16]<F71A_X_T1+F71A_X_T2&&F71A_X_T1+F71A_X_T2<cPPAC_Tsum_up[16])&&!(cPPAC_Tsum_low[17]<F71B_X_T1+F71B_X_T2&&F71B_X_T1+F71B_X_T2<cPPAC_Tsum_up[17])){
      F71_X = F71A_X;
    }else if(!(cPPAC_Tsum_low[16]<F71A_X_T1+F71A_X_T2&&F71A_X_T1+F71A_X_T2<cPPAC_Tsum_up[16])&&(cPPAC_Tsum_low[17]<F71B_X_T1+F71B_X_T2&&F71B_X_T1+F71B_X_T2<cPPAC_Tsum_up[17])){
      F71_X = F71B_X;
    }else{
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[18]<F72A_X_T1+F72A_X_T2&&F72A_X_T1+F72A_X_T2<cPPAC_Tsum_up[18])&&(cPPAC_Tsum_low[19]<F72B_X_T1+F72B_X_T2&&F72B_X_T1+F72B_X_T2<cPPAC_Tsum_up[19])){
      F72_X = (F72A_X + F72B_X)/2.;
    }else if((cPPAC_Tsum_low[18]<F72A_X_T1+F72A_X_T2&&F72A_X_T1+F72A_X_T2<cPPAC_Tsum_up[18])&&!(cPPAC_Tsum_low[19]<F72B_X_T1+F72B_X_T2&&F72B_X_T1+F72B_X_T2<cPPAC_Tsum_up[19])){
      F72_X = F72A_X;
    }else if(!(cPPAC_Tsum_low[18]<F72A_X_T1+F72A_X_T2&&F72A_X_T1+F72A_X_T2<cPPAC_Tsum_up[18])&&(cPPAC_Tsum_low[19]<F72B_X_T1+F72B_X_T2&&F72B_X_T1+F72B_X_T2<cPPAC_Tsum_up[19])){
      F72_X = F72B_X;
    }else{
      BG_flag = 2;
    }
     if((cPPAC_Tsum_low[20]<F71A_Y_T1+F71A_Y_T2&&F71A_Y_T1+F71A_Y_T2<cPPAC_Tsum_up[20])&&(cPPAC_Tsum_low[21]<F71B_Y_T1+F71B_Y_T2&&F71B_Y_T1+F71B_Y_T2<cPPAC_Tsum_up[21])){
      F71_Y = (F71A_Y + F71B_Y)/2.;
    }else if((cPPAC_Tsum_low[20]<F71A_Y_T1+F71A_Y_T2&&F71A_Y_T1+F71A_Y_T2<cPPAC_Tsum_up[20])&&!(cPPAC_Tsum_low[21]<F71B_Y_T1+F71B_Y_T2&&F71B_Y_T1+F71B_Y_T2<cPPAC_Tsum_up[21])){
      F71_Y = F71A_Y;
    }else if(!(cPPAC_Tsum_low[20]<F71A_Y_T1+F71A_Y_T2&&F71A_Y_T1+F71A_Y_T2<cPPAC_Tsum_up[20])&&(cPPAC_Tsum_low[21]<F71B_Y_T1+F71B_Y_T2&&F71B_Y_T1+F71B_Y_T2<cPPAC_Tsum_up[21])){
      F71_Y = F71B_Y;
    }else{
      BG_flag = 2;
    }
    if((cPPAC_Tsum_low[22]<F72A_Y_T1+F72A_Y_T2&&F72A_Y_T1+F72A_Y_T2<cPPAC_Tsum_up[22])&&(cPPAC_Tsum_low[23]<F72B_Y_T1+F72B_Y_T2&&F72B_Y_T1+F72B_Y_T2<cPPAC_Tsum_up[23])){
      F72_Y = (F72A_Y + F72B_Y)/2.;
    }else if((cPPAC_Tsum_low[22]<F72A_Y_T1+F72A_Y_T2&&F72A_Y_T1+F72A_Y_T2<cPPAC_Tsum_up[22])&&!(cPPAC_Tsum_low[23]<F72B_Y_T1+F72B_Y_T2&&F72B_Y_T1+F72B_Y_T2<cPPAC_Tsum_up[23])){
      F72_Y = F72A_Y;
    }else if(!(cPPAC_Tsum_low[22]<F72A_Y_T1+F72A_Y_T2&&F72A_Y_T1+F72A_Y_T2<cPPAC_Tsum_up[22])&&(cPPAC_Tsum_low[23]<F72B_Y_T1+F72B_Y_T2&&F72B_Y_T1+F72B_Y_T2<cPPAC_Tsum_up[23])){
      F72_Y = F72B_Y;
    }else{
      BG_flag = 2;
    }
    
    F3X = (F31_X + F32_X)/2.;
    F3Y = (F31_Y + F32_Y)/2.;
    F3A = 1000.*TMath::ATan((F31_X-F32_X)/DistF3PPAC);
    F3B = 1000.*TMath::ATan((F31_Y-F32_Y)/DistF3PPAC);
   
    F5X = (F51_X + F52_X)/2.;
    F5Y = (F51_Y + F52_Y)/2.;
    F5A = 1000.*TMath::ATan((F51_X-F52_X)/DistF5PPAC);
    F5B = 1000.*TMath::ATan((F51_Y-F52_Y)/DistF5PPAC);
    
    F7X = (F71_X + F72_X)/2.;
    F7Y = (F71_Y + F72_Y)/2.;
    F7A = 1000.*TMath::ATan((F71_X-F72_X)/DistF7PPAC);
    F7B = 1000.*TMath::ATan((F71_Y-F72_Y)/DistF7PPAC);
    
    deltaF3F5 = (F5X - XXF3F5*F3X - XAF3F5*F3A)/XDF3F5; //[%] 
    deltaF5F7 = (F7X - XXF5F7*F5X - XAF5F7*F5A)/XDF5F7;

    brhoF3F5 = Brho0F3F5*(1 + deltaF3F5*0.01);
    brhoF5F7 = Brho0F5F7*(1 + deltaF5F7*0.01);

    aoqF3F5 = brhoF3F5*clight/m_u/betaF3F5/gammaF3F5 + 0.0017*F3X + 0.00005*F5X + 0.000002*F5X*F5X;
    aoqF5F7 = brhoF5F7*clight/m_u/betaF5F7/gammaF5F7 + 0.0006*F5A + 0.00015*F7A;

    aoqBR = aoqF3F5;

    recoF3A = (ADF3F5*F5X - XDF3F5*F5A - (ADF3F5*XXF3F5 - XDF3F5*AXF3F5)*F3X)/(ADF3F5*XAF3F5 - XDF3F5*AAF3F5);
    
    anaF3_Time = F3_Time;
    anaF5_Time = F5_Time;
    anaF7_Time = F7_Time;    

    //===== Cut by graphical cut ==========================================================
    if(!cF3pla->IsInside(F3_TR-F3_TL,log(F3_QL/F3_QR))
       ||!cF7pla->IsInside(F7_TR-F7_TL,log(F7_QL/F7_QR))
       /*||!cF5Qchange->IsInside(aoqF3F5/aoqF5F7,zetBR)*/
       ){
      BG_flag = 1;
    }
    anatreeB->Fill();
  }
  anafile->cd();
  anatreeB->Write();
  anafile->Close();
}
