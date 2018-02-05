#include<stdio.h>
#include<iostream>
#include<math.h>

#include"TTree.h"

using namespace std;

void test(){

  //=====Load input file=================================================
  TFile *infile = TFile::Open("rootfiles/run0056/run0056_BEAM.root");

  //=====in tree variables===============================================
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;

  //===BEAM===
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
  
  double Dist = 44792.; //mm
  double offset = 00.; //nsec
  double ionpair = 4.866; //keV
  double m_e = 511.; //keV
  double clight = 299.792458; //mm/nsec
  double zet_c1 = 0.0398715; //slope
  double zet_c2 = -4.21233; //const.

  //=====Load CUT files==================================================
  TFile *cutfileF3pla = new TFile("cutfiles/cut_F3pla.root");
  TCutG *cF3pla = (TCutG*)cutfileF3pla->Get("F3pla");
  TFile *cutfileF7pla = new TFile("cutfiles/cut_F7pla.root");
  TCutG *cF7pla = (TCutG*)cutfileF7pla->Get("F7pla");
  

  //=====Create output file/tree=========================================
  TFile *outfile = new TFile("test.root","RECREATE");
  TTree *outtree = new TTree("outtree","outtree");

  //=====Define out variables============================================
  Double_t vF3F7;
  Double_t outF31A_X;
  Double_t outF31A_X_T1;
  Double_t outF31A_X_T2;
  Double_t betaF3F7;
  Double_t raw_zet;
  Double_t zet;
  Int_t BG_flag; //flag for background
  outtree->Branch("vF3F7",&vF3F7);
  outtree->Branch("outF31A_X",&outF31A_X);
  outtree->Branch("outF31A_X_T1",&outF31A_X_T1);
  outtree->Branch("outF31A_X_T2",&outF31A_X_T2);
  outtree->Branch("betaF3F7",&betaF3F7);
  outtree->Branch("zet",&zet);
  outtree->Branch("BG_flag",&BG_flag);
  infile->cd();
  
  int nEntry = caltreeB->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    caltreeB->GetEntry(iEntry);

    BG_flag = 1;
    vF3F7 = TMath::Sqrt(-1);
    outF31A_X = TMath::Sqrt(-1);
    outF31A_X_T1 = TMath::Sqrt(-1);
    outF31A_X_T2 = TMath::Sqrt(-1);
    zet = TMath::Sqrt(-1);
    
    vF3F7 = Dist/(F7_Time - F3_Time + offset);

    betaF3F7 = vF3F7/clight;

    raw_zet = vF3F7 * TMath::Sqrt(F7IC_E/(TMath::Log(2*m_e*vF3F7*vF3F7/ionpair)-TMath::Log(1-betaF3F7*betaF3F7)-betaF3F7*betaF3F7));

    zet = raw_zet;
    //zet = zet_c1 * raw_zet + zet_c2;
    
    outF31A_X = F31A_X;
    outF31A_X_T1 = F31A_X_T1;
    outF31A_X_T2 = F31A_X_T2;

    if(!cF3pla->IsInside(F3_TR-F3_TL,log(F3_QL/F3_QR))&&!cF7pla->IsInside(F7_TR-F7_TL,log(F7_QL/F7_QR))){
      BG_flag = -1;
    }
    outtree->Fill();
  }
  outfile->cd();
  outtree->Write();
  outfile->Close();
}
