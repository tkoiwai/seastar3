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
using namespace TMath;

void ana_frag(){

  //===== Load iuput file =====
  TFile *infile_hodo = TFile::Open("rootfiles/run0056/run0056_HODO.root");
  TTree *caltreeH;
  infile_hodo->GetObject("caltreeH",caltreeH);
  /*
  TFile *infile_dc   = TFile::Open("rootfiles/run0056/run0056_DC.root");
  TTree *caltreeDC;
  infile_dc->GetObject("caltreeDC",caltreeDC);
  */
  TFile *infile_beam = TFile::Open("rootfiles/ana_beam.root");
  TTree *anatreeB;
  infile_beam->GetObject("anatreeB",anatreeB);
  
  //===== input tree variables =====
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;
  
  //=== Hodo (includes SBT) ===
  Double_t AllHodo_Charge[24];
  Double_t AllHodo_Time[24];
  Int_t Hodo_ID; // ID of the hodoscope with the highest charge
  Int_t Hodo_Multiplicity;
  Double_t Hodo_QCal, Hodo_QRaw; // Highest charge //  Double_t Hodo_Charge;
  Double_t Hodo_TCal, Hodo_TRaw; // Time of the hodoscope with the highest charge  Double_t Hodo_Time;
  Int_t Hodoi_TURaw[24], Hodoi_TDRaw[24], Hodoi_QURaw[24], Hodoi_QDRaw[24];
  Double_t Hodoi_TUCal[24], Hodoi_TDCal[24], Hodoi_QUCal[24], Hodoi_QDCal[24];
  Double_t Hodoi_TCal[24], Hodoi_TRaw[24], Hodoi_QCal[24], Hodoi_QRaw[24];

  Double_t SBT1_Charge, SBT2_Charge, SBT1_Time, SBT2_Time, SBT1_TimeDiff, SBT2_TimeDiff;
  Double_t SBT1_QL, SBT1_QR, SBT1_TL, SBT1_TR;
  Double_t SBT2_QL, SBT2_QR, SBT2_TL, SBT2_TR;
  
  //=== DC ===
  /*
  Int_t BDC1_TDC[8][16], BDC1_TrailTDC[8][16];
  Int_t BDC1_WireID[8][16];
  Double_t BDC1_WirePosition[8][16], BDC1_WireZPosition[8][16];
  Int_t BDC1_Layer[8][16],BDC1_PlaneID[8][16], BDC1_HitID[8][16];

  Int_t BDC2_TDC[8][16], BDC2_TrailTDC[8][16];
  Int_t BDC2_WireID[8][16];
  Double_t BDC2_WirePosition[8][16], BDC2_WireZPosition[8][16];
  Int_t BDC2_Layer[8][16],BDC2_PlaneID[8][16], BDC2_HitID[8][16];

  Int_t FDC1_TDC[14][32], FDC1_TrailTDC[14][32];
  Int_t FDC1_WireID[14][32];
  Double_t FDC1_WirePosition[14][32], FDC1_WireZPosition[14][32];
  Int_t FDC1_Layer[14][32],FDC1_PlaneID[14][32], FDC1_HitID[14][32];

  Int_t FDC2_TDC[14][112], FDC2_TrailTDC[14][112];
  Int_t FDC2_WireID[14][112];
  Double_t FDC2_WirePosition[14][112], FDC2_WireZPosition[14][112];
  Int_t FDC2_Layer[14][112],FDC2_PlaneID[14][112], FDC2_HitID[14][112];
 
  Int_t NumBDC1Hit, NumBDC2Hit, NumFDC1Hit, NumFDC2Hit;
  */

  //=== BEAM ===
  Double_t anaF3_Time, anaF5_Time, anaF7_Time, vF5F7;
  
  //===== SetBranchAddress =====
  caltreeH->SetBranchAddress("RunNumber",&RunNumber);
  caltreeH->SetBranchAddress("EventNumber",&EventNumber);

  //=== Hodo ===
  for (Int_t i=0;i<24;i++){
    caltreeH->SetBranchAddress(Form("Hodo%d_QCal",i+1),&Hodoi_QCal[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_TCal",i+1),&Hodoi_TCal[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_QRaw",i+1),&Hodoi_QRaw[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_TRaw",i+1),&Hodoi_TRaw[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_TURaw",i+1),&Hodoi_TURaw[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_TDRaw",i+1),&Hodoi_TDRaw[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_QURaw",i+1),&Hodoi_QURaw[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_QDRaw",i+1),&Hodoi_QDRaw[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_TUCal",i+1),&Hodoi_TUCal[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_TDCal",i+1),&Hodoi_TDCal[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_QUCal",i+1),&Hodoi_QUCal[i]);
    caltreeH->SetBranchAddress(Form("Hodo%d_QDCal",i+1),&Hodoi_QDCal[i]);
  }
  caltreeH->SetBranchAddress("Hodo_Multiplicity",&Hodo_Multiplicity);
  caltreeH->SetBranchAddress("Hodo_ID",&Hodo_ID);
  caltreeH->SetBranchAddress("Hodo_QCal",&Hodo_QCal);
  caltreeH->SetBranchAddress("Hodo_TCal",&Hodo_TCal);
  caltreeH->SetBranchAddress("Hodo_QRaw",&Hodo_QRaw);
  caltreeH->SetBranchAddress("Hodo_TRaw",&Hodo_TRaw);
  
  caltreeH->SetBranchAddress("SBT1_Charge",&SBT1_Charge);
  caltreeH->SetBranchAddress("SBT1_Time",&SBT1_Time);
  caltreeH->SetBranchAddress("SBT1_TimeDiff",&SBT1_TimeDiff);
  caltreeH->SetBranchAddress("SBT2_Charge",&SBT2_Charge);
  caltreeH->SetBranchAddress("SBT2_Time",&SBT2_Time);
  caltreeH->SetBranchAddress("SBT2_TimeDiff",&SBT2_TimeDiff);
  caltreeH->SetBranchAddress("SBT1_QL",&SBT1_QL);
  caltreeH->SetBranchAddress("SBT1_QR",&SBT1_QR);
  caltreeH->SetBranchAddress("SBT1_TL",&SBT1_TL);
  caltreeH->SetBranchAddress("SBT1_TR",&SBT1_TR);
  caltreeH->SetBranchAddress("SBT2_QL",&SBT2_QL);
  caltreeH->SetBranchAddress("SBT2_QR",&SBT2_QR);
  caltreeH->SetBranchAddress("SBT2_TL",&SBT2_TL);
  caltreeH->SetBranchAddress("SBT2_TR",&SBT2_TR);
  //=== DC ===
  //=== BEAM ===
  anatreeB->SetBranchAddress("anaF3_Time",&anaF3_Time);
  anatreeB->SetBranchAddress("anaF5_Time",&anaF5_Time);
  anatreeB->SetBranchAddress("anaF7_Time",&anaF7_Time);
  anatreeB->SetBranchAddress("vF5F7",&vF5F7);
  
  //===== Load CUT files =====
  //=== SBT (ln(QL/QR):TR-TL) ===
  TFile *cutfileSBT1 = new TFile("cutfiles/cut_SBT1.root");
  TCutG *cSBT1 = (TCutG*)cutfileSBT1->Get("CUTG");
  cSBT1->SetName("cutSBT1");
  
  TFile *cutfileSBT2 = new TFile("cutfiles/cut_SBT2.root");
  TCutG *cSBT2 = (TCutG*)cutfileSBT2->Get("CUTG");
  cSBT2->SetName("cutSBT2");
  
  //===== Create output file/tree =====
  TFile *anafile_frag = new TFile("rootfiles/ana_frag.root","RECREATE");
  TTree *anatreeF     = new TTree("anatreeF","anatreeF");

  //===== Declear const.s =====
  Double_t Dist_F7F13 = 36617.;
  Double_t Dist_SBTTarget = 2737.67;
  Double_t Offset_F7F13 = 570.;
  Double_t Offset_TargetHOD = 230.;
  
  //===== Define ana variables =====
  Int_t RunNum, EventNum;
  
  Double_t vF7F13;
  Double_t F13_Time;
  Double_t tofF7F13;
  Double_t tofSBTTarget;
  Double_t tofTargetHOD;
  Double_t Hodo_Time;

  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatreeF->Branch("RunNum",&RunNum);
  anatreeF->Branch("EventNum",&EventNum);

  anatreeF->Branch("vF7F13",&vF7F13);
  anatreeF->Branch("F13_Time",&F13_Time);
  anatreeF->Branch("tofF7F13",&tofF7F13);
  anatreeF->Branch("tofSBTTarget",&tofSBTTarget);
  anatreeF->Branch("tofTargetHOD",&tofTargetHOD);
  anatreeF->Branch("Hodo_Time",&Hodo_Time);
  
  anatreeF->Branch("BG_flag",&BG_flag);
  
  //===== Begin LOOP =====
  int nEntry = caltreeH->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){

    
    if(iEntry%100 == 0){
      clog<< iEntry/1000 << "k events treated..." << "\r";
    }

    caltreeH->GetEntry(iEntry);
    //caltreeDC->GetEntry(iEntry);
    anatreeB->GetEntry(iEntry);
    
    RunNum = RunNumber;
    EventNum = EventNumber;
    
    //=== Initialization ===
    vF7F13 = Sqrt(-1);
    F13_Time = Sqrt(-1);
    tofF7F13 = Sqrt(-1);
    tofSBTTarget = Sqrt(-1);
    tofTargetHOD = Sqrt(-1);
    Hodo_Time = Sqrt(-1);
    
    //=== Calculation ===
    F13_Time = (SBT1_Time + SBT2_Time)/2.;
    tofF7F13 = F13_Time - anaF7_Time + Offset_F7F13;
    vF7F13 = Dist_F7F13/tofF7F13;
    tofSBTTarget = Dist_SBTTarget/vF7F13;
    tofTargetHOD = Hodo_TCal - F13_Time - tofSBTTarget  + Offset_TargetHOD;
    Hodo_Time = Hodo_TCal;
    
    //=== Cut by graphical cut ===
    if(  !cSBT1->IsInside(SBT1_TR-SBT1_TL,log(SBT1_QL/SBT1_QR))
       ||!cSBT2->IsInside(SBT2_TR-SBT2_TL,log(SBT2_QL/SBT2_QR))
	 ){
      BG_flag = 1;
    }//graphical cut
    
    anatreeF->Fill();
  }//for(int iEntry=0;iEntry<nEntry;++iEntry i.e. LOOP
  anafile_frag->cd();
  anatreeF->Write();
  anafile_frag->Close();
}//main i.e. ana_frag()
