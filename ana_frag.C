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

void ana_frag(){

  //===== Load iuput file =====
  TFile *infile_hodo = TFile::Open("rootfiles/run0056/run0056_HODO.root");
  TTree *caltreeH;
  infile_hodo->GetObject("caltreeH",caltreeH);

  TFile *infile_dc   = TFile::Open("rootfiles/run0056/run0056_DC.root");
  TTree *caltreeDC;
  infile_dc->GetObject("caltreeDC",caltreeDC);

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
  
  //===== SetBranchAddress =====
  caltreeH->SetBranchAddress("RunNumber",&RunNumber);
  caltreeH->SetBranchAddress("EventNumber",&EventNumber);

  //=== Hodo ===
   for (Int_t i=0;i<24;i++)
      {
	caltreeH->Branch(Form("Hodo%d_QCal",i+1),&Hodoi_QCal[i]);
	caltreeH->Branch(Form("Hodo%d_TCal",i+1),&Hodoi_TCal[i]);
	caltreeH->Branch(Form("Hodo%d_QRaw",i+1),&Hodoi_QRaw[i]);
	caltreeH->Branch(Form("Hodo%d_TRaw",i+1),&Hodoi_TRaw[i]);
	caltreeH->Branch(Form("Hodo%d_TURaw",i+1),&Hodoi_TURaw[i]);
	caltreeH->Branch(Form("Hodo%d_TDRaw",i+1),&Hodoi_TDRaw[i]);
	caltreeH->Branch(Form("Hodo%d_QURaw",i+1),&Hodoi_QURaw[i]);
	caltreeH->Branch(Form("Hodo%d_QDRaw",i+1),&Hodoi_QDRaw[i]);
	caltreeH->Branch(Form("Hodo%d_TUCal",i+1),&Hodoi_TUCal[i]);
	caltreeH->Branch(Form("Hodo%d_TDCal",i+1),&Hodoi_TDCal[i]);
	caltreeH->Branch(Form("Hodo%d_QUCal",i+1),&Hodoi_QUCal[i]);
	caltreeH->Branch(Form("Hodo%d_QDCal",i+1),&Hodoi_QDCal[i]);
      }
    caltreeH->Branch("Hodo_Multiplicity",&Hodo_Multiplicity);
    caltreeH->Branch("Hodo_ID",&Hodo_ID);
    caltreeH->Branch("Hodo_QCal",&Hodo_QCal);
    caltreeH->Branch("Hodo_TCal",&Hodo_TCal);
    caltreeH->Branch("Hodo_QRaw",&Hodo_QRaw);
    caltreeH->Branch("Hodo_TRaw",&Hodo_TRaw);

    caltreeH->Branch("SBT1_Charge",&SBT1_Charge);
    caltreeH->Branch("SBT1_Time",&SBT1_Time);
    caltreeH->Branch("SBT1_TimeDiff",&SBT1_TimeDiff);
    caltreeH->Branch("SBT2_Charge",&SBT2_Charge);
    caltreeH->Branch("SBT2_Time",&SBT2_Time);
    caltreeH->Branch("SBT2_TimeDiff",&SBT2_TimeDiff);
    caltreeH->Branch("SBT1_QL",&SBT1_QL);
    caltreeH->Branch("SBT1_QR",&SBT1_QR);
    caltreeH->Branch("SBT1_TL",&SBT1_TL);
    caltreeH->Branch("SBT1_TR",&SBT1_TR);
    caltreeH->Branch("SBT2_QL",&SBT2_QL);
    caltreeH->Branch("SBT2_QR",&SBT2_QR);
    caltreeH->Branch("SBT2_TL",&SBT2_TL);
    caltreeH->Branch("SBT2_TR",&SBT2_TR);
  //=== DC ===
  
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

  //===== Define ana variables =====

  //===== Create anatree Branch =====

  //===== Begin LOOP =====
  int nEntry = caltreeH->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    caltreeH->GetEntry(iEntry);
    caltreeDC->GetEntry(iEntry);

    //=== Initialization ===


    //=== Calculation ===


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
