#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<bitset>

#include"TROOT.h"
#include"TDirectory.h"
#include"TFile.h"
#include"TH1I.h"
#include"TTree.h"
#include"TCut.h"
#include"TCutG.h"
#include"TMath.h"
#include"TString.h"

using namespace std;
using namaspace TMath;

int main(int argc, char *atgv[]){
  
  Int_t FileNum = TString(atgv[1]).Atoi();

  //===== Load input file =====
  TFile *infile = TFile::Open(Form("/home/koiwai/analysis/strage2018feb20/run%04_ALL.root".FileNum));

  TTree *caltr;
  infile->GetObject("caltr",caltr);

  //===== Input tree variables =====
  Long64_t EventNumber;
  Int_t RunNumber;

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

  //===== SetBranchAddress of input tree valiables =====
  caltr->SetBranchAddress("Hodo_Multiplicity",&Hodo_Multiplicity);
  caltr->SetBranchAddress("Hodo_IDcal",&Hodo_ID);
  caltr->SetBranchAddress("Hodo_QCal",&Hodo_QCal);
  caltr->SetBranchAddress("Hodo_TCal",&Hodo_TCal);

  for(Int_t i=0;i<24;++i){
    caltr->SetBranchAddress(Form("Hodo%02d_QCal",i+1),&Hodoi_QCal[i]);
    caltr->SetBranchAddress(Form("Hodo%02d_TCal",i+1),&Hodoi_TCal[i]);
  }

  caltr->SetBranchAddress("SBT1_Charge",&SBT1_Charge);
  caltr->SetBranchAddress("SBT2_Charge",&SBT2_Charge);
  caltr->SetBranchAddress("SBT1_Time",&SBT1_Time);
  caltr->SetBranchAddress("SBT2_Time",&SBT2_Time);
  caltr->SetBranchAddress("SBT1_TimeDiff",&SBT1_TimeDiff);
  caltr->SetBranchAddress("SBT2_TimeDiff",&SBT2_TimeDiff);

  //===== Create output file =====
  TFile *anafile = new TFile("/home/koiwai/analysis/anafiles/ana_hodo.root","recreate");
  TTree *anatrH = new TTree("anatrH","anatrH");

  //===== Declare const.s =====

  //===== Declare anatree const.s =====

  //===== Declare anatree variables =====
  Double_t Hodo_Q, Hodo_T, Hodo_ID;

  Double_t SBT_Q, SBT_T;
  
  //===== Create anatree Branches =====
  anatrH->Branch("Hodo_Q",&Hodo_Q);
  anatrH->Branch("Hodo_T",&Hodo_T);
  anatrH->Branch("Hodo_ID",&Hodo_ID);
  
  anatrH->Branch("SBT_Q",&SBE_Q);
  anatrH->Branch("SBT_T,&SBT_T");
  
  //===== Begin LOOP =====
  int nEntry = caltr->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
  //for(int iEntry=0;iEntry<nEntry;++iEntry){

    if(iEntry%100==0) clog<<iEntry/1000<<"k enents treated..."<<"\r";

    caltr->GetEntry(iEntry);
    
    RunNum = RunNumber;
    EventNum = EventNumber;

    //=== Initialzation ===
    Hodo_Q = Sqrt(-1);
    Hodo_T = Sqrt(-1);
    Hodo_ID = -9999;

    SBT_Q = Sqrt(-1);
    SBT_T = Sqrt(-1);
    
    //=== Calclation ===
    
    
    
    anatrH->Fill();
  }//LOOP
  anafile->cd();
  anatrH->Write();
  anafile->Close();
}//main
