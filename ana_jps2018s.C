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

using namespace std;
using namespace TMath;

int main(int argc, char *argv[]){

  Int_t FileNumber = TString(argv[1]).Atoi();

  //B===== Load input file =================================================
  TString FileName_beam = Form("/home/koiwai/analysis/anailes/beam/ana_beam%04d.root",FileNumber);
  TFile *infile_beam = TFile::Open(FileName);

  TTree *anatrB;
  infile->GetObject("anatrB",anatrB);

  //B===== beam variables =====
  Double_t ENum, RNum;
  Double_t zetBR, aoqBR;
  Double_t betaF5F7;
  Int_t BG_flagB;

  //B===== SetBranchAddress =====
  anatrB->SetBranchAddress("ENum",&ENum);
  anatrB->SetBranchAddress("RNum",&RNum);
  anatrB->SetBranchAddress("zetBR",&zetBR);
  anatrB->SetBranchAddress("aoqBR",&aoqBR);
  anatrB->SetBranchAddress("betaF5F7",&betaF5F7);
  anatrB->SetBranchAddress("BG_flagB",&BG_flag);

  //D===== Load input file =====
  TString FileName_dali = Form("/home/koiwai/analysis/rootfiles/dali/run%04d_DALI.root",FileNumber);
  TFile *infile_dali = TFile::Open(FileName);

  TTree *caltrD;
  infile->GetObject("caltrD",caltrD);
  
  //D===== dali variables =====
  Double_t DALI_Energy, DALI_CosTheta, DALI_Time, DALI_ID, DALI_Layer;
  Double_t DALI_X, DALI_Y, DALI_Z, DALI_Pos;
  Int_t DALI_Multi;

  //D===== SetBranchAddress =====
  caltrD->SetBranchAddress("DALI_Energy",&DALI_Energy);
  caltrD->SetBranchAddress("DALI_CosTheta",&CosTheta);
  caltrD->SetBranchAddress("DALI_Time",&DALI_Time);
  caltrD->SetBranchAddress("DALI_ID",&DALI_ID);
  caltrD->SetBranchAddress("DALI_Layer",&DALI_Layer);
  caltrD->SetBranchAddress("DALI_X",&DALI_X);
  caltrD->SetBranchAddress("DALI_Y",&DALI_Y);
  caltrD->SetBranchAddress("DALI_Z",&DALI_Z);
  caltrD->SetBranchAddress("DALI_Pos",&DALI_Pos);
  caltrD->SetBranchAddress("DALI_Multi",&DALI_Multi);

  //A===== Load input file =====
  TString FileName_all = Form("/home/koiwai/analysis/rootfiles/run%04d/run%04d.root",FileNumber,FileNumber);
  TFile *infile_all = TFile::Open(FileName);

  TTree *caltr;
  infile->GetObject("caltr",caltr);

  //A===== input variables =====
  Int_t Hodo_ID; // ID of the hodoscope with the highest charge
  Int_t Hodo_Multiplicity;
  Double_t Hodo_QCal, Hodo_TCal;
  Double_t Hodoi_QCal[24], Hodoi_TCal[24];
  Double_t SBT1_Charge, SBT2_Charge, SBT1_Time, SBT2_Time, SBT1_TimeDiff, SBT2_TimeDiff;

  //A===== SetBranchAddress =====
  caltr->SetBranchAddress("Hodo_Multiplicity",&Hodo_Multiplicity);
  caltr->SetBranchAddress("Hodo_ID",&Hodo_ID);
  caltr->SetBranchAddress("Hodo_QCal",&Hodo_QCal);
  caltr->SetBranchAddress("Hodo_TCal",&Hodo_TCal);
  caltr->SetBranchAddress("SBT1_Charge",&SBT1_Charge);
  caltr->SetBranchAddress("SBT1_Time",&SBT1_Time);
  caltr->SetBranchAddress("SBT1_TimeDiff",&SBT1_TimeDiff);
  caltr->SetBranchAddress("SBT2_Charge",&SBT2_Charge);
  caltr->SetBranchAddress("SBT2_Time",&SBT2_Time);
  caltr->SetBranchAddress("SBT2_TimeDiff",&SBT2_TimeDiff);
  for (Int_t i=0;i<24;i++)
    {
      caltr->SetBranchAddress(Form("Hodo%d_QCal",i+1),&Hodoi_QCal[i]);
      caltr->SetBranchAddress(Form("Hodo%d_TCal",i+1),&Hodoi_TCal[i]);
    }

  //===== Load CUT file =====
  TFile *cutfile = TFile::Open("/home/koiwai/analysis/cutfiles/cutBR56Ca.root");
  TCutG *cBR56Ca = (TCutG*)cutfile->Get("CUTG");
  
  //===== Create output file/tree =====
  TString AnaFileName = Form("/home/koiwai/analysis/anafiles/jps/ana_jps%04d.root",FuileNumber);
  TFile *ana_file = new TFile(AnaFileName,"recreate");
  TTree * anatrJ = new TTree("anatrJ","anatrJ");

  //===== Dealare const.s =====
  

  //===== Declare ana variables =====
  Int_t EventNumber, RunNumber;

  Int_t BR56Ca;
  //===== Set Branches ======
  anatrJ->Branch("EventNumber",&EventNumber);
  anatrJ->Branch("RunNumber",&RunNumber);

  anatrJ->Branch("BR56Ca",&BR56Ca);
  //===== LOOP =====
  int nEntry = caltr->GetEntries();
  for(Int_t iEntry=0;iEntry<nEntry;++iEntry){
    //for(Int_t iEntry;iEmtry<100,++iEntry){

    if(iEntry%100==0) clog << iEntry/1000 << "k events treated..." << "\r";

    caltr->GetEntry(iEntry);
    
    RunNumber = RNum;
    EventNumber = ENum;

    //=== Initialization ===
    BR56Ca = 0;


    //=== Calc. ===


    //=== CUT ===
    if(cutfile->IsInside(aoqBR,zetBR)) BR56Ca = 1;

    anatrJ->Fill();
  }//LOOP
  anafile->cd();
  anatrJ->Write();
  anafile->Close();
}//main
