#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>

#include"TROOT.h"
#include"TDirectory.h"
#include"TFile.h"
#include"TH1I.h"
#include"TTree.h"
#include"TCut.h"
#include"TCutG.h"
#include"TMath.h"

using namespace std;

//int main(int argc, char *argv[]){
void ana_dc_tdcdist(){
  //Int_t FileNum = TString(argv[1]).Atoi();
  
  //===== Load input file =====
  //TString infname = Form("/home/koiwai/analysis/rootfiles/all/run%04d_ALL.root",FileNum);
  //TFile *infile = TFile::Open(infname);
  //TFile *infile = TFile::Open("/home/koiwai/analysis/rootfiles/all/run0056_ALL.root");
  TFile *infile = TFile::Open("/home/koiwai/analysis/chain/calch.root");
  TTree *caltr;
  infile->GetObject("caltr",caltr);

  //===== input tree variables =====
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;

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
  caltr->SetBranchAddress("RunNumber",&RunNumber);
  caltr->SetBranchAddress("EventNumber",&EventNumber);

  caltr->SetBranchAddress("BDC1_TDC",BDC1_TDC);
  caltr->SetBranchAddress("BDC1_TrailTDC",BDC1_TrailTDC);
  caltr->SetBranchAddress("BDC1_WireID",BDC1_WireID);
  caltr->SetBranchAddress("BDC1_WirePosition",BDC1_WirePosition);
  caltr->SetBranchAddress("BDC1_WireZPosition",BDC1_WireZPosition);
  caltr->SetBranchAddress("BDC1_Layer",BDC1_Layer);
  caltr->SetBranchAddress("BDC1_PlaneID",BDC1_PlaneID);
  caltr->SetBranchAddress("BDC1_HitID",BDC1_HitID);
  
  caltr->SetBranchAddress("BDC2_TDC",BDC2_TDC);
  caltr->SetBranchAddress("BDC2_TrailTDC",BDC2_TrailTDC);
  caltr->SetBranchAddress("BDC2_WireID",BDC2_WireID);
  caltr->SetBranchAddress("BDC2_WirePosition",BDC2_WirePosition);
  caltr->SetBranchAddress("BDC2_WireZPosition",BDC2_WireZPosition);
  caltr->SetBranchAddress("BDC2_Layer",BDC2_Layer);
  caltr->SetBranchAddress("BDC2_PlaneID",BDC2_PlaneID);
  caltr->SetBranchAddress("BDC2_HitID",BDC2_HitID);
  
  caltr->SetBranchAddress("FDC1_TDC",FDC1_TDC);
  caltr->SetBranchAddress("FDC1_TrailTDC",FDC1_TrailTDC);
  caltr->SetBranchAddress("FDC1_WireID",FDC1_WireID);
  caltr->SetBranchAddress("FDC1_WirePosition",FDC1_WirePosition);
  caltr->SetBranchAddress("FDC1_WireZPosition",FDC1_WireZPosition);
  caltr->SetBranchAddress("FDC1_Layer",FDC1_Layer);
  caltr->SetBranchAddress("FDC1_PlaneID",FDC1_PlaneID);
  caltr->SetBranchAddress("FDC1_HitID",FDC1_HitID);
  
  caltr->SetBranchAddress("FDC2_TDC",FDC2_TDC);
  caltr->SetBranchAddress("FDC2_TrailTDC",FDC2_TrailTDC);
  caltr->SetBranchAddress("FDC2_WireID",FDC2_WireID);
  caltr->SetBranchAddress("FDC2_WirePosition",FDC2_WirePosition);
  caltr->SetBranchAddress("FDC2_WireZPosition",FDC2_WireZPosition);
  caltr->SetBranchAddress("FDC2_Layer",FDC2_Layer);
  caltr->SetBranchAddress("FDC2_PlaneID",FDC2_PlaneID);
  caltr->SetBranchAddress("FDC2_HitID",FDC2_HitID);
  
  //===== Create output file/tree =====
  TFile *anafile_dc_tdcdist = new TFile("/home/koiwai/analysis/rootfiles/tdc_dist/ana_dc_tdcdist.root","RECREATE");

  //===== Create TDC Distributions =====
  TH1I *hbdc1tdc0 = new TH1I("hbdc1tdc0","BDC1 TDC Layer 0",600,1400,2000);
  TH1I *hbdc1tdc1 = new TH1I("hbdc1tdc1","BDC1 TDC Layer 1",600,1400,2000);
  TH1I *hbdc1tdc2 = new TH1I("hbdc1tdc2","BDC1 TDC Layer 2",600,1400,2000);
  TH1I *hbdc1tdc3 = new TH1I("hbdc1tdc3","BDC1 TDC Layer 3",600,1400,2000);
  TH1I *hbdc1tdc4 = new TH1I("hbdc1tdc4","BDC1 TDC Layer 4",600,1400,2000);
  TH1I *hbdc1tdc5 = new TH1I("hbdc1tdc5","BDC1 TDC Layer 5",600,1400,2000);
  TH1I *hbdc1tdc6 = new TH1I("hbdc1tdc6","BDC1 TDC Layer 6",600,1400,2000);
  TH1I *hbdc1tdc7 = new TH1I("hbdc1tdc7","BDC1 TDC Layer 7",600,1400,2000);
  
  TH1I *hbdc2tdc0 = new TH1I("hbdc2tdc0","BDC2 TDC Layer 0",600,1400,2000);
  TH1I *hbdc2tdc1 = new TH1I("hbdc2tdc1","BDC2 TDC Layer 1",600,1400,2000);
  TH1I *hbdc2tdc2 = new TH1I("hbdc2tdc2","BDC2 TDC Layer 2",600,1400,2000);
  TH1I *hbdc2tdc3 = new TH1I("hbdc2tdc3","BDC2 TDC Layer 3",600,1400,2000);
  TH1I *hbdc2tdc4 = new TH1I("hbdc2tdc4","BDC2 TDC Layer 4",600,1400,2000);
  TH1I *hbdc2tdc5 = new TH1I("hbdc2tdc5","BDC2 TDC Layer 5",600,1400,2000);
  TH1I *hbdc2tdc6 = new TH1I("hbdc2tdc6","BDC2 TDC Layer 6",600,1400,2000);
  TH1I *hbdc2tdc7 = new TH1I("hbdc2tdc7","BDC2 TDC Layer 7",600,1400,2000);
  
  TH1I *hfdc1tdc0  = new TH1I("hfdc1tdc0", "FDC1 TDC Layer 0", 600,1400,2000);
  TH1I *hfdc1tdc1  = new TH1I("hfdc1tdc1", "FDC1 TDC Layer 1", 600,1400,2000);
  TH1I *hfdc1tdc2  = new TH1I("hfdc1tdc2", "FDC1 TDC Layer 2", 600,1400,2000);
  TH1I *hfdc1tdc3  = new TH1I("hfdc1tdc3", "FDC1 TDC Layer 3", 600,1400,2000);
  TH1I *hfdc1tdc4  = new TH1I("hfdc1tdc4", "FDC1 TDC Layer 4", 600,1400,2000);
  TH1I *hfdc1tdc5  = new TH1I("hfdc1tdc5", "FDC1 TDC Layer 5", 600,1400,2000);
  TH1I *hfdc1tdc6  = new TH1I("hfdc1tdc6", "FDC1 TDC Layer 6", 600,1400,2000);
  TH1I *hfdc1tdc7  = new TH1I("hfdc1tdc7", "FDC1 TDC Layer 7", 600,1400,2000);
  TH1I *hfdc1tdc8  = new TH1I("hfdc1tdc8", "FDC1 TDC Layer 8", 600,1400,2000);
  TH1I *hfdc1tdc9  = new TH1I("hfdc1tdc9", "FDC1 TDC Layer 9", 600,1400,2000);
  TH1I *hfdc1tdc10 = new TH1I("hfdc1tdc10","FDC1 TDC Layer 10",600,1400,2000);
  TH1I *hfdc1tdc11 = new TH1I("hfdc1tdc11","FDC1 TDC Layer 11",600,1400,2000);
  TH1I *hfdc1tdc12 = new TH1I("hfdc1tdc12","FDC1 TDC Layer 12",600,1400,2000);
  TH1I *hfdc1tdc13 = new TH1I("hfdc1tdc13","FDC1 TDC Layer 13",600,1400,2000);

  TH1I *hfdc2tdc0  = new TH1I("hfdc2tdc0", "FDC2 TDC Layer 0", 1000,1000,2000);
  TH1I *hfdc2tdc1  = new TH1I("hfdc2tdc1", "FDC2 TDC Layer 1", 1000,1000,2000);
  TH1I *hfdc2tdc2  = new TH1I("hfdc2tdc2", "FDC2 TDC Layer 2", 1000,1000,2000);
  TH1I *hfdc2tdc3  = new TH1I("hfdc2tdc3", "FDC2 TDC Layer 3", 1000,1000,2000);
  TH1I *hfdc2tdc4  = new TH1I("hfdc2tdc4", "FDC2 TDC Layer 4", 1000,1000,2000);
  TH1I *hfdc2tdc5  = new TH1I("hfdc2tdc5", "FDC2 TDC Layer 5", 1000,1000,2000);
  TH1I *hfdc2tdc6  = new TH1I("hfdc2tdc6", "FDC2 TDC Layer 6", 1000,1000,2000);
  TH1I *hfdc2tdc7  = new TH1I("hfdc2tdc7", "FDC2 TDC Layer 7", 1000,1000,2000);
  TH1I *hfdc2tdc8  = new TH1I("hfdc2tdc8", "FDC2 TDC Layer 8", 1000,1000,2000);
  TH1I *hfdc2tdc9  = new TH1I("hfdc2tdc9", "FDC2 TDC Layer 9", 1000,1000,2000);
  TH1I *hfdc2tdc10 = new TH1I("hfdc2tdc10","FDC2 TDC Layer 10",1000,1000,2000);
  TH1I *hfdc2tdc11 = new TH1I("hfdc2tdc11","FDC2 TDC Layer 11",1000,1000,2000);
  TH1I *hfdc2tdc12 = new TH1I("hfdc2tdc12","FDC2 TDC Layer 12",1000,1000,2000);
  TH1I *hfdc2tdc13 = new TH1I("hfdc2tdc13","FDC2 TDC Layer 13",1000,1000,2000);
  
  //===== Declear const.s =====
  Int_t BDCNumLayer = 8;
  Int_t FDCNumLayer = 14;

  //===== Begin LOOP =====
  int nEntry = caltr->GetEntries();
  cout << TMath::Ceil(nEntry/1000) << "k events to be treated." << endl;
  
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    caltr->GetEntry(iEntry);

    if(iEntry%100==0) clog << iEntry/1000 << "k events done..." << "\r";
    
    for(int i=0;i<16;++i){
      if(BDC1_TDC[0][i]>0) hbdc1tdc0 -> Fill(BDC1_TDC[0][i]);
      if(BDC1_TDC[1][i]>0) hbdc1tdc1 -> Fill(BDC1_TDC[1][i]);
      if(BDC1_TDC[2][i]>0) hbdc1tdc2 -> Fill(BDC1_TDC[2][i]);
      if(BDC1_TDC[3][i]>0) hbdc1tdc3 -> Fill(BDC1_TDC[3][i]);
      if(BDC1_TDC[4][i]>0) hbdc1tdc4 -> Fill(BDC1_TDC[4][i]);
      if(BDC1_TDC[5][i]>0) hbdc1tdc5 -> Fill(BDC1_TDC[5][i]);
      if(BDC1_TDC[6][i]>0) hbdc1tdc6 -> Fill(BDC1_TDC[6][i]);
      if(BDC1_TDC[7][i]>0) hbdc1tdc7 -> Fill(BDC1_TDC[7][i]);
      
      if(BDC2_TDC[0][i]>0) hbdc2tdc0 -> Fill(BDC2_TDC[0][i]);
      if(BDC2_TDC[1][i]>0) hbdc2tdc1 -> Fill(BDC2_TDC[1][i]);
      if(BDC2_TDC[2][i]>0) hbdc2tdc2 -> Fill(BDC2_TDC[2][i]);
      if(BDC2_TDC[3][i]>0) hbdc2tdc3 -> Fill(BDC2_TDC[3][i]);
      if(BDC2_TDC[4][i]>0) hbdc2tdc4 -> Fill(BDC2_TDC[4][i]);
      if(BDC2_TDC[5][i]>0) hbdc2tdc5 -> Fill(BDC2_TDC[5][i]);
      if(BDC2_TDC[6][i]>0) hbdc2tdc6 -> Fill(BDC2_TDC[6][i]);
      if(BDC2_TDC[7][i]>0) hbdc2tdc7 -> Fill(BDC2_TDC[7][i]);
    }
    for(int i=0;i<32;++i){
      if(FDC1_TDC[0][i]>0)  hfdc1tdc0 ->  Fill(FDC1_TDC[0][i]);
      if(FDC1_TDC[1][i]>0)  hfdc1tdc1 ->  Fill(FDC1_TDC[1][i]);
      if(FDC1_TDC[2][i]>0)  hfdc1tdc2 ->  Fill(FDC1_TDC[2][i]);
      if(FDC1_TDC[3][i]>0)  hfdc1tdc3 ->  Fill(FDC1_TDC[3][i]);
      if(FDC1_TDC[4][i]>0)  hfdc1tdc4 ->  Fill(FDC1_TDC[4][i]);
      if(FDC1_TDC[5][i]>0)  hfdc1tdc5 ->  Fill(FDC1_TDC[5][i]);
      if(FDC1_TDC[6][i]>0)  hfdc1tdc6 ->  Fill(FDC1_TDC[6][i]);
      if(FDC1_TDC[7][i]>0)  hfdc1tdc7 ->  Fill(FDC1_TDC[7][i]);
      if(FDC1_TDC[8][i]>0)  hfdc1tdc8 ->  Fill(FDC1_TDC[8][i]);
      if(FDC1_TDC[9][i]>0)  hfdc1tdc9 ->  Fill(FDC1_TDC[9][i]);
      if(FDC1_TDC[10][i]>0) hfdc1tdc10 -> Fill(FDC1_TDC[10][i]);
      if(FDC1_TDC[11][i]>0) hfdc1tdc11 -> Fill(FDC1_TDC[11][i]);
      if(FDC1_TDC[12][i]>0) hfdc1tdc12 -> Fill(FDC1_TDC[12][i]);
      if(FDC1_TDC[13][i]>0) hfdc1tdc13 -> Fill(FDC1_TDC[13][i]);
    }
    for(int i=0;i<112;++i){
      if(FDC2_TDC[0][i]>0)  hfdc2tdc0 ->  Fill(FDC2_TDC[0][i]);
      if(FDC2_TDC[1][i]>0)  hfdc2tdc1 ->  Fill(FDC2_TDC[1][i]);
      if(FDC2_TDC[2][i]>0)  hfdc2tdc2 ->  Fill(FDC2_TDC[2][i]);
      if(FDC2_TDC[3][i]>0)  hfdc2tdc3 ->  Fill(FDC2_TDC[3][i]);
      if(FDC2_TDC[4][i]>0)  hfdc2tdc4 ->  Fill(FDC2_TDC[4][i]);
      if(FDC2_TDC[5][i]>0)  hfdc2tdc5 ->  Fill(FDC2_TDC[5][i]);
      if(FDC2_TDC[6][i]>0)  hfdc2tdc6 ->  Fill(FDC2_TDC[6][i]);
      if(FDC2_TDC[7][i]>0)  hfdc2tdc7 ->  Fill(FDC2_TDC[7][i]);
      if(FDC2_TDC[8][i]>0)  hfdc2tdc8 ->  Fill(FDC2_TDC[8][i]);
      if(FDC2_TDC[9][i]>0)  hfdc2tdc9 ->  Fill(FDC2_TDC[9][i]);
      if(FDC2_TDC[10][i]>0) hfdc2tdc10 -> Fill(FDC2_TDC[10][i]);
      if(FDC2_TDC[11][i]>0) hfdc2tdc11 -> Fill(FDC2_TDC[11][i]);
      if(FDC2_TDC[12][i]>0) hfdc2tdc12 -> Fill(FDC2_TDC[12][i]);
      if(FDC2_TDC[13][i]>0) hfdc2tdc13 -> Fill(FDC2_TDC[13][i]);
    }
    
  }//for LOOP
  anafile_dc_tdcdist->cd();
  //anatreeDC->Write();
  hbdc1tdc0->Write();
  hbdc1tdc1->Write();
  hbdc1tdc2->Write();
  hbdc1tdc3->Write();
  hbdc1tdc4->Write();
  hbdc1tdc5->Write();
  hbdc1tdc6->Write();
  hbdc1tdc7->Write();
  
  hbdc2tdc0->Write();
  hbdc2tdc1->Write();
  hbdc2tdc2->Write();
  hbdc2tdc3->Write();
  hbdc2tdc4->Write();
  hbdc2tdc5->Write();
  hbdc2tdc6->Write();
  hbdc2tdc7->Write();

  hfdc1tdc0->Write();
  hfdc1tdc1->Write();
  hfdc1tdc2->Write();
  hfdc1tdc3->Write();
  hfdc1tdc4->Write();
  hfdc1tdc5->Write();
  hfdc1tdc6->Write();
  hfdc1tdc7->Write();
  hfdc1tdc8->Write();
  hfdc1tdc9->Write();
  hfdc1tdc10->Write();
  hfdc1tdc11->Write();
  hfdc1tdc12->Write();
  hfdc1tdc13->Write();

  hfdc2tdc0->Write();
  hfdc2tdc1->Write();
  hfdc2tdc2->Write();
  hfdc2tdc3->Write();
  hfdc2tdc4->Write();
  hfdc2tdc5->Write();
  hfdc2tdc6->Write();
  hfdc2tdc7->Write();
  hfdc2tdc8->Write();
  hfdc2tdc9->Write();
  hfdc2tdc10->Write();
  hfdc2tdc11->Write();
  hfdc2tdc12->Write();
  hfdc2tdc13->Write();

  anafile_dc_tdcdist->Close();
}//ana_dc()
