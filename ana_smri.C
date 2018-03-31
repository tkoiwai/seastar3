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

#include"/home/koiwai/analysis/brho_func/Brho_A56Z20_br56Ca_sa56Ca.C"
#include"/home/koiwai/analysis/brho_func/Len_A56Z20_br56Ca_sa56Ca.C"

using namespace std;

int main(int argc, char *argv[]){

  Int_t FileNum = TString(argv[1]).Atoi();
  
  //===== Load input file =====
  TString infname = Form("/home/koiwai/analysis/rootfiles/all/run%04d_ALL.root",FileNum);
  TFile *infile = TFile::Open(infname);
  
  TTree *caltr;
  infile->GetObject("caltr",caltr);

  //===== input tree variables =====
  Long64_t EventNumber_all;
  Int_t RunNumber_all;

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

  //===== SetBranchAddress =====
  caltr->SetBranchAddress("RunNumber",&RunNumber_all);
  caltr->SetBranchAddress("EventNumber",&EventNumber_all);

  for (Int_t i=0;i<24;i++)
      {
	caltr->SetBranchAddress(Form("Hodo%d_QCal",i+1),&Hodoi_QCal[i]);
	caltr->SetBranchAddress(Form("Hodo%d_TCal",i+1),&Hodoi_TCal[i]);
	//caltr->SetBranchAddress(Form("Hodo%d_Charge",i+1),&AllHodo_Charge[i]);
	//caltr->SetBranchAddress(Form("Hodo%d_Time",i+1),&AllHodo_Time[i]);
	caltr->SetBranchAddress(Form("Hodo%d_QRaw",i+1),&Hodoi_QRaw[i]);
	caltr->SetBranchAddress(Form("Hodo%d_TRaw",i+1),&Hodoi_TRaw[i]);
	caltr->SetBranchAddress(Form("Hodo%d_TURaw",i+1),&Hodoi_TURaw[i]);
	caltr->SetBranchAddress(Form("Hodo%d_TDRaw",i+1),&Hodoi_TDRaw[i]);
	caltr->SetBranchAddress(Form("Hodo%d_QURaw",i+1),&Hodoi_QURaw[i]);
	caltr->SetBranchAddress(Form("Hodo%d_QDRaw",i+1),&Hodoi_QDRaw[i]);
	caltr->SetBranchAddress(Form("Hodo%d_TUCal",i+1),&Hodoi_TUCal[i]);
	caltr->SetBranchAddress(Form("Hodo%d_TDCal",i+1),&Hodoi_TDCal[i]);
	caltr->SetBranchAddress(Form("Hodo%d_QUCal",i+1),&Hodoi_QUCal[i]);
	caltr->SetBranchAddress(Form("Hodo%d_QDCal",i+1),&Hodoi_QDCal[i]);
      }
  caltr->SetBranchAddress("Hodo_Multiplicity",&Hodo_Multiplicity);
  caltr->SetBranchAddress("Hodo_ID",&Hodo_ID);
  caltr->SetBranchAddress("Hodo_QCal",&Hodo_QCal);
  caltr->SetBranchAddress("Hodo_TCal",&Hodo_TCal);
  caltr->SetBranchAddress("Hodo_QRaw",&Hodo_QRaw);
  caltr->SetBranchAddress("Hodo_TRaw",&Hodo_TRaw);
  
  caltr->SetBranchAddress("SBT1_Charge",&SBT1_Charge);
  caltr->SetBranchAddress("SBT1_Time",&SBT1_Time);
  caltr->SetBranchAddress("SBT1_TimeDiff",&SBT1_TimeDiff);
  caltr->SetBranchAddress("SBT2_Charge",&SBT2_Charge);
  caltr->SetBranchAddress("SBT2_Time",&SBT2_Time);
  caltr->SetBranchAddress("SBT2_TimeDiff",&SBT2_TimeDiff);
  caltr->SetBranchAddress("SBT1_QL",&SBT1_QL);
  caltr->SetBranchAddress("SBT1_QR",&SBT1_QR);
  caltr->SetBranchAddress("SBT1_TL",&SBT1_TL);
  caltr->SetBranchAddress("SBT1_TR",&SBT1_TR);
  caltr->SetBranchAddress("SBT2_QL",&SBT2_QL);
  caltr->SetBranchAddress("SBT2_QR",&SBT2_QR);
  caltr->SetBranchAddress("SBT2_TL",&SBT2_TL);
  caltr->SetBranchAddress("SBT2_TR",&SBT2_TR);

  //===== Load input DC file =====
  TString infnameDC = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/ana_mwdc%04d.root",FileNum);
  TFile *infileDC = TFile::Open(infnameDC);
  
  TTree *anatrDC;
  infileDC->GetObject("anatrDC",anatrDC);

  //===== Declare DC variables =====
  Int_t EventNumber_dc, RunNumber_dc;

  Double_t BDC1_X, BDC1_Y;
  Double_t BDC2_X, BDC2_Y;
  Double_t BDC_X, BDC_Y, BDC_A, BDC_B;
  Double_t Target_X, Target_Y;

  Double_t FDC1_X, FDC1_Y, FDC1_A, FDC1_B;
  Double_t FDC2_X, FDC2_Y, FDC2_A, FDC2_B;
  
  Int_t BG_flag_dc;

  //===== DC SetBranchAddress =====
  anatrDC->SetBranchAddress("EventNum",&EventNumber_dc);
  anatrDC->SetBranchAddress("RunNum",&RunNumber_dc);

  anatrDC->SetBranchAddress("BDC1_X",&BDC1_X);
  anatrDC->SetBranchAddress("BDC1_Y",&BDC1_Y);
  anatrDC->SetBranchAddress("BDC2_X",&BDC2_X);
  anatrDC->SetBranchAddress("BDC2_Y",&BDC2_Y);
  anatrDC->SetBranchAddress("BDC_X",&BDC_X);
  anatrDC->SetBranchAddress("BDC_Y",&BDC_Y);
  anatrDC->SetBranchAddress("BDC_A",&BDC_A);
  anatrDC->SetBranchAddress("BDC_B",&BDC_B);
  anatrDC->SetBranchAddress("Target_X",&Target_X);
  anatrDC->SetBranchAddress("Target_Y",&Target_Y);
  anatrDC->SetBranchAddress("FDC1_X",&FDC1_X);
  anatrDC->SetBranchAddress("FDC1_Y",&FDC1_Y);
  anatrDC->SetBranchAddress("FDC1_A",&FDC1_A);
  anatrDC->SetBranchAddress("FDC1_B",&FDC1_B);
  anatrDC->SetBranchAddress("FDC2_X",&FDC2_X);
  anatrDC->SetBranchAddress("FDC2_Y",&FDC2_Y);
  anatrDC->SetBranchAddress("FDC2_A",&FDC2_A);
  anatrDC->SetBranchAddress("FDC2_B",&FDC2_B);
  anatrDC->SetBranchAddress("BG_flag",&BG_flag_dc);

  //===== AddFriend =====
  caltr->AddFriend(anatrDC);

  //===== Load .dat files =====
  TEnv *env = new TEnv("/home/koiwai/analysis/db/geometry_psp.dat");
  TEnv *env_hodot = new TEnv("/home/koiwai/analysis/db/hodo_toff.dat");
  TEnv *env_hodoq = new TEnv("/home/koiwai/analysis/db/hodo_qcorr.dat");
  
  //===== Create output file/tree =====
  TString ofname = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNum);
  TFile *anafile_smri = new TFile(ofname,"RECREATE");
  TTree *anatrS  = new TTree("anatrS","anatrS");
  
  /*
  //===== Create TDC Distributions =====
  Int_t BDCNumLayer = 8;
  Int_t FDCNumLayer = 14;
  Int_t DCLayerNum[4] = {8,8,14,14};
  Int_t DCNum = 4;
  Int_t tdcwindow[2] = {1000,2000};
  Double_t tdc2mm[DCNum][14][2000] = {};
  Int_t tdcint[DCNum][14] = {};
  Double_t wire_gap[DCNum] = {2.5, 2.5, 5, 10};
  char DCName[DCNum] = {'b','b','f','f'};
  Int_t DCid[DCNum] = {1,2,1,2};
  
  TFile *RootFile = new TFile("/home/koiwai/analysis/rootfiles/tdc_dist/ana_dc_tdcdist.root","READ");
  if(RootFile){
    gROOT->cd();
    for(Int_t n=0;n<DCNum;++n){
      //TH1I *h = NULL;    
      for(Int_t l=0;l<DCLayerNum[n];++l){
	TH1I *h = NULL;
	h = (TH1I*)RootFile->Get(Form("h%cdc%dtdc%d",DCName[n],DCid[n],l));
	
	if(h){
	  tdcint[n][l] = (Double_t)h->Integral(h->FindBin(tdcwindow[0]),h->FindBin(tdcwindow[1]));
	  for(Int_t i=0;i<2000;++i){
	    if(i>=0&&i<=1000) tdc2mm[n][l][i] = 0;
	    else tdc2mm[n][l][i] = wire_gap[n]*(Double_t)h->Integral(h->FindBin(i),h->FindBin(tdcwindow[1]))/tdcint[n][l];
	  }	  	
	}
      }
    }
  }
  */
  
  //===== Declare const.s =====
  Int_t Dist_BDC1BDC2 = env->GetValue("Dist_BDC1BDC2",998.); //mm
  Int_t Dist_BDC1Target = env->GetValue("Dist_BDC1Target",2231.);
  Int_t Dist_BDC1FDC1 = env->GetValue("Dist_BDC1FDC1",3533.);
  Int_t Dist_SBTTarget = env->GetValue("Dist_SBT_Target",2795.);
  Int_t Width_BDC1 = env->GetValue("BDC1_Width",68.);
  Int_t Width_FDC1 = env->GetValue("FDC1_Width",180.);
  
  //===== Declare anatree const.s =====
  
  //===== Declare valables for calc. =====
  //===== Declare anatree variables =====
  Int_t RunNum, EventNum;
  
 
  Double_t brhoSAMURAI;

  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatrS->Branch("RunNum",&RunNum);
  anatrS->Branch("EventNum",&EventNum);
  anatrS->Branch("brhoSAMURAI",&brhoSAMURAI);

  anatrS->Branch("BG_flag",&BG_flag);
  
  //===== Begin LOOP =====
  int nEntry = caltr->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    //for(int iEntry=0;iEntry<5;++iEntry){
  
    
    if(iEntry%100 == 0){
      clog<< iEntry/1000 << "k events treated..." << "\r";
    }

    RunNum = RunNumber;
    EventNum = EventNumber;
    
    caltr->GetEntry(iEntry);
 
    //@@@ Brho Length Function @@@

    Double_t x[6];

    x[0] = FDC1_X;
    x[1] = FDC1_A;
    x[2] = FDC1_Y;
    x[3] = FDC1_B;
    x[4] = FDC2_X;
    x[5] = FDC2_A;

    brhoSAMURAI = MDF_Brho_A56Z20(x);
    
    //@@@ HODO @@@









    //@@@ HODO end @@@
    
    anatrS->Fill();
  }//for LOOP
  anafile_smri->cd();
  anatrS->Write();
  anafile_smri->Close();
}//main()
