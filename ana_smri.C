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
  Long64_t EventNumber;
  Int_t RunNumber;

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
  
  //===== Create output file/tree =====
  TString ofname = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNum);
  TFile *anafile_smri = new TFile(ofname,"RECREATE");
  TTree *anatrS  = new TTree("anatrS","anatrS");
  
  
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
  
  TFile *RootFile = new TFile("/home/koiwai/analysis/rootfiles/tdc_dist/ana_dc_tdcdist.root","READ");
  if(RootFile){
    gROOT->cd();
    for(Int_t n=0;n<DCNum;++n){
      TH1I *h = NULL;    
      for(Int_t l=0;l<DCLayerNum[n];++l){
	h = (TH1I*)RootFile->Get(Form("h%cdc%dtdc%d",DCName[n],n+1,l));
	
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
  
  //===== Declare const.s =====
  Int_t Dist_BDC1BDC2 = 998.; //mm
  Int_t Dist_BDC1Target = 2231.;
  Int_t Dist_BDC1FDC1 = 3533.;
  Int_t Dist_SBTTarget = 2795.;
  Int_t Width_BDC1 = 68.;
  Int_t Width_FDC1 = 180.;
  
  //===== Declare anatree const.s =====
  Int_t nohit_layerX, nohit_layerY, nohit_layerU, nohit_layerV;
  Int_t tmphitnumX, tmphitnumY, tmphitnumU, tmphitnumV;
  Int_t tr;

  Double_t FDC1_trackX_z[6] = {-65,-55,-5,5,55,65};
  Double_t FDC1_trackU_z[4] = {-45,-35,15,25};
  Double_t FDC1_trackV_z[4] = {-25,15,35,45};

  Double_t FDC2_trackX_z[6] = {-308.66, -291.34, -8.66, 8.66, 291.34, 308.66};
  Double_t FDC2_trackU_z[4] = {-208.66, -191.34, 91.34, 108.66};
  Double_t FDC2_trackV_z[4] = {-108.66, -91.34, 191.34, 208.66};


  //===== Declare valables for calc. =====
  Int_t BDC1_hit_tdc[8][16];
  Double_t BDC1_hit_pos[8][16];
  Double_t BDC1_trackX_pos[16][4], BDC1_trackY_pos[16][4];
  Double_t BDC1_trackX_mm[16][4], BDC1_trackY_mm[16][4];
  Int_t BDC1_allhitnumX;
  Int_t BDC1_allhitnumY;

  Int_t BDC2_hit_tdc[8][16];
  Double_t BDC2_hit_pos[8][16];
  Double_t BDC2_trackX_pos[16][4], BDC2_trackY_pos[16][4];
  Double_t BDC2_trackX_mm[16][4], BDC2_trackY_mm[16][4];
  Int_t BDC2_allhitnumX;
  Int_t BDC2_allhitnumY;

  Int_t FDC1_hit_tdc[14][32];
  Double_t FDC1_hit_pos[14][32];
  Int_t FDC1_hit_num[14];
  Double_t FDC1_trackX_pos[16][6], FDC1_trackU_pos[16][4], FDC1_trackV_pos[16][4];
  Double_t FDC1_trackX_mm[16][6], FDC1_trackU_mm[16][4], FDC1_trackV_mm[16][4];
  Double_t FDC1_trackX[64][6], FDC1_trackU[64][4], FDC1_trackV[64][4];
  Int_t FDC1_allhitnumX;
  Int_t FDC1_allhitnumU;
  Int_t FDC1_allhitnumV;  

  Int_t FDC2_hit_tdc[14][112];
  Double_t FDC2_hit_pos[14][112];
  Int_t FDC2_hit_num[14];
  Double_t FDC2_trackX_pos[16][6], FDC2_trackU_pos[16][4], FDC2_trackV_pos[16][4];
  Double_t FDC2_trackX_mm[16][6], FDC2_trackU_mm[16][4], FDC2_trackV_mm[16][4];
  Double_t FDC2_trackX[64][6], FDC2_trackU[64][4], FDC2_trackV[64][4];
  

  //===== Declare anatree variables =====
  Int_t RunNum, EventNum;
  
  Int_t BDC1_hit_num[8];
  Double_t BDC1_X, BDC1_Y, BDC1_Chi2X, BDC1_Chi2Y;

  Int_t BDC2_hit_num[8];
  Double_t BDC2_X, BDC2_Y, BDC2_Chi2X, BDC2_Chi2Y;

  Double_t BDC_X, BDC_Y;
  Double_t BDC_A, BDC_B;
  Double_t Target_X, Target_Y;

  Double_t FDC1_Xpos, FDC1_Upos, FDC1_Vpos, FDC1_Chi2X, FDC1_Chi2U, FDC1_Chi2V;
  Double_t FDC1_X, FDC1_Y, FDC1_A, FDC1_B;
  
  Double_t FDC2_Xpos, FDC2_Upos, FDC2_Vpos, FDC2_Chi2X, FDC2_Chi2U, FDC2_Chi2V;
  Double_t FDC2_X, FDC2_Y, FDC2_A, FDC2_B;

  Double_t brhoSAMURAI;

  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatrS->Branch("RunNum",&RunNum);
  anatrS->Branch("EventNum",&EventNum);
  
  anatrS->Branch("BDC1_hit_num",BDC1_hit_num,"BDC1_hit_num[8]/I");
  anatrS->Branch("BDC1_X",&BDC1_X);
  anatrS->Branch("BDC1_Y",&BDC1_Y);
  anatrS->Branch("BDC1_Chi2X",&BDC1_Chi2X);
  anatrS->Branch("BDC1_Chi2Y",&BDC1_Chi2Y);

  anatrS->Branch("BDC2_hit_num",BDC2_hit_num,"BDC2_hit_num[8]/I");
  anatrS->Branch("BDC2_X",&BDC2_X);
  anatrS->Branch("BDC2_Y",&BDC2_Y);
  anatrS->Branch("BDC2_Chi2X",&BDC2_Chi2X);
  anatrS->Branch("BDC2_Chi2Y",&BDC2_Chi2Y);

  anatrS->Branch("BDC_X",&BDC_X);
  anatrS->Branch("BDC_Y",&BDC_Y);
  anatrS->Branch("BDC_A",&BDC_A);
  anatrS->Branch("BDC_B",&BDC_B);

  anatrS->Branch("Target_X",&Target_X);
  anatrS->Branch("Target_Y",&Target_Y);


  anatrS->Branch("FDC1_Xpos",&FDC1_Xpos);
  anatrS->Branch("FDC1_Upos",&FDC1_Upos);
  anatrS->Branch("FDC1_Vpos",&FDC1_Vpos);
  anatrS->Branch("FDC1_Chi2X",&FDC1_Chi2X);
  anatrS->Branch("FDC1_Chi2U",&FDC1_Chi2U);
  anatrS->Branch("FDC1_Chi2V",&FDC1_Chi2V);

  anatrS->Branch("FDC1_X",&FDC1_X);
  anatrS->Branch("FDC1_Y",&FDC1_Y);
  anatrS->Branch("FDC1_A",&FDC1_A);
  anatrS->Branch("FDC1_B",&FDC1_B);

  anatrS->Branch("FDC2_Xpos",&FDC2_Xpos);
  anatrS->Branch("FDC2_Upos",&FDC2_Upos);
  anatrS->Branch("FDC2_Vpos",&FDC2_Vpos);
  anatrS->Branch("FDC2_Chi2X",&FDC2_Chi2X);
  anatrS->Branch("FDC2_Chi2U",&FDC2_Chi2U);
  anatrS->Branch("FDC2_Chi2V",&FDC2_Chi2V);

  anatrS->Branch("FDC2_X",&FDC2_X);
  anatrS->Branch("FDC2_Y",&FDC2_Y);
  anatrS->Branch("FDC2_A",&FDC2_A);
  anatrS->Branch("FDC2_B",&FDC2_B);

  anatrS->Branch("brhoSAMURAI",&brhoSAMURAI);

  anatrS->Branch("BG_flag",&BG_flag);
  
  //===== Begin LOOP =====
  int nEntry = caltr->GetEntries();
  //for(int iEntry=0;iEntry<100;++iEntry){
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    
    if(iEntry%100 == 0){
      clog<< iEntry/1000 << "k events treated..." << "\r";
    }

    RunNum = RunNumber;
    EventNum = EventNumber;
    
    caltr->GetEntry(iEntry);

    //@@@ BDC1 @@@
    //=== Initialization ===    
    for(Int_t l=0;l<8;++l){
      for(Int_t i=0;i<16;++i){	
	BDC1_hit_tdc[l][i] = TMath::Sqrt(-1);
	BDC1_hit_pos[l][i] = TMath::Sqrt(-1);
      }
    }
    
    for(Int_t l=0;l<8;++l){
      BDC1_hit_num[l] = -9999;
    }

    for(Int_t i=0;i<16;++i){
      for(Int_t j=0;j<4;++j){
	BDC1_trackX_pos[i][j] = TMath::Sqrt(-1);
	BDC1_trackX_mm[i][j] = TMath::Sqrt(-1);
      }
    }

    BDC1_allhitnumX = 0;
    BDC1_allhitnumY = 0;
    
    nohit_layerX = 0;
    nohit_layerY = 0;
    BDC1_X = TMath::Sqrt(-1);
    BDC1_Y = TMath::Sqrt(-1);
    BDC1_Chi2X = 9999.;
    BDC1_Chi2Y = 9999.;

    BG_flag = 0;
    
    //=== Calclation ===

    for(Int_t l=0;l<8;++l){
      Int_t h = 0;
      for(Int_t i=0;i<16;++i){
	if(BDC1_TDC[l][i]>0&&BDC1_WirePosition[l][i]>-100.){
	  BDC1_hit_tdc[l][h] = BDC1_TDC[l][i];
	  BDC1_hit_pos[l][h] = BDC1_WirePosition[l][i];
	  h++;
	}else continue;
      }
      BDC1_hit_num[l] = h;
      if(h==0){
	if(l%4==0||l%4==1){nohit_layerX++;}
	else if(l%4==2||l%4==3){nohit_layerY++;}
      }
    }
    
    tmphitnumX = 0;
    tmphitnumY = 0;

    tmphitnumX = BDC1_hit_num[0]*BDC1_hit_num[1]*BDC1_hit_num[4]*BDC1_hit_num[5];
    tmphitnumY = BDC1_hit_num[2]*BDC1_hit_num[3]*BDC1_hit_num[6]*BDC1_hit_num[7];
    if(tmphitnumX>BDC1_allhitnumX) BDC1_allhitnumX = tmphitnumX;
    if(tmphitnumY>BDC1_allhitnumY) BDC1_allhitnumY = tmphitnumY;
    
    //cout << "BDC1_hit_ num = {" << BDC1_hit_num[0] << BDC1_hit_num[1] << BDC1_hit_num[2] << BDC1_hit_num[3] << "}" <<endl;
    
    if(nohit_layerX>=2||nohit_layerY>=2) BG_flag = 1; //if(nohit_layer<=1) BG_flag = 0
    
    tr = 0;
    for(Int_t hl0=0;hl0<BDC1_hit_num[0]&&tr<16;++hl0){
      for(Int_t hl1=0;hl1<BDC1_hit_num[1]&&tr<16;++hl1){
	for(Int_t hl4=0;hl4<BDC1_hit_num[4]&&tr<16;++hl4){
	  for(Int_t hl5=0;hl5<BDC1_hit_num[5]&&tr<16;++hl5){
	    BDC1_trackX_pos[tr][0] = BDC1_hit_pos[0][hl0];
	    BDC1_trackX_pos[tr][1] = BDC1_hit_pos[1][hl1];
	    BDC1_trackX_pos[tr][2] = BDC1_hit_pos[4][hl4];
	    BDC1_trackX_pos[tr][3] = BDC1_hit_pos[5][hl5];
	    BDC1_trackX_mm[tr][0] = tdc2mm[0][0][BDC1_hit_tdc[0][hl0]];
	    BDC1_trackX_mm[tr][1] = tdc2mm[0][1][BDC1_hit_tdc[1][hl1]];
	    BDC1_trackX_mm[tr][2] = tdc2mm[0][4][BDC1_hit_tdc[4][hl4]];
	    BDC1_trackX_mm[tr][3] = tdc2mm[0][5][BDC1_hit_tdc[5][hl5]];
	    tr++;
	  }
	}
      }
    }
    
    tr = 0;
    for(Int_t hl2=0;hl2<BDC1_hit_num[2]&&tr<16;++hl2){
      for(Int_t hl3=0;hl3<BDC1_hit_num[3]&&tr<16;++hl3){
	for(Int_t hl6=0;hl6<BDC1_hit_num[6]&&tr<16;++hl6){
	  for(Int_t hl7=0;hl7<BDC1_hit_num[7]&&tr<16;++hl7){
	    BDC1_trackY_pos[tr][0] = BDC1_hit_pos[2][hl2];
	    BDC1_trackY_pos[tr][1] = BDC1_hit_pos[3][hl3];
	    BDC1_trackY_pos[tr][2] = BDC1_hit_pos[6][hl6];
	    BDC1_trackY_pos[tr][3] = BDC1_hit_pos[7][hl7];
	    BDC1_trackY_mm[tr][0] = tdc2mm[0][2][BDC1_hit_tdc[2][hl2]];
	    BDC1_trackY_mm[tr][1] = tdc2mm[0][3][BDC1_hit_tdc[3][hl3]];
	    BDC1_trackY_mm[tr][2] = tdc2mm[0][6][BDC1_hit_tdc[6][hl6]];
	    BDC1_trackY_mm[tr][3] = tdc2mm[0][7][BDC1_hit_tdc[7][hl7]];
	    tr++;
	  }
	}
      }
    }
    

    
    for(tr=0;tr<16;++tr){      if(!(BDC1_trackX_pos[tr][0]>-100&&BDC1_trackX_pos[tr][0]<100)) continue;
      for(Int_t i=0;i<16;++i){
	Double_t tempX = TMath::Sqrt(-1);
	Double_t tempY = TMath::Sqrt(-1);
	Double_t tempChi2X = 0.;
	Double_t tempChi2Y = 0.;
	bitset<4>bit_i(i);
	Int_t sign[4] = {0,0,0,0};
	for(Int_t j=0;j<4;++j){
	  if(bit_i[j]){
	    sign[j] = 1;
	  }else{
	    sign[j] = 0;
	  }
	}
	//cout << i << "bit_i = {" << bit_i[0] << bit_i[1] << bit_i[2] << bit_i[3] << "}, sign = {" << pow(-1,sign[0]) << sign[1] << sign[2] << sign[3] << "}" << endl;
	
	tempX = (BDC1_trackX_pos[tr][0] + pow(-1,sign[0])*BDC1_trackX_mm[tr][0] +
		 BDC1_trackX_pos[tr][1] + pow(-1,sign[1])*BDC1_trackX_mm[tr][1] +
		 BDC1_trackX_pos[tr][2] + pow(-1,sign[2])*BDC1_trackX_mm[tr][2] +
		 BDC1_trackX_pos[tr][3] + pow(-1,sign[3])*BDC1_trackX_mm[tr][3])/4.;
	tempY = (BDC1_trackY_pos[tr][0] + pow(-1,sign[0])*BDC1_trackY_mm[tr][0] +
		 BDC1_trackY_pos[tr][1] + pow(-1,sign[1])*BDC1_trackY_mm[tr][1] +
		 BDC1_trackY_pos[tr][2] + pow(-1,sign[2])*BDC1_trackY_mm[tr][2] +
		 BDC1_trackY_pos[tr][3] + pow(-1,sign[3])*BDC1_trackY_mm[tr][3])/4.;
	for(Int_t j=0;j<4;++j){
	  tempChi2X += pow(tempX - (BDC1_trackX_pos[tr][j] + pow(-1,sign[j])*BDC1_trackX_mm[tr][j]),2);
	  tempChi2Y += pow(tempY - (BDC1_trackY_pos[tr][j] + pow(-1,sign[j])*BDC1_trackY_mm[tr][j]),2);
	}
	//cout << "tempX = " << tempX << ", tempChi2 = " << tempChi2 << endl;
	if(tempChi2X<BDC1_Chi2X){
	  BDC1_Chi2X = tempChi2X;
	  BDC1_X = tempX;
	}
	if(tempChi2Y<BDC1_Chi2Y){
	  BDC1_Chi2Y = tempChi2Y;
	  BDC1_Y = tempY;
	}
      }
    }
    //@@@ BDC1 end @@@

    //@@@ BDC2 @@@
    //=== Initialization ===    
    for(Int_t l=0;l<8;++l){
      for(Int_t i=0;i<16;++i){
	BDC2_hit_tdc[l][i] = TMath::Sqrt(-1);
	BDC2_hit_pos[l][i] = TMath::Sqrt(-1);
      }
    }
    
    for(Int_t l=0;l<8;++l){
      BDC2_hit_num[l] = -9999;
    }

    for(Int_t i=0;i<16;++i){
      for(Int_t j=0;j<4;++j){
	BDC2_trackX_pos[i][j] = TMath::Sqrt(-1);
	BDC2_trackX_mm[i][j] = TMath::Sqrt(-1);
      }
    }

    BDC2_allhitnumX = 0;
    BDC2_allhitnumY = 0;
    
    nohit_layerX = 0;
    nohit_layerY = 0;
    BDC2_X = TMath::Sqrt(-1);
    BDC2_Y = TMath::Sqrt(-1);
    BDC2_Chi2X = 9999.;
    BDC2_Chi2Y = 9999.;

    BG_flag = 0;
    
    //=== Calclation ===

    for(Int_t l=0;l<8;++l){
      Int_t h = 0;
      for(Int_t i=0;i<16;++i){
	if(BDC2_TDC[l][i]>0&&BDC2_WirePosition[l][i]>-100.){
	  BDC2_hit_tdc[l][h] = BDC2_TDC[l][i];
	  BDC2_hit_pos[l][h] = BDC2_WirePosition[l][i];
	  h++;
	}else continue;
      }
      BDC2_hit_num[l] = h;
      if(h==0){
	if(l%4==0||l%4==1){nohit_layerX++;}
	else if(l%4==2||l%4==3){nohit_layerY++;}
      }
    }
    
    tmphitnumX = 0;
    tmphitnumY = 0;

    tmphitnumX = BDC2_hit_num[0]*BDC2_hit_num[1]*BDC2_hit_num[4]*BDC2_hit_num[5];
    tmphitnumY = BDC2_hit_num[2]*BDC2_hit_num[3]*BDC2_hit_num[6]*BDC2_hit_num[7];
    if(tmphitnumX>BDC2_allhitnumX) BDC2_allhitnumX = tmphitnumX;
    if(tmphitnumY>BDC2_allhitnumY) BDC2_allhitnumY = tmphitnumY;
    
    //cout << "BDC2_hit_ num = {" << BDC2_hit_num[0] << BDC2_hit_num[1] << BDC2_hit_num[2] << BDC2_hit_num[3] << "}" <<endl;
    
    if(nohit_layerX>=2||nohit_layerY>=2) BG_flag = 1; //if(nohit_layer<=1) BG_flag = 0
    
    tr = 0;
    for(Int_t hl0=0;hl0<BDC2_hit_num[0]&&tr<16;++hl0){
      for(Int_t hl1=0;hl1<BDC2_hit_num[1]&&tr<16;++hl1){
	for(Int_t hl4=0;hl4<BDC2_hit_num[4]&&tr<16;++hl4){
	  for(Int_t hl5=0;hl5<BDC2_hit_num[5]&&tr<16;++hl5){
	    BDC2_trackX_pos[tr][0] = BDC2_hit_pos[0][hl0];
	    BDC2_trackX_pos[tr][1] = BDC2_hit_pos[1][hl1];
	    BDC2_trackX_pos[tr][2] = BDC2_hit_pos[4][hl4];
	    BDC2_trackX_pos[tr][3] = BDC2_hit_pos[5][hl5];
	    BDC2_trackX_mm[tr][0] = tdc2mm[0][0][BDC2_hit_tdc[0][hl0]];
	    BDC2_trackX_mm[tr][1] = tdc2mm[0][1][BDC2_hit_tdc[1][hl1]];
	    BDC2_trackX_mm[tr][2] = tdc2mm[0][4][BDC2_hit_tdc[4][hl4]];
	    BDC2_trackX_mm[tr][3] = tdc2mm[0][5][BDC2_hit_tdc[5][hl5]];
	    tr++;
	  }
	}
      }
    }
    
    tr = 0;
    for(Int_t hl2=0;hl2<BDC2_hit_num[2]&&tr<16;++hl2){
      for(Int_t hl3=0;hl3<BDC2_hit_num[3]&&tr<16;++hl3){
	for(Int_t hl6=0;hl6<BDC2_hit_num[6]&&tr<16;++hl6){
	  for(Int_t hl7=0;hl7<BDC2_hit_num[7]&&tr<16;++hl7){
	    BDC2_trackY_pos[tr][0] = BDC2_hit_pos[2][hl2];
	    BDC2_trackY_pos[tr][1] = BDC2_hit_pos[3][hl3];
	    BDC2_trackY_pos[tr][2] = BDC2_hit_pos[6][hl6];
	    BDC2_trackY_pos[tr][3] = BDC2_hit_pos[7][hl7];
	    BDC2_trackY_mm[tr][0] = tdc2mm[0][2][BDC2_hit_tdc[2][hl2]];
	    BDC2_trackY_mm[tr][1] = tdc2mm[0][3][BDC2_hit_tdc[3][hl3]];
	    BDC2_trackY_mm[tr][2] = tdc2mm[0][6][BDC2_hit_tdc[6][hl6]];
	    BDC2_trackY_mm[tr][3] = tdc2mm[0][7][BDC2_hit_tdc[7][hl7]];
	    tr++;
	  }
	}
      }
    }
    

    
    for(tr=0;tr<16;++tr){
      if(!(BDC2_trackX_pos[tr][0]>-100&&BDC2_trackX_pos[tr][0]<100)) continue;
      for(Int_t i=0;i<16;++i){
	Double_t tempX = TMath::Sqrt(-1);
	Double_t tempY = TMath::Sqrt(-1);
	Double_t tempChi2X = 0.;
	Double_t tempChi2Y = 0.;
	bitset<4>bit_i(i);
	Int_t sign[4] = {0,0,0,0};
	for(Int_t j=0;j<4;++j){
	  if(bit_i[j]){
	    sign[j] = 1;
	  }else{
	    sign[j] = 0;
	  }
	}
	//cout << i << "bit_i = {" << bit_i[0] << bit_i[1] << bit_i[2] << bit_i[3] << "}, sign = {" << pow(-1,sign[0]) << sign[1] << sign[2] << sign[3] << "}" << endl;
	
	tempX = (BDC2_trackX_pos[tr][0] + pow(-1,sign[0])*BDC2_trackX_mm[tr][0] +
		 BDC2_trackX_pos[tr][1] + pow(-1,sign[1])*BDC2_trackX_mm[tr][1] +
		 BDC2_trackX_pos[tr][2] + pow(-1,sign[2])*BDC2_trackX_mm[tr][2] +
		 BDC2_trackX_pos[tr][3] + pow(-1,sign[3])*BDC2_trackX_mm[tr][3])/4.;
	tempY = (BDC2_trackY_pos[tr][0] + pow(-1,sign[0])*BDC2_trackY_mm[tr][0] +
		 BDC2_trackY_pos[tr][1] + pow(-1,sign[1])*BDC2_trackY_mm[tr][1] +
		 BDC2_trackY_pos[tr][2] + pow(-1,sign[2])*BDC2_trackY_mm[tr][2] +
		 BDC2_trackY_pos[tr][3] + pow(-1,sign[3])*BDC2_trackY_mm[tr][3])/4.;
	for(Int_t j=0;j<4;++j){
	  tempChi2X += pow(tempX - (BDC2_trackX_pos[tr][j] + pow(-1,sign[j])*BDC2_trackX_mm[tr][j]),2);
	  tempChi2Y += pow(tempY - (BDC2_trackY_pos[tr][j] + pow(-1,sign[j])*BDC2_trackY_mm[tr][j]),2);
	}
	//cout << "tempX = " << tempX << ", tempChi2 = " << tempChi2 << endl;
	if(tempChi2X<BDC2_Chi2X){
	  BDC2_Chi2X = tempChi2X;
	  BDC2_X = tempX;
	}
	if(tempChi2Y<BDC2_Chi2Y){
	  BDC2_Chi2Y = tempChi2Y;
	  BDC2_Y = tempY;
	}
      }
    }
    //@@@ BDC2 end @@@

    //@@@ BDC @@@
    //===== Initialization =====
    BDC_X = TMath::Sqrt(-1);
    BDC_Y = TMath::Sqrt(-1);
    BDC_A = TMath::Sqrt(-1);
    BDC_B = TMath::Sqrt(-1);
    Target_X = TMath::Sqrt(-1);
    Target_Y = TMath::Sqrt(-1);
    
    //===== Calc. =====
    BDC_X = (BDC1_X + BDC2_X)/2.;
    BDC_Y = (BDC1_Y + BDC2_Y)/2.;

    BDC_A = 1000.*TMath::ATan((BDC2_X - BDC1_X)/Dist_BDC1BDC2); //[mrad]
    BDC_B = 1000.*TMath::ATan((BDC2_Y - BDC1_Y)/Dist_BDC1BDC2); //[mrad]

    Target_X = BDC1_X + Dist_BDC1Target*TMath::Tan(BDC_A/1000.);
    Target_Y = BDC1_Y + Dist_BDC1Target*TMath::Tan(BDC_B/1000.);
    //@@@ BDC end @@@

    //@@@ FDC1 @@@
    //=== Initialization ===    
    for(Int_t l=0;l<14;++l){
      for(Int_t i=0;i<32;++i){
	FDC1_hit_tdc[l][i] = -9999;
	FDC1_hit_pos[l][i] = TMath::Sqrt(-1);
      }
    }
    
    for(Int_t l=0;l<14;++l){
      FDC1_hit_num[l] = -9999;
    }

    for(Int_t i=0;i<16;++i){
      for(Int_t j=0;j<6;++j){
	FDC1_trackX_pos[i][j] = TMath::Sqrt(-1);
	FDC1_trackX_mm[i][j] = TMath::Sqrt(-1);
      }
      for(Int_t j=0;j<4;++j){
	FDC1_trackU_pos[i][j] = TMath::Sqrt(-1);
	FDC1_trackU_mm[i][j] = TMath::Sqrt(-1);
	FDC1_trackV_pos[i][j] = TMath::Sqrt(-1);
	FDC1_trackV_mm[i][j] = TMath::Sqrt(-1);
      }
    }

    FDC1_allhitnumX = 0;
    FDC1_allhitnumU = 0;
    FDC1_allhitnumV = 0;
    
    nohit_layerX = 0;
    nohit_layerU = 0;
    nohit_layerV = 0;
    FDC1_Xpos = TMath::Sqrt(-1);
    FDC1_Upos = TMath::Sqrt(-1);
    FDC1_Vpos = TMath::Sqrt(-1);
    FDC1_X = TMath::Sqrt(-1);
    FDC1_Y = TMath::Sqrt(-1);
    FDC1_A = TMath::Sqrt(-1);
    FDC1_B = TMath::Sqrt(-1);
    FDC1_Chi2X = 9999.;
    FDC1_Chi2U = 9999.;
    FDC1_Chi2V = 9999.;

    BG_flag = 0;
    
    //=== Calclation ===

    for(Int_t l=0;l<14;++l){
      Int_t h = 0;
      for(Int_t i=0;i<32;++i){
	if(FDC1_TDC[l][i]>0&&FDC1_WirePosition[l][i]>-200.){
	  FDC1_hit_tdc[l][h] = FDC1_TDC[l][i];
	  FDC1_hit_pos[l][h] = FDC1_WirePosition[l][i];
	  h++;
	}else continue;
      }
      FDC1_hit_num[l] = h;
      if(h==0){
	Int_t l6 = l%6;
	switch(l6){
	case 0: nohit_layerX++; break;
	case 1: nohit_layerX++; break;
	case 2: nohit_layerU++; break;
	case 3: nohit_layerU++; break;
	case 4: nohit_layerV++; break;
	case 5: nohit_layerV++; break;
	}
	//if(l%4==0||l%4==1){nohit_layerX++;}
	//else if(l%4==2||l%4==3){nohit_layerY++;}
      }
    }
    
    tmphitnumX = 0;
    tmphitnumU = 0;
    tmphitnumV = 0;

    tmphitnumX = FDC1_hit_num[0]*FDC1_hit_num[1]*FDC1_hit_num[6]*FDC1_hit_num[7]*FDC1_hit_num[12]*FDC1_hit_num[13];
    tmphitnumU = FDC1_hit_num[2]*FDC1_hit_num[3]*FDC1_hit_num[8]*FDC1_hit_num[9];
    tmphitnumV = FDC1_hit_num[4]*FDC1_hit_num[5]*FDC1_hit_num[10]*FDC1_hit_num[11];
    if(tmphitnumX>FDC1_allhitnumX) FDC1_allhitnumX = tmphitnumX;
    if(tmphitnumU>FDC1_allhitnumU) FDC1_allhitnumU = tmphitnumU;
    if(tmphitnumV>FDC1_allhitnumV) FDC1_allhitnumV = tmphitnumV;
    
    //cout << "BDC1_hit_ num = {" << BDC1_hit_num[0] << BDC1_hit_num[1] << BDC1_hit_num[2] << BDC1_hit_num[3] << "}" <<endl;
    
    if(nohit_layerX>=3||nohit_layerU>=2||nohit_layerV>=2) BG_flag = 1; //if(nohit_layer<=1) BG_flag = 0

    //X
    tr = 0;
    for(Int_t hl0=0;hl0<FDC1_hit_num[0]&&tr<16;++hl0){
      for(Int_t hl1=0;hl1<FDC1_hit_num[1]&&tr<16;++hl1){
	for(Int_t hl6=0;hl6<FDC1_hit_num[6]&&tr<16;++hl6){
	  for(Int_t hl7=0;hl7<FDC1_hit_num[7]&&tr<16;++hl7){
	    for(Int_t hl12=0;hl12<FDC1_hit_num[12]&&tr<16;++hl12){
	      for(Int_t hl13=0;hl13<FDC1_hit_num[13]&&tr<16;++hl13){
		FDC1_trackX_pos[tr][0] = FDC1_hit_pos[0][hl0];
		FDC1_trackX_pos[tr][1] = FDC1_hit_pos[1][hl1];
		FDC1_trackX_pos[tr][2] = FDC1_hit_pos[6][hl6];
		FDC1_trackX_pos[tr][3] = FDC1_hit_pos[7][hl7];
		FDC1_trackX_pos[tr][4] = FDC1_hit_pos[12][hl12];
		FDC1_trackX_pos[tr][5] = FDC1_hit_pos[13][hl13];
		FDC1_trackX_mm[tr][0] = tdc2mm[2][0][FDC1_hit_tdc[0][hl0]];
		FDC1_trackX_mm[tr][1] = tdc2mm[2][1][FDC1_hit_tdc[1][hl1]];
		FDC1_trackX_mm[tr][2] = tdc2mm[2][6][FDC1_hit_tdc[6][hl6]];
		FDC1_trackX_mm[tr][3] = tdc2mm[2][7][FDC1_hit_tdc[7][hl7]];
		FDC1_trackX_mm[tr][4] = tdc2mm[2][12][FDC1_hit_tdc[12][hl12]];
		FDC1_trackX_mm[tr][5] = tdc2mm[2][13][FDC1_hit_tdc[13][hl13]];
		tr++;
	      }
	    }
	  }
	}
      }
    }
    
    //U
    tr = 0;
    for(Int_t hl2=0;hl2<FDC1_hit_num[2]&&tr<16;++hl2){
      for(Int_t hl3=0;hl3<FDC1_hit_num[3]&&tr<16;++hl3){
	for(Int_t hl8=0;hl8<FDC1_hit_num[8]&&tr<16;++hl8){
	  //for(Int_t hl9=0;hl9<FDC1_hit_num[9]+1&&tr<16;++hl9){
	    FDC1_trackU_pos[tr][0] = FDC1_hit_pos[2][hl2];
	    FDC1_trackU_pos[tr][1] = FDC1_hit_pos[3][hl3];
	    FDC1_trackU_pos[tr][2] = FDC1_hit_pos[8][hl8];
	    //FDC1_trackU_pos[tr][3] = FDC1_hit_pos[9][hl9];
	    FDC1_trackU_mm[tr][0] = tdc2mm[2][2][FDC1_hit_tdc[2][hl2]];
	    FDC1_trackU_mm[tr][1] = tdc2mm[2][3][FDC1_hit_tdc[3][hl3]];
	    FDC1_trackU_mm[tr][2] = tdc2mm[2][8][FDC1_hit_tdc[8][hl8]];
	    //FDC1_trackU_mm[tr][3] = tdc2mm[2][9][FDC1_hit_tdc[9][hl9]];

	    //cout << "U " << tr << " " << FDC1_trackU_pos[tr][0] << endl;

	    tr++;
	    //}
	}
      }
    }

    //V
    tr = 0;
    for(Int_t hl4=0;hl4<FDC1_hit_num[4]&&tr<16;++hl4){
      for(Int_t hl5=0;hl5<FDC1_hit_num[5]&&tr<16;++hl5){
	for(Int_t hl10=0;hl10<FDC1_hit_num[10]&&tr<16;++hl10){
	  for(Int_t hl11=0;hl11<FDC1_hit_num[11]&&tr<16;++hl11){
	    FDC1_trackV_pos[tr][0] = FDC1_hit_pos[4][hl4];
	    FDC1_trackV_pos[tr][1] = FDC1_hit_pos[5][hl5];
	    FDC1_trackV_pos[tr][2] = FDC1_hit_pos[10][hl10];
	    FDC1_trackV_pos[tr][3] = FDC1_hit_pos[11][hl11];
	    FDC1_trackV_mm[tr][0] = tdc2mm[2][4][FDC1_hit_tdc[4][hl4]];
	    FDC1_trackV_mm[tr][1] = tdc2mm[2][5][FDC1_hit_tdc[5][hl5]];
	    FDC1_trackV_mm[tr][2] = tdc2mm[2][10][FDC1_hit_tdc[10][hl10]];
	    FDC1_trackV_mm[tr][3] = tdc2mm[2][11][FDC1_hit_tdc[11][hl11]];

	    //cout << "V " << tr << " " << FDC1_trackV_pos[tr][0] << endl;
	    
	    tr++;
	  }
	}
      }
    }
    
    Double_t FDC1_slopeX = TMath::Sqrt(-1);
    Double_t FDC1_interseptX = TMath::Sqrt(-1);
    Double_t FDC1_slopeU = TMath::Sqrt(-1);
    Double_t FDC1_interseptU = TMath::Sqrt(-1);
    Double_t FDC1_slopeV = TMath::Sqrt(-1);
    Double_t FDC1_interseptV = TMath::Sqrt(-1);

    
    for(tr=0;tr<16;++tr){

      for(Int_t i=0;i<64;++i){
	for(Int_t l=0;l<6;++l) FDC1_trackX[i][l] = TMath::Sqrt(-1);
	for(Int_t l=0;l<4;++l){
	  FDC1_trackU[i][l] = TMath::Sqrt(-1);
	  FDC1_trackV[i][l] = TMath::Sqrt(-1);      
	}
      }

 
      //if(!(FDC1_trackX_pos[tr][0]>-200&&FDC1_trackX_pos[tr][0]<200)) continue;
      //X
      for(Int_t i=0;i<64;++i){
	Double_t tempX = TMath::Sqrt(-1);
	Double_t tempChi2X = 0.;
	bitset<6>bitX(i);
	Int_t signX[6] = {0,0,0,0,0,0};
	for(Int_t j=0;j<6;++j){
	  if(bitX[j]){
	    signX[j] = 1;
	  }else{
	    signX[j] = 0;
	  }
	}

	for(Int_t l=0;l<6;++l){
	  FDC1_trackX[i][l] = FDC1_trackX_pos[tr][l] + pow(-1,signX[l])*FDC1_trackX_mm[tr][l];
	}

	if(FDC1_trackX[i][0]>-1200&&FDC1_trackX[i][0]<1200){
	
	//===== Least square =====
	Double_t aX[2] = {-9999,-9999};
	Double_t sX[3] = {0};
	Double_t tX[2] = {0};
	
	for(Int_t l=0;l<6;++l){
	  sX[0] += 1.;
	  sX[1] += FDC1_trackX_z[l];
	  sX[2] += pow(FDC1_trackX_z[l],2);
	  tX[0] += FDC1_trackX[i][l];
	  tX[1] += FDC1_trackX[i][l]*FDC1_trackX_z[l];
	}
	aX[0] = (sX[2]*tX[0]-sX[1]*tX[1])/(sX[0]*sX[2]-sX[1]*sX[1]);
	aX[1] = (sX[0]*tX[1]-sX[1]*tX[0])/(sX[0]*sX[2]-sX[1]*sX[1]);

	tempX = aX[0];

	//cout << aX[0] << " " << aX[1] << endl;
	
	for(Int_t l=0;l<6;++l){
	  tempChi2X += pow(FDC1_trackX[i][l] - (aX[1]*FDC1_trackX_z[l] + aX[0]),2);
	}

	//cout << tempChi2X << endl;
	
	if(tempChi2X<FDC1_Chi2X){
	  FDC1_Chi2X = tempChi2X;
	  FDC1_Xpos = tempX;
	  FDC1_slopeX = aX[1];
	  FDC1_interseptX = aX[0];
	}

	//cout << "Chi2 " << FDC2_Chi2X << "X " << FDC2_Xpos <<  endl;
	//cout << "slope " << slope << " const " << intersept << endl;
	
	
	}else continue;

      }
	
	//===== Least square end =====
	

      /*
	
	//cout << i << "bitX = {" << bitX[0] << bitX[1] << bitX[2] << bitX[3] << bitX[4] << bitX[5] <<  "}, signX = {" << pow(-1,signX[0]) << signX[1] << signX[2] << signX[3] << signX[4] << signX[5] << "}" << endl;
	
	tempX = (FDC1_trackX_pos[tr][0] + pow(-1,signX[0])*FDC1_trackX_mm[tr][0] +
		 FDC1_trackX_pos[tr][1] + pow(-1,signX[1])*FDC1_trackX_mm[tr][1] +
		 FDC1_trackX_pos[tr][2] + pow(-1,signX[2])*FDC1_trackX_mm[tr][2] +
		 FDC1_trackX_pos[tr][3] + pow(-1,signX[3])*FDC1_trackX_mm[tr][3] +
		 FDC1_trackX_pos[tr][4] + pow(-1,signX[4])*FDC1_trackX_mm[tr][4] +
		 FDC1_trackX_pos[tr][5] + pow(-1,signX[5])*FDC1_trackX_mm[tr][5])/6.;
	//cout << "tempX = " << tempX << endl;
	for(Int_t j=0;j<6;++j){
	  tempChi2X += pow(tempX - (FDC1_trackX_pos[tr][j] + pow(-1,signX[j])*FDC1_trackX_mm[tr][j]),2);
	}
	//cout << "tempX = " << tempX << ", tempChi2X = " << tempChi2X << endl;
	if(tempChi2X<FDC1_Chi2X){
	  FDC1_Chi2X = tempChi2X;
	  FDC1_Xpos = tempX;
	}
      }
      */

      
      //U,V
      for(Int_t i=0;i<16;++i){
	Double_t tempU = TMath::Sqrt(-1);
	Double_t tempV = TMath::Sqrt(-1);
	Double_t tempChi2U = 0.;
	Double_t tempChi2V = 0.;
	bitset<4>bitUV(i);
	Int_t signUV[4] = {0,0,0,0};
	for(Int_t j=0;j<4;++j){
	  if(bitUV[j]){
	    signUV[j] = 1;
	  }else{
	    signUV[j] = 0;
	  }
	}

	for(Int_t l=0;l<4;++l){
	  FDC1_trackU[i][l] = FDC1_trackU_pos[tr][l] + pow(-1,signUV[l])*FDC1_trackU_mm[tr][l];
	  FDC1_trackV[i][l] = FDC1_trackV_pos[tr][l] + pow(-1,signUV[l])*FDC1_trackV_mm[tr][l];
	}

	if(FDC2_trackU[i][0]>-1200&&FDC2_trackU[i][0]<1200){
	  
	  //===== Least square =====
	  Double_t aU[2] = {-9999,-9999};
	  Double_t sU[3] = {0};
	  Double_t tU[2] = {0};
	  
	  for(Int_t l=0;l<3;++l){
	    sU[0] += 1.;
	    sU[1] += FDC1_trackU_z[l];
	    sU[2] += pow(FDC1_trackU_z[l],2);
	    tU[0] += FDC1_trackU[i][l];
	    tU[1] += FDC1_trackU[i][l]*FDC1_trackU_z[l];
	  }
	  aU[0] = (sU[2]*tU[0]-sU[1]*tU[1])/(sU[0]*sU[2]-sU[1]*sU[1]);
	  aU[1] = (sU[0]*tU[1]-sU[1]*tU[0])/(sU[0]*sU[2]-sU[1]*sU[1]);
	  
	  tempU = aU[0];
	  
	  //cout << aU[0] << " " << aU[1] << endl;
	  
	  for(Int_t l=0;l<3;++l){
	    tempChi2U += pow(FDC1_trackU[i][l] - (aU[1]*FDC1_trackU_z[l] + aU[0]),2);
	  }
	  
	  //cout << tempChi2U << endl;
	  
	  if(tempChi2U<FDC1_Chi2U){
	    FDC1_Chi2U = tempChi2U;
	    FDC1_Upos = tempU;
	    FDC1_slopeU = aU[1];
	    FDC1_interseptU = aU[0];
	  }
	  
	  //cout << "Chi2 " << FDC2_Chi2U << "U " << FDC2_Upos <<  endl;
	  //cout << "slope " << slope << " const " << intersept << endl;
	  
	  
	}else continue;

	if(FDC1_trackV[i][0]>-1200&&FDC1_trackV[i][0]<1200){
	  
	  //===== Least square =====
	  Double_t aV[2] = {-9999,-9999};
	  Double_t sV[3] = {0};
	  Double_t tV[2] = {0};
	  
	  for(Int_t l=0;l<4;++l){
	    sV[0] += 1.;
	    sV[1] += FDC1_trackV_z[l];
	    sV[2] += pow(FDC1_trackV_z[l],2);
	    tV[0] += FDC1_trackV[i][l];
	  tV[1] += FDC1_trackV[i][l]*FDC1_trackV_z[l];
	  }
	  aV[0] = (sV[2]*tV[0]-sV[1]*tV[1])/(sV[0]*sV[2]-sV[1]*sV[1]);
	  aV[1] = (sV[0]*tV[1]-sV[1]*tV[0])/(sV[0]*sV[2]-sV[1]*sV[1]);
	  
	  tempV = aV[0];
	  
	//cout << aV[0] << " " << aV[1] << endl;
	  
	  for(Int_t l=0;l<4;++l){
	    tempChi2V += pow(FDC1_trackV[i][l] - (aV[1]*FDC1_trackV_z[l] + aV[0]),2);
	  }
	  
	  //cout << tempChi2V << endl;
	  
	  if(tempChi2V<FDC1_Chi2V){
	    FDC1_Chi2V = tempChi2V;
	    FDC1_Vpos = tempV;
	    FDC1_slopeV = aV[1];
	    FDC1_interseptV = aV[0];
	  }
	  
	  //cout << "Chi2 " << FDC1_Chi2V << "V " << FDC1_Vpos <<  endl;
	  //cout << "slope " << slope << " const " << intersept << endl;
	  
	  
	}else continue;

      }//i

    }//tr
    
    FDC1_X = (FDC1_Xpos + (FDC1_Vpos + FDC1_Upos)/2.)/2.;
    FDC1_Y = TMath::Sqrt(3)/2.*(FDC1_Vpos-FDC1_Upos);
    FDC1_A = (6*TMath::ATan(FDC1_slopeX) + 2/TMath::Sqrt(3)*4*TMath::ATan(FDC1_slopeU) + 2/TMath::Sqrt(3)*4*TMath::ATan(FDC1_slopeV))/14.;
    FDC1_B = (-TMath::ATan(FDC1_slopeU) + TMath::ATan(FDC1_slopeV))/2.;
    /*
	//if(FDC1_trackU_pos[tr][abc]>-200&&FDC1_trackU_pos[tr][abc]<200){
	tempU = (FDC1_trackU_pos[tr][0] + pow(-1,signUV[0])*FDC1_trackU_mm[tr][0] +
		 FDC1_trackU_pos[tr][1] + pow(-1,signUV[1])*FDC1_trackU_mm[tr][1] +
		 FDC1_trackU_pos[tr][2] + pow(-1,signUV[2])*FDC1_trackU_mm[tr][2]
		 FDC1_trackU_pos[tr][3] + pow(-1,signUV[3])*FDC1_trackU_mm[tr][3])/3.;	    
	//tempU += FDC1_trackU_pos[tr][abc] + pow(-1,signUV[abc])*FDC1_trackU_mm[tr][abc];
	//}else continue;
	//}
	//tempU = tempU/de;

	tempV = (FDC1_trackV_pos[tr][0] + pow(-1,signUV[0])*FDC1_trackV_mm[tr][0] +
		 FDC1_trackV_pos[tr][1] + pow(-1,signUV[1])*FDC1_trackV_mm[tr][1] +
		 FDC1_trackV_pos[tr][2] + pow(-1,signUV[2])*FDC1_trackV_mm[tr][2] +
		 FDC1_trackV_pos[tr][3] + pow(-1,signUV[3])*FDC1_trackV_mm[tr][3])/4.;
	for(Int_t j=0;j<3;++j){
	  //if(FDC1_trackU_pos[tr][abc]>-200&&FDC1_trackU_pos[tr][abc]<200){
	    tempChi2U += pow(tempU - (FDC1_trackU_pos[tr][j] + pow(-1,signUV[j])*FDC1_trackU_mm[tr][j]),2);
	    //}
	}
	for(Int_t j=0;j<4;++j){
	    tempChi2V += pow(tempV - (FDC1_trackV_pos[tr][j] + pow(-1,signUV[j])*FDC1_trackV_mm[tr][j]),2);	  
	}
	//cout << "tempU = " << tempU << ", tempChi2U = " << tempChi2U << endl;
	if(tempChi2U<FDC1_Chi2U){
	  FDC1_Chi2U = tempChi2U;
	  FDC1_Upos = tempU;
	}
	if(tempChi2V<FDC1_Chi2V){
	  FDC1_Chi2V = tempChi2V;
	  FDC1_Vpos = tempV;
	}
      }
    }
    
    FDC1_X = (FDC1_Xpos + (FDC1_Vpos + FDC1_Upos)/2.)/2.;
    FDC1_Y = TMath::Sqrt(3)/2.*(FDC1_Vpos-FDC1_Upos);

    */
    //@@@ FDC1 end @@@

    //@@@ FDC2 @@@
    //=== Initialization ===    
    for(Int_t l=0;l<14;++l){
      for(Int_t i=0;i<112;++i){
	FDC2_hit_tdc[l][i] = -9999;
	FDC2_hit_pos[l][i] = TMath::Sqrt(-1);
      }
    }
    
    for(Int_t l=0;l<14;++l){
      FDC2_hit_num[l] = -9999;
    }

    for(Int_t i=0;i<16;++i){
      for(Int_t j=0;j<6;++j){
	FDC2_trackX_pos[i][j] = TMath::Sqrt(-1);
	FDC2_trackX_mm[i][j] = TMath::Sqrt(-1);
      }
      for(Int_t j=0;j<4;++j){
	FDC2_trackU_pos[i][j] = TMath::Sqrt(-1);
	FDC2_trackU_mm[i][j] = TMath::Sqrt(-1);
	FDC2_trackV_pos[i][j] = TMath::Sqrt(-1);
	FDC2_trackV_mm[i][j] = TMath::Sqrt(-1);
      }
    }
 
    nohit_layerX = 0;
    nohit_layerU = 0;
    nohit_layerV = 0;
    FDC2_Xpos = TMath::Sqrt(-1);
    FDC2_Upos = TMath::Sqrt(-1);
    FDC2_Vpos = TMath::Sqrt(-1);
    FDC2_X = TMath::Sqrt(-1);
    FDC2_Y = TMath::Sqrt(-1);
    FDC2_A = TMath::Sqrt(-1);
    FDC2_B = TMath::Sqrt(-1);
    FDC2_Chi2X = 9999.;
    FDC2_Chi2U = 9999.;
    FDC2_Chi2V = 9999.;

    BG_flag = 0;
    
    //=== Calclation ===
    
    for(Int_t l=0;l<14;++l){
      Int_t h = 0;
      for(Int_t i=0;i<112;++i){
	if(FDC2_TDC[l][i]>0&&FDC2_WirePosition[l][i]>-1200.){
	  FDC2_hit_tdc[l][h] = FDC2_TDC[l][i];
	  FDC2_hit_pos[l][h] = FDC2_WirePosition[l][i];
	  h++;
	}else continue;
      }
      FDC2_hit_num[l] = h;
      if(h==0){
	Int_t l6 = l%6;
	switch(l6){
	case 0: nohit_layerX++; break;
	case 1: nohit_layerX++; break;
	case 2: nohit_layerU++; break;
	case 3: nohit_layerU++; break;
	case 4: nohit_layerV++; break;
	case 5: nohit_layerV++; break;
	}
      }
    }
   
    if(nohit_layerX>=3||nohit_layerU>=2||nohit_layerV>=2) BG_flag = 1; //if(nohit_layer<=1) BG_flag = 0

    //X
    tr = 0;
    for(Int_t hl0=0;hl0<FDC2_hit_num[0]&&tr<16;++hl0){
      for(Int_t hl1=0;hl1<FDC2_hit_num[1]&&tr<16;++hl1){
	for(Int_t hl6=0;hl6<FDC2_hit_num[6]&&tr<16;++hl6){
	  for(Int_t hl7=0;hl7<FDC2_hit_num[7]&&tr<16;++hl7){
	    for(Int_t hl12=0;hl12<FDC2_hit_num[12]&&tr<16;++hl12){
	      for(Int_t hl13=0;hl13<FDC2_hit_num[13]&&tr<16;++hl13){
		FDC2_trackX_pos[tr][0] = FDC2_hit_pos[0][hl0];
		FDC2_trackX_pos[tr][1] = FDC2_hit_pos[1][hl1];
		FDC2_trackX_pos[tr][2] = FDC2_hit_pos[6][hl6];
		FDC2_trackX_pos[tr][3] = FDC2_hit_pos[7][hl7];
		FDC2_trackX_pos[tr][4] = FDC2_hit_pos[12][hl12];
		FDC2_trackX_pos[tr][5] = FDC2_hit_pos[13][hl13];
		FDC2_trackX_mm[tr][0] = tdc2mm[3][0][FDC2_hit_tdc[0][hl0]];
		FDC2_trackX_mm[tr][1] = tdc2mm[3][1][FDC2_hit_tdc[1][hl1]];
		FDC2_trackX_mm[tr][2] = tdc2mm[3][6][FDC2_hit_tdc[6][hl6]];
		FDC2_trackX_mm[tr][3] = tdc2mm[3][7][FDC2_hit_tdc[7][hl7]];
		FDC2_trackX_mm[tr][4] = tdc2mm[3][12][FDC2_hit_tdc[12][hl12]];
		FDC2_trackX_mm[tr][5] = tdc2mm[3][13][FDC2_hit_tdc[13][hl13]];	
		tr++;
	      }
	    }
	  }
	}
      }
    }
    
    //U
    tr = 0;
    for(Int_t hl2=0;hl2<FDC2_hit_num[2]&&tr<16;++hl2){
      for(Int_t hl3=0;hl3<FDC2_hit_num[3]&&tr<16;++hl3){
	for(Int_t hl8=0;hl8<FDC2_hit_num[8]&&tr<16;++hl8){
	  for(Int_t hl9=0;hl9<FDC2_hit_num[9]+1&&tr<16;++hl9){
	    FDC2_trackU_pos[tr][0] = FDC2_hit_pos[2][hl2];
	    FDC2_trackU_pos[tr][1] = FDC2_hit_pos[3][hl3];
	    FDC2_trackU_pos[tr][2] = FDC2_hit_pos[8][hl8];
	    FDC2_trackU_pos[tr][3] = FDC2_hit_pos[9][hl9];
	    FDC2_trackU_mm[tr][0] = tdc2mm[3][2][FDC2_hit_tdc[2][hl2]];
	    FDC2_trackU_mm[tr][1] = tdc2mm[3][3][FDC2_hit_tdc[3][hl3]];
	    FDC2_trackU_mm[tr][2] = tdc2mm[3][8][FDC2_hit_tdc[8][hl8]];
	    FDC2_trackU_mm[tr][3] = tdc2mm[3][9][FDC2_hit_tdc[9][hl9]];
	    tr++;
	  }
	}
      }
    }

    //V
    tr = 0;
    for(Int_t hl4=0;hl4<FDC2_hit_num[4]&&tr<16;++hl4){
      for(Int_t hl5=0;hl5<FDC2_hit_num[5]&&tr<16;++hl5){
	for(Int_t hl10=0;hl10<FDC2_hit_num[10]&&tr<16;++hl10){
	  for(Int_t hl11=0;hl11<FDC2_hit_num[11]&&tr<16;++hl11){
	    FDC2_trackV_pos[tr][0] = FDC2_hit_pos[4][hl4];
	    FDC2_trackV_pos[tr][1] = FDC2_hit_pos[5][hl5];
	    FDC2_trackV_pos[tr][2] = FDC2_hit_pos[10][hl10];
	    FDC2_trackV_pos[tr][3] = FDC2_hit_pos[11][hl11];
	    FDC2_trackV_mm[tr][0] = tdc2mm[3][4][FDC2_hit_tdc[4][hl4]];
	    FDC2_trackV_mm[tr][1] = tdc2mm[3][5][FDC2_hit_tdc[5][hl5]];
	    FDC2_trackV_mm[tr][2] = tdc2mm[3][10][FDC2_hit_tdc[10][hl10]];
	    FDC2_trackV_mm[tr][3] = tdc2mm[3][11][FDC2_hit_tdc[11][hl11]];
	    
	    //cout << "V " << tr << " " << FDC2_trackV_pos[tr][0] << endl;
	    
	    tr++;
	  }
	}
      }
    }

    Double_t slopeX = TMath::Sqrt(-1);
    Double_t interseptX = TMath::Sqrt(-1);
    Double_t slopeU = TMath::Sqrt(-1);
    Double_t interseptU = TMath::Sqrt(-1);
    Double_t slopeV = TMath::Sqrt(-1);
    Double_t interseptV = TMath::Sqrt(-1);

      
    
    for(tr=0;tr<16;++tr){

      for(Int_t i=0;i<64;++i){
	for(Int_t l=0;l<6;++l) FDC2_trackX[i][l] = TMath::Sqrt(-1);
	for(Int_t l=0;l<4;++l){
	  FDC2_trackU[i][l] = TMath::Sqrt(-1);
	  FDC2_trackV[i][l] = TMath::Sqrt(-1);      
      }
      }

      
      
      //if(!(FDC1_trackX_pos[tr][0]>-200&&FDC1_trackX_pos[tr][0]<200)) continue;
      //X
      for(Int_t i=0;i<64;++i){
	Double_t tempX = TMath::Sqrt(-1);
	Double_t tempChi2X = 0.;
	bitset<6>bitX(i);
	Int_t signX[6] = {0,0,0,0,0,0};
	for(Int_t j=0;j<6;++j){
	  if(bitX[j]){
	    signX[j] = 1;
	  }else{
	    signX[j] = 0;
	  }
	}
	
	for(Int_t l=0;l<6;++l){
	  FDC2_trackX[i][l] = FDC2_trackX_pos[tr][l] + pow(-1,signX[l])*FDC2_trackX_mm[tr][l];
	}

	if(FDC2_trackX[i][0]>-1200&&FDC2_trackX[i][0]<1200){
	
	//===== Least square =====
	Double_t aX[2] = {-9999,-9999};
	Double_t sX[3] = {0};
	Double_t tX[2] = {0};
	
	for(Int_t l=0;l<6;++l){
	  sX[0] += 1.;
	  sX[1] += FDC2_trackX_z[l];
	  sX[2] += pow(FDC2_trackX_z[l],2);
	  tX[0] += FDC2_trackX[i][l];
	  tX[1] += FDC2_trackX[i][l]*FDC2_trackX_z[l];
	}
	aX[0] = (sX[2]*tX[0]-sX[1]*tX[1])/(sX[0]*sX[2]-sX[1]*sX[1]);
	aX[1] = (sX[0]*tX[1]-sX[1]*tX[0])/(sX[0]*sX[2]-sX[1]*sX[1]);

	tempX = aX[0];

	//cout << aX[0] << " " << aX[1] << endl;
	
	for(Int_t l=0;l<6;++l){
	  tempChi2X += pow(FDC2_trackX[i][l] - (aX[1]*FDC2_trackX_z[l] + aX[0]),2);
	}

	//cout << tempChi2X << endl;
	
	if(tempChi2X<FDC2_Chi2X){
	  FDC2_Chi2X = tempChi2X;
	  FDC2_Xpos = tempX;
	  slopeX = aX[1];
	  interseptX = aX[0];
	}

	//cout << "Chi2 " << FDC2_Chi2X << "X " << FDC2_Xpos <<  endl;
	//cout << "slope " << slope << " const " << intersept << endl;
	
	
	}else continue;

      }
	
	//===== Least square end =====
	
      
      //U,V
      for(Int_t i=0;i<16;++i){
	Double_t tempU = TMath::Sqrt(-1);
	Double_t tempV = TMath::Sqrt(-1);
	Double_t tempChi2U = 0.;
	Double_t tempChi2V = 0.;
	bitset<4>bitUV(i);
	Int_t signUV[4] = {0,0,0,0};
	for(Int_t j=0;j<4;++j){
	  if(bitUV[j]){
	    signUV[j] = 1;
	  }else{
	    signUV[j] = 0;
	  }
	}
	
	for(Int_t l=0;l<4;++l){
	  FDC2_trackU[i][l] = FDC2_trackU_pos[tr][l] + pow(-1,signUV[l])*FDC2_trackU_mm[tr][l];
	  FDC2_trackV[i][l] = FDC2_trackV_pos[tr][l] + pow(-1,signUV[l])*FDC2_trackV_mm[tr][l];
	}

	if(FDC2_trackU[i][0]>-1200&&FDC2_trackU[i][0]<1200){
	  
	  //===== Least square =====
	  Double_t aU[2] = {-9999,-9999};
	  Double_t sU[3] = {0};
	  Double_t tU[2] = {0};
	  
	  for(Int_t l=0;l<4;++l){
	    sU[0] += 1.;
	    sU[1] += FDC2_trackU_z[l];
	    sU[2] += pow(FDC2_trackU_z[l],2);
	    tU[0] += FDC2_trackU[i][l];
	    tU[1] += FDC2_trackU[i][l]*FDC2_trackU_z[l];
	  }
	  aU[0] = (sU[2]*tU[0]-sU[1]*tU[1])/(sU[0]*sU[2]-sU[1]*sU[1]);
	  aU[1] = (sU[0]*tU[1]-sU[1]*tU[0])/(sU[0]*sU[2]-sU[1]*sU[1]);
	  
	  tempU = aU[0];
	  
	  //cout << aU[0] << " " << aU[1] << endl;
	  
	  for(Int_t l=0;l<4;++l){
	    tempChi2U += pow(FDC2_trackU[i][l] - (aU[1]*FDC2_trackU_z[l] + aU[0]),2);
	  }
	  
	  //cout << tempChi2U << endl;
	  
	  if(tempChi2U<FDC2_Chi2U){
	    FDC2_Chi2U = tempChi2U;
	    FDC2_Upos = tempU;
	    slopeU = aU[1];
	    interseptU = aU[0];
	  }
	  
	  //cout << "Chi2 " << FDC2_Chi2U << "U " << FDC2_Upos <<  endl;
	  //cout << "slope " << slope << " const " << intersept << endl;
	  
	  
	}else continue;

	
	if(FDC2_trackV[i][0]>-1200&&FDC2_trackV[i][0]<1200){
	  
	  //===== Least square =====
	  Double_t aV[2] = {-9999,-9999};
	  Double_t sV[3] = {0};
	  Double_t tV[2] = {0};
	  
	  for(Int_t l=0;l<4;++l){
	    sV[0] += 1.;
	    sV[1] += FDC2_trackV_z[l];
	    sV[2] += pow(FDC2_trackV_z[l],2);
	    tV[0] += FDC2_trackV[i][l];
	  tV[1] += FDC2_trackV[i][l]*FDC2_trackV_z[l];
	  }
	  aV[0] = (sV[2]*tV[0]-sV[1]*tV[1])/(sV[0]*sV[2]-sV[1]*sV[1]);
	  aV[1] = (sV[0]*tV[1]-sV[1]*tV[0])/(sV[0]*sV[2]-sV[1]*sV[1]);
	  
	  tempV = aV[0];
	  
	//cout << aV[0] << " " << aV[1] << endl;
	  
	  for(Int_t l=0;l<4;++l){
	    tempChi2V += pow(FDC2_trackV[i][l] - (aV[1]*FDC2_trackV_z[l] + aV[0]),2);
	  }
	  
	  //cout << tempChi2V << endl;
	  
	  if(tempChi2V<FDC2_Chi2V){
	    FDC2_Chi2V = tempChi2V;
	    FDC2_Vpos = tempV;
	    slopeV = aV[1];
	    interseptV = aV[0];
	  }
	  
	  //cout << "Chi2 " << FDC2_Chi2V << "V " << FDC2_Vpos <<  endl;
	  //cout << "slope " << slope << " const " << intersept << endl;
	  
	  
	}else continue;
	
	
	
      }//i = right left 2^4 UV

      //cout << "slope " << slope << " intersept " << intersept << endl;

    }//tr
    FDC2_X = (FDC2_Xpos + (FDC2_Vpos + FDC2_Upos)/2.)/2.;
    FDC2_Y = TMath::Sqrt(3)/2.*(FDC2_Vpos-FDC2_Upos);
    FDC2_A = (6*TMath::ATan(slopeX) + 2/TMath::Sqrt(3)*4*TMath::ATan(slopeU) + 2/TMath::Sqrt(3)*4*TMath::ATan(slopeV))/14.;
    FDC2_B = (-TMath::ATan(slopeU) + TMath::ATan(slopeV))/2.;
    //@@@ FDC2 end @@@

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
}//ana_smiri()
