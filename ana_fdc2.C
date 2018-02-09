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

using namespace std;

void ana_fdc2(){

  //===== Load input file =====
  TFile *infile_dc   = TFile::Open("rootfiles/run0056/run0056_DC.root");
  
  TTree *caltreeDC;
  infile_dc->GetObject("caltreeDC",caltreeDC);

  //===== input tree variables =====
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;
  /*
  Int_t BDC1_TDC[8][16], BDC1_TrailTDC[8][16];
  Int_t BDC1_WireID[8][16];
  Double_t BDC1_WirePosition[8][16], BDC1_WireZPosition[8][16];
  Int_t BDC1_Layer[8][16],BDC1_PlaneID[8][16], BDC1_HitID[8][16];
  
  Int_t BDC2_TDC[8][16], BDC2_TrailTDC[8][16];
  Int_t BDC2_WireID[8][16];
  Double_t BDC2_WirePosition[8][16], BDC2_WireZPosition[8][16];
  Int_t BDC2_Layer[8][16],BDC2_PlaneID[8][16], BDC2_HitID[8][16];
  /*
  Int_t FDC1_TDC[14][32], FDC1_TrailTDC[14][32];
  Int_t FDC1_WireID[14][32];
  Double_t FDC1_WirePosition[14][32], FDC1_WireZPosition[14][32];
  Int_t FDC1_Layer[14][32],FDC1_PlaneID[14][32], FDC1_HitID[14][32];
  */
  Int_t FDC2_TDC[14][112], FDC2_TrailTDC[14][112];
  Int_t FDC2_WireID[14][112];
  Double_t FDC2_WirePosition[14][112], FDC2_WireZPosition[14][112];
  Int_t FDC2_Layer[14][112],FDC2_PlaneID[14][112], FDC2_HitID[14][112];
  
  Int_t NumBDC1Hit, NumBDC2Hit, NumFDC1Hit, NumFDC2Hit;
  
  //===== SetBranchAddress =====
  caltreeDC->SetBranchAddress("RunNumber",&RunNumber);
  caltreeDC->SetBranchAddress("EventNumber",&EventNumber);
  /*
  caltreeDC->SetBranchAddress("BDC1_TDC",BDC1_TDC);
  caltreeDC->SetBranchAddress("BDC1_TrailTDC",BDC1_TrailTDC);
  caltreeDC->SetBranchAddress("BDC1_WireID",BDC1_WireID);
  caltreeDC->SetBranchAddress("BDC1_WirePosition",BDC1_WirePosition);
  caltreeDC->SetBranchAddress("BDC1_WireZPosition",BDC1_WireZPosition);
  caltreeDC->SetBranchAddress("BDC1_Layer",BDC1_Layer);
  caltreeDC->SetBranchAddress("BDC1_PlaneID",BDC1_PlaneID);
  caltreeDC->SetBranchAddress("BDC1_HitID",BDC1_HitID);
  
  caltreeDC->SetBranchAddress("BDC2_TDC",BDC2_TDC);
  caltreeDC->SetBranchAddress("BDC2_TrailTDC",BDC2_TrailTDC);
  caltreeDC->SetBranchAddress("BDC2_WireID",BDC2_WireID);
  caltreeDC->SetBranchAddress("BDC2_WirePosition",BDC2_WirePosition);
  caltreeDC->SetBranchAddress("BDC2_WireZPosition",BDC2_WireZPosition);
  caltreeDC->SetBranchAddress("BDC2_Layer",BDC2_Layer);
  caltreeDC->SetBranchAddress("BDC2_PlaneID",BDC2_PlaneID);
  caltreeDC->SetBranchAddress("BDC2_HitID",BDC2_HitID);
  
  caltreeDC->SetBranchAddress("FDC1_TDC",FDC1_TDC);
  caltreeDC->SetBranchAddress("FDC1_TrailTDC",FDC1_TrailTDC);
  caltreeDC->SetBranchAddress("FDC1_WireID",FDC1_WireID);
  caltreeDC->SetBranchAddress("FDC1_WirePosition",FDC1_WirePosition);
  caltreeDC->SetBranchAddress("FDC1_WireZPosition",FDC1_WireZPosition);
  caltreeDC->SetBranchAddress("FDC1_Layer",FDC1_Layer);
  caltreeDC->SetBranchAddress("FDC1_PlaneID",FDC1_PlaneID);
  caltreeDC->SetBranchAddress("FDC1_HitID",FDC1_HitID);
  */
  caltreeDC->SetBranchAddress("FDC2_TDC",FDC2_TDC);
  caltreeDC->SetBranchAddress("FDC2_TrailTDC",FDC2_TrailTDC);
  caltreeDC->SetBranchAddress("FDC2_WireID",FDC2_WireID);
  caltreeDC->SetBranchAddress("FDC2_WirePosition",FDC2_WirePosition);
  caltreeDC->SetBranchAddress("FDC2_WireZPosition",FDC2_WireZPosition);
  caltreeDC->SetBranchAddress("FDC2_Layer",FDC2_Layer);
  caltreeDC->SetBranchAddress("FDC2_PlaneID",FDC2_PlaneID);
  caltreeDC->SetBranchAddress("FDC2_HitID",FDC2_HitID);
  
  //===== Load CUT files =====
  
  //===== Create output file/tree =====
  TFile *anafile_fdc2 = new TFile("rootfiles/ana_fdc2.root","RECREATE");
  TTree *anatreeFDC2  = new TTree("anatreeFDC2","anatreeFDC2");
  
  
  //===== Create TDC Distributions =====
  Int_t BDCNumLayer = 8;
  Int_t FDCNumLayer = 14;
  Int_t DCNum = 4; //BDC1, BDC2, FDC1 and FDC2
  Int_t tdcwindow[2] = {1000,2000};
  Double_t tdc2mm[DCNum][14][2000] = {};
  Int_t tdcint[DCNum][14] = {};
  Double_t wire_gap[DCNum] = {2.5,2.5,5,10};
  
  TFile *RootFile = new TFile("/home/koiwai/analysis/rootfiles/ana_dc_tdcdist.root","READ");
  if(RootFile){
    gROOT->cd();
    for(Int_t n=0;n<DCNum;++n){
      TH1I *h = NULL;
      if(n==0||n==1){
	for(Int_t l=0;l<BDCNumLayer;++l){
	  h = (TH1I*)RootFile->Get(Form("hbdc%dtdc%d",n+1,l));
	  
	  if(h){
	    tdcint[n][l] = (Double_t)h->Integral(h->FindBin(tdcwindow[0]),h->FindBin(tdcwindow[1]));
	    for(Int_t i=0;i<2000;++i){
	      if(i>=0&&i<=1000) tdc2mm[n][l][i] = 0;
	      else tdc2mm[n][l][i] = wire_gap[n]*(Double_t)h->Integral(h->FindBin(i),h->FindBin(tdcwindow[1]))/tdcint[n][l]; //2.5 only for BDCs & FDC1?	  
	    }	  	
	  }
	}
      }else if(n==2||n==3)
      for(Int_t l=0;l<FDCNumLayer;++l){
	h = (TH1I*)RootFile->Get(Form("hfdc%dtdc%d",n-1,l));
	
	if(h){
	  tdcint[n][l] = (Double_t)h->Integral(h->FindBin(tdcwindow[0]),h->FindBin(tdcwindow[1]));
	  for(Int_t i=0;i<2000;++i){
	    if(i>=0&&i<=1000) tdc2mm[n][l][i] = 0;
	    else tdc2mm[n][l][i] = wire_gap[n]*(Double_t)h->Integral(h->FindBin(i),h->FindBin(tdcwindow[1]))/tdcint[n][l]; // CHECK the WIRE_GAP!!!	  
	  }	  	
	}
      }
    }
  }
  
  //===== Declear const.s =====
  /*Int_t Dist_BDC1BDC2 = 998.; //mm
  Int_t Dist_BDC1Target = 2231.;
  Int_t Dist_BDC1FDC1 = 3533.;
  Int_t Dist_SBTTarget = 2795.;
  Int_t Width_BDC1 = 68.;*/
  Int_t Width_FDC1 = 180.;

    //===== Declrear anatree const.s =====
  Int_t nohit_layerX, nohit_layerU, nohit_layerV;
  Int_t tmphitnumX, tmphitnumU, tmphitnumV;
  Int_t tr;

  //===== Define anatree variables =====
  Int_t RunNum, EventNum;
  /*
  Double_t BDC1_WirePos[8][16];
  Int_t BDC1_WireTDC[8][16];

  Int_t BDC1_hit_tdc[8][16];
  Double_t BDC1_hit_pos[8][16];
  Int_t BDC1_hit_num[8];
  Double_t BDC1_trackX_pos[16][4], BDC1_trackY_pos[16][4];
  Double_t BDC1_trackX_mm[16][4], BDC1_trackY_mm[16][4];

  Double_t BDC1_X, BDC1_Y, BDC1_Chi2X, BDC1_Chi2Y;

  Int_t BDC1_allhitnumX;
  Int_t BDC1_allhitnumY;

  Double_t BDC2_WirePos[8][16];
  Int_t BDC2_WireTDC[8][16];

  Int_t BDC2_hit_tdc[8][16];
  Double_t BDC2_hit_pos[8][16];
  Int_t BDC2_hit_num[8];
  Double_t BDC2_trackX_pos[16][4], BDC2_trackY_pos[16][4];
  Double_t BDC2_trackX_mm[16][4], BDC2_trackY_mm[16][4];

  Double_t BDC2_X, BDC2_Y, BDC2_Chi2X, BDC2_Chi2Y;

  Int_t BDC2_allhitnumX;
  Int_t BDC2_allhitnumY;

  Double_t BDC_X, BDC_Y;
  Double_t BDC_A, BDC_B;
  Double_t Target_X, Target_Y;
  
  Double_t FDC1_WirePos[14][32];
  Int_t FDC1_WireTDC[14][32];

  Int_t FDC1_hit_tdc[14][32];
  Double_t FDC1_hit_pos[14][32];
  Int_t FDC1_hit_num[14];
  Double_t FDC1_trackX_pos[16][6], FDC1_trackU_pos[16][4], FDC1_trackV_pos[16][4];
  Double_t FDC1_trackX_mm[16][6], FDC1_trackU_mm[16][4], FDC1_trackV_mm[16][4];

  Double_t FDC1_Xpos, FDC1_Upos, FDC1_Vpos, FDC1_Chi2X, FDC1_Chi2U, FDC1_Chi2V;

  Double_t FDC1_X, FDC1_Y;

  Int_t FDC1_allhitnumX;
  Int_t FDC1_allhitnumU;
  Int_t FDC1_allhitnumV;  
  */
  Double_t FDC2_WirePos[14][112];
  Int_t FDC2_WireTDC[14][112];

  Int_t FDC2_hit_tdc[14][112];
  Double_t FDC2_hit_pos[14][112];
  Int_t FDC2_hit_num[14];
  Double_t FDC2_trackX_pos[16][6], FDC2_trackU_pos[16][4], FDC2_trackV_pos[16][4];
  Double_t FDC2_trackX_mm[16][6], FDC2_trackU_mm[16][4], FDC2_trackV_mm[16][4];
  Double_t FDC2_trackX[64][6], FDC2_trackU[64][4], FDC2_trackV[64][4];
  Double_t FDC2_trackX_z[6] = {-308.66, -291.34, -8.66, 8.66, 291.34, 308.66};
  Double_t FDC2_trackU_z[6] = {-208.66, -191.34, 91.34, 108.66};
  Double_t FDC2_trackV_z[6] = {-108.66, -91.34, 191.34, 208.66};

  Double_t FDC2_Xpos, FDC2_Upos, FDC2_Vpos, FDC2_Chi2X, FDC2_Chi2U, FDC2_Chi2V;

  Double_t FDC2_X, FDC2_Y, FDC2_A, FDC2_B;

  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatreeFDC2->Branch("RunNum",&RunNum);
  anatreeFDC2->Branch("EventNum",&EventNum);
  /*
  anatreeFDC1->Branch("BDC1_WirePos",BDC1_WirePos,"BDC1_WirePos[8][16]/D");
  anatreeFDC1->Branch("BDC1_WireTDC",BDC1_WireTDC,"BDC1_WireTDC[8][16]/I");

  anatreeFDC1->Branch("BDC1_hit_tdc",BDC1_hit_tdc,"BDC1_hit_tdc[8][16]/I");
  anatreeFDC1->Branch("BDC1_hit_pos",BDC1_hit_pos,"BDC1_hit_pos[8][16]/D");
  anatreeFDC1->Branch("BDC1_hit_num",BDC1_hit_num,"BDC1_hit_num[8]/I");
  
  anatreeFDC1->Branch("BDC1_trackX_pos",BDC1_trackX_pos,"BDC1_trackX_pos[16][4]/D");
  anatreeFDC1->Branch("BDC1_trackX_mm",BDC1_trackX_mm,"BDC1_trackX_mm[16][4]/D");

  anatreeFDC1->Branch("BDC1_X",&BDC1_X);
  anatreeFDC1->Branch("BDC1_Y",&BDC1_Y);
  anatreeFDC1->Branch("BDC1_Chi2X",&BDC1_Chi2X);
  anatreeFDC1->Branch("BDC1_Chi2Y",&BDC1_Chi2Y);

  anatreeFDC1->Branch("BDC1_allhitnumX",&BDC1_allhitnumX);
  anatreeFDC1->Branch("BDC1_allhitnumY",&BDC1_allhitnumY);

  anatreeFDC1->Branch("BDC2_WirePos",BDC2_WirePos,"BDC2_WirePos[8][16]/D");
  anatreeFDC1->Branch("BDC2_WireTDC",BDC2_WireTDC,"BDC2_WireTDC[8][16]/I");

  anatreeFDC1->Branch("BDC2_hit_tdc",BDC2_hit_tdc,"BDC2_hit_tdc[8][16]/I");
  anatreeFDC1->Branch("BDC2_hit_pos",BDC2_hit_pos,"BDC2_hit_pos[8][16]/D");
  anatreeFDC1->Branch("BDC2_hit_num",BDC2_hit_num,"BDC2_hit_num[8]/I");
  
  anatreeFDC1->Branch("BDC2_trackX_pos",BDC2_trackX_pos,"BDC2_trackX_pos[16][4]/D");
  anatreeFDC1->Branch("BDC2_trackX_mm",BDC2_trackX_mm,"BDC2_trackX_mm[16][4]/D");

  anatreeFDC1->Branch("BDC2_X",&BDC2_X);
  anatreeFDC1->Branch("BDC2_Y",&BDC2_Y);
  anatreeFDC1->Branch("BDC2_Chi2X",&BDC2_Chi2X);
  anatreeFDC1->Branch("BDC2_Chi2Y",&BDC2_Chi2Y);

  anatreeFDC1->Branch("BDC2_allhitnumX",&BDC2_allhitnumX);
  anatreeFDC1->Branch("BDC2_allhitnumY",&BDC2_allhitnumY);

  anatreeFDC1->Branch("BDC_X",&BDC_X);
  anatreeFDC1->Branch("BDC_Y",&BDC_Y);
  anatreeFDC1->Branch("BDC_A",&BDC_A);
  anatreeFDC1->Branch("BDC_B",&BDC_B);

  anatreeFDC1->Branch("Target_X",&Target_X);
  anatreeFDC1->Branch("Target_Y",&Target_Y);
  

  anatreeFDC1->Branch("FDC1_WirePos",FDC1_WirePos,"FDC1_WirePos[14][32]/D");
  anatreeFDC1->Branch("FDC1_WireTDC",FDC1_WireTDC,"FDC1_WireTDC[14][32]/I");

  anatreeFDC1->Branch("FDC1_hit_tdc",FDC1_hit_tdc,"FDC1_hit_tdc[14][32]/I");
  anatreeFDC1->Branch("FDC1_hit_pos",FDC1_hit_pos,"FDC1_hit_pos[14][32]/D");
  anatreeFDC1->Branch("FDC1_hit_num",FDC1_hit_num,"FDC1_hit_num[14]/I");
  
  anatreeFDC1->Branch("FDC1_trackX_pos",FDC1_trackX_pos,"FDC1_trackX_pos[16][6]/D");
  anatreeFDC1->Branch("FDC1_trackX_mm",FDC1_trackX_mm,"FDC1_trackX_mm[16][6]/D");
  anatreeFDC1->Branch("FDC1_trackU_pos",FDC1_trackU_pos,"FDC1_trackU_pos[16][4]/D");
  anatreeFDC1->Branch("FDC1_trackU_mm",FDC1_trackU_mm,"FDC1_trackU_mm[16][4]/D");
  anatreeFDC1->Branch("FDC1_trackV_pos",FDC1_trackV_pos,"FDC1_trackV_pos[16][4]/D");
  anatreeFDC1->Branch("FDC1_trackV_mm",FDC1_trackV_mm,"FDC1_trackV_mm[16][4]/D");

  
  anatreeFDC1->Branch("FDC1_Xpos",&FDC1_Xpos);
  anatreeFDC1->Branch("FDC1_Upos",&FDC1_Upos);
  anatreeFDC1->Branch("FDC1_Vpos",&FDC1_Vpos);
  anatreeFDC1->Branch("FDC1_Chi2X",&FDC1_Chi2X);
  anatreeFDC1->Branch("FDC1_Chi2U",&FDC1_Chi2U);
  anatreeFDC1->Branch("FDC1_Chi2V",&FDC1_Chi2V);

  anatreeFDC1->Branch("FDC1_X",&FDC1_X);
  anatreeFDC1->Branch("FDC1_Y",&FDC1_Y);

  anatreeFDC1->Branch("FDC1_allhitnumX",&FDC1_allhitnumX);
  anatreeFDC1->Branch("FDC1_allhitnumU",&FDC1_allhitnumU);
  anatreeFDC1->Branch("FDC1_allhitnumV",&FDC1_allhitnumV);
  */
  anatreeFDC2->Branch("FDC2_WirePos",FDC2_WirePos,"FDC2_WirePos[14][112]/D");
  anatreeFDC2->Branch("FDC2_WireTDC",FDC2_WireTDC,"FDC2_WireTDC[14][112]/I");

  anatreeFDC2->Branch("FDC2_hit_tdc",FDC2_hit_tdc,"FDC2_hit_tdc[14][112]/I");
  anatreeFDC2->Branch("FDC2_hit_pos",FDC2_hit_pos,"FDC2_hit_pos[14][112]/D");
  anatreeFDC2->Branch("FDC2_hit_num",FDC2_hit_num,"FDC2_hit_num[14]/I");
  
  anatreeFDC2->Branch("FDC2_trackX_pos",FDC2_trackX_pos,"FDC2_trackX_pos[16][6]/D");
  anatreeFDC2->Branch("FDC2_trackX_mm",FDC2_trackX_mm,"FDC2_trackX_mm[16][6/D");
  anatreeFDC2->Branch("FDC2_trackX_z",FDC2_trackX_z,"FDC2_trackX_z[6]/D");
  anatreeFDC2->Branch("FDC2_trackX",FDC2_trackX,"FDC2_trackX[64][6]/D");

  anatreeFDC2->Branch("FDC2_trackU_pos",FDC2_trackU_pos,"FDC2_trackU_pos[16][4]/D");
  anatreeFDC2->Branch("FDC2_trackU_mm",FDC2_trackU_mm,"FDC2_trackU_mm[16][4]/D");
  anatreeFDC2->Branch("FDC2_trackU_z",FDC2_trackU_z,"FDC2_trackU_z[6]/D");
  anatreeFDC2->Branch("FDC2_trackU",FDC2_trackU,"FDC2_trackU[64][6]/D");

  anatreeFDC2->Branch("FDC2_trackV_pos",FDC2_trackV_pos,"FDC2_trackV_pos[16][4]/D");
  anatreeFDC2->Branch("FDC2_trackV_mm",FDC2_trackV_mm,"FDC2_trackV_mm[16][4]/D");
  anatreeFDC2->Branch("FDC2_trackV_z",FDC2_trackV_z,"FDC2_trackV_z[6]/D");
  anatreeFDC2->Branch("FDC2_trackV",FDC2_trackV,"FDC2_trackV[64][6]/D");

  anatreeFDC2->Branch("FDC2_Xpos",&FDC2_Xpos);
  anatreeFDC2->Branch("FDC2_Upos",&FDC2_Upos);
  anatreeFDC2->Branch("FDC2_Vpos",&FDC2_Vpos);
  anatreeFDC2->Branch("FDC2_Chi2X",&FDC2_Chi2X);
  anatreeFDC2->Branch("FDC2_Chi2U",&FDC2_Chi2U);
  anatreeFDC2->Branch("FDC2_Chi2V",&FDC2_Chi2V);

  anatreeFDC2->Branch("FDC2_X",&FDC2_X);
  anatreeFDC2->Branch("FDC2_Y",&FDC2_Y);
  anatreeFDC2->Branch("FDC2_A",&FDC2_A);
  anatreeFDC2->Branch("FDC2_B",&FDC2_B);

  anatreeFDC2->Branch("BG_flag",&BG_flag);
  
  //===== Begin LOOP =====
  int nEntry = caltreeDC->GetEntries();
  //for(int iEntry=0;iEntry<100;++iEntry){
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    
    if(iEntry%100 == 0){
      clog<< iEntry/1000 << "k events treated..." << "\r";
    }

    RunNum = RunNumber;
    EventNum = EventNumber;
    
    caltreeDC->GetEntry(iEntry);

    //@@@ FDC1 @@@
    //=== Initialization ===    
    for(Int_t l=0;l<14;++l){
      for(Int_t i=0;i<112;++i){
	FDC2_WirePos[l][i] = TMath::Sqrt(-1);
	FDC2_WirePos[l][i] = FDC2_WirePosition[l][i];
	FDC2_WireTDC[l][i] = -9999;
	FDC2_WireTDC[l][i] = FDC2_TDC[l][i];

	//FDC2_hit_tdc[l][i] = TMath::Sqrt(-1);
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
    
    for(tr=0;tr<16;++tr){

      for(Int_t i=0;i<64;++i){
	for(Int_t l=0;l<6;++l) FDC2_trackX[i][l] = TMath::Sqrt(-1);
	for(Int_t l=0;l<4;++l){
	  FDC2_trackU[i][l] = TMath::Sqrt(-1);
	  FDC2_trackV[i][l] = TMath::Sqrt(-1);      
      }
      }

      Double_t slope, intersept;

      
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
	  slope = aX[1];
	  intersept = aX[0];
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
	    slope = aU[1];
	    intersept = aU[0];
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
	    slope = aV[1];
	    intersept = aV[0];
	  }
	  
	  //cout << "Chi2 " << FDC2_Chi2V << "V " << FDC2_Vpos <<  endl;
	  //cout << "slope " << slope << " const " << intersept << endl;
	  
	  
	}else continue;
	
	
	
      }//i = right left 2^4 UV

      //cout << "slope " << slope << " intersept " << intersept << endl;

    }//tr
    //FDC2_X = (FDC2_Xpos + (FDC2_Vpos + FDC2_Upos)/2.)/2.;
    //FDC2_Y = TMath::Sqrt(3)/2.*(FDC2_Vpos-FDC2_Upos);
    //@@@ FDC2 end @@@

    
    
    anatreeFDC2->Fill();
  }//for LOOP
  anafile_fdc2->cd();
  anatreeFDC2->Write();
  anafile_fdc2->Close();
}//ana_dc()
