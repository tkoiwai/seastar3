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

void ana_fdc1(){

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
  */
  Int_t FDC1_TDC[14][32], FDC1_TrailTDC[14][32];
  Int_t FDC1_WireID[14][32];
  Double_t FDC1_WirePosition[14][32], FDC1_WireZPosition[14][32];
  Int_t FDC1_Layer[14][32],FDC1_PlaneID[14][32], FDC1_HitID[14][32];
  /*
  Int_t FDC2_TDC[14][112], FDC2_TrailTDC[14][112];
  Int_t FDC2_WireID[14][112];
  Double_t FDC2_WirePosition[14][112], FDC2_WireZPosition[14][112];
  Int_t FDC2_Layer[14][112],FDC2_PlaneID[14][112], FDC2_HitID[14][112];
  */
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
  */
  caltreeDC->SetBranchAddress("FDC1_TDC",FDC1_TDC);
  caltreeDC->SetBranchAddress("FDC1_TrailTDC",FDC1_TrailTDC);
  caltreeDC->SetBranchAddress("FDC1_WireID",FDC1_WireID);
  caltreeDC->SetBranchAddress("FDC1_WirePosition",FDC1_WirePosition);
  caltreeDC->SetBranchAddress("FDC1_WireZPosition",FDC1_WireZPosition);
  caltreeDC->SetBranchAddress("FDC1_Layer",FDC1_Layer);
  caltreeDC->SetBranchAddress("FDC1_PlaneID",FDC1_PlaneID);
  caltreeDC->SetBranchAddress("FDC1_HitID",FDC1_HitID);
  /*
  caltreeDC->SetBranchAddress("FDC2_TDC",FDC2_TDC);
  caltreeDC->SetBranchAddress("FDC2_TrailTDC",FDC2_TrailTDC);
  caltreeDC->SetBranchAddress("FDC2_WireID",FDC2_WireID);
  caltreeDC->SetBranchAddress("FDC2_WirePosition",FDC2_WirePosition);
  caltreeDC->SetBranchAddress("FDC2_WireZPosition",FDC2_WireZPosition);
  caltreeDC->SetBranchAddress("FDC2_Layer",FDC2_Layer);
  caltreeDC->SetBranchAddress("FDC2_PlaneID",FDC2_PlaneID);
  caltreeDC->SetBranchAddress("FDC2_HitID",FDC2_HitID);
  */
  //===== Load CUT files =====
  
  //===== Create output file/tree =====
  TFile *anafile_fdc1 = new TFile("rootfiles/ana_fdc1.root","RECREATE");
  TTree *anatreeFDC1  = new TTree("anatreeFDC1","anatreeFDC1");
  
  
  //===== Create TDC Distributions =====
  Int_t BDCNumLayer = 8;
  Int_t FDCNumLayer = 14;
  Int_t DCNum = 3; //BDC1, BDC2 and FDC1 currentry
  Int_t tdcwindow[2] = {1000,2000};
  Double_t tdc2mm[DCNum][14][2000] = {};
  Int_t tdcint[DCNum][14] = {};
  Double_t wire_gap[DCNum] = {2.5,2.5,5};
  
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
  */
  Double_t FDC1_WirePos[14][32];
  Int_t FDC1_WireTDC[14][32];

  Int_t FDC1_hit_tdc[14][32];
  Double_t FDC1_hit_pos[14][32];
  Int_t FDC1_hit_num[14];
  Double_t FDC1_trackX_pos[16][6], FDC1_trackU_pos[16][4], FDC1_trackV_pos[16][4];
  Double_t FDC1_trackX_mm[16][6], FDC1_trackU_mm[16][4], FDC1_trackV_mm[16][4];

  Double_t FDC1_X, FDC1_U, FDC1_V, FDC1_Chi2X, FDC1_Chi2U, FDC1_Chi2V;

  Int_t FDC1_allhitnumX;
  Int_t FDC1_allhitnumU;
  Int_t FDC1_allhitnumV;  
  /*
  Double_t FDC2_WirePos[14][32];
  Int_t FDC2_WireTDC[14][32];

  Int_t FDC2_hit_tdc[14][32];
  Double_t FDC2_hit_pos[14][32];
  Int_t FDC2_hit_num[14];
  Double_t FDC2_trackX_pos[32][4], FDC2_trackY_pos[32][4];
  Double_t FDC2_trackX_mm[32][4], FDC2_trackY_mm[32][4];

  Double_t FDC2_X, FDC2_Y, FDC2_Chi2X, FDC2_Chi2Y;

  Int_t FDC2_allhitnumX;
  Int_t FDC2_allhitnumY;

  Double_t FDC_X, FDC_Y;
  Double_t FDC_A, FDC_B;
  */
  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatreeFDC1->Branch("RunNum",&RunNum);
  anatreeFDC1->Branch("EventNum",&EventNum);
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
  */

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

  
  anatreeFDC1->Branch("FDC1_X",&FDC1_X);
  anatreeFDC1->Branch("FDC1_U",&FDC1_U);
  anatreeFDC1->Branch("FDC1_V",&FDC1_V);
  anatreeFDC1->Branch("FDC1_Chi2X",&FDC1_Chi2X);
  anatreeFDC1->Branch("FDC1_Chi2U",&FDC1_Chi2U);
  anatreeFDC1->Branch("FDC1_Chi2V",&FDC1_Chi2V);

  anatreeFDC1->Branch("FDC1_allhitnumX",&FDC1_allhitnumX);
  anatreeFDC1->Branch("FDC1_allhitnumU",&FDC1_allhitnumU);
  anatreeFDC1->Branch("FDC1_allhitnumV",&FDC1_allhitnumV);
  /*
  anatreeFDC1->Branch("FDC2_WirePos",FDC2_WirePos,"FDC2_WirePos[14][32]/D");
  anatreeFDC1->Branch("FDC2_WireTDC",FDC2_WireTDC,"FDC2_WireTDC[14][32]/I");

  anatreeFDC1->Branch("FDC2_hit_tdc",FDC2_hit_tdc,"FDC2_hit_tdc[14][32]/I");
  anatreeFDC1->Branch("FDC2_hit_pos",FDC2_hit_pos,"FDC2_hit_pos[14][32]/D");
  anatreeFDC1->Branch("FDC2_hit_num",FDC2_hit_num,"FDC2_hit_num[14]/I");
  
  anatreeFDC1->Branch("FDC2_trackX_pos",FDC2_trackX_pos,"FDC2_trackX_pos[32][4]/D");
  anatreeFDC1->Branch("FDC2_trackX_mm",FDC2_trackX_mm,"FDC2_trackX_mm[32][4]/D");

  anatreeFDC1->Branch("FDC2_X",&FDC2_X);
  anatreeFDC1->Branch("FDC2_Y",&FDC2_Y);
  anatreeFDC1->Branch("FDC2_Chi2X",&FDC2_Chi2X);
  anatreeFDC1->Branch("FDC2_Chi2Y",&FDC2_Chi2Y);

  anatreeFDC1->Branch("FDC2_allhitnumX",&FDC2_allhitnumX);
  anatreeFDC1->Branch("FDC2_allhitnumY",&FDC2_allhitnumY);

  anatreeFDC1->Branch("FDC_X",&FDC_X);
  anatreeFDC1->Branch("FDC_Y",&FDC_Y);
  anatreeFDC1->Branch("FDC_A",&FDC_A);
  anatreeFDC1->Branch("FDC_B",&FDC_B);
  */
  
  anatreeFDC1->Branch("BG_flag",&BG_flag);
  
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
      for(Int_t i=0;i<32;++i){
	FDC1_WirePos[l][i] = TMath::Sqrt(-1);
	FDC1_WirePos[l][i] = FDC1_WirePosition[l][i];
	FDC1_WireTDC[l][i] = -9999;
	FDC1_WireTDC[l][i] = FDC1_TDC[l][i];

	//FDC1_hit_tdc[l][i] = TMath::Sqrt(-1);
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
    FDC1_X = TMath::Sqrt(-1);
    FDC1_U = TMath::Sqrt(-1);
    FDC1_V = TMath::Sqrt(-1);
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
    
    for(tr=0;tr<16;++tr){
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
	  FDC1_X = tempX;
	}
      }

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
	//cout << i << " " << "bitUV = {" << bitUV[0] << bitUV[1] << bitUV[2] << bitUV[3] << "}, signUV = {" << signUV[0] << signUV[1] << signUV[2] << signUV[3] << "}" << endl;
	//Int_t de = 0;
	//Int_t abc;
	//for(abc=0;abc<4;++abc){
	//if(FDC1_trackU_pos[tr][abc]>-200&&FDC1_trackU_pos[tr][abc]<200){
	tempU = (FDC1_trackU_pos[tr][0] + pow(-1,signUV[0])*FDC1_trackU_mm[tr][0] +
		 FDC1_trackU_pos[tr][1] + pow(-1,signUV[1])*FDC1_trackU_mm[tr][1] +
		 FDC1_trackU_pos[tr][2] + pow(-1,signUV[2])*FDC1_trackU_mm[tr][2]
		 /*FDC1_trackU_pos[tr][3] + pow(-1,signUV[3])*FDC1_trackU_mm[tr][3]*/)/3.;	    
	//tempU += FDC1_trackU_pos[tr][abc] + pow(-1,signUV[abc])*FDC1_trackU_mm[tr][abc];
	//cout << pow(-1,signUV[abc])*FDC1_trackU_mm[tr][abc] << " ";
	//++de;
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
	  FDC1_U = tempU;
	}
	if(tempChi2V<FDC1_Chi2V){
	  FDC1_Chi2V = tempChi2V;
	  FDC1_V = tempV;
	}
      }
    }
    //cout << FDC1_U << " ";
      //@@@ FDC1 end @@@
    anatreeFDC1->Fill();
  }//for LOOP
  anafile_fdc1->cd();
  anatreeFDC1->Write();
  anafile_fdc1->Close();
}//ana_dc()
