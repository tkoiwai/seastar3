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

void ana_smri(int argc, char *argv[]){

  Int_t RunNumber = TString(argv[optind]).Atoi();
  
  //===== Load input file =====
  TFile *infile_dc   = TFile::Open(Form("rootfiles/run00%d/run00%d_ALL.root"));
  
  TTree *caltreeDC;
  infile_dc->GetObject("caltreeDC",caltreeDC);

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
  caltreeDC->SetBranchAddress("RunNumber",&RunNumber);
  caltreeDC->SetBranchAddress("EventNumber",&EventNumber);

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
  
  caltreeDC->SetBranchAddress("FDC2_TDC",FDC2_TDC);
  caltreeDC->SetBranchAddress("FDC2_TrailTDC",FDC2_TrailTDC);
  caltreeDC->SetBranchAddress("FDC2_WireID",FDC2_WireID);
  caltreeDC->SetBranchAddress("FDC2_WirePosition",FDC2_WirePosition);
  caltreeDC->SetBranchAddress("FDC2_WireZPosition",FDC2_WireZPosition);
  caltreeDC->SetBranchAddress("FDC2_Layer",FDC2_Layer);
  caltreeDC->SetBranchAddress("FDC2_PlaneID",FDC2_PlaneID);
  caltreeDC->SetBranchAddress("FDC2_HitID",FDC2_HitID);
  
  //===== Create output file/tree =====
  TFile *anafile_smri = new TFile("rootfiles/ana_smri.root","RECREATE");
  TTree *anatrS  = new TTree("anatrS","anatrS");
  
  
  //===== Create TDC Distributions =====
  Int_t BDCNumLayer = 8;
  Int_t FDCNumLayer = 14;
  Int_t DCLayerNum = {8,8,14,14};
  Int_t DCNum = 4;
  Int_t tdcwindow[2] = {1000,2000};
  Double_t tdc2mm[DCNum][14][2000] = {};
  Int_t tdcint[DCNum][14] = {};
  Double_t wire_gap[DCNum] = {2.5, 2.5, 5, 10};
  char DCName[DCNum] = {b,b,f,f};
  
  TFile *RootFile = new TFile("/home/koiwai/analysis/rootfiles/ana_dc_tdcdist.root","READ");
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
  
  //===== Declear const.s =====
  Int_t Dist_BDC1BDC2 = 998.; //mm
  Int_t Dist_BDC1Target = 2231.;
  Int_t Dist_BDC1FDC1 = 3533.;
  Int_t Dist_SBTTarget = 2795.;
  Int_t Width_BDC1 = 68.;
  Int_t Width_FDC1 = 180.;

    //===== Declrear anatree const.s =====
  Int_t nohit_layerX, nohit_layerY;
  Int_t tmphitnumX, tmphitnumY;
  Int_t tr;

  //===== Define anatree variables =====
  Int_t RunNum, EventNum;
  
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
  
  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatreeDC->Branch("RunNum",&RunNum);
  anatreeDC->Branch("EventNum",&EventNum);
  
  anatreeDC->Branch("BDC1_WirePos",BDC1_WirePos,"BDC1_WirePos[8][16]/D");
  anatreeDC->Branch("BDC1_WireTDC",BDC1_WireTDC,"BDC1_WireTDC[8][16]/I");

  anatreeDC->Branch("BDC1_hit_tdc",BDC1_hit_tdc,"BDC1_hit_tdc[8][16]/I");
  anatreeDC->Branch("BDC1_hit_pos",BDC1_hit_pos,"BDC1_hit_pos[8][16]/D");
  anatreeDC->Branch("BDC1_hit_num",BDC1_hit_num,"BDC1_hit_num[8]/I");
  
  anatreeDC->Branch("BDC1_trackX_pos",BDC1_trackX_pos,"BDC1_trackX_pos[16][4]/D");
  anatreeDC->Branch("BDC1_trackX_mm",BDC1_trackX_mm,"BDC1_trackX_mm[16][4]/D");

  anatreeDC->Branch("BDC1_X",&BDC1_X);
  anatreeDC->Branch("BDC1_Y",&BDC1_Y);
  anatreeDC->Branch("BDC1_Chi2X",&BDC1_Chi2X);
  anatreeDC->Branch("BDC1_Chi2Y",&BDC1_Chi2Y);

  anatreeDC->Branch("BDC1_allhitnumX",&BDC1_allhitnumX);
  anatreeDC->Branch("BDC1_allhitnumY",&BDC1_allhitnumY);

  anatreeDC->Branch("BDC2_WirePos",BDC2_WirePos,"BDC2_WirePos[8][16]/D");
  anatreeDC->Branch("BDC2_WireTDC",BDC2_WireTDC,"BDC2_WireTDC[8][16]/I");

  anatreeDC->Branch("BDC2_hit_tdc",BDC2_hit_tdc,"BDC2_hit_tdc[8][16]/I");
  anatreeDC->Branch("BDC2_hit_pos",BDC2_hit_pos,"BDC2_hit_pos[8][16]/D");
  anatreeDC->Branch("BDC2_hit_num",BDC2_hit_num,"BDC2_hit_num[8]/I");
  
  anatreeDC->Branch("BDC2_trackX_pos",BDC2_trackX_pos,"BDC2_trackX_pos[16][4]/D");
  anatreeDC->Branch("BDC2_trackX_mm",BDC2_trackX_mm,"BDC2_trackX_mm[16][4]/D");

  anatreeDC->Branch("BDC2_X",&BDC2_X);
  anatreeDC->Branch("BDC2_Y",&BDC2_Y);
  anatreeDC->Branch("BDC2_Chi2X",&BDC2_Chi2X);
  anatreeDC->Branch("BDC2_Chi2Y",&BDC2_Chi2Y);

  anatreeDC->Branch("BDC2_allhitnumX",&BDC2_allhitnumX);
  anatreeDC->Branch("BDC2_allhitnumY",&BDC2_allhitnumY);

  anatreeDC->Branch("BDC_X",&BDC_X);
  anatreeDC->Branch("BDC_Y",&BDC_Y);
  anatreeDC->Branch("BDC_A",&BDC_A);
  anatreeDC->Branch("BDC_B",&BDC_B);

  anatreeDC->Branch("Target_X",&Target_X);
  anatreeDC->Branch("Target_Y",&Target_Y);

  anatreeDC->Branch("BG_flag",&BG_flag);
  
  //===== Begin LOOP =====
  int nEntry = caltreeDC->GetEntries();
  //for(int iEntry=0;iEntry<10;++iEntry){
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    
    if(iEntry%100 == 0){
      clog<< iEntry/1000 << "k events treated..." << "\r";
    }

    RunNum = RunNumber;
    EventNum = EventNumber;
    
    caltreeDC->GetEntry(iEntry);

    //@@@ BDC1 @@@
    //=== Initialization ===    
    for(Int_t l=0;l<8;++l){
      for(Int_t i=0;i<16;++i){
	BDC1_WirePos[l][i] = TMath::Sqrt(-1);
	BDC1_WirePos[l][i] = BDC1_WirePosition[l][i];
	BDC1_WireTDC[l][i] = -9999;
	BDC1_WireTDC[l][i] = BDC1_TDC[l][i];

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
	BDC2_WirePos[l][i] = TMath::Sqrt(-1);
	BDC2_WirePos[l][i] = BDC2_WirePosition[l][i];
	BDC2_WireTDC[l][i] = -9999;
	BDC2_WireTDC[l][i] = BDC2_TDC[l][i];

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
    
    anatrS->Fill();
  }//for LOOP
  anafile_smri->cd();
  anatrS->Write();
  anafile_smri->Close();
}//ana_smiri()
