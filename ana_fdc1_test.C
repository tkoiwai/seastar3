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

void ana_fdc1_test(){

  //===== Load input file =====
  TFile *infile_dc   = TFile::Open("rootfiles/run0056/run0056_DC.root");
  
  TTree *caltreeDC;
  infile_dc->GetObject("caltreeDC",caltreeDC);

  //===== input tree variables =====
  Long64_t EventNumber = 0;
  Int_t RunNumber = -1;
  Int_t CoincidenceTrigger;

  Int_t FDC1_TDC[14][32], FDC1_TrailTDC[14][32];
  Int_t FDC1_WireID[14][32];
  Double_t FDC1_WirePosition[14][32], FDC1_WireZPosition[14][32];
  Int_t FDC1_Layer[14][32],FDC1_PlaneID[14][32], FDC1_HitID[14][32];
  
  Int_t NumFDC1Hit;
  
  //===== SetBranchAddress =====
  caltreeDC->SetBranchAddress("RunNumber",&RunNumber);
  caltreeDC->SetBranchAddress("EventNumber",&EventNumber);
  
  caltreeDC->SetBranchAddress("FDC1_TDC",FDC1_TDC);
  caltreeDC->SetBranchAddress("FDC1_TrailTDC",FDC1_TrailTDC);
  caltreeDC->SetBranchAddress("FDC1_WireID",FDC1_WireID);
  caltreeDC->SetBranchAddress("FDC1_WirePosition",FDC1_WirePosition);
  caltreeDC->SetBranchAddress("FDC1_WireZPosition",FDC1_WireZPosition);
  caltreeDC->SetBranchAddress("FDC1_Layer",FDC1_Layer);
  caltreeDC->SetBranchAddress("FDC1_PlaneID",FDC1_PlaneID);
  caltreeDC->SetBranchAddress("FDC1_HitID",FDC1_HitID);
  
  //===== Load CUT files =====
  
  //===== Create output file/tree =====
  TFile *anafile_fdc1_test = new TFile("rootfiles/ana_fdc1_test.root","RECREATE");
  TTree *atr  = new TTree("atr","atr");
  
  
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
  Int_t Width_FDC1 = 180.;

    //===== Declrear anatree const.s =====
  Int_t nohit_layerX, nohit_layerU, nohit_layerV;
  Int_t tmphitnumX, tmphitnumU, tmphitnumV;
  Int_t tr;

  //===== Define anatree variables =====
  Int_t RunNum, EventNum;
  
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
 
  Int_t BG_flag;
  
  //===== Create anatree Branch =====
  anatreeFDC1->Branch("RunNum",&RunNum);
  anatreeFDC1->Branch("EventNum",&EventNum);
  
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
  
  anatreeFDC1->Branch("BG_flag",&BG_flag);
  
  //===== Begin LOOP =====
  int nEntry = caltreeDC->GetEntries();
  for(int iEntry=0;iEntry<100;++iEntry){
    //for(int iEntry=0;iEntry<nEntry;++iEntry){
    
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

    vector<Double_t> hit_posU;
    vector<Double_t> hit_tdcU;

    for(Int_t l=0;l<14;++l){
      Int_t nhit = 0;
      for(Int_t i=0;i<32;++i){	
	if(FDC1_TDC[l][i]>0&&FDC1_WirePosition[l][i]>-200.){
	  hit_posU.push_back(FDC1_WirePosition[l][i]);
	  hit_tdcU.push_back(FDC1_TDC[l][i]);
	  ++nhit;
	}else continue;
      }//i
      FDC1_hit_num[l] = nhit;
      if(h==0){
	Int_t l6 = l%6;
	switch(l6){
	case 2: nohit_layerU++; break;
	case 3: nohit_layerU++; break;
	default: break;
	}//hctiws
      }//fi
    }//l
    /*
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
    */
    
    if(nohit_layerX>=3||nohit_layerU>=2||nohit_layerV>=2) BG_flag = 1; //if(nohit_layer<=1) BG_flag = 0

    Int_t lU = 4;
    Int_t liveU = lU - nohit_layerU;

    
    //U
    tr = 0;
    for(Int_t hl2=0;hl2<FDC1_hit_num[2]&&tr<16;++hl2){
      for(Int_t hl3=0;hl3<FDC1_hit_num[3]&&tr<16;++hl3){
	for(Int_t hl8=0;hl8<FDC1_hit_num[8]&&tr<16;++hl8){
	  for(Int_t hl9=0;hl9<FDC1_hit_num[9]+1&&tr<16;++hl9){
	    FDC1_trackU_pos[tr][0] = FDC1_hit_pos[2][hl2];
	    FDC1_trackU_pos[tr][1] = FDC1_hit_pos[3][hl3];
	    FDC1_trackU_pos[tr][2] = FDC1_hit_pos[8][hl8];
	    FDC1_trackU_pos[tr][3] = FDC1_hit_pos[9][hl9];
	    FDC1_trackU_mm[tr][0] = tdc2mm[2][2][FDC1_hit_tdc[2][hl2]];
	    FDC1_trackU_mm[tr][1] = tdc2mm[2][3][FDC1_hit_tdc[3][hl3]];
	    FDC1_trackU_mm[tr][2] = tdc2mm[2][8][FDC1_hit_tdc[8][hl8]];
	    FDC1_trackU_mm[tr][3] = tdc2mm[2][9][FDC1_hit_tdc[9][hl9]];

	    //cout << "U " << tr << " " << FDC1_trackU_pos[tr][0] << endl;

	    tr++;
	  }
	}
      }
    }

   
    for(tr=0;tr<16;++tr){
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
	Int_t de = 0;
	Int_t abc;
	for(abc=0;abc<4;++abc){
	  if(FDC1_trackU_pos[tr][abc]>-200&&FDC1_trackU_pos[tr][abc]<200){
	    /*tempU = (FDC1_trackU_pos[tr][0] + pow(-1,signUV[0])*FDC1_trackU_mm[tr][0] +
		     FDC1_trackU_pos[tr][1] + pow(-1,signUV[1])*FDC1_trackU_mm[tr][1] +
		     FDC1_trackU_pos[tr][2] + pow(-1,signUV[2])*FDC1_trackU_mm[tr][2] +
		     FDC1_trackU_pos[tr][3] + pow(-1,signUV[3])*FDC1_trackU_mm[tr][3])/4.;*/	    
	    tempU += FDC1_trackU_pos[tr][abc] + pow(-1,signUV[abc])*FDC1_trackU_mm[tr][abc];
	    cout << pow(-1,signUV[abc])*FDC1_trackU_mm[tr][abc] << " ";
	    ++de;
	  }else continue;
	}
	tempU = tempU/de;

	tempV = (FDC1_trackV_pos[tr][0] + pow(-1,signUV[0])*FDC1_trackV_mm[tr][0] +
		 FDC1_trackV_pos[tr][1] + pow(-1,signUV[1])*FDC1_trackV_mm[tr][1] +
		 FDC1_trackV_pos[tr][2] + pow(-1,signUV[2])*FDC1_trackV_mm[tr][2] +
		 FDC1_trackV_pos[tr][3] + pow(-1,signUV[3])*FDC1_trackV_mm[tr][3])/4.;
	for(Int_t j=0;j<4;++j){
	  if(FDC1_trackU_pos[tr][abc]>-200&&FDC1_trackU_pos[tr][abc]<200){
	    tempChi2U += pow(tempU - (FDC1_trackU_pos[tr][j] + pow(-1,signUV[j])*FDC1_trackU_mm[tr][j]),2);
	  }
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
