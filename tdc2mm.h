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

void LoadDCTDCDist(char *FileName){

  TFile *RootFile = new TFile(FileName,"READ");
  if(RootFile){
    gROOT->cd();
    TH1I *h = NULL;
    Int_t BDCNumLayer = 8;
    Int_t FDCNumLayer = 14;
    Int_t tdcwindow[2];
    tdcwindow[0] = 1000;
    tdcwindow[1] = 2000;
    Double_t tdc2mm[14][1000];
    Int_t tdcint[14];
    
    for(Int_t l=0;l<BDCNumLayer;++l){
      h = (TH1I*)RootFile->Get(Form("hbdc1tdc%d",l));
      
      if(h){
	tdcint[l] = (Double_t)h->Integral(h->FindBin(tdcwindow[0]),h->FindBin(tdcwindow[1]));
	for(Int_t i=0;i<1000;++i){
	  tdc2mm[l][i] = 2.5*(Double_t)h->Integral(h->FindBin(i),h->FindBin(tdcwindow[1]))/tdcint[l];
	}	  	
      }else{
	cout << "\e[35m" << "No histogram found." << "\e[37m" << endl;
      }
    }else{
      cout << "\e[35m" << "No RootFile: " << FileName << "found." << "\e[37m" << endl;
    }
    
  }
}
