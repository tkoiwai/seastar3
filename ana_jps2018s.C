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
#include"TVector3.h"
#include"TClass.h"
#include"TTreeIndex.h"
#include"TEnv.h"

using namespace std;
using namespace TMath;

int main(int argc, char *argv[]){

  TEnv *env = new TEnv("/home/koiwai/analysis/db/geometry_psp17.dat");

  Int_t FileNumber = TString(argv[1]).Atoi();
  
  //B===== Load input file =================================================
  TString FileName_beam = Form("/home/koiwai/analysis/anafiles/beam/ana_beam%04d.root",FileNumber);
  TFile *infile_beam = TFile::Open(FileName_beam);

  TTree *anatrB;
  infile_beam->GetObject("anatrB",anatrB);

  //B===== beam variables =====
  Int_t EventNumber, RunNumber;
  Double_t zetBR, aoqBR;
  Double_t betaF5F7;
  Double_t anaF7_Time;
  

  //B===== SetBranchAddress =====
  anatrB->SetBranchAddress("EventNumber",&EventNumber);
  anatrB->SetBranchAddress("RunNumber",&RunNumber);
  anatrB->SetBranchAddress("zetBR",&zetBR);
  anatrB->SetBranchAddress("aoqBR",&aoqBR);
  anatrB->SetBranchAddress("betaF5F7",&betaF5F7);
  anatrB->SetBranchAddress("anaF7_Time",&anaF7_Time);
  /*
  //D===== Load input file =====
  TString FileName_dali = Form("/home/koiwai/analysis/rootdali/run%04d_DALI.root",FileNumber);
  TFile *infile_dali = TFile::Open(FileName_dali);

  TTree *caltrD;
  infile_dali->GetObject("caltrD",caltrD);
  
  //D===== dali variables =====
  vector<Double_t> *DALI_Energy;
  //TBranch *bDALI_Energy;
  vector<Double_t> *DALI_CosTheta;
  //TBranch *bDALI_CosTheta;
  vector<Double_t> *DALI_Time;
  //TBranch *bDALI_Time;
  //vector<Int_t> *DALI_ID;
  //vector<Double_t> *DALI_Layer;
  //vector<Double_t> *DALI_X;
  //vector<Double_t> *DALI_Y;
  //vector<Double_t> *DALI_Z;
  //vector<TVector3> DALI_Pos;
  Int_t DALI_Multi;

  //D===== SetBranchAddress =====
  caltrD->SetBranchAddress("DALI_Energy",&DALI_Energy);
  caltrD->SetBranchAddress("DALI_CosTheta",&DALI_CosTheta);
  caltrD->SetBranchAddress("DALI_Time",&DALI_Time);
  //caltrD->SetBranchAddress("DALI_ID",&DALI_ID);
  //caltrD->SetBranchAddress("DALI_Layer",&DALI_Layer);
  //caltrD->SetBranchAddress("DALI_X",&DALI_X);
  //caltrD->SetBranchAddress("DALI_Y",&DALI_Y);
  //caltrD->SetBranchAddress("DALI_Z",&DALI_Z);
  //caltrD->SetBranchAddress("DALI_Pos",&DALI_Pos);
  caltrD->SetBranchAddress("DALI_Multi",&DALI_Multi);
  */
  //A===== Load input file =====
  TString FileName_all = Form("/home/koiwai/analysis/rootfiles/run%04d/run%04d_ALL.root",FileNumber,FileNumber);
  TFile *infile_all = TFile::Open(FileName_all);

  TTree *caltr;
  infile_all->GetObject("caltr",caltr);

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
  TString AnaFileName = Form("/home/koiwai/analysis/anafiles/jps/ana_jps%04d.root",FileNumber);
  TFile *anafile = new TFile(AnaFileName,"recreate");
  TTree *anatrJ = new TTree("anatrJ","anatrJ");

  //===== Dealare const.s =====
  Double_t distSBTMinos = env->GetValue("Dist_SBTTarget",0.0);
  Double_t distF7SBT = env->GetValue("Dist_F7SBT",0.0);
  Double_t offsetF7SBT = env->GetValue("offsetF7SBT",0.0);
  Double_t offsethodo = env->GetValue("offsethodo",0.0); //each bar

  //===== Declare ana variables =====
  Int_t ENum, RNum;
  /*
  vector<Double_t> dali_e;
  vector<Double_t> dali_cos;
  vector<Double_t> dali_t;
  */

  Double_t hodo_q, hodo_t, hodoi_q[24], hodoi_t[24];
  Int_t hodo_id;
  Double_t SBT_t;
  Double_t tofF7SBT, tofSBThodo, tofSBTMinos, tofMinoshodo;

  Double_t hodoi_Z[24];
  
  Int_t BR56Ca;
  //===== Set Branches ======
  anatrJ->Branch("ENum",&ENum);
  anatrJ->Branch("RNum",&RNum);
  /*
  anatrJ->Branch("dali_e",&dali_e);
  anatrJ->Branch("dali_cos",&dali_cos);
  anatrJ->Branch("dali_t",&dali_t);
  */

  anatrJ->Branch("zetBR",&zetBR);
  anatrJ->Branch("aoqBR",&aoqBR);
    
  anatrJ->Branch("hodo_q",&hodo_q);
  anatrJ->Branch("hodo_t",&hodo_t);
  for(Int_t i=0;i<24;++i){
    anatrJ->Branch(Form("hodo%02d_q",i+1),&hodoi_q[i]);
    anatrJ->Branch(Form("hodo%02d_t",i+1),&hodoi_t[i]);
  }
  anatrJ->Branch("hodo_id",&hodo_id);
  anatrJ->Branch("SBT_t",&SBT_t);
  anatrJ->Branch("tofF7SBT",&tofF7SBT);
  anatrJ->Branch("tofSBThodo",&tofSBThodo);
  anatrJ->Branch("tofSBTMinos",&tofSBTMinos);
  anatrJ->Branch("tofMinoshodo",&tofMinoshodo);
  for(Int_t i=0;i<24;++i){
    anatrJ->Branch(Form("hodo%02d_Z",i+1),&hodoi_Z[i]);
  }
  
  anatrJ->Branch("BR56Ca",&BR56Ca);

  //===== LOOP =====
  int nEntry = caltr->GetEntries();
  for(Int_t iEntry=0;iEntry<nEntry;++iEntry){
    //for(Int_t iEntry;iEntry<100,++iEntry){
    
    if(iEntry%100==0){
      clog << iEntry/1000 << "k events treated..." << "\r";
    }

    anatrB->GetEntry(iEntry);
    //caltrD->GetEntry(iEntry);
    caltr->GetEntry(iEntry);
    
    RNum = RunNumber;
    ENum = EventNumber;
    
    //=== Initialization ===
    /*
    dali_e.clear();
    dali_cos.clear();
    dali_t.clear();
    */
    hodo_q = Sqrt(-1);
    hodo_t = Sqrt(-1);
    for(Int_t i=0;i<24;++i){
      hodoi_q[i] = Sqrt(-1);
      hodoi_t[i] = Sqrt(-1);
      hodoi_Z[i] = Sqrt(-1);
    }
    hodo_id = 0;
    SBT_t = Sqrt(-1);
    tofF7SBT = Sqrt(-1);
    tofSBThodo = Sqrt(-1);
    tofSBTMinos = Sqrt(-1);
    tofMinoshodo = Sqrt(-1);
    
    BR56Ca = 0;
    
    //=== Calc. ===
    /*
    for(Int_t i=0;i<DALI_Energy->size();++i){
      dali_e.push_back(DALI_Energy->at(i));
      dali_cos.push_back(DALI_CosTheta->at(i));
      dali_t.push_back(DALI_Time->at(i));
    }
    */
    hodo_q = Hodo_QCal;
    hodo_t = Hodo_TCal;
    for(Int_t i=0;i<24;++i){
      hodoi_q[i] = Hodoi_QCal[i];
      hodoi_t[i] = Hodoi_TCal[i];
    }
    hodo_id = Hodo_ID;
    
    SBT_t = (SBT1_Time + SBT2_Time)/2.;
    tofF7SBT = (SBT_t - anaF7_Time + offsetF7SBT);
    tofSBThodo = (hodo_t - SBT_t);
    tofSBTMinos = distSBTMinos/(distF7SBT/(SBT_t - anaF7_Time + offsetF7SBT));
    tofMinoshodo = tofSBThodo - tofSBTMinos + offsethodo;

    Double_t zraw[24] = {0};
    if(hodo_id==12){
      zraw[11] = hodo_q - (58.3888*tofMinoshodo - 3570.69);
      hodoi_Z[11] = 18.9417 + 0.00533762*zraw[11] - 8.10244*pow(10,-7)*tofMinoshodo*tofMinoshodo;
    }
    //=== CUT ===
    if(cBR56Ca->IsInside(aoqBR,zetBR)) BR56Ca = 1;
    
    anatrJ->Fill();
  }//LOOP
  anafile->cd();
  anatrJ->Write();
  anafile->Close();
}//main
