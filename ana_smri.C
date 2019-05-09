#include<stdio.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<bitset>

#include<iomanip>
#include<sys/time.h>
#include<signal.h>

#include"TROOT.h"
#include"TDirectory.h"
#include"TFile.h"
#include"TH1I.h"
#include"TTree.h"
#include"TCut.h"
#include"TCutG.h"
#include"TMath.h"
#include"TString.h"
#include"TEnv.h"

#include"/home/koiwai/analysis/brho_func/Brho_A56Z20_br56Ca_sa56Ca.C"
#include"/home/koiwai/analysis/brho_func/Len_A56Z20_br56Ca_sa56Ca.C"
#include"/home/koiwai/analysis/brho_func/Brho_A54Z20_br54Ca_sa54Ca.C"
#include"/home/koiwai/analysis/brho_func/Len_A54Z20_br54Ca_sa54Ca.C"

using namespace std;
using namespace TMath;

bool signal_recieved = false;
void signalhandler(int sig){
  if(sig==SIGINT){
    signal_recieved = true;
  }
}

double get_time(){
  struct timeval t;
  gettimeofday(&t,NULL);
  double d = t.tv_sec + (double)t.tv_usec/1000000;
  return d;
}

int main(int argc, char *argv[]){
  double time_start = get_time();
  signal(SIGINT,signalhandler);

  
  time_t start, stop;
  time(&start);

  Int_t FileNum = TString(argv[1]).Atoi();

  printf("\n%s %d %s \n\n","=== Exexute ana_smri for RUN",FileNum,"===");
  
  //===== Load input file =====
  TString infname = Form("/home/koiwai/analysis/rootfiles/unpacked/run%04d.root",FileNum);
  TFile *infile = TFile::Open(infname);
  
  TTree *caltr;
  infile->GetObject("caltr",caltr);

  printf("%-20s %s \n","Input data file:",infname.Data());

  //===== input tree variables =====
  Long64_t EventNumber_all = 0;
  Int_t RunNumber_all = -1;

  //Double_t AllHodo_Charge[24];
  //Double_t AllHodo_Time[24];
  Int_t Hodo_ID; // ID of the hodoscope with the highest charge
  Int_t Hodo_Multiplicity;
  Double_t Hodo_QCal, Hodo_QRaw; // Highest charge //  Double_t Hodo_Charge;
  Double_t Hodo_TCal, Hodo_TRaw; // Time of the hodoscope with the highest charge  Double_t Hodo_Time;
  Int_t Hodoi_TURaw[24], Hodoi_TDRaw[24], Hodoi_QURaw[24], Hodoi_QDRaw[24];
  Double_t Hodoi_TUCal[24], Hodoi_TDCal[24], Hodoi_QUCal[24], Hodoi_QDCal[24];
  Double_t Hodoi_TCal[24], Hodoi_TRaw[24], Hodoi_QCal[24], Hodoi_QRaw[24];

  Double_t sbt1_Q, sbt2_Q, sbt1_T, sbt2_T, sbt1_dT, sbt2_dT;
  Double_t sbt1_QL, sbt1_QR, sbt1_TL, sbt1_TR;
  Double_t sbt2_QL, sbt2_QR, sbt2_TL, sbt2_TR;
  Double_t sbt1_Tslew, sbt2_Tslew;

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
  
  caltr->SetBranchAddress("sbt1_Q",&sbt1_Q);
  caltr->SetBranchAddress("sbt1_T",&sbt1_T);
  caltr->SetBranchAddress("sbt1_dT",&sbt1_dT);
  caltr->SetBranchAddress("sbt2_Q",&sbt2_Q);
  caltr->SetBranchAddress("sbt2_T",&sbt2_T);
  caltr->SetBranchAddress("sbt2_dT",&sbt2_dT);
  caltr->SetBranchAddress("sbt1_QL",&sbt1_QL);
  caltr->SetBranchAddress("sbt1_QR",&sbt1_QR);
  caltr->SetBranchAddress("sbt1_TL",&sbt1_TL);
  caltr->SetBranchAddress("sbt1_TR",&sbt1_TR);
  caltr->SetBranchAddress("sbt2_QL",&sbt2_QL);
  caltr->SetBranchAddress("sbt2_QR",&sbt2_QR);
  caltr->SetBranchAddress("sbt2_TL",&sbt2_TL);
  caltr->SetBranchAddress("sbt2_TR",&sbt2_TR);
  caltr->SetBranchAddress("sbt1_Tslew",&sbt1_Tslew);
  caltr->SetBranchAddress("sbt2_Tslew",&sbt2_Tslew);

  //===== Load input DC file =====
  TString infnameDC = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/ana_mwdc%04d.root",FileNum);
  TFile *infileDC = TFile::Open(infnameDC);
  
  TTree *anatrDC;
  infileDC->GetObject("anatrDC",anatrDC);

  printf("%-20s %s \n","Input MWDC file:",infnameDC.Data());

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
  
  //===== Load input Beam file =====
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNum);
  TFile *infileB = TFile::Open(infnameB);
  
  TTree *anatrB;
  infileB->GetObject("anatrB",anatrB);

  printf("%-20s %s \n","Input beam file:",infnameB.Data());

  //===== Declare Beam Variable =====
  Int_t EventNumber_beam, RunNumber_beam;

  Double_t betaF7F13, betaF3F13;

  Double_t zetBR, aoqBR;

  Bool_t BG_flag_beam;
  //Int_t BR56Sc;

  //=== PID ===
  Bool_t br64v; 
  Bool_t br63v; 
  Bool_t br62v; 
  Bool_t br61v; 
  Bool_t br60v; 
  Bool_t br62ti;
  Bool_t br61ti;
  Bool_t br60ti;
  Bool_t br59ti;
  Bool_t br58ti;
  Bool_t br57ti;
  Bool_t br59sc;
  Bool_t br58sc;
  Bool_t br57sc;
  Bool_t br56sc;
  Bool_t br55sc;
  Bool_t br54sc;
  Bool_t br56ca;
  Bool_t br55ca;
  Bool_t br54ca;
  Bool_t br53ca;
  Bool_t br52ca;
  Bool_t br54k; 
  Bool_t br53k; 
  Bool_t br52k; 
  Bool_t br51k; 
  Bool_t br50k; 
  Bool_t br49k; 
  Bool_t br51ar;
  Bool_t br50ar;
  Bool_t br49ar;
  Bool_t br48ar;
  
  //===== Beam SetBranchAddress =====
  anatrB->SetBranchAddress("EventNumber",&EventNumber_beam);
  anatrB->SetBranchAddress("RunNumber",&RunNumber_beam);

  anatrB->SetBranchAddress("betaF7F13",&betaF7F13); //
  anatrB->SetBranchAddress("betaF3F13",&betaF3F13); // 
  
  anatrB->SetBranchAddress("zetBR",&zetBR);
  anatrB->SetBranchAddress("aoqBR",&aoqBR);

  anatrB->SetBranchAddress("BG_flag",&BG_flag_beam);
  //anatrB->SetBranchAddress("BR56Sc",&BR56Sc);

  anatrB->SetBranchAddress("br64v",&br64v); 
  anatrB->SetBranchAddress("br63v",&br63v); 
  anatrB->SetBranchAddress("br62v",&br62v); 
  anatrB->SetBranchAddress("br61v",&br61v); 
  anatrB->SetBranchAddress("br60v",&br60v); 
  anatrB->SetBranchAddress("br62ti",&br62ti);
  anatrB->SetBranchAddress("br61ti",&br61ti);
  anatrB->SetBranchAddress("br60ti",&br60ti);
  anatrB->SetBranchAddress("br59ti",&br59ti);
  anatrB->SetBranchAddress("br58ti",&br58ti);
  anatrB->SetBranchAddress("br57ti",&br57ti);
  anatrB->SetBranchAddress("br59sc",&br59sc);
  anatrB->SetBranchAddress("br58sc",&br58sc);
  anatrB->SetBranchAddress("br57sc",&br57sc);
  anatrB->SetBranchAddress("br56sc",&br56sc);
  anatrB->SetBranchAddress("br55sc",&br55sc);
  anatrB->SetBranchAddress("br54sc",&br54sc);
  anatrB->SetBranchAddress("br56ca",&br56ca);
  anatrB->SetBranchAddress("br55ca",&br55ca);
  anatrB->SetBranchAddress("br54ca",&br54ca);
  anatrB->SetBranchAddress("br53ca",&br53ca);
  anatrB->SetBranchAddress("br52ca",&br52ca);
  anatrB->SetBranchAddress("br54k",&br54k); 
  anatrB->SetBranchAddress("br53k",&br53k); 
  anatrB->SetBranchAddress("br52k",&br52k); 
  anatrB->SetBranchAddress("br51k",&br51k); 
  anatrB->SetBranchAddress("br50k",&br50k); 
  anatrB->SetBranchAddress("br49k",&br49k); 
  anatrB->SetBranchAddress("br51ar",&br51ar);
  anatrB->SetBranchAddress("br50ar",&br50ar);
  anatrB->SetBranchAddress("br49ar",&br49ar);
  anatrB->SetBranchAddress("br48ar",&br48ar);
  
  //===== AddFriend =====
  caltr->AddFriend(anatrDC);
  caltr->AddFriend(anatrB);

  //===== Load cut files =====
  //TFile *cutfilesbt1 = new TFile("/home/koiwai/analysis/cutfiles/cutsbt1.root");
  //TCutG *csbt1 = (TCutG*)cutfilesbt1->Get("sbt1");

  //TFile *SApid_temp = new TFile("/home/koiwai/analysis/macros/SApid_temp.root");
  //TCutG *cSA56Sc_temp = (TCutG*)SApid_temp->Get("SA56Sc_temp");

  
  //===== Load .dat files =====
  TEnv *env              = new TEnv("/home/koiwai/analysis/db/geometry_psp.dat");
  TEnv *env_hodot        = new TEnv("/home/koiwai/analysis/db/hodo_toff.dat");
  TEnv *env_hodoq        = new TEnv("/home/koiwai/analysis/db/hodo_qcor.dat");
  TEnv *env_hodoq2z      = new TEnv("/home/koiwai/analysis/db/hodo_q2z.dat");
  TEnv *env_hodozraw2z   = new TEnv("/home/koiwai/analysis/db/hodo_zraw2z.dat");
  /*
  //TEnv *env_hodo02tofcor = new TEnv("/home/koiwai/analysis/db/hodo02_tofcor.dat");
  TEnv *env_hodo03tofcor = new TEnv("/home/koiwai/analysis/db/hodo03_tofcor.dat");
  TEnv *env_hodo04tofcor = new TEnv("/home/koiwai/analysis/db/hodo04_tofcor.dat");
  TEnv *env_hodo05tofcor = new TEnv("/home/koiwai/analysis/db/hodo05_tofcor.dat");
  TEnv *env_hodo06tofcor = new TEnv("/home/koiwai/analysis/db/hodo06_tofcor.dat");
  TEnv *env_hodo07tofcor = new TEnv("/home/koiwai/analysis/db/hodo07_tofcor.dat");
  TEnv *env_hodo08tofcor = new TEnv("/home/koiwai/analysis/db/hodo08_tofcor.dat");
  TEnv *env_hodo09tofcor = new TEnv("/home/koiwai/analysis/db/hodo09_tofcor.dat");
  TEnv *env_hodo10tofcor = new TEnv("/home/koiwai/analysis/db/hodo10_tofcor.dat");
  TEnv *env_hodo11tofcor = new TEnv("/home/koiwai/analysis/db/hodo11_tofcor.dat");
  TEnv *env_hodo12tofcor = new TEnv("/home/koiwai/analysis/db/hodo12_tofcor.dat");
  TEnv *env_hodo13tofcor = new TEnv("/home/koiwai/analysis/db/hodo13_tofcor.dat");
  TEnv *env_hodo14tofcor = new TEnv("/home/koiwai/analysis/db/hodo14_tofcor.dat");
  TEnv *env_hodo15tofcor = new TEnv("/home/koiwai/analysis/db/hodo15_tofcor.dat");
  TEnv *env_hodo16tofcor = new TEnv("/home/koiwai/analysis/db/hodo16_tofcor.dat");
  TEnv *env_hodo17tofcor = new TEnv("/home/koiwai/analysis/db/hodo17_tofcor.dat");
  TEnv *env_hodo18tofcor = new TEnv("/home/koiwai/analysis/db/hodo18_tofcor.dat");
  TEnv *env_hodo19tofcor = new TEnv("/home/koiwai/analysis/db/hodo19_tofcor.dat");
  TEnv *env_hodo20tofcor = new TEnv("/home/koiwai/analysis/db/hodo20_tofcor.dat");
  TEnv *env_hodo21tofcor = new TEnv("/home/koiwai/analysis/db/hodo21_tofcor.dat");
  TEnv *env_hodo22tofcor = new TEnv("/home/koiwai/analysis/db/hodo22_tofcor.dat");
  TEnv *env_hodo23tofcor = new TEnv("/home/koiwai/analysis/db/hodo23_tofcor.dat");
  //TEnv *env_hodo24tofcor = new TEnv("/home/koiwai/analysis/db/hodo24_tofcor.dat");
  */
  TEnv *env_hodotofcor[24];
  for(int id=0;id<24;++id){
    env_hodotofcor[id] = new TEnv(Form("/home/koiwai/analysis/db/hodo%02d_tofcor.dat",id+1));
  }
  
  TEnv *env_hodoaoqcor   = new TEnv("/home/koiwai/analysis/db/hodo_aoqcor.dat");
  
  //===== Create output file/tree =====
  TString ofname = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNum);
  TFile *anafile_smri = new TFile(ofname,"RECREATE");
  TTree *anatrS  = new TTree("anatrS","anatrS");

  printf("\n%-20s %s \n\n","Output file:",ofname.Data());
  
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
  //Double_t toff_hodo = env->GetValue("toff_hodo",250.);
  Double_t toff_hodo = 227.;
  Double_t clight = 299.79258; //[mm/nsec]
  Double_t mu = 931.49432; //[MeV]
  Double_t me = 511.; //[keV]
  Double_t ionpair = 4.866; //[keV

  
  //===== Declare anatree const.s =====
  Double_t hodo_toff[24], hodo_qcor[24];
  //Double_t hodo_t2q0[24], hodo_t2q1[24];
  //Double_t hodo_zraw2z0[24], hodo_zraw2z1[24], hodo_zraw2z2[24];
  Double_t hodo_235T_zraw2z_p0[24], hodo_235T_zraw2z_p1[24];
  Double_t hodo_270T_zraw2z_p0[24], hodo_270T_zraw2z_p1[24];
  Double_t hodo_zraw2z[2];
  for(Int_t id=0;id<24;id++){
    hodo_toff[id] = env_hodot->GetValue(Form("hodo_toff_%02d",id+1),0.0);
    hodo_qcor[id] = env_hodoq->GetValue(Form("hodo_qcor_%02d",id+1),0.0);
    //hodo_t2q0[id] = env_hodoq2z->GetValue(Form("hodo%02dt2q0",id+1),0.0);
    //hodo_t2q1[id] = env_hodoq2z->GetValue(Form("hodo%02dt2q1",id+1),0.0);
    //hodo_zraw2z0[id] = env_hodoq2z->GetValue(Form("hodo%02dzraw2z0",id+1),0.0);
    //hodo_zraw2z1[id] = env_hodoq2z->GetValue(Form("hodo%02dzraw2z1",id+1),0.0);
    //hodo_zraw2z2[id] = env_hodoq2z->GetValue(Form("hodo%02dzraw2z2",id+1),0.0);
    hodo_235T_zraw2z_p0[id] = env_hodozraw2z->GetValue(Form("235T_%02d_p0",id+1),0.0);
    hodo_235T_zraw2z_p1[id] = env_hodozraw2z->GetValue(Form("235T_%02d_p1",id+1),0.0);
    hodo_270T_zraw2z_p0[id] = env_hodozraw2z->GetValue(Form("270T_%02d_p0",id+1),0.0);
    hodo_270T_zraw2z_p1[id] = env_hodozraw2z->GetValue(Form("270T_%02d_p1",id+1),0.0);
  }
  hodo_zraw2z[0] = env_hodozraw2z->GetValue("zraw2z_p0",0.0);
  hodo_zraw2z[1] = env_hodozraw2z->GetValue("zraw2z_p1",0.0);

  Double_t hodo_aoqcor[24][2];
  for(Int_t i=0;i<24;i++){
    hodo_aoqcor[i][0] = env_hodoaoqcor->GetValue(Form("%02dp0",i+1),0.0);
    hodo_aoqcor[i][1] = env_hodoaoqcor->GetValue(Form("%02dp1",i+1),1.0);
    //cout << hodo_aoqcor[i][0] << endl;
  }
  /*
  Double_t hodo02_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo02_tofcor[i] = env_hodo02tofcor->GetValue(n,0.0);
  }
  Double_t hodo03_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo03_tofcor[i] = env_hodo03tofcor->GetValue(n,0.0);
  }
  Double_t hodo04_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo04_tofcor[i] = env_hodo04tofcor->GetValue(n,0.0);
  }
  Double_t hodo05_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo05_tofcor[i] = env_hodo05tofcor->GetValue(n,0.0);
  }
  Double_t hodo06_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo06_tofcor[i] = env_hodo06tofcor->GetValue(n,0.0);
  }
  Double_t hodo07_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo07_tofcor[i] = env_hodo07tofcor->GetValue(n,0.0);
  }
  Double_t hodo08_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo08_tofcor[i] = env_hodo08tofcor->GetValue(n,0.0);
  }
  Double_t hodo09_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo09_tofcor[i] = env_hodo09tofcor->GetValue(n,0.0);
  }
  Double_t hodo10_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo10_tofcor[i] = env_hodo10tofcor->GetValue(n,0.0);
  }
  Double_t hodo11_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo11_tofcor[i] = env_hodo11tofcor->GetValue(n,0.0);
  }
  Double_t hodo12_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo12_tofcor[i] = env_hodo12tofcor->GetValue(n,0.0);
  }
  Double_t hodo13_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo13_tofcor[i] = env_hodo13tofcor->GetValue(n,0.0);
  }
  Double_t hodo14_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo14_tofcor[i] = env_hodo14tofcor->GetValue(n,0.0);
  }
  Double_t hodo15_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo15_tofcor[i] = env_hodo15tofcor->GetValue(n,0.0);
  }
  Double_t hodo16_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo16_tofcor[i] = env_hodo16tofcor->GetValue(n,0.0);
  }
  Double_t hodo17_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo17_tofcor[i] = env_hodo17tofcor->GetValue(n,0.0);
  }
  Double_t hodo18_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo18_tofcor[i] = env_hodo18tofcor->GetValue(n,0.0);
  }
  Double_t hodo19_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo19_tofcor[i] = env_hodo19tofcor->GetValue(n,0.0);
  }
  Double_t hodo20_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo20_tofcor[i] = env_hodo20tofcor->GetValue(n,0.0);
  }
  Double_t hodo21_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo21_tofcor[i] = env_hodo21tofcor->GetValue(n,0.0);
  }
  Double_t hodo22_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo22_tofcor[i] = env_hodo22tofcor->GetValue(n,0.0);
  }
  Double_t hodo23_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo23_tofcor[i] = env_hodo23tofcor->GetValue(n,0.0);
  }
  Double_t hodo24_tofcor[231];
  for(Int_t i=0;i<231;i++){
    const char *n = Form("%d",FileNum);
    hodo24_tofcor[i] = env_hodo24tofcor->GetValue(n,0.0);
  }
  */

  Double_t hodo_tofcor[24];
  const char *nn = Form("%d",FileNum);
  for(int id=0;id<24;++id){
    //TEnv *env_hodoIDtofcor = (TEnv*)Form("env_hodo%02dtofcor",id+1);
    hodo_tofcor[id] = env_hodotofcor[id]->GetValue(nn,0.0);
    //hodo_tofcor[id] = env_hodoIDtofcor->GetValue(nn,0.0);
    //cout << "hodo_tofcor" << id << " " << hodo_tofcor[id] << endl;
  }
  /* //try some day...
  Double_t hodo_tofcor[24][231];
  //for(Int_t id=0;id<24;id++){
  TEnv *env_hodotofcor = (TEnv*)Form("env_hodo%02dtofcor",12);
    //if(!Form("env_hodo%02dtofcor",id+1)) continue;
    //if(!env_hodotofcor) continue;
    for(Int_t runnum=0;runnum<231;runnum++){
      const char *n = Form("%d",FileNum);
      //hodo_tofcor[11][runnum] = Form("env_hodo%02dtofcor",12)->GetValue(n,0.0);
      hodo_tofcor[11][runnum] = env_hodotofcor->GetValue(n,0.0);
    }
    //}
    */
  //===== Declare valables for calc. =====
  Double_t dev;

  
  //===== Declare anatree variables =====
  Int_t RunNum, EventNum;
  

  Int_t hodo_id, hodo_multi;
  Double_t hodo_q, hodo_t;

  Double_t t_minoshodo_notcor, t_minoshodo, v_minoshodo, beta_minoshodo, gamma_minoshodo;


  
  Double_t brhoSA, lengSA;
  //Double_t brhoSA_tan, lengSA_tan;
  Double_t brhoSA_rad, lengSA_rad;

  
  Double_t zraw;

  Double_t zetSA235, zetSA270;
  
  Double_t aoqSA, zetSA;

  Double_t aoqSA_notcor, aoqSA_tmpcor;
  
  Bool_t BG_flag;

  //Int_t SA56Sc_temp;

  //double sbt_t;

  double hodo09_tofcor;
  
  //===== Create anatree Branch =====
  anatrS->Branch("RunNumber",&RunNum);
  anatrS->Branch("EventNumber",&EventNum);

  anatrS->Branch("t_minoshodo_notcor",&t_minoshodo_notcor);
  anatrS->Branch("t_minoshodo",&t_minoshodo);
  anatrS->Branch("v_minoshodo",&v_minoshodo);
  anatrS->Branch("beta_minoshodo",&beta_minoshodo);
  anatrS->Branch("gamma_minoshodo",&gamma_minoshodo);

  anatrS->Branch("hodo_id",&hodo_id);
  anatrS->Branch("hodo_multi",&hodo_multi);
  anatrS->Branch("hodo_q",&hodo_q);
  anatrS->Branch("hodo_t",&hodo_t);

  anatrS->Branch("aoqSA",&aoqSA);
  anatrS->Branch("aoqSA_notcor",&aoqSA_notcor);
  anatrS->Branch("aoqSA_tmpcor",&aoqSA_tmpcor);
  anatrS->Branch("zetSA",&zetSA);

  anatrS->Branch("zetSA235",&zetSA235);
  anatrS->Branch("zetSA270",&zetSA270);

  anatrS->Branch("zraw",&zraw);
  
  anatrS->Branch("brhoSA",&brhoSA);
  anatrS->Branch("lengSA",&lengSA);
  //anatrS->Branch("brhoSA_tan",&brhoSA_tan);
  //anatrS->Branch("lengSA_tan",&lengSA_tan);
  anatrS->Branch("brhoSA_rad",&brhoSA_rad);
  anatrS->Branch("lengSA_rad",&lengSA_rad);

  anatrS->Branch("BG_flag",&BG_flag,"BG_flag/O");
  anatrS->Branch("BG_flag_beam",&BG_flag_beam);

  anatrS->Branch("br64v",&br64v); 
  anatrS->Branch("br63v",&br63v); 
  anatrS->Branch("br62v",&br62v); 
  anatrS->Branch("br61v",&br61v); 
  anatrS->Branch("br60v",&br60v); 
  anatrS->Branch("br62ti",&br62ti);
  anatrS->Branch("br61ti",&br61ti);
  anatrS->Branch("br60ti",&br60ti);
  anatrS->Branch("br59ti",&br59ti);
  anatrS->Branch("br58ti",&br58ti);
  anatrS->Branch("br57ti",&br57ti);
  anatrS->Branch("br59sc",&br59sc);
  anatrS->Branch("br58sc",&br58sc);
  anatrS->Branch("br57sc",&br57sc);
  anatrS->Branch("br56sc",&br56sc);
  anatrS->Branch("br55sc",&br55sc);
  anatrS->Branch("br54sc",&br54sc);
  anatrS->Branch("br56ca",&br56ca);
  anatrS->Branch("br55ca",&br55ca);
  anatrS->Branch("br54ca",&br54ca);
  anatrS->Branch("br53ca",&br53ca);
  anatrS->Branch("br52ca",&br52ca);
  anatrS->Branch("br54k",&br54k); 
  anatrS->Branch("br53k",&br53k); 
  anatrS->Branch("br52k",&br52k); 
  anatrS->Branch("br51k",&br51k); 
  anatrS->Branch("br50k",&br50k); 
  anatrS->Branch("br49k",&br49k); 
  anatrS->Branch("br51ar",&br51ar);
  anatrS->Branch("br50ar",&br50ar);
  anatrS->Branch("br49ar",&br49ar);
  anatrS->Branch("br48ar",&br48ar);

  anatrS->Branch("hodo09_tofcor",&hodo09_tofcor);
  
  //anatrS->Branch("BR56Sc",&BR56Sc);
  //anatrS->Branch("SA56Sc_temp",&SA56Sc_temp);

  //anatrS->Branch("sbt_t",&sbt_t);
  
  //===== Begin LOOP =====
  double time_prev = 0.;
  int nEntry = caltr->GetEntries();
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    //for(int iEntry=0;iEntry<1;++iEntry){
    
    if(iEntry%1000 == 0){
      //clog << iEntry/1000 << "k events treated..." << endl;
      //time(&t1);
      //cout <<
      double time_end = get_time();
      cout << (100.*iEntry)/nEntry << " % (" << iEntry << " events) done\t" << iEntry/(time_end - time_start) << " events/s\t" << (nEntry - iEntry)*(time_end - time_start)/(double)iEntry << " s to go\r " ;
      if(iEntry!=1000) cout << "current speed: " << 1000./(time_end - time_prev) << " events/s\t" << endl;
      else cout << endl;
      time_prev = get_time();
    }

    caltr->GetEntry(iEntry);

    RunNum = RunNumber_all;
    EventNum = EventNumber_all;
    
 
    //@@@ Brho Length Function @@@
    //=== Initialize ===
    brhoSA = Sqrt(-1);
    lengSA = Sqrt(-1);
    //brhoSA_tan = Sqrt(-1);
    //lengSA_tan = Sqrt(-1);
    brhoSA_rad = Sqrt(-1);
    lengSA_rad = Sqrt(-1);

    BG_flag = 0;
    //SA56Sc_temp = 0;
    
    Double_t x[6];

    x[0] = FDC1_X;
    x[1] = FDC1_A;
    x[2] = FDC1_Y;
    x[3] = FDC1_B;
    x[4] = FDC2_X;
    x[5] = FDC2_A;

    /*
    Double_t tan[6];

    tan[0] = FDC1_X;
    tan[1] = Tan(FDC1_A)*1000;
    tan[2] = FDC1_Y;
    tan[3] = Tan(FDC1_B)*1000;
    tan[4] = FDC2_X;
    tan[5] = Tan(FDC2_A)*1000;
    */
    
    Double_t rad[6];

    rad[0] = FDC1_X;
    rad[1] = FDC1_A*1000;
    rad[2] = FDC1_Y;
    rad[3] = FDC1_B*1000;
    rad[4] = FDC2_X;
    rad[5] = FDC2_A*1000;

    brhoSA = MDF_Brho_A56Z20(x);
    lengSA = MDF_Len_A56Z20(x);
    //brhoSA = MDF_Brho_A54Z20(x);
    //lengSA = MDF_Len_A54Z20(x);

    //brhoSA_tan = MDF_Brho_A56Z20(tan);
    //lengSA_tan = MDF_Len_A56Z20(tan);

    brhoSA_rad = MDF_Brho_A56Z20(rad);
    lengSA_rad = MDF_Len_A56Z20(rad);
    //brhoSA_rad = MDF_Brho_A54Z20(rad);
    //lengSA_rad = MDF_Len_A54Z20(rad);
    
    
    //@@@ HODO @@@
    //=== Initialize ===
    t_minoshodo_notcor = Sqrt(-1);
    t_minoshodo     = Sqrt(-1);
    v_minoshodo     = Sqrt(-1);
    beta_minoshodo  = Sqrt(-1);
    gamma_minoshodo = Sqrt(-1);

    hodo_q     = 0.;
    hodo_t     = Sqrt(-1);
    hodo_id    = 0;
    hodo_multi = 0;

    aoqSA        = Sqrt(-1);
    aoqSA_notcor = Sqrt(-1);
    aoqSA_tmpcor = Sqrt(-1);

    zetSA    = Sqrt(-1);
    zetSA235 = Sqrt(-1);
    zetSA270 = Sqrt(-1);
    zraw     = Sqrt(-1);
    dev      = Sqrt(-1);

    Double_t allHodo_Q[24];
    Double_t allHodo_T[24];
    for(Int_t i=0;i<24;i++){
      allHodo_Q[i] = Sqrt(-1);
      allHodo_T[i] = Sqrt(-1);
    }

    hodo09_tofcor = Sqrt(-1);

    //=== Calc. ===
    for(Int_t i=0;i<24;i++){
      allHodo_Q[i] = Hodoi_QCal[i]*hodo_qcor[i];
      allHodo_T[i] = Hodoi_TCal[i]+hodo_toff[i];
      if(allHodo_Q[i]>hodo_q){
	hodo_q  = allHodo_Q[i];
	hodo_t  = allHodo_T[i];
	hodo_id = i+1;
      }
    }
    hodo_multi = Hodo_Multiplicity;

    //@@@ HODO end @@@

    //@@@ SBT1 slew correction @@@

    //SBT1_Time += -35./sqrt(SBT1_QL)+1.095 -25./sqrt(SBT1_QR)+1.179;
    //SBT2_Time += -28./sqrt(SBT2_QL)+1.230 -40./sqrt(SBT2_QR)+2.118;

    //sbt_t = SBT1_Time;
    
    //t_minoshodo_notcor = hodo_t - sbt1_Tslew - (Dist_SBTTarget/betaF7F13/clight) + toff_hodo;
    t_minoshodo_notcor = hodo_t - sbt1_Tslew - (Dist_SBTTarget/betaF3F13/clight) + toff_hodo;

    /*
    switch(hodo_id){
    case  2: t_minoshodo = t_minoshodo_notcor + hodo02_tofcor[RunNum]; break;
    case  3: t_minoshodo = t_minoshodo_notcor + hodo03_tofcor[RunNum]; break;
    case  4: t_minoshodo = t_minoshodo_notcor + hodo04_tofcor[RunNum]; break;
    case  5: t_minoshodo = t_minoshodo_notcor + hodo05_tofcor[RunNum]; break;
    case  6: t_minoshodo = t_minoshodo_notcor + hodo06_tofcor[RunNum]; break;
    case  7: t_minoshodo = t_minoshodo_notcor + hodo07_tofcor[RunNum]; break;
    case  8: t_minoshodo = t_minoshodo_notcor + hodo08_tofcor[RunNum]; break;
    case  9: t_minoshodo = t_minoshodo_notcor + hodo09_tofcor[RunNum]; break;
    case 10: t_minoshodo = t_minoshodo_notcor + hodo10_tofcor[RunNum]; break;
    case 11: t_minoshodo = t_minoshodo_notcor + hodo11_tofcor[RunNum]; break;
    case 12: t_minoshodo = t_minoshodo_notcor + hodo12_tofcor[RunNum]; break;
    case 13: t_minoshodo = t_minoshodo_notcor + hodo13_tofcor[RunNum]; break;
    case 14: t_minoshodo = t_minoshodo_notcor + hodo14_tofcor[RunNum]; break;
    case 15: t_minoshodo = t_minoshodo_notcor + hodo15_tofcor[RunNum]; break;
    case 16: t_minoshodo = t_minoshodo_notcor + hodo16_tofcor[RunNum]; break;
    case 17: t_minoshodo = t_minoshodo_notcor + hodo17_tofcor[RunNum]; break;
    case 18: t_minoshodo = t_minoshodo_notcor + hodo18_tofcor[RunNum]; break;
    case 19: t_minoshodo = t_minoshodo_notcor + hodo19_tofcor[RunNum]; break;
    case 20: t_minoshodo = t_minoshodo_notcor + hodo20_tofcor[RunNum]; break;
    case 21: t_minoshodo = t_minoshodo_notcor + hodo21_tofcor[RunNum]; break;
    case 22: t_minoshodo = t_minoshodo_notcor + hodo22_tofcor[RunNum]; break;
    case 23: t_minoshodo = t_minoshodo_notcor + hodo23_tofcor[RunNum]; break;
    case 24: t_minoshodo = t_minoshodo_notcor + hodo24_tofcor[RunNum]; break;					 
    }
    */
    //t_minoshodo = t_minoshodo_notcor + hodo_tofcor[hodo_id-1];
    t_minoshodo = hodo_t - sbt1_Tslew - (Dist_SBTTarget/betaF3F13/clight) + toff_hodo + hodo_tofcor[hodo_id-1];

    hodo09_tofcor = hodo_tofcor[9-1];
    
    v_minoshodo = lengSA_rad/t_minoshodo;
    beta_minoshodo  = v_minoshodo/clight;
    gamma_minoshodo = 1./Sqrt(1.-beta_minoshodo*beta_minoshodo);

    dev = Log(2*me*beta_minoshodo*beta_minoshodo/ionpair) - Log(1 - beta_minoshodo*beta_minoshodo) - beta_minoshodo*beta_minoshodo; 
    
    zraw = v_minoshodo*Sqrt(hodo_q/dev);
    
    //zetSA235 = hodo_235T_zraw2z_p0[hodo_id-1] + hodo_235T_zraw2z_p1[hodo_id-1]*zraw;
    //zetSA270 = hodo_270T_zraw2z_p0[hodo_id-1] + hodo_270T_zraw2z_p1[hodo_id-1]*zraw;
    if(hodo_id==3||hodo_id==4||hodo_id==5)
      zetSA = hodo_235T_zraw2z_p0[hodo_id-1] + hodo_235T_zraw2z_p1[hodo_id-1]*zraw;
    else
      zetSA = hodo_270T_zraw2z_p0[hodo_id-1] + hodo_270T_zraw2z_p1[hodo_id-1]*zraw;

    zetSA = hodo_zraw2z[0] + hodo_zraw2z[1]*zetSA;
    
    aoqSA_notcor = brhoSA_rad/beta_minoshodo/gamma_minoshodo*clight/mu;
    aoqSA_tmpcor = 0.545 + 0.7551*aoqSA_notcor;
    //aoqSA = hodo_aoqcor[hodo_id-1][0] + hodo_aoqcor[hodo_id-1][1]*aoqSA_tmpcor;
    aoqSA = hodo_aoqcor[hodo_id-1][0] + hodo_aoqcor[hodo_id-1][1]*aoqSA_notcor;
        
    //zraw = hodo_q - (hodo_t2q0[hodo_id+1] + hodo_t2q1[hodo_id+1]*t_minoshodo);
    //zetSA = hodo_zraw2z0[hodo_id+1] + hodo_zraw2z1[hodo_id+1]*zraw + hodo_zraw2z2[hodo_id+1]*zraw*zraw;

    //cout << "z " << zetSA <<endl;
    



    //===== BG cut =====
    //=== cut by CUTG ===
    //if(!csbt1->IsInside(SBT1_TR-SBT1_TL,log(SBT1_QL/SBT1_QR))) BG_flag = 1;
    //if(cSA56Sc_temp->IsInside(aoqSA,zetSA)) SA56Sc_temp = 1;
    
    //=== cut by Hodo Time ===
    for(Int_t i=0;i<24;i++){
      if(Hodoi_TURaw[i]==-1000||Hodoi_TDRaw[i]==-1000) BG_flag = 1;
    }
    anatrS->Fill();
  }//for LOOP
  anafile_smri->cd();
  anatrS->Write();
  anafile_smri->Close();

  time(&stop);
  printf("Elapsed time: %.1f seconds\n",difftime(stop,start));

  double time_end = get_time();
  cout << endl << "Program Run Time " << time_end - time_start << " s." << endl;
  cout << nEntry/(time_end - time_start) << " events/s." << endl;



  
  //printf("%d k events have been treated.\n",nEntry/1000);
  cout << "RUN " << FileNum << ":" << "Conversion finished! : smri" << FileNum << "ok" << endl;
}//main()
