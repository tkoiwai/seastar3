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

//#include"/home/koiwai/analysis/brho_func/Brho_A56Z20_br56Ca_sa56Ca.C"
#include"/home/koiwai/analysis/macros/brho_func/Brho_A56Z20_br57Sc_sa56Ca.C"
//#include"/home/koiwai/analysis/brho_func/Len_A56Z20_br56Ca_sa56Ca.C"
#include"/home/koiwai/analysis/macros/brho_func/Len_A56Z20_br57Sc_sa56Ca.C"
#include"/home/koiwai/analysis/brho_func/Brho_A54Z20_br54Ca_sa54Ca.C"
#include"/home/koiwai/analysis/brho_func/Len_A54Z20_br54Ca_sa54Ca.C"

#include"/home/koiwai/analysis/include/smridef.h"

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

  printf("\n%s %d %s \n\n","=== Execute ana_smri for RUN",FileNum,"===");
  
  //===== Load input file =======================================================================
  TString infname = Form("/home/koiwai/analysis/rootfiles/unpacked/run%04d.root",FileNum);
  TFile *infile = TFile::Open(infname);
  
  TTree *caltr;
  infile->GetObject("caltr",caltr);

  printf("%-20s %s \n","Input data file:",infname.Data());

  Get_Branch_unpacked(caltr);
  
  //===== Load input DC file ====================================================================
  TString infnameDC = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/ana_mwdc%04d.root",FileNum);
  TFile *infileDC = TFile::Open(infnameDC);
  
  TTree *anatrDC;
  infileDC->GetObject("anatrDC",anatrDC);

  printf("%-20s %s \n","Input MWDC file:",infnameDC.Data());

  Get_Branch_mwdc(anatrDC);
  
  //===== Load input Beam file ==================================================================
  TString infnameB = Form("/home/koiwai/analysis/rootfiles/ana/beam/ana_beam%04d.root",FileNum);
  TFile *infileB = TFile::Open(infnameB);
  
  TTree *anatrB;
  infileB->GetObject("anatrB",anatrB);

  printf("%-20s %s \n","Input beam file:",infnameB.Data());

  Get_Branch_beam(anatrB);
  
  //===== AddFriend =============================================================================
  caltr->AddFriend(anatrDC);
  caltr->AddFriend(anatrB);

  //===== Load cut files ========================================================================
  
  //===== Load .dat files =======================================================================
  TEnv *env              = new TEnv("/home/koiwai/analysis/db/geometry_psp.dat");
  TEnv *env_hodot        = new TEnv("/home/koiwai/analysis/db/hodo_toff.dat"); //unpacked -> hodo_t
  TEnv *env_hodoq        = new TEnv("/home/koiwai/analysis/db/hodo_qcor.dat"); //unpacked -> hodo_q
  //TEnv *env_hodoq2z      = new TEnv("/home/koiwai/analysis/db/hodo_q2z.dat");  
  TEnv *env_hodozraw2z   = new TEnv("/home/koiwai/analysis/db/hodo_zraw2z.dat"); // zraw -> z
  TEnv *env_hodotofcor[24]; // t_minoshodo_notcor + hodo_tofcor = t_minodhodo
  for(int id=0;id<24;++id)
    env_hodotofcor[id] = new TEnv(Form("/home/koiwai/analysis/db/hodo%02d_tofcor.dat",id+1)); //tofcor over runs
  TEnv *env_hodoaoqcor   = new TEnv("/home/koiwai/analysis/db/hodo_aoqcor.dat"); //aoqcor for each bar
  
  //===== Create output file/tree ===============================================================
  TString ofname = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",FileNum);
  TFile *anafile_smri = new TFile(ofname,"RECREATE");
  TTree *anatrS  = new TTree("anatrS","anatrS");

  printf("\n%-20s %s \n\n","Output file:",ofname.Data());
  
  //===== Declare const.s =======================================================================
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

  
  //===== Declare anatree const.s ===============================================================
  Double_t hodo_toff[24], hodo_qcor[24];
  //Double_t hodo_t2q0[24], hodo_t2q1[24];
  Double_t hodo_zraw2z_p0[24], hodo_zraw2z_p1[24];// hodo_zraw2z2[24];
  //Double_t hodo_235T_zraw2z_p0[24], hodo_235T_zraw2z_p1[24];
  //Double_t hodo_270T_zraw2z_p0[24], hodo_270T_zraw2z_p1[24];
  Double_t hodo_zraw2z[2];
  for(Int_t id=0;id<24;id++){
    hodo_toff[id] = env_hodot->GetValue(Form("hodo_toff_%02d",id+1),0.0);
    hodo_qcor[id] = env_hodoq->GetValue(Form("hodo_qcor_%02d",id+1),0.0);
    //hodo_235T_zraw2z_p0[id] = env_hodozraw2z->GetValue(Form("235T_%02d_p0",id+1),0.0);
    //hodo_235T_zraw2z_p1[id] = env_hodozraw2z->GetValue(Form("235T_%02d_p1",id+1),0.0);
    //hodo_270T_zraw2z_p0[id] = env_hodozraw2z->GetValue(Form("270T_%02d_p0",id+1),0.0);
    //hodo_270T_zraw2z_p1[id] = env_hodozraw2z->GetValue(Form("270T_%02d_p1",id+1),0.0);
    hodo_zraw2z_p0[id] = env_hodozraw2z->GetValue(Form("zraw2z_%02d_p0",id+1),0.0);
    hodo_zraw2z_p1[id] = env_hodozraw2z->GetValue(Form("zraw2z_%02d_p1",id+1),1.0);
  }
  //hodo_zraw2z[0] = env_hodozraw2z->GetValue("zraw2z_p0",0.0);
  //hodo_zraw2z[1] = env_hodozraw2z->GetValue("zraw2z_p1",0.0);

  Double_t hodo_aoqcor[24][2];
  for(Int_t i=0;i<24;i++){
    hodo_aoqcor[i][0] = env_hodoaoqcor->GetValue(Form("%02dp0",i+1),0.0);
    hodo_aoqcor[i][1] = env_hodoaoqcor->GetValue(Form("%02dp1",i+1),1.0);
  }

  Double_t hodo_tofcor[24];
  const char *nn = Form("%d",FileNum);
  for(int id=0;id<24;++id){
    //TEnv *env_hodoIDtofcor = (TEnv*)Form("env_hodo%02dtofcor",id+1);
    hodo_tofcor[id] = env_hodotofcor[id]->GetValue(nn,0.0);
    //hodo_tofcor[id] = env_hodoIDtofcor->GetValue(nn,0.0);
    //cout << "hodo_tofcor" << id << " " << hodo_tofcor[id] << endl;
  }

  //===== Declare valables for calc. =============================================================
  Double_t dev, dev54;

  Set_Branch_smri(anatrS);
 
  //===== Begin LOOP =============================================================================
  double time_prev = 0.;
  int AllEntry = caltr->GetEntries();
  int tmpEntry = 8000000;
  int nEntry = 0;
  if(AllEntry>tmpEntry)
    nEntry = tmpEntry;
  else
    nEntry = AllEntry;

  cout << "Number of events to treat: " << nEntry << endl;

  double time_startloop = get_time();
  
  for(int iEntry=0;iEntry<nEntry;++iEntry){
  //for(int iEntry=0;iEntry<100000;++iEntry){
    
    if(iEntry%1000 == 0){
      //clog << iEntry/1000 << "k events treated..." << endl;
      //time(&t1);
      //cout <<
      double time_end = get_time();

      double t_diff_a = time_end - time_start;
      int t_hour_a    = t_diff_a/3600;
      int t_min_a     = (t_diff_a - 3600*t_hour_a)/60;
      double t_sec_a  = t_diff_a -3600*t_hour_a - 60*t_min_a;
      
      cout << "\r"
	   << t_hour_a <<"h"<< t_min_a <<"m"<< t_sec_a <<"s elapsed:  "
	   << (100.*iEntry)/nEntry << " % (" << iEntry << " events) done:  "
	   << iEntry/(time_end - time_startloop) << " events/s:  "
	   << (nEntry - iEntry)*(time_end - time_startloop)/(double)iEntry << " s to go:  " ;
      if(iEntry!=1000) cout << "current speed: " << 1000./(time_end - time_prev) << " events/s     " << flush;
      cout << endl << flush;
      //else cout << endl;
      time_prev = get_time();
    }

    caltr->GetEntry(iEntry);

    RunNum = RunNumber_all;
    EventNum = EventNumber_all;
    
 
    //@@@ Brho Length Function @@@
    //=== Initialize ===-------------------------------------------------------------------------
    brhoSA     = Sqrt(-1);
    lengSA     = Sqrt(-1);
    brhoSA_tan = Sqrt(-1);
    lengSA_tan = Sqrt(-1);
    brhoSA_rad = Sqrt(-1);
    lengSA_rad = Sqrt(-1);
    brho54     = Sqrt(-1);
    leng54     = Sqrt(-1);

    BG_flag = 0;

    goodEvt      = false;
    goodEvt_beam = true;
    goodEvt_smri = true;
    goodEvt_mwdc = true;

    //20191111
    BDC2_X -= 0.5096;
    BDC2_Y -= 0.4481;

    Target_X = BDC1_X +Dist_BDC1Target / Dist_BDC1BDC2 * (BDC2_X - BDC1_X);
    Target_Y = BDC1_Y +Dist_BDC1Target / Dist_BDC1BDC2 * (BDC2_Y - BDC1_Y);
    //Target_A = (BDC2_X - BDC1_X) / Dist_BDC1BDC2;
    //Target_B = (BDC2_Y - BDC1_Y) / Dist_BDC1BDC2;

    if(TMath::Abs(Target_X)<100 && TMath::Abs(Target_Y)<100 &&
       TMath::Abs(FDC1_X)<5000 && TMath::Abs(FDC1_Y) < 5000){
      FDC1_A = (FDC1_X - Target_X) / (Dist_BDC1FDC1-Dist_BDC1Target);
      FDC1_B = (FDC1_Y - Target_Y) / (Dist_BDC1FDC1-Dist_BDC1Target);
    }


    
    Double_t x[6];

    x[0] = FDC1_X;
    x[1] = FDC1_A;
    x[2] = FDC1_Y;
    x[3] = FDC1_B;
    x[4] = FDC2_X;
    x[5] = FDC2_A;

    
    Double_t tan[6];

    tan[0] = FDC1_X;
    tan[1] = Tan(FDC1_A)*1000;
    tan[2] = FDC1_Y;
    tan[3] = Tan(FDC1_B)*1000;
    tan[4] = FDC2_X;
    tan[5] = Tan(FDC2_A)*1000;
    
    
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

    brhoSA_tan = MDF_Brho_A56Z20(tan);
    lengSA_tan = MDF_Len_A56Z20(tan);

    brhoSA_rad = MDF_Brho_A56Z20(rad);
    lengSA_rad = MDF_Len_A56Z20(rad);
    //brhoSA_rad = MDF_Brho_A54Z20(rad);
    //lengSA_rad = MDF_Len_A54Z20(rad);


    brho54 = MDF_Brho_A54Z20(tan);
    leng54 = MDF_Len_A54Z20(tan);
    
    
    //@@@ HODO @@@
    //=== Initialize ===---------------------------------------------------------------------------
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

    zetSA    = Sqrt(-1);
    zraw     = Sqrt(-1);
    dev      = Sqrt(-1);
    dev54    = Sqrt(-1);

    Double_t allHodo_Q[24];
    Double_t allHodo_T[24];
    for(Int_t i=0;i<24;i++){
      allHodo_Q[i] = Sqrt(-1);
      allHodo_T[i] = Sqrt(-1);
    }

    //hodo09_tofcor = Sqrt(-1);

    tof13T   = Sqrt(-1);
    tof13H   = Sqrt(-1);
    tofTH_nc = Sqrt(-1);
    tofTH    = Sqrt(-1);

    betaTH = Sqrt(-1);
    gammaTH = Sqrt(-1);

    betaTH54  = Sqrt(-1);
    gammaTH54 = Sqrt(-1);
    zraw54    = Sqrt(-1);
    aoqSA_nc54 = Sqrt(-1);
    
    //Initialize_smri();
    //if(EventNum%1000==0) init_test = kTRUE;

    //=== Calc. ===--------------------------------------------------------------------------------
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

    tof13T   = Dist_SBTTarget/betaF7F13/clight;
    tof13H   = hodo_t - sbt1_Tslew;
    tofTH_nc = tof13H - tof13T + toff_hodo;
    tofTH    = tofTH_nc + hodo_tofcor[hodo_id-1];

    betaTH  = lengSA_tan/tofTH/clight;
    gammaTH = 1./Sqrt(1.-betaTH*betaTH);

    betaTH54  = leng54/tofTH/clight;
    gammaTH54 = 1./Sqrt(1.-betaTH54*betaTH54);

    
    //t_minoshodo_notcor = hodo_t - sbt1_Tslew - (Dist_SBTTarget/betaF7F13/clight) + toff_hodo;
    t_minoshodo_notcor = hodo_t - sbt1_Tslew - (Dist_SBTTarget/betaF3F13/clight) + toff_hodo;
    //t_minoshodo = t_minoshodo_notcor + hodo_tofcor[hodo_id-1];
    t_minoshodo = hodo_t - sbt1_Tslew - (Dist_SBTTarget/betaF3F13/clight) + toff_hodo + hodo_tofcor[hodo_id-1];
    v_minoshodo = lengSA_tan/t_minoshodo;
    beta_minoshodo  = v_minoshodo/clight;
    gamma_minoshodo = 1./Sqrt(1.-beta_minoshodo*beta_minoshodo);
    //dev = Log(2*me*beta_minoshodo*beta_minoshodo/ionpair) - Log(1 - beta_minoshodo*beta_minoshodo) - beta_minoshodo*beta_minoshodo;

    
    dev = Log(2*me*betaTH*betaTH/ionpair) - Log(1 - betaTH*betaTH) - betaTH*betaTH; 
    //zraw = v_minoshodo*Sqrt(hodo_q/dev);
    zraw = betaTH*clight*Sqrt(hodo_q/dev);

    dev54 = Log(2*me*betaTH54*betaTH54/ionpair) - Log(1 - betaTH54*betaTH54) - betaTH54*betaTH54; 
    zraw54 = betaTH54*clight*Sqrt(hodo_q/dev54);
    
    zetSA = hodo_zraw2z_p0[hodo_id-1] + hodo_zraw2z_p1[hodo_id-1]*zraw;
    //if(hodo_id==9)
    //  zetSA = 0.0060259114 * zraw - 3.8686351311;
    
    //aoqSA_notcor = brhoSA_tan/beta_minoshodo/gamma_minoshodo*clight/mu;
    aoqSA_notcor = brhoSA_tan/betaTH/gammaTH*clight/mu;

    aoqSA_nc54 = brho54/betaTH54/gammaTH54*clight/mu;
    
    aoqSA = hodo_aoqcor[hodo_id-1][0] + hodo_aoqcor[hodo_id-1][1]*aoqSA_notcor;
        
    //zraw = hodo_q - (hodo_t2q0[hodo_id+1] + hodo_t2q1[hodo_id+1]*t_minoshodo);
    //zetSA = hodo_zraw2z0[hodo_id+1] + hodo_zraw2z1[hodo_id+1]*zraw + hodo_zraw2z2[hodo_id+1]*zraw*zraw;

    //===== BG cut =================================================================================
    //=== cut by CUTG ===---------------------------------------------------------------------------
    //if(!csbt1->IsInside(SBT1_TR-SBT1_TL,log(SBT1_QL/SBT1_QR))) BG_flag = 1;
    //if(cSA56Sc_temp->IsInside(aoqSA,zetSA)) SA56Sc_temp = 1;
    
    //=== cut by Hodo Time ===----------------------------------------------------------------------
    for(Int_t i=0;i<24;i++){
      if(Hodoi_TURaw[i]==-1000||Hodoi_TDRaw[i]==-1000) goodEvt_smri = false;
    }
    if(!BG_flag_beam) goodEvt_beam = false;
    if(!BG_flag_dc)   goodEvt_mwdc = false;

    if(goodEvt_smri&&goodEvt_beam) goodEvt = true;
    
    anatrS->Fill();
  }//for LOOP
  anafile_smri->cd();
  anatrS->Write();
  anafile_smri->Close();

  time(&stop);
  cout << endl;
  int t_hour = (int)difftime(stop,start)/3600;
  int t_min  = (int)(difftime(stop,start) - t_hour*3600)/60;
  double t_sec  = difftime(stop,start) - t_hour*3600 - t_min*60;
  printf("Elapsed time: %dh %dm %.1f seconds\n",t_hour,t_min,t_sec);

  double time_end = get_time();
  cout << endl << "Program Run Time " << time_end - time_start << " s." << endl;
  cout << nEntry/(time_end - time_start) << " events/s." << endl;



  
  //printf("%d k events have been treated.\n",nEntry/1000);
  cout << "RUN " << FileNum << ":" << "Conversion finished! : smri" << FileNum << "ok" << endl;
}//main()
