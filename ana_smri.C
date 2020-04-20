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
//#include"/home/koiwai/analysis/brho_func/Brho_A54Z20_br54Ca_sa54Ca.C"
//#include"/home/koiwai/analysis/brho_func/Len_A54Z20_br54Ca_sa54Ca.C"

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

  int tmpEntry = 2000000;

  if (argc < 2){
    printf("Usage: ./ana_smri_test RUNNUMBER\nOR\n        ./ana_smri_test RUNNUMBER MAXEVENTS\n");
    exit(EXIT_FAILURE); 
  }
  printf("=======================================\n");
  if (argc == 3) {
    tmpEntry = TString(argv[2]).Atoi();    
    printf(" You will process %d events\n",tmpEntry);
  }
  
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
  TString infnameDC = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/anaDC%04d.root",FileNum);
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
  Double_t ionpair = 15796.0;

  
  //===== Declare anatree const.s ===============================================================
  Double_t hodo_toff[24], hodo_qcor[24];
  Double_t hodo_zraw2z_p0[24], hodo_zraw2z_p1[24];// hodo_zraw2z2[24];
  Double_t hodo_zraw2z[2];
  for(Int_t id=0;id<24;id++){
    hodo_toff[id] = env_hodot->GetValue(Form("hodo_toff_%02d",id+1),0.0);
    hodo_qcor[id] = env_hodoq->GetValue(Form("hodo_qcor_%02d",id+1),0.0);
    hodo_zraw2z_p0[id] = env_hodozraw2z->GetValue(Form("zraw2z_%02d_p0",id+1),0.0);
    hodo_zraw2z_p1[id] = env_hodozraw2z->GetValue(Form("zraw2z_%02d_p1",id+1),1.0);
  }
  Double_t hodo_aoqcor[24][2];
  for(Int_t i=0;i<24;i++){
    hodo_aoqcor[i][0] = env_hodoaoqcor->GetValue(Form("%02dp0",i+1),0.0);
    hodo_aoqcor[i][1] = env_hodoaoqcor->GetValue(Form("%02dp1",i+1),1.0);
  }

  Double_t hodo_tofcor[24];
  const char *nn = Form("%d",FileNum);
  for(int id=0;id<24;++id){
    hodo_tofcor[id] = env_hodotofcor[id]->GetValue(nn,0.0);
  }
  Double_t hodo_qoff[24] = {
       0,-500,-500,-500,-500,-780,
    -500,-500,-700,-500,-500,-600,
    -600,-600,-500,-750,-500,-400,
    -250,-350,-200,-150,-200,   0
  };

  Double_t hodo_slew[24][3] = {
           {0,0,0}, {-50,-50,2.48}, {-50,-50,2.38}, {-50,-50,2.29}, {-55,-55,2.56}, {-80,-80,3.30},
    {-50,-50,2.31}, {-50,-50,2.29}, {-55,-55,2.47}, {-50,-50,2.24}, {-55,-55,2.50}, {-60,-60,2.70},
    {-50,-50,2.22}, {-55,-55,2.40}, {-60,-60,2.59}, {-60,-60,2.64}, {-70,-70,2.95}, {-60,-60,2.49},
    {-60,-60,2.44}, {-60,-60,2.42}, {-60,-60,2.40}, {-65,-65,2.57}, {-90,-90,3.58}, {0,0,0}
  };

  //===== Declare valables for calc. =============================================================
  Double_t dev;

  Set_Branch_smri(anatrS);
 
  //===== Begin LOOP =============================================================================
  double time_prev = 0.;
  int AllEntry = caltr->GetEntries();
  int nEntry = 0;
  if(AllEntry>tmpEntry)
    nEntry = tmpEntry;
  else
    nEntry = AllEntry;

  printf("Number of events to treat: %d\n",nEntry);

  double time_startloop = get_time();
  
  for(int iEntry=0;iEntry<nEntry;++iEntry){
    
    if(iEntry%1000 == 0){

      double time_end = get_time();

      double t_diff_a = time_end - time_start;
      int t_hour_a    = t_diff_a/3600;
      int t_min_a     = (t_diff_a - 3600*t_hour_a)/60;
      double t_sec_a  = t_diff_a -3600*t_hour_a - 60*t_min_a;
      
      printf("%dh %dm %.2fs elapsed: %.1f%% (%dk events) done: %.2f events/s: %.2fs to go: current speed: %.2f events/s \n",t_hour_a,t_min_a,t_sec_a,(100.*iEntry)/nEntry,(int)(iEntry/1000),iEntry/(time_end - time_startloop),(nEntry - iEntry)*(time_end - time_startloop)/(double)iEntry,1000./(time_end - time_prev));
      
      time_prev = get_time();
    }

    caltr->GetEntry(iEntry);

    RunNum = RunNumber_all;
    EventNum = EventNumber_all;    
 
    //=== Initialize ===-------------------------------------------------------------------------
    hodo_q     = 0.;
    hodo_t     = Sqrt(-1);
    hodo_id    = 0;
    hodo_multi = 0;
    
    aoqSA        = Sqrt(-1);
    zetSA    = Sqrt(-1);
    zraw     = Sqrt(-1);
    dev      = Sqrt(-1);

    Double_t allHodo_Q[24];
    Double_t allHodo_T[24];
     for(Int_t i=0;i<24;i++){
      allHodo_Q[i] = Sqrt(-1);
      allHodo_T[i] = Sqrt(-1);
    }

    tof13T   = Sqrt(-1);
    tof13H   = Sqrt(-1);
    tofTH_nc = Sqrt(-1);
    tofTH    = Sqrt(-1);

    betaTH = Sqrt(-1);
    gammaTH = Sqrt(-1);

    brhoSA     = Sqrt(-1);
    lengSA     = Sqrt(-1);

    BG_flag = 0;
    goodEvt      = false;
    goodEvt_beam = true;
    goodEvt_smri = true;
    goodEvt_mwdc = true;

    //=== Calc. ===--------------------------------------------------------------------------------
    //@@@ Brho func. @@@
    
    Double_t x[6];
    x[0] = FDC1_X;
    x[1] = FDC1_A*1000.;
    x[2] = FDC1_Y;
    x[3] = FDC1_B*1000.;
    x[4] = FDC2_X;
    x[5] = Tan(FDC2_A)*1000.;

    brhoSA = MDF_Brho_A56Z20(x);
    lengSA = MDF_Len_A56Z20(x);

    //@@@ Hodo @@@
    
    for(Int_t i=0;i<24;i++){
      allHodo_Q[i] = Hodoi_QCal[i]*hodo_qcor[i];
      allHodo_T[i] = Hodoi_TCal[i]+hodo_toff[i];
      if(allHodo_Q[i]>hodo_q){
	hodo_q  = allHodo_Q[i];
	hodo_t  = allHodo_T[i];
	hodo_t += hodo_slew[i][0]/sqrt(Hodoi_QUCal[i]-300.)
	        + hodo_slew[i][1]/sqrt(Hodoi_QDCal[i]-300.) + hodo_slew[i][2];
	hodo_id = i+1;
      }
    }
    hodo_multi = Hodo_Multiplicity;
    
    //@@@ TOF @@@
    tof13T   = Dist_SBTTarget/betaF7F13/clight;
    tof13H   = hodo_t - sbt1_Tslew;
    tofTH_nc = tof13H - tof13T + toff_hodo;
    tofTH    = tofTH_nc + hodo_tofcor[hodo_id-1];

    betaTH  = lengSA/tofTH/clight;
    gammaTH = 1./Sqrt(1.-betaTH*betaTH);
    
    dev = Log(ionpair*betaTH*betaTH) - Log(1 - betaTH*betaTH) - betaTH*betaTH; 
    zraw = betaTH*Sqrt((hodo_q+hodo_qoff[hodo_id-1])/dev);
    zetSA = hodo_zraw2z_p0[hodo_id-1] + hodo_zraw2z_p1[hodo_id-1]*zraw;

    
    aoqSA = brhoSA/betaTH/gammaTH*clight/mu;

    //===== BG cut =================================================================================
    //=== cut by CUTG ===---------------------------------------------------------------------------
    //=== cut by Hodo Time ===----------------------------------------------------------------------
    for(Int_t i=0;i<24;i++){
      if(Hodoi_TURaw[i]==-1000||Hodoi_TDRaw[i]==-1000) goodEvt_smri = false;
    }
    if(!BG_flag_beam) goodEvt_beam = false;
    //if(!BG_flag_dc)   goodEvt_mwdc = false;

    if(goodEvt_smri&&goodEvt_beam) goodEvt = true;
    
    anatrS->Fill();
  }//for LOOP
  anafile_smri->cd();
  anatrS->Write();
  anafile_smri->Close();

  time(&stop);

  int t_hour = (int)difftime(stop,start)/3600;
  int t_min  = (int)(difftime(stop,start) - t_hour*3600)/60;
  double t_sec  = difftime(stop,start) - t_hour*3600 - t_min*60;
  printf("Elapsed time: %dh %dm %.1f seconds\n",t_hour,t_min,t_sec);

  double time_end = get_time();
  printf("Average process speed: %f events/s\n",nEntry/(time_end - time_start));
  printf("RUN%d: Conversion finished!: smriok%d\n",FileNum,FileNum);
}//main()
