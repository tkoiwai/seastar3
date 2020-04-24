#include "../include/piddef.h"

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

//=====main Function==========================================================
int main(int argc, char *argv[]){

  double time_start = get_time();
  signal(SIGINT,signalhandler);
  
  time_t start, stop;
  time(&start);
  
  gInterpreter->GenerateDictionary("vector<TVector3>","TVector3.h");

  Int_t FileNumber = TString(argv[1]).Atoi();
  
  if(FileNumber==0){
    std::cerr <<  " You should provide either a runnumber" << endl;
  }
  if (argc < 2 || argc > 5){
    printf("Usage: ./hist_PID RUNNUMBER\nOR     ./hist_PID RUNNUMBER MAXEVENTS\nOR     ./hist_PID RUNNUMBER MAXEVENTS TEST\n");
    exit(EXIT_FAILURE); 
  }

  int MaxEventNumber = 0;
  
  if (argc > 2) {
    MaxEventNumber = TString(argv[2]).Atoi();
    printf(" You will process %d events\n",MaxEventNumber);
  }

  //===== Load input files =====
  TString infnameB = Form("rootfiles/ana/beam/ana_beam%04d.root",FileNumber);
  TFile   *infileB = TFile::Open(infnameB);
  TTree   *intrB   = (TTree*)infileB->Get("anatrB");
  PID_Get_Branch_beam(intrB);

  TString infnameS = Form("rootfiles/ana/smri/ana_smri%04d.root",FileNumber);
  TFile   *infileS = TFile::Open(infnameS);
  TTree   *intrS   = (TTree*)infileS->Get("anatrS");
  PID_Get_Branch_smri(intrS);

  TString infnameDC = Form("rootfiles/ana/mwdc_new/anaDC%04d.root",FileNumber);
  TFile   *infileDC = TFile::Open(infnameDC);
  TTree   *intrDC   = (TTree*)infileDC->Get("anatrDC");
  PID_Get_Branch_mwdc(intrDC);

  //TString infnameV = Form("rootfiles/minos/vertex/vertex%04d.root",FileNumber);
  //TFile   *infileV = TFile::Open(infnameV);
  //TTree   *intrV   = (TTree*)infileV->Get("tr");
  //PID_Get_Branch_vertex(intrV);

  TString infnameM = Form("rootfiles/minos/cal/cal_minos%04d.root",FileNumber);
  TFile   *infileM = TFile::Open(infnameM);
  TTree   *intrM   = (TTree*)infileM->Get("caltrM");
  PID_Get_Branch_minos(intrM);
  
  //=== AddFriend ===
  intrB->AddFriend(intrS);
  intrB->AddFriend(intrDC);
  intrB->AddFriend(intrM);

  //=====ROOT file setting==========================================================

  TString ofname;

  if(argc < 4)
    ofname = Form("/home/koiwai/analysis/rootfiles/pid_hist/hist_pid%04d.root",FileNumber);
  else if(argc == 4)
    ofname = Form("/home/koiwai/analysis/macros/testhist_pid%04d.root",FileNumber);

  TFile *outfile = new TFile(ofname,"RECREATE");

  //=====Define variables===================================================

  //===== LOAD CUTS =====================================================================

  TFile *fcutSA_K  = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_K.root");
  TCutG *csa55k    = (TCutG*)fcutSA_K->Get("sa55k");

  TFile *fcutSA_Ca = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca.root");
  TCutG *csa57ca   = (TCutG*)fcutSA_Ca->Get("sa57ca");
  TCutG *csa55ca   = (TCutG*)fcutSA_Ca->Get("sa55ca");
  TCutG *csa54ca   = (TCutG*)fcutSA_Ca->Get("sa54ca");
  TCutG *csa53ca   = (TCutG*)fcutSA_Ca->Get("sa53ca");
  
  //===== DEFINE HIST ====================================================================
  
  char* cnames[10]    = {(char*)"br54ca_saAll",(char*)"br54ca_sa53ca",(char*)"sa55k" ,(char*)"sa55ca",(char*)"sa57ca",(char*)"sa57ca",(char*)"",(char*)"",(char*)"",(char*)""};
  char* cnamesgate[10] = {(char*)"all",(char*)"f5x",(char*)"targetR",(char*)"both",(char*)"",(char*)"",(char*)"",(char*)"",(char*)"",(char*)""};
  char* cnamesminos[2] = {(char*)"wominos",(char*)"wminos"};
  
  TH2F *hpid[20];
 
  for(int i=0;i<1;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<2;k++){
	hpid[(i+0)*10+j*2+k] = new TH2F(Form("h_BR_%s_%s_%s",cnames[i],cnamesgate[j],cnamesminos[k]),Form("h_BR_%s_%s_%s",cnames[i],cnamesgate[j],cnamesminos[k]),500,2.55,2.85,500,16,25); 
	hpid[(i+1)*10+j*2+k] = new TH2F(Form("h_SA_%s_%s_%s",cnames[i+1],cnamesgate[j],cnamesminos[k]),Form("h_SA_%s_%s_%s",cnames[i+1],cnamesgate[j],cnamesminos[k]),500,2.55,2.85,500,16,25); 
      }
    }
  }
  
  
  //=== To check cal_minos PID gates ===

  TH2F *hminosBR = new TH2F("hminosBR","hminosBR",1000,2,3.2,1000,14,27);
  TH2F *hminosSA = new TH2F("hminosSA","hminosSA",1000,2,3.2,1000,14,27);

  //=== whole PID plots ===

  TH2F *hBR = new TH2F("hBR","hBR",1000,2.55,2.85,1000,16,25); 
  TH2F *hSA = new TH2F("hSA","hSA",1000,2.,3.,1000,14,26); 

  //===== LOOP =========================================================================

  Int_t nEntry = intrB->GetEntries();
  int iEntry = 0;
  int AllEntry;
  if(argc > 2 && MaxEventNumber < nEntry)
    AllEntry = MaxEventNumber;
  else
    AllEntry = nEntry;
  
  double time_prev = 0;
  double time_startloop = get_time();
  
  for(Int_t iEntry=0;iEntry<AllEntry;iEntry++){
    
    intrB->GetEntry(iEntry);

    const int showstat = 1000;
    if(iEntry%showstat == 0){	
      double time_end = get_time();	
      double t_diff_a = time_end - time_start;
      int t_hour_a    = t_diff_a/3600;
      int t_min_a     = (t_diff_a - 3600*t_hour_a)/60;
      double t_sec_a  = t_diff_a -3600*t_hour_a - 60*t_min_a;
      
      printf("%dh %dm %.2fs elapsed: %dk events done: %.2f events/s: about %.2fs to go: current speed: %.2f events/s \n",
	     t_hour_a,
	     t_min_a,
	     t_sec_a,
	     (int)(iEntry/1000),
	     iEntry/(time_end - time_startloop),
	     (AllEntry - iEntry)*(time_end - time_startloop)/(double)iEntry,
	     showstat/(time_end - time_prev));	
      time_prev = get_time();
    }
 
    //===== GATES =================================================================

    //===== PID gate ==============================================================

    //===== INIT ==================================================================

    //Double_t vertexZ_cor = vertexZ + MINOSoffsetZ;
    
    //===== FILL HIST =============================================================

    if(br54ca){ //0x
      hpid[0]->Fill(aoqBR,zetBR);
      if(-90<F5X&&F5X<20)
	hpid[2]->Fill(aoqBR,zetBR);
      if(sqrt(Target_X*Target_X+Target_Y*Target_Y)<15)
	hpid[4]->Fill(aoqBR,zetBR);
      if(-90<F5X&&F5X<20&&sqrt(Target_X*Target_X+Target_Y*Target_Y)<15)
	hpid[6]->Fill(aoqBR,zetBR);

      //if(-10<vertexZ_cor&&vertexZ_cor<160){
      if( NumberTracks > 0 ){
	hpid[1]->Fill(aoqBR,zetBR);
	if(-90<F5X&&F5X<20)
	  hpid[3]->Fill(aoqBR,zetBR);
	if(sqrt(Target_X*Target_X+Target_Y*Target_Y)<15)
	  hpid[5]->Fill(aoqBR,zetBR);
	if(-90<F5X&&F5X<20&&sqrt(Target_X*Target_X+Target_Y*Target_Y)<15)
	  hpid[7]->Fill(aoqBR,zetBR);
      }
    }

    if(br54ca&&csa53ca->IsInside(aoqSA,zetSA)){ //1x
      hpid[10]->Fill(aoqSA,zetSA);
      if(-90<F5X&&F5X<20)
	hpid[12]->Fill(aoqSA,zetSA);
      if(sqrt( Target_X*Target_X + Target_Y*Target_Y ) < 15)
	hpid[14]->Fill(aoqSA,zetSA);
      if( -90<F5X && F5X<20 && sqrt( Target_X*Target_X + Target_Y*Target_Y ) < 15)
	hpid[16]->Fill(aoqSA,zetSA);

      //if(-10<vertexZ_cor&&vertexZ_cor<160){
      if( NumberTracks > 0 ){
	hpid[11]->Fill(aoqSA,zetSA);
	if(-90<F5X&&F5X<20)
	  hpid[13]->Fill(aoqSA,zetSA);
	if(sqrt(Target_X*Target_X+Target_Y*Target_Y)<15)
	  hpid[15]->Fill(aoqSA,zetSA);
	if(-90<F5X&&F5X<20&&sqrt(Target_X*Target_X+Target_Y*Target_Y)<15)
	  hpid[17]->Fill(aoqSA,zetSA);
      }
    }

    if(NumberTracks>0){
      hminosBR->Fill(aoqBR,zetBR);
      hminosSA->Fill(aoqSA,zetSA);
    }

    hBR->Fill(aoqBR,zetBR);
    hSA->Fill(aoqSA,zetSA);
    
  }//while loop
  std::clog << std::endl;


  outfile->cd();

  for(int i=0;i<18;i++){
    if(i==8||i==9) continue;
    hpid[i]->Write();
  }
  hminosBR->Write();
  hminosSA->Write();

  hBR->Write();
  hSA->Write();
  
  outfile->Write();
  outfile->Close("R");

  time(&stop);

  int t_hour = (int)difftime(stop,start)/3600;
  int t_min  = (int)(difftime(stop,start) - t_hour*3600)/60;
  double t_sec  = difftime(stop,start) - t_hour*3600 - t_min*60;
  printf("Elapsed time: %dh %dm %.1f seconds\n",t_hour,t_min,t_sec);

  double time_end = get_time();
  printf("Average process speed: %f events/s\n",AllEntry/(time_end - time_start));
  printf("RUN%d: Conversion finished!: histpidok%d\n",FileNumber,FileNumber);

  return 0;
}//main()


inline bool exists_test (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
inline bool exists_test (const TString& name) {
  return ( access( name.Data(), F_OK ) != -1 );
}
