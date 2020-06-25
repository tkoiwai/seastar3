#include "../include/piddef.h"

//&=====main Function==========================================================
int main(int argc, char* argv[]) {
  initiate_timer_tk();

  gInterpreter->GenerateDictionary("vector<TVector3>", "TVector3.h");

  bool TestMode       = false;
  bool ENum_flag      = false;
  int  FileNumber     = -1;
  int  MaxEventNumber = 0;

  struct option longopts[] = {
      {"eventnumber", required_argument, NULL, 'e'},
      {"runnumber", required_argument, NULL, 'r'},
      {"testmode", no_argument, NULL, 't'},
      {0, 0, 0, 0},
  };

  int opt;
  int longindex;

  while((opt = getopt_long(argc, argv, "tr:e:", longopts, &longindex)) != -1) {
    switch(opt) {
      case 'r':
        FileNumber = atoi(optarg);
        break;
      case 'e':
        ENum_flag      = true;
        MaxEventNumber = atoi(optarg);
        break;
      case 't':
        TestMode = true;
        break;

      default:
        break;
    }
  }

  if(argc == 1) {
    printf("\nUsage: ./hist_dali -r <run number> -e <max event number to treat> -t (to activate test mode)\n \n");
    exit(EXIT_FAILURE);
  }

  if(FileNumber == -1) {
    std::cerr << " You should provide a runnumber" << endl;
    exit(EXIT_FAILURE);
  }

  //+===== Load input files =====
  TString infnameB = Form("rootfiles/ana/beam/ana_beam%04d.root", FileNumber);
  TFile*  infileB  = TFile::Open(infnameB);
  TTree*  intrB    = (TTree*)infileB->Get("anatrB");
  PID_Get_Branch_beam(intrB);

  TString infnameS = Form("rootfiles/ana/smri/ana_smri%04d.root", FileNumber);
  TFile*  infileS  = TFile::Open(infnameS);
  TTree*  intrS    = (TTree*)infileS->Get("anatrS");
  PID_Get_Branch_smri(intrS);

  TString infnameDC = Form("rootfiles/ana/mwdc_new/anaDC%04d.root", FileNumber);
  TFile*  infileDC  = TFile::Open(infnameDC);
  TTree*  intrDC    = (TTree*)infileDC->Get("anatrDC");
  PID_Get_Branch_mwdc(intrDC);

  TString infnameM = Form("rootfiles/minos/cal_new/Tracks_run_%04d.root", FileNumber);
  TFile*  infileM  = TFile::Open(infnameM);
  TTree*  intrM    = (TTree*)infileM->Get("tree");
  PID_Get_Branch_minos(intrM);

  intrB->AddFriend(intrS);
  intrB->AddFriend(intrDC);
  intrB->AddFriend(intrM);

  //+=====output ROOT file setting===============================================

  TString ofname;

  if(!TestMode)
    ofname = Form("/home/koiwai/analysis/rootfiles/pid_hist/hist_pid%04d.root", FileNumber);
  else
    ofname = Form("/home/koiwai/analysis/macros/testhist_pid%04d.root", FileNumber);

  TFile* outfile = new TFile(ofname, "RECREATE");

  //+=====Define variables===================================================

  //+===== LOAD/Define CUTS ========================================================

  TFile* fcutSA_K = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_K.root");
  TCutG* csa55k   = (TCutG*)fcutSA_K->Get("sa55k");

  TFile* fcutSA_Ca = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca.root");
  TCutG* csa57ca   = (TCutG*)fcutSA_Ca->Get("sa57ca");
  TCutG* csa55ca   = (TCutG*)fcutSA_Ca->Get("sa55ca");
  TCutG* csa54ca   = (TCutG*)fcutSA_Ca->Get("sa54ca");
  TCutG* csa53ca   = (TCutG*)fcutSA_Ca->Get("sa53ca");

  TFile* fcutSA_Ca_wMINOS      = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_Ca_wMINOS.root");
  TCutG* csa53ca_br54ca_wMINOS = (TCutG*)fcutSA_Ca_wMINOS->Get("csa53ca_wminos");

  TFile* fcutSA_50Ar = TFile::Open("/home/koiwai/analysis/cutfiles/cutSA_50Ar.root");
  TCutG* csa50ar     = (TCutG*)fcutSA_50Ar->Get("sa50ar");

  TCut f5x     = "-90 < F5X && F5X < 20";
  TCut targetR = "sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15.";

  //+===== DEFINE HIST ======================================================

  //char* cnames[10]     = {(char*)"br54ca_saAll", (char*)"br54ca_sa53ca", (char*)"sa55k", (char*)"sa55ca", (char*)"sa57ca", (char*)"sa57ca", (char*)"", (char*)"", (char*)"", (char*)""};
  char* cnames[15] = {
      (char*)"br54ca_saAll",
      (char*)"br51k_saAll",
      (char*)"br56ca_saAll",
      (char*)"br56sc_saAll",
      (char*)"br58sc_saAll",
      (char*)"br59sc_saAll",
      (char*)"br59ti_saAll",
      (char*)"br60ti_saAll",
      (char*)"brAll_sa53ca",
      (char*)"brAll_sa50ar",
      (char*)"brAll_sa55ca",
      (char*)"brAll_sa55k",
      (char*)"brAll_sa57ca",
      (char*)"",
      (char*)""};
  char* cnamesreaction[10] = {
      (char*)"br54ca_sa53ca",
      (char*)"br51k_sa50ar",
      (char*)"br56ca_sa55k",
      (char*)"br56ca_sa55ca",
      (char*)"br56sc_sa55ca",
      (char*)"br58sc_sa57ca",
      (char*)"br59sc_sa57ca",
      (char*)"br59ti_sa57ca",
      (char*)"br60ti_sa57ca",
      (char*)""};
  char* cnamesgate[5]  = {(char*)"all", (char*)"f5x", (char*)"targetR", (char*)"both", (char*)""};
  char* cnamesminos[2] = {(char*)"wominos", (char*)"wminos"};

  TH2F* hpidBR[30];
  TH2F* hpidSA[30];
  TH2F* hpid[30];

  for(int ch = 0; ch < 13; ch++) {
    hpidBR[ch] = new TH2F(
        Form("h_BR_%s", cnames[ch]),
        Form("h_BR_%s", cnames[ch]),
        500, 2.55, 3., 500, 16, 25);
    hpidSA[ch] = new TH2F(
        Form("h_SA_%s", cnames[ch]),
        Form("h_SA_%s", cnames[ch]),
        500, 2., 3., 500, 16, 25);
  }

  for(int i = 0; i < 1; i++) {
    for(int j = 0; j < 4; j++) {
      for(int k = 0; k < 2; k++) {
        hpid[(i + 0) * 10 + j * 2 + k] = new TH2F(
            Form("h_BR_%s_%s_%s", cnames[i], cnamesgate[j], cnamesminos[k]),
            Form("h_BR_%s_%s_%s", cnames[i], cnamesgate[j], cnamesminos[k]),
            500, 2.55, 2.85, 500, 16, 25);
        hpid[(i + 1) * 10 + j * 2 + k] = new TH2F(
            Form("h_SA_%s_%s_%s", cnames[i + 1], cnamesgate[j], cnamesminos[k]),
            Form("h_SA_%s_%s_%s", cnames[i + 1], cnamesgate[j], cnamesminos[k]),
            500, 2.55, 2.85, 500, 16, 25);
        hpid[(i + 2) * 10 + j * 2 + k] = new TH2F(
            Form("h_SA_%s_%s_%s", cnames[i], cnamesgate[j], cnamesminos[k]),
            Form("h_SA_%s_%s_%s", cnames[i], cnamesgate[j], cnamesminos[k]),
            500, 2.55, 2.85, 500, 16, 25);
      }
    }
  }

  //
  TH2F* hminostrack =
      new TH2F("hminostrack", "hminostrack", 230, 0, 230, 10, 0, 10);

  //+=== To check MINOS efficiency ===
  TH2F* hminoseff[6];
  hminoseff[0] = new TH2F(Form("h_minoseff_wminos0_%s", cnamesreaction[0]), Form("for MINOS eff (woMINOS) (%s)", cnamesreaction[0]), 500, 2., 3., 500, 16, 25);
  hminoseff[1] = new TH2F(Form("h_minoseff_wminos1_%s", cnamesreaction[0]), Form("for MINOS eff (NumTrack=1) (%s)", cnamesreaction[0]), 500, 2., 3., 500, 16, 25);
  hminoseff[2] = new TH2F(Form("h_minoseff_wminosall_%s", cnamesreaction[0]), Form("for MINOS eff (NumTrack>=1) (%s)", cnamesreaction[0]), 500, 2., 3., 500, 16, 25);
  hminoseff[3] = new TH2F(Form("h_minoseff_wminos0_%s", cnamesreaction[1]), Form("for MINOS eff (woMINOS) (%s)", cnamesreaction[1]), 500, 2., 3., 500, 16, 25);
  hminoseff[4] = new TH2F(Form("h_minoseff_wminos1_%s", cnamesreaction[1]), Form("for MINOS eff (NumTrack=1or2) (%s)", cnamesreaction[1]), 500, 2., 3., 500, 16, 25);
  hminoseff[5] = new TH2F(Form("h_minoseff_wminosall_%s", cnamesreaction[1]), Form("for MINOS eff (NumTrack>=1) (%s)", cnamesreaction[1]), 500, 2., 3., 500, 16, 25);

  //+=== whole PID plots ===

  TH2F* hBR = new TH2F("hBR", "hBR", 1000, 2.55, 2.85, 1000, 16, 25);
  TH2F* hSA = new TH2F("hSA", "hSA", 1000, 2., 3., 1000, 14, 26);

  //&===== LOOP =================================================================

  Int_t nEntry = intrB->GetEntries();
  int   iEntry = 0;
  int   AllEntry;

  if(ENum_flag && MaxEventNumber < nEntry)
    AllEntry = MaxEventNumber;
  else
    AllEntry = nEntry;

  prepare_timer_tk();

  for(Int_t iEntry = 0; iEntry < AllEntry; iEntry++) {
    intrB->GetEntry(iEntry);

    start_timer_tk(iEntry, AllEntry, 1000);

    //+===== GATES ================================================================

    //+===== PID gate ===========================================================

    //+===== INIT ==================================================================

    //+===== FILL HIST ==========================================================

    bool PIDgates[13] = {
        br54ca,
        br51k,
        br56ca,
        br56sc,
        br58sc,
        br59sc,
        br59ti,
        br60ti,
        csa53ca->IsInside(aoqSA, zetSA),
        csa50ar->IsInside(aoqSA, zetSA),
        csa55ca->IsInside(aoqSA, zetSA),
        csa55k->IsInside(aoqSA, zetSA),
        csa57ca->IsInside(aoqSA, zetSA)};

    for(int i = 0; i < 13; i++) {
      if(PIDgates[i]) {
        hpidBR[i]->Fill(aoqBR, zetBR);
        hpidSA[i]->Fill(aoqSA, zetSA);
      }
    }

    if(br54ca) {  // 0x
      hpid[0]->Fill(aoqBR, zetBR);
      if(-90 < F5X && F5X < 20) hpid[2]->Fill(aoqBR, zetBR);
      if(sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
        hpid[4]->Fill(aoqBR, zetBR);
      if(-90 < F5X && F5X < 20 &&
         sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
        hpid[6]->Fill(aoqBR, zetBR);

      // if(-10<vertexZ_cor&&vertexZ_cor<160){
      if(NumberTracks > 0) {
        hpid[1]->Fill(aoqBR, zetBR);
        if(-90 < F5X && F5X < 20) hpid[3]->Fill(aoqBR, zetBR);
        if(sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
          hpid[5]->Fill(aoqBR, zetBR);
        if(-90 < F5X && F5X < 20 &&
           sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
          hpid[7]->Fill(aoqBR, zetBR);
      }
    }

    if(br54ca && csa53ca_br54ca_wMINOS->IsInside(aoqSA, zetSA)) {  // 1x
      hpid[10]->Fill(aoqSA, zetSA);
      if(-90 < F5X && F5X < 20) hpid[12]->Fill(aoqSA, zetSA);
      if(sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
        hpid[14]->Fill(aoqSA, zetSA);
      if(-90 < F5X && F5X < 20 &&
         sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
        hpid[16]->Fill(aoqSA, zetSA);

      // if(-10<vertexZ_cor&&vertexZ_cor<160){
      if(NumberTracks > 0) {
        hpid[11]->Fill(aoqSA, zetSA);
        if(-90 < F5X && F5X < 20) hpid[13]->Fill(aoqSA, zetSA);
        if(sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
          hpid[15]->Fill(aoqSA, zetSA);
        if(-90 < F5X && F5X < 20 &&
           sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
          hpid[17]->Fill(aoqSA, zetSA);
      }
    }

    // BR54Ca && SAall

    if(br54ca) {  // 1x
      hpid[20]->Fill(aoqSA, zetSA);
      if(-90 < F5X && F5X < 20) hpid[22]->Fill(aoqSA, zetSA);
      if(sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
        hpid[24]->Fill(aoqSA, zetSA);
      if(-90 < F5X && F5X < 20 &&
         sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
        hpid[26]->Fill(aoqSA, zetSA);
      // if(-10<vertexZ_cor&&vertexZ_cor<160){
      if(NumberTracks > 0) {
        hpid[21]->Fill(aoqSA, zetSA);
        if(-90 < F5X && F5X < 20) hpid[23]->Fill(aoqSA, zetSA);
        if(sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
          hpid[25]->Fill(aoqSA, zetSA);
        if(-90 < F5X && F5X < 20 &&
           sqrt(Target_X * Target_X + Target_Y * Target_Y) < 15)
          hpid[27]->Fill(aoqSA, zetSA);
      }
    }

    hBR->Fill(aoqBR, zetBR);
    hSA->Fill(aoqSA, zetSA);

    hminostrack->Fill(FileNumber, NumberTracks);

    if(br54ca && csa53ca->IsInside(aoqSA, zetSA)) {
      hminoseff[0]->Fill(aoqSA, zetSA);
      if(NumberTracks == 1)
        hminoseff[1]->Fill(aoqSA, zetSA);
      if(NumberTracks > 0)
        hminoseff[2]->Fill(aoqSA, zetSA);
    }
    if(br51k && csa50ar->IsInside(aoqSA, zetSA)) {
      hminoseff[3]->Fill(aoqSA, zetSA);
      if(NumberTracks == 1 || NumberTracks == 2)
        hminoseff[4]->Fill(aoqSA, zetSA);
      if(NumberTracks > 0)
        hminoseff[5]->Fill(aoqSA, zetSA);
    }

  }  //- while loop
  std::clog << std::endl;

  outfile->cd();

  for(int ii = 0; ii < 13; ii++) {
    hpidBR[ii]->Write();
    hpidSA[ii]->Write();
  }

  for(int i = 0; i < 28; i++) {
    if(i == 8 || i == 9 || i == 18 || i == 19) continue;
    hpid[i]->Write();
  }

  for(int i = 0; i < 6; i++) {
    hminoseff[i]->Write();
  }

  hBR->Write();
  hSA->Write();

  hminostrack->Write();

  outfile->Write();
  outfile->Close("R");

  stop_timer_tk(FileNumber, AllEntry);

  return 0;
}  // main()

inline bool exists_test(const std::string& name) {
  return (access(name.c_str(), F_OK) != -1);
}
inline bool exists_test(const TString& name) {
  return (access(name.Data(), F_OK) != -1);
}
