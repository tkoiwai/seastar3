#include <TSystem.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/stat.h> //get the status of files. "st_"
#include <sys/types.h>
#include <unistd.h> //UNIx STanDard Header file
#include <climits> //char limits

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TEnv.h"


void makechain(){
  Int_t a;
  vector<Int_t> run;

  ifstream ifs("/home/koiwai/analysis/physics_runs.txt");
  while(ifs >> a){
    run.push_back(a);
    //cout << run.back() << endl;
  }
  //ifs.Close();
  
  TFile *f = new TFile("calch.root","recreate");
  TChain *ch = new TChain("caltr","caltree");
  /*
  for(Int_t i=0;i<run.size();++i){
    ch->Add("Form(./rootfiles/run%04d_ALL.root)",run[i]);
  }
  */
  ch->Add("/home/koiwai/analysis/strage2018feb20/*.root");
  
  f->cd();
  ch->Write();
  f->Close("R");
}
