////////////////////////////////////////////////////////////////////////////////
//
// Program used for unpacking and SAMURAI (DC) data from SEASTAR2017
//  NB: this will be basically identical to caltree.C, just much neater 
// Author: Frank Browne (frank@ribf.riken.jp)
//
//  2019Nov22: modified by TK.
//
////////////////////////////////////////////////////////////////////////////////
#define ana_DC_cxx
#include "/home/koiwai/analysis/macros/include/dctree.h"

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
  
  Long64_t MaxEventNumber = LLONG_MAX; 
  if (argc < 2){
    printf("Usage: ./dctree_tk RUNNUMBER\nOR\n        ./dctree_tk RUNNUMBER MAXEVENTS\n");
    exit(EXIT_FAILURE); 
  }
  printf("=======================================\n");
  if (argc == 3) {
    MaxEventNumber = TString(argv[2]).Atoi();
    printf(" You will process %lld events\n",MaxEventNumber);
  }

  Int_t FileNumber = TString(argv[1]).Atoi();
  TString RIDFFile = Form("/home/koiwai/analysis/ridf/sdaq02/run%04d.ridf.gz",FileNumber);
  if ( !exists_test(RIDFFile) ){
    cerr << " ERROR - '" << RIDFFile << "' does not exist '" << endl;
    //printf("ERROR - \' %s \' does not exist",RIDFFile);
    exit(EXIT_FAILURE); 
  }
  TArtEventStore *EventStore = new TArtEventStore;
  if (!EventStore->Open(RIDFFile)) {
    cerr << " ERROR - EventStore couldn't open: " << RIDFFile << endl;
    exit(EXIT_FAILURE); 
  }
  
  //cout << " Starting processing of file " << RIDFFile << endl;
  //cout << "--------------------------------------------" << endl;
  printf("Starting processing of file %s\n",RIDFFile.Data());
  printf("----------------------------------------------\n");
  
  // Load parameters
  TArtSAMURAIParameters *SamuraiPara = new TArtSAMURAIParameters(); 
  TArtStoreManager *sman = TArtStoreManager::Instance();  
  SamuraiPara->LoadParameter("/home/koiwai/analysis/db/SAMURAIBDC1.xml");
  SamuraiPara->LoadParameter((char*)"/home/koiwai/analysis/db/SAMURAIBDC2.xml");
  SamuraiPara->LoadParameter((char*)"/home/koiwai/analysis/db/SAMURAIFDC1.xml");
  SamuraiPara->LoadParameter((char*)"/home/koiwai/analysis/db/SAMURAIFDC2.xml");


  TArtCalibBDC1Hit   *CalibBDC1Hit   = new TArtCalibBDC1Hit();
  TArtCalibBDC1Track *CalibBDC1Track = new TArtCalibBDC1Track();
  TArtCalibBDC2Hit   *CalibBDC2Hit   = new TArtCalibBDC2Hit();
  TArtCalibBDC2Track *CalibBDC2Track = new TArtCalibBDC2Track();

  TArtCalibFDC1Hit   *CalibFDC1Hit   = new TArtCalibFDC1Hit();  
  TArtCalibFDC1Track *CalibFDC1Track = new TArtCalibFDC1Track();
  TArtCalibFDC2Hit   *CalibFDC2Hit   = new TArtCalibFDC2Hit();  
  //TArtCalibFDC2Track *CalibFDC2Track = new TArtCalibFDC2Track();
  TArtCalibFDC2Track *CalibFDC2Track = new TArtCalibFDC2Track();

  // Creating output
  //TString ROOTFile = Form("/home/koiwai/analysis/macros/testbench/anaDC%04d.root",FileNumber);
  //TString ROOTFile = Form("/home/koiwai/analysis/rootfiles/ana/mwdc/anaDC%04d.root",FileNumber);
  TString ROOTFile = Form("/home/koiwai/analysis/macros/testanaDC%04d.root",FileNumber);
  TFile * outf = new TFile(ROOTFile,"recreate");
  TTree * tree_DC = new TTree("anatrDC","Unpacked DC tree");

  SetBranches(tree_DC);
  //  TString DCTDCdistfile = Form("/data1/SEASTAR/rootfiles/dcdist/run%d_dcdist.root",FileNumber);
  char DCTDCdistfile[200];
  sprintf(DCTDCdistfile,"/home/koiwai/analysis/rootfiles/tdc_dist/ana_dc_tdcdist.root");

  LoadDCTDCDistributionFull( (char*)DCTDCdistfile,CalibBDC1Track,CalibBDC2Track,CalibFDC1Track,CalibFDC2Track );

  int iEntry;

  double time_prev = 0;
  double time_startloop = get_time();
  
  while (EventStore->GetNextEvent() && EventNumber < MaxEventNumber) {
    iEntry = EventNumber;
    EventNumber++;
    //if (EventNumber%100000 == 0)
    //clog << " ANALYSIS-Info: " << EventNumber << " events processed"<< endl;
    
    if(iEntry%1000 == 0){
	
	double time_end = get_time();
	
	double t_diff_a = time_end - time_start;
	int t_hour_a    = t_diff_a/3600;
	int t_min_a     = (t_diff_a - 3600*t_hour_a)/60;
	double t_sec_a  = t_diff_a -3600*t_hour_a - 60*t_min_a;
	

	printf("%dh %dm %.2fs elapsed: %dk events done: %.2f events/s: about %.2fs to go: current speed: %.2f events/s \n",t_hour_a,t_min_a,t_sec_a,(int)(iEntry/1000),iEntry/(time_end - time_startloop),(100000 - iEntry)*(time_end - time_startloop)/(double)iEntry,1000./(time_end - time_prev));
	
	time_prev = get_time();
    }
      

    CalibBDC1Hit->ClearData();
    CalibBDC1Track->ClearData();
    CalibBDC2Hit->ClearData();
    CalibBDC2Track->ClearData();

    CalibFDC1Hit->ClearData();
    CalibFDC1Track->ClearData();
    CalibFDC2Hit->ClearData();
    CalibFDC2Track->ClearData();

    // Reading event segments 
    TArtRawEventObject *fEvent = (TArtRawEventObject*)sman->FindDataContainer("RawEvent");
    
    for (int i=0; i<fEvent->GetNumSeg(); i++){
      TArtRawSegmentObject * seg = fEvent->GetSegment(i);
      Int_t device   = seg->GetDevice();
      Int_t detector = seg->GetDetector();

      if(device == SAMURAI){
      	switch(detector){
	case       BDC: CalibBDC1Hit->LoadData(seg); break;
	case      FDC1: CalibFDC1Hit->LoadData(seg); break;
	case      FDC2: CalibFDC2Hit->LoadData(seg); break;
      	default: break;
      	}
      }
      TClonesArray *EventInfo = (TClonesArray*)sman->FindDataContainer("EventInfo");
      CoincidenceTrigger = ((TArtEventInfo *)EventInfo->At(0))->GetTriggerBit();
      RunNumber = ((TArtEventInfo *)EventInfo->At(0))->GetRunNumber();
      
    }
    ProcBDC(sman, CalibBDC1Hit, CalibBDC1Track, CalibBDC2Hit, CalibBDC2Track);
    ProcFDC(sman, CalibFDC1Hit, CalibFDC1Track, CalibFDC2Hit, CalibFDC2Track);    
    tree_DC->Fill();
  }
  
  outf->cd();
  tree_DC->Write();
  outf->Close();
  time(&stop);
  
  int t_hour = (int)difftime(stop,start)/3600;
  int t_min  = (int)(difftime(stop,start) - t_hour*3600)/60;
  double t_sec  = difftime(stop,start) - t_hour*3600 - t_min*60;
  printf("Elapsed time: %dh %dm %.1f seconds\n",t_hour,t_min,t_sec);
  
  double time_end = get_time();
  printf("Average process speed: %f events/s\n",100000/(time_end - time_start));
  //cout << endl << "Program Run Time " << time_end - time_start << " s." << endl;
  //cout << nEntry/(time_end - time_start) << " events/s." << endl;

  //cout << "End of beamtree for Run: " <<  RunNumber << endl;
  printf("RUN%d: Conversion finished!: mwdcok%d\n",FileNumber,FileNumber);
  printf("==================================================\n");
  cout << "============================================" << endl;

  
  return 0; 
}
inline bool exists_test (const std::string& name) {
  return ( access( name.c_str(), F_OK ) != -1 );
}
inline bool exists_test (const TString& name) {
  return ( access( name.Data(), F_OK ) != -1 );
}

void ProcBDC(TArtStoreManager * sman, TArtCalibBDC1Hit * CalibBDC1Hit, TArtCalibBDC1Track * CalibBDC1Track,
	                              TArtCalibBDC2Hit * CalibBDC2Hit, TArtCalibBDC2Track * CalibBDC2Track)
{
  //BDC1
  BDC1_X = NAN;
  BDC1_Y = NAN;
  BDC1_ThetaX = NAN;
  BDC1_ThetaY = NAN;
  CalibBDC1Hit->ReconstructData();
  CalibBDC1Track->ReconstructData();

  TClonesArray *BDC1Tracks = (TClonesArray *)sman->FindDataContainer("SAMURAIBDC1Track");
  if (BDC1Tracks) {
    Int_t BDC1NumberOfTracks = BDC1Tracks->GetEntries();
    Double_t TempXPosition, TempYPosition, TempChi2, MinChi2x = 1e6, MinChi2y = 1e6;

    if (BDC1NumberOfTracks > 0){
      TArtDCTrack *TrackBDC1;

      for(Int_t i=0; i<BDC1NumberOfTracks; i++){
	TrackBDC1 = (TArtDCTrack *)BDC1Tracks->At(i);

	if(TrackBDC1){
	  TempXPosition = TrackBDC1->GetPosition(0);
	  TempYPosition = TrackBDC1->GetPosition(1);
	  TempChi2 = TrackBDC1->GetChi2() / (Double_t)TrackBDC1->GetNDF();

	  if(TempChi2 > 0 ){

	    if(TMath::Abs(TempXPosition) < 5000 && TempChi2<MinChi2x){
	      BDC1_X = TempXPosition;
	      BDC1_ThetaX = TMath::ATan(TrackBDC1->GetAngle(0));
	      MinChi2x = TempChi2;
	    }

	    if (TMath::Abs(TempYPosition) < 5000 && TempChi2<MinChi2y){
	      BDC1_Y = TempYPosition;
	      BDC1_ThetaY = TMath::ATan(TrackBDC1->GetAngle(1));
	      MinChi2y = TempChi2;
	    }
	  }
	}
      }
    }
  }

  //BDC2
  BDC2_X = NAN;
  BDC2_Y = NAN;
  BDC2_ThetaX = NAN;
  BDC2_ThetaY = NAN;

  CalibBDC2Hit->ReconstructData();
  CalibBDC2Track->ReconstructData();

  TClonesArray *BDC2Tracks = (TClonesArray *)sman->FindDataContainer("SAMURAIBDC2Track");
  if (BDC2Tracks) {
    Int_t BDC2NumberOfTracks = BDC2Tracks->GetEntries();
    Double_t TempXPosition, TempYPosition, TempChi2, MinChi2x=1e6, MinChi2y=1e6;

    if (BDC2NumberOfTracks > 0){
      TArtDCTrack *TrackBDC2;

      for(Int_t i=0; i<BDC2NumberOfTracks; i++){
	TrackBDC2 = (TArtDCTrack *)BDC2Tracks->At(i);

	if(TrackBDC2){
	  TempXPosition = TrackBDC2->GetPosition(0);
	  TempYPosition = TrackBDC2->GetPosition(1);
	  TempChi2 = TrackBDC2->GetChi2() / (Double_t)TrackBDC2->GetNDF();

	  if(TempChi2 > 0 ){
	    if(TMath::Abs(TempXPosition) < 5000 && TempChi2<MinChi2x){
	      BDC2_X = TempXPosition;
	      BDC2_ThetaX = TMath::ATan(TrackBDC2->GetAngle(0));
	      MinChi2x = TempChi2;
	    }
	    if (TMath::Abs(TempYPosition) < 5000 && TempChi2<MinChi2y){
	      BDC2_Y = TempYPosition;
	      BDC2_ThetaY = TMath::ATan(TrackBDC2->GetAngle(1));
	      MinChi2y = TempChi2;
	    }
	  }
	}
      }
    }
  }

  BDC2_X -= 0.5096;
  BDC2_Y -= 0.4481;
  
  //TARGET
  Target_X = NAN;
  Target_Y = NAN;
  Target_ThetaX = NAN;
  Target_ThetaY = NAN;

  Target_X = BDC1_X + Dist_BDC1Target / Dist_BDC1BDC2 * (BDC2_X - BDC1_X);
  Target_Y = BDC1_Y + Dist_BDC1Target / Dist_BDC1BDC2 * (BDC2_Y - BDC1_Y);
  Target_ThetaX = TMath::ATan( (BDC2_X - BDC1_X) / Dist_BDC1BDC2 );
  Target_ThetaY = TMath::ATan( (BDC2_Y - BDC1_Y) / Dist_BDC1BDC2 );

  //  tree_DC ->Fill();
} // END BDC ANALYSIS

void ProcFDC(TArtStoreManager * sman, TArtCalibFDC1Hit * CalibFDC1Hit, TArtCalibFDC1Track * CalibFDC1Track,
	     TArtCalibFDC2Hit * CalibFDC2Hit, TArtCalibFDC2Track * CalibFDC2Track)
{
  // ============== FDC1 ================
  FDC1_X = NAN;
  FDC1_Y = NAN;
  FDC1_ThetaX = NAN;
  FDC1_ThetaY = NAN;

  CalibFDC1Hit->ReconstructData();
  CalibFDC1Track->ReconstructData();

  TClonesArray *FDC1Tracks = (TClonesArray *)sman->FindDataContainer((char*)"SAMURAIFDC1Track");
  Int_t FDC1NumberOfTracks = -1;
  Double_t Chi2_FDC1 = -1.;

  if(FDC1Tracks) {
    FDC1NumberOfTracks = FDC1Tracks->GetEntries();
    if(FDC1NumberOfTracks > 0) { //Tracks are already sorted by chi2
      TArtDCTrack *FDC1Track = (TArtDCTrack *)FDC1Tracks->At(0);
      if(FDC1Track) {
	FDC1_X = FDC1Track->GetPosition(0);
	FDC1_Y = FDC1Track->GetPosition(1);
	Chi2_FDC1 = FDC1Track->GetChi2() / (Double_t)FDC1Track->GetNDF();
      }

      if( FDC1NumberOfTracks > 0 &&
	  TMath::Abs(Target_X)<100 && // 100mm for safety (can be used in most of the experiment ! to be ajusted ?
	  TMath::Abs(Target_Y)<100 && // 100mm for safety (can be used in most of the experiment ! to be ajusted ?
	  TMath::Abs(FDC1_X) < 5000 && TMath::Abs(FDC1_Y) < 5000 ){
	FDC1_ThetaX = TMath::ATan((FDC1_X - Target_X) / (Dist_BDC1FDC1-Dist_BDC1Target));
	FDC1_ThetaY = TMath::ATan((FDC1_Y - Target_Y) / (Dist_BDC1FDC1-Dist_BDC1Target));
      }
    }
  }

  // ============== FDC2 ================
  FDC2_X = NAN;
  FDC2_Y = NAN;
  FDC2_ThetaX = NAN;
  FDC2_ThetaY = NAN;

  CalibFDC2Hit->ReconstructData();
  CalibFDC2Track->ReconstructData();

  TClonesArray *FDC2Tracks = (TClonesArray *)sman->FindDataContainer((char*)"SAMURAIFDC2Track");
  Int_t FDC2NumberOfTracks = -1;
  if(FDC2Tracks) {
    FDC2NumberOfTracks = FDC2Tracks->GetEntries();
    if(FDC2NumberOfTracks > 0) { //Tracks are already sorted by chi2
      TArtDCTrack *FDC2Track = (TArtDCTrack *)FDC2Tracks->At(0);
      if(FDC2Track) {
	FDC2_X = FDC2Track->GetPosition(0);
	FDC2_Y = FDC2Track->GetPosition(1);
	FDC2_ThetaX = TMath::ATan(FDC2Track->GetAngle(0));
	FDC2_ThetaY = TMath::ATan(FDC2Track->GetAngle(1));
      }
    }
  }

} // END FDCs ANALYSIS


void LoadDCTDCDistributionFull(char *FileName, TArtCalibBDC1Track *CalibBDC1Track, TArtCalibBDC2Track *CalibBDC2Track,
			       TArtCalibFDC1Track *CalibFDC1Track, TArtCalibFDC2Track *CalibFDC2Track) {

  TFile *RootFile = new TFile(FileName,"READ");

  if(RootFile) {
    gROOT->cd();
    TH1I *Hist1D = NULL;
    TH1F *Hist2D = NULL;
    Int_t BDCNumberOfLayers = 8;
    Int_t FDCNumberOfLayers = 14;

    for(Int_t i=0; i<BDCNumberOfLayers; i++) {
      Hist1D = (TH1I*) RootFile->Get(Form("hbdc1tdc%d",i));
     
      if(Hist1D) {
	CalibBDC1Track->SetTDCDistribution(Hist1D,i);
	delete Hist1D;
	Hist1D = NULL;
	delete Hist2D;
	Hist2D = NULL;
      }
      else
	cout << "\e[35m " << "Warning LoadDCTDCDistribution :: Could not find the following histogram " << Form("BDC1TDCDistLayer%d",i) << "\e[37m" << endl;
    }

    for(Int_t i=0; i<BDCNumberOfLayers; i++) {
      Hist1D = (TH1I*) RootFile->Get(Form("hbdc2tdc%d",i));

      if(Hist1D) {
	CalibBDC2Track->SetTDCDistribution(Hist1D,i);
	delete Hist1D;
	Hist1D = NULL;
	delete Hist2D;
	Hist2D = NULL;
	
      } else
	cout << "\e[35m " << "Warning LoadDCTDCDistribution :: Could not find the following histogram " << Form("BDC2TDCDistLayer%d",i) << "\e[37m" << endl;
    }

    for(Int_t i=0; i<FDCNumberOfLayers; i++) {
      Hist1D = (TH1I*) RootFile->Get(Form("hfdc1tdc%d",i));

      if(Hist1D) {
     	CalibFDC1Track->SetTDCDistribution(Hist1D,i);

     	delete Hist1D;
     	Hist1D = NULL;
     	delete Hist2D;
     	Hist2D = NULL;

      } else
     	cout << "\e[35m " << "Warning LoadDCTDCDistribution :: Could not find the following histogram " << Form("FDC1TDCDistLayer%d",i) << "\e[37m" << endl;
    }

    for(Int_t i=0; i<FDCNumberOfLayers; i++) {
      Hist1D = (TH1I*) RootFile->Get(Form("hfdc2tdc%d",i));
      if(Hist1D) {
	CalibFDC2Track->SetTDCDistribution(Hist1D,i);
	delete Hist1D;
	Hist1D = NULL;
	delete Hist2D;
	Hist2D = NULL;
      }
      else
	cout << "\e[35m " << "Warning LoadDCTDCDistribution :: Could not find the following histogram " << Form("FDC2TDCDistLayer%d",i) << "\e[37m" << endl;
    }
  } else
    cout << "\e[35m " << "Could not find the following file : " << FileName << "\e[37m" << endl;

  delete RootFile;
}

void SetBranches(TTree *tree_DC)
{

  tree_DC->Branch("EventNumber",&EventNumber);
  tree_DC->Branch("RunNumber",&RunNumber);

  // BDC
  tree_DC->Branch("BDC1_X",&BDC1_X);
  tree_DC->Branch("BDC1_Y",&BDC1_Y);
  tree_DC->Branch("BDC1_ThetaX",&BDC1_ThetaX);
  tree_DC->Branch("BDC1_ThetaY",&BDC1_ThetaY);
  tree_DC->Branch("BDC2_X",&BDC2_X);
  tree_DC->Branch("BDC2_Y",&BDC2_Y);
  tree_DC->Branch("BDC2_ThetaX",&BDC2_ThetaX);
  tree_DC->Branch("BDC2_ThetaY",&BDC2_ThetaY);
  tree_DC->Branch("Target_X",&Target_X);
  tree_DC->Branch("Target_Y",&Target_Y);
  tree_DC->Branch("Target_ThetaX",&Target_ThetaX);
  tree_DC->Branch("Target_ThetaY",&Target_ThetaY);

  // FDC
  tree_DC->Branch("FDC1_X",&FDC1_X);
  tree_DC->Branch("FDC1_Y",&FDC1_Y);
  tree_DC->Branch("FDC1_ThetaX",&FDC1_ThetaX);
  tree_DC->Branch("FDC1_ThetaY",&FDC1_ThetaY);
  tree_DC->Branch("FDC2_X",&FDC2_X);
  tree_DC->Branch("FDC2_Y",&FDC2_Y);
  tree_DC->Branch("FDC2_ThetaX",&FDC2_ThetaX);
  tree_DC->Branch("FDC2_ThetaY",&FDC2_ThetaY);
  
}

