{
  TChain * ch = new TChain("AnalysisTree");
  const int nruns = 181;
  //                 1               5                  10     
  int runs[nruns] ={32, 36, 37, 38, 40, 41, 42, 43, 44, 45, 
		    46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 
		    56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 
		    67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 
		    78, 79, 80, 81, 82, 83, 84, 85, 87, //50
		    88, 90, 91, 92, 93, 94, 95, 96 ,97, 
                    98, 99, 100,101,102,103,104,105,106,107, //70
                    108,109,110,111,112,113,114,115,116,117,
                    119,120,121,122,123,124,125,126,127,128,
                    129,130,131,132,133,134,135,136,137,138, //100
                    140,141,142,143,144,145,146,147,148,149,
		    150,151,152,153,154,155,156,157,158,159,
		    160,161,162,163,164,165,166,168,169,170,
		    174,175,179,180,181,182,183,184,185,186,
                    190,191,192,193,194,195,196,197,198,199, //150
		    200,201,202,203,204,205,206,207,208,209,
                    210,211,212,213,214,215,216,217,218,219,
                    220,221,222,223,224,225,226,227,228,229,
                    230};

  int truns = nruns;

  for(int i=0;i<truns;i++)
    ch->Add(Form("/mnt/raid/SEASTAR/home/exp/exp1705_psp17/anaroot/users/wimmer/rootA/run%04d_analysed.root",runs[i]));

  /*
  TFile *brcuts = new TFile("brcuts.root","");
  TCutG* brcut[50];
  int j=0;
  brcut[j] = (TCutG*)brcuts->Get("in64V"); brcut[j++]->SetName("in64V");
  brcut[j] = (TCutG*)brcuts->Get("in63V"); brcut[j++]->SetName("in63V");
  brcut[j] = (TCutG*)brcuts->Get("in62V"); brcut[j++]->SetName("in62V");
  brcut[j] = (TCutG*)brcuts->Get("in61V"); brcut[j++]->SetName("in61V");

  brcut[j] = (TCutG*)brcuts->Get("in62Ti");brcut[j++]->SetName("in62Ti");
  brcut[j] = (TCutG*)brcuts->Get("in61Ti");brcut[j++]->SetName("in61Ti");
  brcut[j] = (TCutG*)brcuts->Get("in60Ti");brcut[j++]->SetName("in60Ti");
  brcut[j] = (TCutG*)brcuts->Get("in59Ti");brcut[j++]->SetName("in59Ti");
  brcut[j] = (TCutG*)brcuts->Get("in58Ti");brcut[j++]->SetName("in58Ti");

  brcut[j] = (TCutG*)brcuts->Get("in59Sc");brcut[j++]->SetName("in59Sc");
  brcut[j] = (TCutG*)brcuts->Get("in58Sc");brcut[j++]->SetName("in58Sc");
  brcut[j] = (TCutG*)brcuts->Get("in57Sc");brcut[j++]->SetName("in57Sc");
  brcut[j] = (TCutG*)brcuts->Get("in56Sc");brcut[j++]->SetName("in56Sc");
  brcut[j] = (TCutG*)brcuts->Get("in55Sc");brcut[j++]->SetName("in55Sc");
  //brcut[j] = (TCutG*)brcuts->Get("in54Sc");brcut[j++]->SetName("in54Sc");

  brcut[j] = (TCutG*)brcuts->Get("in57Ca");brcut[j++]->SetName("in57Ca");
  brcut[j] = (TCutG*)brcuts->Get("in56Ca");brcut[j++]->SetName("in56Ca");
  brcut[j] = (TCutG*)brcuts->Get("in55Ca");brcut[j++]->SetName("in55Ca");
  brcut[j] = (TCutG*)brcuts->Get("in54Ca");brcut[j++]->SetName("in54Ca");
  brcut[j] = (TCutG*)brcuts->Get("in53Ca");brcut[j++]->SetName("in53Ca");
  brcut[j] = (TCutG*)brcuts->Get("in52Ca");brcut[j++]->SetName("in52Ca");

  brcut[j] = (TCutG*)brcuts->Get("in53K"); brcut[j++]->SetName("in53K");
  brcut[j] = (TCutG*)brcuts->Get("in52K"); brcut[j++]->SetName("in52K");
  brcut[j] = (TCutG*)brcuts->Get("in51K"); brcut[j++]->SetName("in51K");

  brcut[j] = (TCutG*)brcuts->Get("in51Ar");brcut[j++]->SetName("in51Ar");
  brcut[j] = (TCutG*)brcuts->Get("in50Ar");brcut[j++]->SetName("in50Ar");
  brcut[j] = (TCutG*)brcuts->Get("in49Ar");brcut[j++]->SetName("in49Ar");

  brcut[j] = (TCutG*)brcuts->Get("in49Cl");brcut[j++]->SetName("in49Cl");
  brcut[j] = (TCutG*)brcuts->Get("in48Cl");brcut[j++]->SetName("in48Cl");

  cout << j << " incoming BR cuts found " << endl;
  int nbr = j;
  TFile *sacuts = new TFile("sacuts.root","");
  TCutG* sacut[50];
  j=0;
  sacut[j] = (TCutG*)sacuts->Get("out62Ti");sacut[j++]->SetName("out62Ti");
  sacut[j] = (TCutG*)sacuts->Get("out61Ti");sacut[j++]->SetName("out61Ti");
  sacut[j] = (TCutG*)sacuts->Get("out60Ti");sacut[j++]->SetName("out60Ti");
  sacut[j] = (TCutG*)sacuts->Get("out59Ti");sacut[j++]->SetName("out59Ti");
  sacut[j] = (TCutG*)sacuts->Get("out58Ti");sacut[j++]->SetName("out58Ti");
  sacut[j] = (TCutG*)sacuts->Get("out57Ti");sacut[j++]->SetName("out57Ti");

  sacut[j] = (TCutG*)sacuts->Get("out61Sc");sacut[j++]->SetName("out61Sc");
  sacut[j] = (TCutG*)sacuts->Get("out60Sc");sacut[j++]->SetName("out60Sc");
  sacut[j] = (TCutG*)sacuts->Get("out59Sc");sacut[j++]->SetName("out59Sc");
  sacut[j] = (TCutG*)sacuts->Get("out58Sc");sacut[j++]->SetName("out58Sc");
  sacut[j] = (TCutG*)sacuts->Get("out57Sc");sacut[j++]->SetName("out57Sc");
  sacut[j] = (TCutG*)sacuts->Get("out56Sc");sacut[j++]->SetName("out56Sc");
  sacut[j] = (TCutG*)sacuts->Get("out55Sc");sacut[j++]->SetName("out55Sc");
  sacut[j] = (TCutG*)sacuts->Get("out54Sc");sacut[j++]->SetName("out54Sc");

  sacut[j] = (TCutG*)sacuts->Get("out58Ca");sacut[j++]->SetName("out58Ca");
  sacut[j] = (TCutG*)sacuts->Get("out57Ca");sacut[j++]->SetName("out57Ca");
  sacut[j] = (TCutG*)sacuts->Get("out56Ca");sacut[j++]->SetName("out56Ca");
  sacut[j] = (TCutG*)sacuts->Get("out55Ca");sacut[j++]->SetName("out55Ca");
  sacut[j] = (TCutG*)sacuts->Get("out54Ca");sacut[j++]->SetName("out54Ca");
  sacut[j] = (TCutG*)sacuts->Get("out53Ca");sacut[j++]->SetName("out53Ca");
  sacut[j] = (TCutG*)sacuts->Get("out52Ca");sacut[j++]->SetName("out52Ca");
  sacut[j] = (TCutG*)sacuts->Get("out51Ca");sacut[j++]->SetName("out51Ca");

  sacut[j] = (TCutG*)sacuts->Get("out55K"); sacut[j++]->SetName("out55K"); 
  sacut[j] = (TCutG*)sacuts->Get("out54K"); sacut[j++]->SetName("out54K"); 
  sacut[j] = (TCutG*)sacuts->Get("out53K"); sacut[j++]->SetName("out53K"); 
  sacut[j] = (TCutG*)sacuts->Get("out52K"); sacut[j++]->SetName("out52K"); 
  sacut[j] = (TCutG*)sacuts->Get("out51K"); sacut[j++]->SetName("out51K"); 
  sacut[j] = (TCutG*)sacuts->Get("out50K"); sacut[j++]->SetName("out50K"); 

  sacut[j] = (TCutG*)sacuts->Get("out53Ar");sacut[j++]->SetName("out53Ar");
  sacut[j] = (TCutG*)sacuts->Get("out52Ar");sacut[j++]->SetName("out52Ar");
  sacut[j] = (TCutG*)sacuts->Get("out51Ar");sacut[j++]->SetName("out51Ar");
  sacut[j] = (TCutG*)sacuts->Get("out50Ar");sacut[j++]->SetName("out50Ar");
  sacut[j] = (TCutG*)sacuts->Get("out49Ar");sacut[j++]->SetName("out49Ar");

  sacut[j] = (TCutG*)sacuts->Get("out50Cl");sacut[j++]->SetName("out50Cl");
  sacut[j] = (TCutG*)sacuts->Get("out49Cl");sacut[j++]->SetName("out49Cl");
  sacut[j] = (TCutG*)sacuts->Get("out48Cl");sacut[j++]->SetName("out48Cl");

  sacut[j] = (TCutG*)sacuts->Get("out48S"); sacut[j++]->SetName("out48S"); 
  sacut[j] = (TCutG*)sacuts->Get("out47S"); sacut[j++]->SetName("out47S"); 
  sacut[j] = (TCutG*)sacuts->Get("out46S"); sacut[j++]->SetName("out46S"); 


  cout << j << " outgoming SA cuts found " << endl;
 
  int nsa = j;
  
  TCanvas * cc = new TCanvas("cc","cc",1200,600);
  cc->Divide(2,1);
  cc->cd(1);
  ch->Draw("BigRIPS_Zic:BigRIPS_AoQ>>hbr(1000,2.3,3,1000,15,25)","","colz");
  for(int i=0;i<nbr;i++)
    brcut[i]->Draw("same");


  cc->cd(2);
  ch->Draw("SAMURAI_Zho:SAMURAI_AoQ2>>hsa(1000,2.3,3.2,1000,15,25)","","colz");
  for(int i=0;i<nsa;i++)
    sacut[i]->Draw("same");

  window();

  */
  /*  
  TH2F* hbr = (TH2F*)gROOT->FindObject("hbr");hbr->SetName("hbr");
  TH2F* hsa = (TH2F*)gROOT->FindObject("hsa");hsa->SetName("hsa");

  TH1F* hg[12];
  j=0;
  ch->Draw("DALI_Edop_vertex_AB>>g0(300,0,6000)","in54Ca&&out53Ca","");    hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g1(300,0,6000)","in54Ca&&out53K","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g2(300,0,6000)","in52Ca&&out51Ca","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g3(300,0,6000)","in52Ca&&out51K","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g4(300,0,6000)","in53Ca&&out52Ca","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g5(300,0,6000)","in55Ca&&out54Ca","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));

  ch->Draw("DALI_Edop_vertex_AB>>g6(300,0,6000)","in56Ca&&out55Ca","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g7(300,0,6000)","in56Sc&&out55Ca","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));

  ch->Draw("DALI_Edop_vertex_AB>>g8(300,0,6000)","in58Ti&&out57Ti","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g9(300,0,6000)","in58Sc&&out57Ca","");	   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));

  ch->Draw("DALI_Edop_vertex_AB>>g10(300,0,6000)","in59Sc&&out58Ca","");   hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++));
  ch->Draw("DALI_Edop_vertex_AB>>g11(300,0,6000)","out58Ca","");           hg[j]->SetName(Form("g%d",j));	hg[j] = (TH1F*)gROOT->FindObject(Form("g%d",j++)); 
  TFile *fout = new TFile("histos_32_90.root","RECREATE");
  hbr->Write();
  hsa->Write();
  for(int i=0;i<12;i++){
    cout << i << endl;
    hg[i]->Write();
  }
  fout->ls();
  fout->Close();
  */
}
