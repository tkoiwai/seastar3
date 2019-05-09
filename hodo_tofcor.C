//#include "loadsmri.C"
#include <stdio.h>
#include <typeinfo.h>
void hodo_tofcor(Int_t runnum, int hodo_id){

  gROOT->ProcessLine( "gErrorIgnoreLevel = 2001;");
  
  //loadsmri(runnum);
  TString FileNameS = Form("/home/koiwai/analysis/rootfiles/ana/smri/ana_smri%04d.root",runnum);
  TFile *f = TFile::Open(FileNameS);

    
  switch(hodo_id){
    //cout << "in switch" << endl;
  case 1:
    cout << "case 1" << endl;
    break;
  case 2:
    cout << "case 2" << endl;
    //anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==2");
    //TF1 *func = new TF1("func","gaus");
    //h->Fit("func");
    //Double_t tofcor02 = 59.1884 - func->GetParameter(1);
    //if(abs(tofcor02)>5.) break;
    //ofstream fout("/home/koiwai/analysis/db/hodo02_tofcor.dat",ios::app);
    //fout << runnum << ":" << tofcor02 << endl;
    break;
  case 3:
    cout << "case 3" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==3");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor03 = 57.2172 - func->GetParameter(1);
    if(abs(tofcor03)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo03_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor03 << endl;
    break;
  case 4:
    cout << "case 4" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==04");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor04 = 57.0676 - func->GetParameter(1);
    if(abs(tofcor04)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo04_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor04 << endl;
    break;
  case 5:
    cout << "case 5" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==5");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor05 = 56.9419 - func->GetParameter(1);
    if(abs(tofcor05)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo05_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor05 << endl;
    break;
  case 6:
    cout << "case 6" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==6");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor06 = 56.8162 - func->GetParameter(1);
    if(abs(tofcor06)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo06_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor06 << endl;
    break;
  case 7:
    cout << "case 7" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==7&&br50ar");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor07 = 58.0926 - func->GetParameter(1);
    if(abs(tofcor07)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo07_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor07 << endl;
    break;
  case 8:
    cout << "case 8" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,57,62)","BG_flag_beam&&hodo_id==8&&br50ar");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor08 = 58.0190 - func->GetParameter(1);
    if(abs(tofcor08)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo08_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor08 << endl;
    break;
  case 9:
    cout << "case 9" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==9&&br51k");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor09 = 56.4931 - func->GetParameter(1);
    if(abs(tofcor09)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo09_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor09 << endl;
    break;
  case 10:
    cout << "case 10" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==10&&br51k");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor10 = 56.4213 - func->GetParameter(1);
    if(abs(tofcor10)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo10_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor10 << endl;
    break;
  case 11:
    cout << "case 11" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==11&&br51k");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor11 = 56.3674 - func->GetParameter(1);
    if(abs(tofcor11)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo11_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor11 << endl;
    break;
  case 12:
    cout << "case 12" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==12&&br51k");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor12 = 56.3255 - func->GetParameter(1);
    if(abs(tofcor12)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo12_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor12 << endl;
    break;
  case 13:
    cout << "case 13" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,65)","BG_flag_beam&&hodo_id==13&&br56sc");
    if(h) cout << "histogram created" << endl;
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor13 = 57.3169 - func->GetParameter(1);
    if(abs(tofcor13)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo13_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor13 << endl;
    break;
  case 14:
    cout << "case 14" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==14&&br58ti");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor14 = 57.6192 - func->GetParameter(1);
    if(abs(tofcor14)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo14_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor14 << endl;
    break;
  case 15:
    cout << "case 15" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==15&&br61v");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor15 = 60.9090 - func->GetParameter(1);
    if(abs(tofcor15)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo15_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor15 << endl;
    break;
  case 16:
    cout << "case 16" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==16&&br62v");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor16 = 62.3208 - func->GetParameter(1);
    if(abs(tofcor16)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo16_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor16 << endl;
    break;
  case 17:
    cout << "case 17" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==17&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor17 = 57.2255 - func->GetParameter(1);
    if(abs(tofcor17)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo17_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor17 << endl;
    break;
  case 18:
    cout << "case 18" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==18&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor18 = 57.2255 - func->GetParameter(1);
    if(abs(tofcor18)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo18_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor18 << endl;
    break;
  case 19:
    cout << "case 19" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==19&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor19 = 57.2377 - func->GetParameter(1);
    if(abs(tofcor19)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo19_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor19 << endl;
    break;
  case 20:
    cout << "case 20" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==20&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor20 = 57.2621 - func->GetParameter(1);
    if(abs(tofcor20)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo20_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor20 << endl;
    break;
  case 21:
    cout << "case 21" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==21&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor21 = 57.2987 - func->GetParameter(1);
    if(abs(tofcor21)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo21_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor21 << endl;
    break;
  case 22:
    cout << "case 22" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==22&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor22 = 57.3413 - func->GetParameter(1);
    if(abs(tofcor22)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo22_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor22 << endl;
    break;
  case 23:
    cout << "case 23" << endl;
    anatrS->Draw("t_minoshodo_notcor>>h(200,55,67)","BG_flag_beam&&hodo_id==23&&br56sc");
    TF1 *func = new TF1("func","gaus");
    h->Fit("func");
    Double_t tofcor23 = 57.3840 - func->GetParameter(1);
    if(abs(tofcor23)>10.) break;
    ofstream fout("/home/koiwai/analysis/db/hodo23_tofcor.dat",ios::app);
    fout << runnum << ":" << tofcor23 << endl;
    break;
  case 24:
    cout << "case 24" << endl;
    //anatrS->Draw("t_minoshodo_notcor>>h(100,57,77)","BG_flag_beam&&hodo_id==24");
    //TF1 *func = new TF1("func","gaus");
    //h->Fit("func");
    //Double_t tofcor24 = 66.3343 - func->GetParameter(1);
    //if(abs(tofcor24)>5.) break;
    //ofstream fout("/home/koiwai/analysis/db/hodo24_tofcor.dat",ios::app);
    //fout << runnum << ":" << tofcor24 << endl;
    break;
  default:
    //cout << id << endl;
    cout << "default called" << endl;
  }
    
  f->Close();
  
    
}
