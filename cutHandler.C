char *cutfname =NULL;
bool filecalled = false;
void setFile(char* fname){
  cout << "setting filename to " << fname << endl;
  cutfname = fname;
  filecalled = true;
}
void cutHandler(){
  cout << "================================================================================" << endl;
  cout << "Load / Draw and save 2D cuts in ROOT" << endl;
  cout << "     " << endl;
  cout << "--------------------------------------------------------------------------------" << endl;  
  cout << "- load cuts by calling \"LoadCut(\"SETTINGFILE\")\" " << endl;
  cout << "- set output filename by calling \"setFile(\"YOURFILE.root\")\" "<<endl;
  cout << "  default name is \"cutfile.root\" " <<endl;
  cout << "- start cut by calling \"cut(\"YOURNAME\")\" " << endl;
  cout << "================================================================================" << endl;
 
}
void cut(char* name){  
  TCutG *cut;
  cout << "You can draw a cut by clicking View->Toolbar->Graphical Cut" << endl;
  cout << "Double-click to close cut" << endl;
  cut = (TCutG*)gPad->WaitPrimitive("CUTG");
  cut->SetName(name); 
  cout << "Cut " << name << " with "<< cut->GetN() << " points "<< endl;
  vector<double> x;
  vector<double> y;
  x.resize(cut->GetN());
  y.resize(cut->GetN());
  for(int n=0;n<cut->GetN();n++){
    cut->GetPoint(n,x[n],y[n]);
    cout << x[n] << "\t" << y[n] << endl;
  }
  TFile cutf*;
  if(filecalled==false||cutfname ==NULL)
    cutf = new TFile("cutfile.root","update");
  else
    cutf = new TFile(cutfname,"update");

  cout << "Cut " << name << " will be saved to file " << cutf->GetName() << endl;
  cutf->cd();
  cut->Write(name,TObject::kOverwrite);
  cutf->ls();
  cutf->Close();
}
void loadCut(char* setfilename){
  cout << "setting file is set to " << setfilename << endl;
  TEnv *setting = new TEnv(setfilename);
  int cutnum = setting->GetValue("NumberOfCuts",0);
  cout << "Number of cuts: " << cutnum << endl;
  TCutG *cuts[10];
  char *cutnames[10] = {};
  for(int i=0;i<cutnum;i++)
    cutnames[i] = setting->GetValue(Form("CutName%d",i),"");

  char *cutfilename = setting->GetValue("CutFileName","");
  TFile *fcuts = new TFile(cutfilename);

  for(int i=0;i<cutnum;i++){
    cuts[i] =(TCutG*)fcuts->Get(cutnames[i]);
    if(cuts[i] =(TCutG*)fcuts->Get(cutnames[i]))
      cout << "Load " << cutnames[i] << endl;
  }

  fcuts->Close();
}
