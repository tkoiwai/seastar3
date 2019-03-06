void 1Dfit2Dhisto(TTree *tr,TString title,int dim){
  TGraph *gr = new TGraph(tr->GetSelectedRows(),tr->GetV2(),tr->GetV1());
  gr->SetTitle( title);
  gr->Draw("ap");
  gr->Fit(Form("pol%d",dim));
  //Form("pol%d",dim)->Draw("same");
}
