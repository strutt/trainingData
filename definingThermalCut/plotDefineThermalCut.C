void plotDefineThermalCut(const char* fileName){
  TFile* f = TFile::Open(fileName);


  TH1D* hFishSig = dynamic_cast<TH1D*>(f->Get("hFishWais"));
  TH1D* hFishBkg = dynamic_cast<TH1D*>(f->Get("hFishQuietMinBias"));

  auto c1 = new TCanvas();
  hFishBkg->SetLineColor(kRed);
  hFishBkg->Draw();


  hFishSig->SetLineColor(kBlue);
  hFishSig->Draw("same");
  auto l1 = c1->BuildLegend();

  c1->SetLogy(1);


  auto c2 = new TCanvas();
  
  TH2D* h2DWais = dynamic_cast<TH2D*>(f->Get("hImagePeakHilbertPeakWais"));
  TH2D* h2DQuietMinBias = dynamic_cast<TH2D*>(f->Get("hImagePeakHilbertPeakQuietMinBias"));

  const int rebinXFactor = 8;
  const int rebinYFactor = 8;  
  h2DWais->RebinX(rebinXFactor);
  h2DQuietMinBias->RebinX(rebinXFactor);
  h2DWais->RebinY(rebinYFactor);
  h2DQuietMinBias->RebinY(rebinYFactor);  
  
  h2DWais->Draw("box");
  h2DWais->SetFillColor(0);
  h2DWais->SetLineColor(kBlue);    

  
  h2DQuietMinBias->Draw("box same");
  h2DQuietMinBias->SetFillColor(0);
  h2DQuietMinBias->SetLineColor(kRed);

  auto l2 = c2->BuildLegend();

  TF1* fCutLine = dynamic_cast<TF1*>(f->Get("fCutLine"));
  new TCanvas();
  fCutLine->Draw();
}
