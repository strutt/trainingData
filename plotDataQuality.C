TGraph* getLimitGraph(TH2D* hMaxAbsSecondDeriv_0){

  const double deltaY = hMaxAbsSecondDeriv_0->GetYaxis()->GetBinLowEdge(2) - hMaxAbsSecondDeriv_0->GetYaxis()->GetBinLowEdge(1);

  const double integral = 3.52904e+06;
  const double frac = 1 - 1e-4;
  std::cout << "the frac = " << frac << std::endl;
  TGraph* gr0 = new TGraph();
  for(int binx=1; binx<=hMaxAbsSecondDeriv_0->GetNbinsX(); binx++){
    double theIntegral = 0;
    for(int biny=1; biny<=hMaxAbsSecondDeriv_0->GetNbinsY(); biny++){
      theIntegral += hMaxAbsSecondDeriv_0->GetBinContent(binx, biny);
      if(theIntegral >= integral*frac){
	const double theCutVal = hMaxAbsSecondDeriv_0->GetYaxis()->GetBinLowEdge(biny) + deltaY;
	// cerr << binx << "\t" << theCutVal << endl;
	gr0->SetPoint(gr0->GetN(), binx-0.5, theCutVal);
	break;
      }
    }
    // std::cerr << binx << "\t" << theIntegral << std::endl;
  }

  gr0->SetLineWidth(2);
  gr0->SetLineColor(kMagenta);

  std::cout << "limits[" << gr0->GetN() << "] = {";
  for(int ant=0; ant < gr0->GetN(); ant++){
    std::cout << gr0->GetY()[ant];
    if(ant < gr0->GetN()-1){
      std::cout << ", ";
    }
  }
  std::cout << "};" << std::endl;
  return gr0;

}



void plotDataQuality(){

  // TFile* f = TFile::Open("plotDataQualityPlots_130_434_2016-04-25_21-28-44.root");
  // TFile* f = TFile::Open("plotDataQualityPlots_130_434_2016-04-25_22-09-32.root");
  //  TFile* f = TFile::Open("plotDataQualityPlots_130_434_2016-04-26_18-16-09.root");
  TFile* f = TFile::Open("plotDataQualityPlots_130_434_2016-05-06_10-56-20.root");
  
  TH2D* hMaxAbsSecondDeriv_0 = (TH2D*) f->Get("hMaxAbsSecondDeriv_0");
  TH2D* hMaxAbsSecondDeriv_1 = (TH2D*) f->Get("hMaxAbsSecondDeriv_1");

  auto c1 = new TCanvas();
  hMaxAbsSecondDeriv_0->Draw("colz");
  c1->SetLogz(1);
  TGraph* gr0 = getLimitGraph(hMaxAbsSecondDeriv_0);
  gr0->Draw("lsame");
  
  auto c2 = new TCanvas();
  hMaxAbsSecondDeriv_1->Draw("colz");  
  c2->SetLogz(1);
  TGraph* gr1 = getLimitGraph(hMaxAbsSecondDeriv_1);
  gr1->Draw("lsame");  
  // std::cout << "Double_t xMaxLimits[NUM_POL][NUM_SEAVEYS] = {";  
}

  
