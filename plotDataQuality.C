#include "TH1D.h"
#include "OutputConvention.h"
#include "RootTools.h"

void plotDataQuality(){


  TFile* fWais = OutputConvention::getFile("plotWaisDataQualityPlots_*.root");
  TFile* fDec = OutputConvention::getFile("plotDecimatedDataQualityPlots_*.root");

  std::cout << fWais << "\t" <<  fDec << std::endl;

  TH1D* hMaxWais = (TH1D*) fWais->Get("hMaxMaxVolts");
  TH1D* hMaxDec = (TH1D*) fDec->Get("hMaxMaxVolts");


  hMaxWais->Rebin(8);
  hMaxDec->Rebin(8);

  hMaxWais->Scale(1./hMaxWais->Integral());
  hMaxDec->Scale(1./hMaxDec->Integral());


  auto l1 = new TLegend(0.8, 0.8, 1, 1);
  l1->AddEntry(hMaxWais, "WAIS pulses", "l");
  l1->AddEntry(hMaxDec, "All decimated", "l");

  auto c1 = new TCanvas();
  hMaxDec->Draw("");
  hMaxWais->Draw("same");
  l1->Draw();

  hMaxWais->SetLineColor(kRed);
  hMaxDec->SetLineColor(kBlue);
  c1->SetLogy(1);



  TH1D* hMinWais = (TH1D*) fWais->Get("hMinMinVolts");
  TH1D* hMinDec = (TH1D*) fDec->Get("hMinMinVolts");

  hMinWais->Scale(1./hMinWais->Integral());
  hMinDec->Scale(1./hMinDec->Integral());

  hMinWais->Rebin(8);
  hMinDec->Rebin(8);


  auto c2 = new TCanvas();

  hMinDec->Draw("");
  hMinWais->Draw("same");


  hMinWais->SetLineColor(kRed);
  hMinDec->SetLineColor(kBlue);
  c2->SetLogy(1);


  auto l2 = new TLegend(0.8, 0.8, 1, 1);
  l2->AddEntry(hMaxWais, "WAIS pulses", "l");
  l2->AddEntry(hMaxDec, "All decimated", "l");
  l2->Draw();


  TH1D* hDiffWais = (TH1D*) fWais->Get("hMaxMinVolts");
  TH1D* hDiffDec = (TH1D*) fDec->Get("hMaxMinVolts");

  hDiffWais->Scale(1./hDiffWais->Integral());
  hDiffDec->Scale(1./hDiffDec->Integral());

  auto c3 = new TCanvas();

  hDiffDec->Draw("");
  hDiffWais->Draw("same");

  hDiffWais->Rebin(8);
  hDiffDec->Rebin(8);


  hDiffWais->SetLineColor(kRed);
  hDiffDec->SetLineColor(kBlue);


  auto l3 = new TLegend(0.8, 0.8, 1, 1);
  l3->AddEntry(hMaxWais, "WAIS pulses", "l");
  l3->AddEntry(hMaxDec, "All decimated", "l");
  l3->Draw();

  c3->SetLogy(1);




  auto c4 = new TCanvas();

  TFile* fWais2 = OutputConvention::getFile("cutFlow4/plotReconstructedWaisPlots_*.root");
  TFile* fDec2 = OutputConvention::getFile("fullTimeWithCuts/plotReconstructedDecimatedPlots_*.root");

  TH2D* h2Wais = (TH2D*) fWais2->Get("hPeakRatio");
  TH2D* h2Dec = (TH2D*) fDec2->Get("hPeakRatio");

  cout << fWais2 << "\t" << fDec2 << endl;
  cout << h2Wais << "\t" << h2Dec << endl;

  auto c5 = new TCanvas();
  h2Wais->Draw("colz");

  // return;

  TH1D* h2WaisY = h2Wais->ProjectionY("h2WaisY");
  TH1D* h2DecY = h2Dec->ProjectionY("h2DecY");

  h2WaisY->SetLineColor(kRed);
  h2DecY->SetLineColor(kBlue);

  h2DecY->Rebin(2);

  h2DecY->Rebin(2);
  h2WaisY->Rebin(2);

  auto l4 = new TLegend(0.8, 0.8, 1, 1);
  l4->AddEntry(h2WaisY, "WAIS Pulses", "l");
  l4->AddEntry(h2DecY, "Decimated (All Time)", "l");

  h2WaisY->Scale(1./h2WaisY->Integral());
  h2DecY->Scale(1./h2DecY->Integral());

  h2DecY->Draw("");
  h2WaisY->Draw("same");

  h2DecY->SetMaximum(0.3);
  c5->SetLogy(1);

  h2DecY->SetTitle("Ratio of map peaks P_{2} / P_{1}; Peak ratio; Fraction of events per bin");
  l4->Draw();

  const double decInt = h2DecY->Integral(0, h2DecY->GetNbinsX()+1);
  const double waisInt = h2WaisY->Integral(0, h2WaisY->GetNbinsX()+1);
  std::cout << "before dec: " << h2DecY->Integral(0, h2DecY->GetNbinsX()+1) << std::endl;
  std::cout << "before wais: " << h2WaisY->Integral(0, h2WaisY->GetNbinsX()+1) << std::endl;

  const double cutVal = 0.5;
  double n1 = 0;
  for(int i=1; i <= h2DecY->GetNbinsX(); i++){
    if(h2DecY->GetXaxis()->GetBinLowEdge(i) < cutVal){
      n1 += h2DecY->GetBinContent(i);
    }
  }
  double n2 = 0;
  for(int i=1; i <= h2WaisY->GetNbinsX(); i++){
    if(h2WaisY->GetXaxis()->GetBinLowEdge(i) < cutVal){
      n2 += h2WaisY->GetBinContent(i);
    }
  }
  std::cerr << "Dec efficiency " << 100*n1/decInt << std::endl;
  std::cerr << "Wais efficiency " << 100*n2/waisInt << std::endl;

}
