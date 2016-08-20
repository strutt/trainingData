#include "OutputConvention.h"
// #include "RootTools.h"

TH1D* getThePeakYVals(TH2D* h2){

  const int nY = h2->GetNbinsY();
  const double minY = h2->GetYaxis()->GetBinLowEdge(0);
  const double maxY = h2->GetYaxis()->GetBinLowEdge(nY+1);

  TString name = TString::Format("%s_peakYVals", h2->GetName());
  
  TH1D* h1 = new TH1D(name, name, nY, minY, maxY);

  for(int bx=1; bx <= h2->GetNbinsX(); bx++){
    double maxBin = -1;
    for(int by=1; by <= nY; by++){
      if(h2->GetBinContent(bx, by) > 0){
	maxBin = by;
      }
    }
    double theMax = h2->GetYaxis()->GetBinLowEdge(maxBin+1);
    // std::cerr << theMax << std::endl;
    h1->Fill(theMax);
  }  
  return h1;
}


TH2D* zeroBinsWithTimesGreaterThanCutVal(TH2D* h2, double cutVal){

  const int nY = h2->GetNbinsY();
  const double minY = h2->GetYaxis()->GetBinLowEdge(0);
  const double maxY = h2->GetYaxis()->GetBinLowEdge(nY+1);

  TString name = TString::Format("%s_clone", h2->GetName());
  auto h2_2 = (TH2D*) h2->Clone(name);
  
  for(int bx=1; bx <= h2->GetNbinsX(); bx++){
    double maxBin = -1;
    for(int by=1; by <= nY; by++){
      if(h2->GetBinContent(bx, by) > 0){
	maxBin = by;
      }
    }

    if(h2->GetYaxis()->GetBinLowEdge(maxBin+1) > cutVal){
      for(int by=1; by <= nY; by++){
	h2_2->SetBinContent(bx, by, 0);
      }
    }

    // std::cerr << theMax << std::endl;
  }  
  return h2_2;
}


double findCutVal(TH1D* h){

  const int nX = h->GetNbinsX();
  // const double minX = h->GetXaxis()->GetBinLowEdge(0);
  // const double maxX = h->GetXaxis()->GetBinLowEdge(nX+1);

  double firstX=0;
  double peakX=0;
  double peakVal=0;    
  for(int bx=1; bx <= nX; bx++){

    double val = h->GetBinContent(bx);
    double X = h->GetBinLowEdge(bx);    

    if(firstX==0 && val > 0){
      firstX = X;
    }
    if(val > peakVal){
      peakVal = val;
      peakX  = X;
    }    
  }

  return 1.5*(peakX + (peakX - firstX));
}


void plotDiscriminatingVariablesVsTime(){
  
  const int numFiles = 1;
  TString fileNames[numFiles] = {"plotReconstructedDecimatedPlots_352_2016-08-19_16-57-43.root"};

  // TString fileNames[numFiles] = {"plotReconstructedDecimatedPlots_130_434_2016-08-19_13-38-17.root"};
  // TString fileNames[numFiles] = {"plotReconstructedDecimatedPlots_130_434_2016-08-19_13-55-04.root"};
  // TString fileNames[numFiles] = {"plotReconstructedDecimatedPlots_130_434_2016-08-17_16-23-51.root"};

  // // this one for first pass through two notches
  // TString fileNames[numFiles] = {"plotReconstructedMinBiasPlots_130*08-16_11-03-11.root"};

  // TString fileNames[numFiles] = {"plotReconstructedMinBiasPlots_130_434_2016-08-17_15-05-33.root"}; //plotReconstructedMinBiasPlots_130_434_2016-08-17_14-54-48.root"};

  // const int numFiles = 2;
  // TString fileNames[numFiles] = {"plotReconstructedMinBiasPlots_352_*.root",
  // 				 "plotReconstructedMinBiasPlots_130_439_2016-08-16_11-03*.root"};

  // const int numFiles = 3;  
  // TString fileNames[numFiles] = {"plotReconstructedMinBiasPlots_130*.root",
  // 				 "plotReconstructedDecimatedPlots_130*.root",
  // 				 "plotReconstructedWaisPlots_130*.root"};

  const int theLogz = 0;
  for(int fileInd=0; fileInd < numFiles; fileInd++){
    TFile* f = OutputConvention::getFile(fileNames[fileInd]);

    std::cerr << f->GetName() << std::endl;
    // gDirectory->ls();
    
    // TH2D* h1 = (TH2D*) f->Get("hPeakDirWrtNorth_HPol");
    auto c1 = new TCanvas();
    TH2D* h1 = (TH2D*) f->Get("hImagePeakTime_HPol");
    if(h1 == NULL){
      h1 = (TH2D*) f->Get("hImagePeakTimeHPol");
    }
    h1->GetXaxis()->SetTimeDisplay(1);    
    h1->Draw("colz");
    c1->SetLogz(theLogz);

    auto c2 = new TCanvas();    
    TH2D* h2 = (TH2D*) f->Get("hImagePeakTime_VPol");
    if(h2 == NULL){
      h2 = (TH2D*) f->Get("hImagePeakTimeVPol");
    }
    h2->Draw("colz");
    h2->GetXaxis()->SetTimeDisplay(1);
    c2->SetLogz(theLogz);


    const double maxHilbert = 512;
    auto c3 = new TCanvas();
    TH2D* h3 = (TH2D*) f->Get("hHilbertPeakTime_HPol");
    if(h3 == NULL){
      h3 = (TH2D*) f->Get("hHilbertPeakTimeHPol");
    }
    h3->GetXaxis()->SetTimeDisplay(1);
    h3->GetYaxis()->SetRangeUser(0, maxHilbert);
    h3->Draw("colz");
    c3->SetLogz(theLogz);
    
    auto c4 = new TCanvas();    
    TH2D* h4 = (TH2D*) f->Get("hHilbertPeakTime_VPol");
    if(h4 == NULL){
      h4 = (TH2D*) f->Get("hHilbertPeakTimeVPol");
    }
    h4->GetXaxis()->SetTimeDisplay(1);
    h4->GetYaxis()->SetRangeUser(0, maxHilbert);
    h4->Draw("colz");
    c4->SetLogz(theLogz);



    auto c2y = new TCanvas();
    c2y->SetLogy(1);    
    TH1D* h2y = h2->ProjectionY(); //getThePeakYVals(h4);
    // auto c2y = new TCanvas();
    h2y->Draw();
    h2y->SetLineColor(kBlue);        
    TH1D* h1y = h1->ProjectionY(); //getThePeakYVals(h4);
    h2y->SetTitle("Min Bias Image Peak Distributions");
    h2y->GetYaxis()->SetTitle("Events per bin");    
    h1y->Draw("same");
    h1y->SetLineColor(kRed);
    auto l2y = new TLegend(0.8, 0.8, 1, 1);
    l2y->AddEntry(h2y, "VPol", "l");
    l2y->AddEntry(h1y, "HPol", "l");
    l2y->Draw();

    // double cutVal1 = findCutVal(h1y);
    // double cutVal2 = findCutVal(h2y);
    
    auto c4y = new TCanvas();
    c4y->SetLogy(1);
    TH1D* h4y = h4->ProjectionY(); //getThePeakYVals(h4);
    // auto c4y = new TCanvas();
    h4y->Draw();
    h4y->SetLineColor(kBlue);        
    TH1D* h3y = h3->ProjectionY(); //getThePeakYVals(h4);
    h4y->GetYaxis()->SetTitle("Events per bin");
    h4y->SetTitle("Min Bias Hilbert Peak Distributions");

    h3y->Draw("same");
    h3y->SetLineColor(kRed);
    auto l4y = new TLegend(0.8, 0.8, 1, 1);
    l4y->AddEntry(h4y, "VPol", "l");
    l4y->AddEntry(h3y, "HPol", "l");
    l4y->Draw();

    // double cutVal3 = findCutVal(h3y);
    // double cutVal4 = findCutVal(h4y);

    // std::cout << cutVal1 << "\t" << cutVal2 << std::endl;
    // std::cout << cutVal3 << "\t" << cutVal4 << std::endl;


    
    // TCanvas* c1z = new TCanvas();
    // auto h1z = zeroBinsWithTimesGreaterThanCutVal(h1, cutVal1);
    // h1z->Draw("colz");
    // TCanvas* c2z = new TCanvas();
    // auto h2z = zeroBinsWithTimesGreaterThanCutVal(h2, cutVal2);
    // h2z->Draw("colz");    

    // TCanvas* c3z = new TCanvas();
    // auto h3z = zeroBinsWithTimesGreaterThanCutVal(h3, cutVal3);
    // h3z->Draw("colz");
    // TCanvas* c4z = new TCanvas();
    // auto h4z = zeroBinsWithTimesGreaterThanCutVal(h4, cutVal4);
    // h4z->Draw("colz");

    // auto c2zy = new TCanvas();
    // TH1D* h2zy = h2z->ProjectionY(); //getThePeakZYVals(h4);
    // // auto c2zy = new TCanvas();
    // h2zy->Draw();
    // h2zy->SetLineColor(kBlue);        
    // TH1D* h1zy = h1z->ProjectionY(); //getThePeakZYVals(h4);

    // h1zy->Draw("same");
    // h1zy->SetLineColor(kRed);
    
    // auto c4zy = new TCanvas();
    // TH1D* h4zy = h4z->ProjectionY(); //getThePeakZYVals(h4);
    // // auto c4zy = new TCanvas();
    // h4zy->Draw();
    // h4zy->SetLineColor(kBlue);        
    // TH1D* h3zy = h3z->ProjectionY(); //getThePeakZYVals(h4);

    // h3zy->Draw("same");
    // h3zy->SetLineColor(kRed);


    
    auto c5 = new TCanvas();
    auto h5 = (TH2D*) f->Get("hPeakHeadingHPol");
    h5->SetTitle("Peak Heading vs. Time HPol; Time; Peak Elevation (Degrees); Events per bin");
    h5->Draw("colz");
    c5->SetLogz(theLogz);
    h5->GetXaxis()->SetTimeDisplay(1);



    auto c6 = new TCanvas();
    auto h6 = (TH2D*) f->Get("hPeakHeadingVPol");
    h6->SetTitle("Peak Heading vs. Time VPol; Time; Peak Heading (Degrees); Events per bin");    
    h6->Draw("colz");
    c6->SetLogz(theLogz);
    h6->GetXaxis()->SetTimeDisplay(1);


    auto c7 = new TCanvas();
    auto h7 = (TH2D*) f->Get("hPeakElevationHPol");
    h7->Draw("colz");
    h7->SetTitle("Peak Elevation vs. Time HPol; Time; Peak Elevation (Degrees); Events per bin");
    c7->SetLogz(theLogz);
    h7->GetXaxis()->SetTimeDisplay(1);
    


    auto c8 = new TCanvas();
    auto h8 = (TH2D*) f->Get("hPeakElevationVPol");
    h8->Draw("colz");
    h8->SetTitle("Peak Elevation vs. Time VPol; Time; Peak Elevation (Degrees); Events per bin");
    c8->SetLogz(theLogz);
    h8->GetXaxis()->SetTimeDisplay(1);



    auto c5x = new TCanvas();
    c5x->SetLogy(1);
    auto h5x = h5->ProjectionX();    

    // auto c6x = new TCanvas();
    auto h6x = h6->ProjectionX();    

    h5x->Rebin(4);
    h6x->Rebin(4);    
    
    h5x->SetLineColor(kRed);
    h6x->SetLineColor(kBlue);
    h5x->Draw();
    h6x->Draw("same");
    auto l5x = new TLegend(0.8, 0.8, 1, 1);
    l5x->AddEntry(h6x, "VPol", "l");
    l5x->AddEntry(h5x, "HPol", "l");
    l5x->Draw();
    
  }
  

  return;
}


