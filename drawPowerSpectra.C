void drawOneEvent(TFile* f){
  auto canName = "canPow";
  TCanvas* c1 = new TCanvas(canName, canName, 1200, 900);
  c1->Divide(3, 2);

  for(int pol=0; pol < AnitaPol::kNotAPol; pol++){
    for(int ring=0; ring<3; ring++){
      c1->cd(pol*3 + ring + 1);
      for(int phi=0; phi<NUM_PHI; phi++){
	int ant = phi + NUM_PHI*ring;
      
	auto name = TString::Format("grAvePs_%d_%d_60871452", pol, ant);
	auto gr = (TGraph*) f->Get(name);
	auto opt = phi==0 ? "" : "lsame";
	
	gr->SetLineColor(gStyle->GetColorPalette(phi*255./NUM_PHI));
	gr->Draw(opt);
	gr->GetXaxis()->SetNoExponent(1);
	gr->SetMaximum(70);
      }
    }
  }
}



void drawPowerSpectra(){

  auto f = TFile::Open("powerSpectraPlots_352_352_2015-11-13_12-17-03.root");

  auto gr = (TGraph*) f->Get("grAvePs_1_32_60881408");
  auto c1 = new TCanvas();
  gr->Draw();

  const Double_t lowerLimit = 0.2;
  const Double_t upperLimit = 1.2;  
  TF1* f1 = new TF1("fitty", "[0]*x*exp(-(x*x)/(2*[1]*[1]))", lowerLimit, upperLimit);
  c1->SetLogy(1);


  TGraph* grChiSqs = new TGraph();
  for(int i=0; i<20; i++){
    
    f1->SetParameter(0, 1e6);
    f1->SetParameter(1, 0.3);
    // f1->Draw("same");
    gr->Fit(f1, "", "", lowerLimit, upperLimit);
    
    Double_t maxSqDiff = -DBL_MAX;
    Int_t maxFreqInd = -1;
    for(Int_t freqInd=0; freqInd < gr->GetN(); freqInd++){
      Double_t freq = gr->GetX()[freqInd];
      if(freq >= lowerLimit && freq < upperLimit){
	Double_t grY = gr->GetY()[freqInd];
	Double_t fitY = f1->Eval(gr->GetX()[freqInd]);

	Double_t sqDiff = pow(grY - fitY, 2);
	if(sqDiff > maxSqDiff){
	  maxSqDiff = sqDiff;
	  maxFreqInd = freqInd;
	}
      }
    }

    grChiSqs->SetPoint(grChiSqs->GetN(), i, f1->GetChisquare());
    std::cout << gr->GetX()[maxFreqInd] << "\t" << f1->GetChisquare() << std::endl;

    new TCanvas();
    TGraph* grTemp = (TGraph*) gr->Clone();
    grTemp->Draw();
    TF1* f1Temp = (TF1*) f1->Clone();
    f1Temp->Draw("same");
    
    gr->RemovePoint(maxFreqInd);
  }

  auto c2 = new TCanvas();
  grChiSqs->Draw();
  c2->SetLogy(1);
}

