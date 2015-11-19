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

  auto f = TFile::Open("powerSpectraPlots_352_2015-11-18_18-34-08.root");

  const Int_t numSamples = 256;
  const Double_t deltaT = 1./2.6;
  Double_t* freqs = FancyFFTs::getFreqArray(numSamples, deltaT);
  const Int_t numFreqs = FancyFFTs::getNumFreqs(numSamples);

  Int_t ant = 32;
  Int_t pol = AnitaPol::kVertical; //AnitaPol::kHorizontal;

  std::vector<TString> baseNames;
  baseNames.push_back("hSouthPs");
  baseNames.push_back("hNorthPs");
  baseNames.push_back("hAvePowSpecs"); 

  std::vector<TString> baseNames2;
  baseNames2.push_back("grSouthAvePs");
  baseNames2.push_back("grNorthAvePs");
  baseNames2.push_back("grAvePs"); 
  
  

  for(UInt_t baseNameInd=0; baseNameInd < baseNames.size(); baseNameInd++){
    TString& baseName = baseNames.at(baseNameInd);
    TH2D* h2 = nullptr;
  
    Int_t count1 = 0;
    TGraph* grChiSquarePerNDF = new TGraph();
    TCanvas* c1 = new TCanvas();
    for(Int_t freqInd=0; freqInd<numFreqs; freqInd++){
      // if(freqs[freqInd] > 0 && freqs[freqInd] < 1.3){
      if(true){      
	TString histName = TString::Format("%s_%d_%d_%d_%d", baseName.Data(), pol, ant, freqInd, numFreqs);
	auto h = (TH1D*) f->Get(histName);
	// h->Rebin(4);
	Int_t nx = h->GetNbinsX();
      
	if(h2==NULL){
	  TString histName2 = TString::Format("%s2_%d_%d_%d_%d", baseName.Data(), pol, ant, freqInd, numFreqs);	  
	  h2 = new TH2D(histName2,
			"Sqrt(PSD) as a function of frequency; Frequency (MHz); Sqrt of PSD (confusing units)",
			numFreqs, freqs[0], numFreqs*(freqs[1]-freqs[0]),
			nx, 1e3*h->GetBinLowEdge(1), 1e3*h->GetBinLowEdge(nx+1));
		      
	}
	for(int binx=1; binx<=nx; binx++){
	  h2->SetBinContent(freqInd+1, binx, h->GetBinContent(binx));
	}
      
	if(h->Integral() > 0){
	  // auto ct = new TCanvas();

	  TString fitName = "fit_" + histName;
	
	  Double_t lowEdge = h->GetBinLowEdge(1);
	  Double_t highEdge = h->GetBinLowEdge(h->GetNbinsX());
	  // cout << lowEdge << "\t" << highEdge << endl;
	  TF1* fit = new TF1(fitName, "[0]*x*exp(-x*x/(2*[1]*[1]))", lowEdge, highEdge);
	  fit->SetParameters(1, h->GetMean());
	  // auto h = (TH1D*) f->Get(TString::Format("hAvePowSpecs_0_23_%d_%d", freqInd, numFreqs));
	  TString opt = count1 == 0 ? "" : "same";
	  h->Draw(opt);
	  h->Fit(fitName, "Q");
	  Double_t chiSq = fit->GetChisquare();
	  chiSq /= fit->GetNDF();
	  if(freqInd > 0 && freqInd < numFreqs-1){
	    grChiSquarePerNDF->SetPoint(count1, freqs[freqInd]*1e3, chiSq);
	    count1++;
	  }
	  // h->Draw(opt);
	  // ct->SetLogy(1);
	}
      }
    }

    new TCanvas();
    // auto gr = (TGraph*) f->Get("grNorthAvePs_0_23_352_352");
    // auto gr = (TGraph*) f->Get("grAvePs_0_23_352_352");
    TString grName = TString::Format("%s_%d_%d_352_352", baseNames2.at(baseNameInd).Data(),
				     pol, ant);    
    auto gr = (TGraph*) f->Get(grName);    
    // std::cout << gr << std::endl;
    gr->Draw();
    grChiSquarePerNDF->SetLineColor(kRed);
    grChiSquarePerNDF->Draw("lsame");

    auto c3 = new TCanvas();
    h2->Draw("colz");
    c3->SetLogz(1);

  }
  
  // auto c2 = new TCanvas();
  // Int_t count2 = 0;
  // for(Int_t freqInd=0; freqInd<numFreqs; freqInd++){
  //   if(freqs[freqInd] > 0.3 && freqs[freqInd] < 1.3){ 
  //     auto h = (TH1D*) f->Get(TString::Format("hSouthPs_0_23_%d_%d", freqInd, numFreqs));
  //     TString opt = count2 == 0 ? "" : "same";
  //     h->Draw(opt);
  //     count2++;
  //   }
  // }
  // c2->SetLogy(1);
}

void drawPowerSpectraOld(){

  auto f = TFile::Open("powerSpectraPlots_352_352_2015-11-13_12-17-03.root");

  auto gr = (TGraph*) f->Get("grAvePs_1_32_60881408");
  auto c1 = new TCanvas();
  gr->Draw();

  const Double_t lowerLimit = 0.2;
  const Double_t upperLimit = 1.2;  
  TF1* f1 = new TF1("fitty", "[0]*x*exp(-(x*x)/(2*[1]*[1]))", lowerLimit, upperLimit);
  c1->SetLogy(1);


  TGraph* grChiSqs = new TGraph();
  auto grFreqsRemoved  = new TGraph();
  const Int_t numPointsToRemove = gr->GetN() - 10;
  for(int i=0; i<numPointsToRemove; i++){
    
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

    if(maxFreqInd<0){
      break;
    }

    grChiSqs->SetPoint(grChiSqs->GetN(), i, f1->GetChisquare()/gr->GetN());
    grFreqsRemoved->SetPoint(grFreqsRemoved->GetN(), i, gr->GetX()[maxFreqInd]);    
    std::cout << gr->GetX()[maxFreqInd] << "\t" << f1->GetChisquare() << std::endl;

    // new TCanvas();
    // TGraph* grTemp = (TGraph*) gr->Clone();
    // grTemp->Draw();
    // TF1* f1Temp = (TF1*) f1->Clone();
    // f1Temp->Draw("same");
    gr->RemovePoint(maxFreqInd);
  }

  auto c2 = new TCanvas();
  grChiSqs->Draw();
  auto c3 = new TCanvas();  
  grFreqsRemoved->SetMarkerStyle(2);
  grFreqsRemoved->Draw("ap");  
  c2->SetLogy(1);
}

