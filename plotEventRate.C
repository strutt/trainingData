const double TrueScaleLat=71;
const double CentralMeridian=0;
const double RadiusOfEarth=6378.1e3; //Metres
const double xOffest=375;
const double yOffset=312.5;
const double scale=271.5/2.19496e+06;
const double xSize=750;
const double ySize=625;

void getRelXYFromLatLongOrig(double latitude, double longitude,
			     double &x, double &y)
{
    //Negative longitude is west
 //    //All latitudes assumed south
    double absLat=TMath::Abs(latitude);
    double r=RadiusOfEarth*TMath::Cos((90.-TrueScaleLat)*TMath::DegToRad())*TMath::Tan((90-absLat)*TMath::DegToRad());

    y=r*TMath::Cos(longitude*TMath::DegToRad());
    x=r*TMath::Sin(longitude*TMath::DegToRad());

    y*=scale;
    y+=yOffset;
    y/=ySize;
    x*=scale;
    x+=xOffest;
    x/=xSize;

    
}

void getRelXYFromLatLong(TH2D* h, double latitude, double longitude,
			 double &x, double &y)
{
    //Negative longitude is west
 //    //All latitudes assumed south
    double absLat=TMath::Abs(latitude);
    double r=RadiusOfEarth*TMath::Cos((90.-TrueScaleLat)*TMath::DegToRad())*TMath::Tan((90-absLat)*TMath::DegToRad());
    y=r*TMath::Cos(longitude*TMath::DegToRad());
    x=r*TMath::Sin(longitude*TMath::DegToRad());

    // std::cout << r << "\t" << x << "\t" << y << std::endl;

    y*=scale;
    y+=yOffset;
    y/=ySize;
    x*=scale;
    x+=xOffest;
    x/=xSize;

    Double_t xMin = h->GetXaxis()->GetBinLowEdge(1);
    Double_t xMax = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    Double_t yMin = h->GetYaxis()->GetBinLowEdge(1);
    Double_t yMax = h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1);

    Double_t histXRange = xMax - xMin;
    Double_t histYRange = yMax - yMin;

    x = xMin + x*histXRange;
    y = yMin + y*histYRange;    
    
    // x = (x+0.5)*histXRange/2;
    // y = (y+0.5)*histYRange/2;

    // x = (-1 + 2*x)*TrueScaleLat;
    // y = (-1 + 2*y)*TrueScaleLat;

}

void plotEventRate(){



  const Int_t firstRun = 127;
  const Int_t lastRun = 439;
  // const Int_t lastRun = 140;

  auto headChain = new TChain("headTree");
  auto gpsChain = new TChain("adu5PatTree");

  for(Int_t run=firstRun; run<=lastRun; run++){

    if(run>=257 && run <= 263) continue;
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
  }

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);

  headChain->BuildIndex("realTime");
  // gpsChain->BuildIndex("realTime");


  TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
  if(!canMap)
    canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
  canMap->Clear();
  // canMap->SetLogz();
  canMap->SetTopMargin(0.03);
  canMap->SetBottomMargin(0.03);
  canMap->SetLeftMargin(0.03);
  canMap->SetRightMargin(0.15);
  TImage *img = TImage::Open("~/Repositories/Install/share/anitaMap/antarcticaIceMapBW.png");
  if (!img) {
    printf("Could not create an image... exit\n");
    return;
  }
  img->SetConstRatio(kFALSE);

  const Int_t nBinsX = 512;
  const Int_t nBinsY = 512;
  auto hTime = new TH2D("hTime", "Event rate", nBinsX, 0, 1, nBinsY, 0, 1);
  auto hDeltaEventNumber = new TH2D("hDeltaEventNumber", "", nBinsX, 0, 1, nBinsY, 0, 1);
  hDeltaEventNumber->GetZaxis()->SetTitle("Event rate (Hz)");
  Double_t x, y;

  // headChain->Show(0);
  // gpsChain->Show(0);
  // gpsChain->Show(1);
  const Long64_t numGps = gpsChain->GetEntries();
  // gpsChain->Show(numGps-1);
  // gpsChain->Show(numGps-2);
  // gpsChain->Show(numGps-3);    
  
  cout << numGps << "\t" << headChain->GetEntries() << endl;
  Long64_t numEntries = headChain->GetEntries();
  // for(Long64_t entry=0; entry<numEntries; entry++){

  UInt_t lastEventNumber = 0;
  UInt_t lastRealTime = 0;
  ProgressBar p(numGps);
  
  for(Long64_t entry=0; entry<numGps; entry++){
  // for(Long64_t entry=0; entry<5; entry++){    
    gpsChain->GetEntry(entry);
    if(headChain->GetEntryWithIndex(pat->realTime) > 0){

      if(lastEventNumber>0){
	Int_t deltaEventNumber = Int_t(header->eventNumber) - Int_t(lastEventNumber);
	// if(deltaEventNumber > 1e5 || entry < 5){
	//   std::cout << endl << entry << "\t" << deltaEventNumber << "\t" << lastEventNumber << "\t" << header->eventNumber << endl;
	// }

	getRelXYFromLatLong(hTime, pat->latitude, pat->longitude, x, y);
	hTime->Fill(x, y);

	getRelXYFromLatLong(hDeltaEventNumber, pat->latitude, pat->longitude, x, y);
	hDeltaEventNumber->Fill(x, y, deltaEventNumber);
      }
    
      lastRealTime = pat->realTime;
      lastEventNumber = header->eventNumber;
    }
    p++;
  }


  hDeltaEventNumber->Divide(hTime);
  
  // hTime->Draw();
  // img->Draw("same");
  // hTime->Draw("colzsame");
  // hTime->GetXaxis()->SetLabelOffset(-99);
  // hTime->GetYaxis()->SetLabelOffset(-99);

  hDeltaEventNumber->Draw();  
  img->Draw("same");
  hDeltaEventNumber->Draw("colzsame");
  hDeltaEventNumber->GetXaxis()->SetLabelOffset(-99);
  hDeltaEventNumber->GetYaxis()->SetLabelOffset(-99);
  
  // gPad->RedrawAxis();
  // gPad->SetGrid();
}
