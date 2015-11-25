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

    // std::cout << histXRange << "\t" << histYRange << std::endl;

    x = xMin + x*histXRange;
    y = yMin + y*histYRange;    
    
    // x = (x+0.5)*histXRange/2;
    // y = (y+0.5)*histYRange/2;

    // x = (-1 + 2*x)*TrueScaleLat;
    // y = (-1 + 2*y)*TrueScaleLat;

}

void plotLocation(){

  const Int_t firstRun = 127;
  const Int_t lastRun = 439;;

  TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
  if(!canMap)
    canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
  canMap->Clear();
  // canMap->SetLogz();
  canMap->SetTopMargin(0);
  canMap->SetBottomMargin(0);
  canMap->SetLeftMargin(0);
  canMap->SetRightMargin(0);
  TImage *img = TImage::Open("~/Repositories/Install/share/anitaMap/antarcticaIceMapBW.png");
  if (!img) {
    printf("Could not create an image... exit\n");
    return;
  }
  img->SetConstRatio(kFALSE);

  const Int_t nBinsX = 512;
  const Int_t nBinsY = 512;
  TString title = "";
  auto hDeltaEventNumber = new TH2D("hDeltaEventNumber", title, nBinsX, 0, 1, nBinsY, 0, 1);
  hDeltaEventNumber->GetZaxis()->SetTitle("Event rate (Hz)");
  
  std::vector<Double_t> xs;
  std::vector<Double_t> ys;

  ProgressBar p(lastRun - firstRun);
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    
    if(run>=257 && run <= 263) continue;

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    auto f = TFile::Open(fileName);
    auto gpsTree = (TTree*) f->Get("adu5PatTree");
  
    Adu5Pat* pat = NULL;
    gpsTree->SetBranchAddress("pat", &pat);


    const Long64_t numGps = gpsTree->GetEntries();
    for(Long64_t entry=0; entry<numGps; entry++){
      gpsTree->GetEntry(entry);

      if((entry%60)==0){
	Double_t x, y;
	getRelXYFromLatLongOrig(pat->latitude, pat->longitude, x, y);      
	// std::cout << std::endl << x << "\t" << y << std::endl;
	
	xs.push_back(x);
	ys.push_back(y);
	
      }
      
    }

    canMap->Clear();
    TString title = TString::Format("Run %d", run);
    hDeltaEventNumber->SetTitle(title);
    hDeltaEventNumber->Draw();
    hDeltaEventNumber->GetXaxis()->SetTitleOffset(-99);
    hDeltaEventNumber->GetYaxis()->SetTitleOffset(-99);
    hDeltaEventNumber->GetXaxis()->SetLabelOffset(-99);
    hDeltaEventNumber->GetYaxis()->SetLabelOffset(-99);
    img->Draw("same");
    TGraph* gr = new TGraph(xs.size(), &xs[0], &ys[0]);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(2);	    
    gr->Draw("lsame");
    TGraph* gr2 = new TGraph(1, &xs[xs.size()-1], &ys[ys.size()-1]);
    gr2->SetMarkerColor(kYellow);
    gr2->SetMarkerStyle(8);
    gr2->SetMarkerSize(0.5);	    
    gr2->Draw("psame");
    canMap->Update();
    
    RootTools::saveCanvas(canMap, TString::Format("iceMap%d", run));

    // delete gr;
    f->Close();
    p++;
    
  }

}
