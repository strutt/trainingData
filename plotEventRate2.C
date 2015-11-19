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

    cout << scale << endl;
    
    
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

    std::cout << histXRange << "\t" << xMin << "\t" << histYRange << "\t" << yMin << "\t" << endl;

    x = xMin + x*histXRange;
    y = yMin + y*histYRange;

    cout << scale << endl;
    
    // x = (x+0.5)*histXRange/2;
    // y = (y+0.5)*histYRange/2;

    // x = (-1 + 2*x)*TrueScaleLat;
    // y = (-1 + 2*y)*TrueScaleLat;

}

void plotEventRate2(){




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

  const Int_t nBinsX = 128;
  const Int_t nBinsY = 128;
  auto hDeltaEventNumber = new TH2D("hDeltaEventNumber", "", nBinsX, 0, 1, nBinsY, 0, 1);
  hDeltaEventNumber->GetZaxis()->SetTitle("Event rate (Hz)");

  Double_t x,y;
  getRelXYFromLatLong(hDeltaEventNumber, -90, 90, x, y);
  hDeltaEventNumber->Fill(x, y);
  getRelXYFromLatLong(hDeltaEventNumber, -80, 90, x, y);
  hDeltaEventNumber->Fill(x, y);
  getRelXYFromLatLong(hDeltaEventNumber, -70, 180, x, y);
  hDeltaEventNumber->Fill(x, y);

  hDeltaEventNumber->Draw();  
  img->Draw("same");
  hDeltaEventNumber->Draw("colzsame");
  hDeltaEventNumber->GetXaxis()->SetLabelOffset(-99);
  hDeltaEventNumber->GetYaxis()->SetLabelOffset(-99);
  
  // gPad->RedrawAxis();
  // gPad->SetGrid();
}
