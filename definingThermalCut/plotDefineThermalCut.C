#include "OutputConvention.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
// #include "RootTools.h"
#include "TXMLEngine.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"
#include "TExec.h"
#include "TROOT.h"
#include "TPaletteAxis.h"

TH1D* makeNormalizedIntegral(TH1D* h, bool invert=false);
std::vector<Float_t> getWeightsFromXml(TString xmlWeightFileName);

void Pal1();
void Pal2();



void plotDefineThermalCut(){


  // const double desiredBackgroundAcceptance = 2.5e-8; //1e-7;
  const double desiredBackgroundAcceptance = 1e-7;

  TFile* f = OutputConvention::getFile("defineThermalCutPlots*.root");
  TTree* t = (TTree*) f->Get("TestTree");
  t->Show(0);
  const double minFischer = -3;
  const double maxFischer = 6;
  const int nBinsCoarse = 512;
  const int nBinsFine = 1024;

  TString isSignal = "classID==0";
  TString isBackground = "classID==1";



  TH1D* hSignal = new TH1D("hSignal", "Signal and Background Acceptance as a function of w_{0}; w_{0} (no units); Fraction of population above w_{0}", nBinsCoarse, minFischer, maxFischer);

  hSignal->Sumw2();
  t->Draw("Fisher>>hSignal", isSignal, "goff");

  TH1D* hBackground = new TH1D("hBackground", "Background", nBinsFine, minFischer, maxFischer);
  hBackground->Sumw2();
  t->Draw("Fisher>>hBackground", isBackground, "goff");


  TH2D* hSignal2 = new TH2D("hSignal2", "Signal", nBinsCoarse, 0, 1, nBinsCoarse, 0, 1024);
  t->Draw("hilbertPeak:imagePeak>>hSignal2", isSignal, "goff");
  TH2D* hBackground2 = new TH2D("hBackground2", "Background", nBinsFine, 0, 1, nBinsFine, 0, 1024);
  t->Draw("hilbertPeak:imagePeak>>hBackground2", isBackground, "goff");


  hSignal->SetLineColor(kRed);
  hSignal2->SetLineColor(kRed);
  hSignal->SetMarkerStyle(0);
  hSignal2->SetMarkerStyle(0);

  hBackground->SetLineColor(kBlue);
  hBackground2->SetLineColor(kBlue);

  // yFi(i) = W0 + W1*x1 + W2*x2;
  std::vector<Float_t> weights = getWeightsFromXml("weights/discriminationTraining_Fisher.weights.xml");
  std::cout << weights[0] << "\t" << weights[1] << "\t" << weights[2] << std::endl;

  auto hSignalInt = makeNormalizedIntegral(hSignal, true);
  auto hBackgroundInt = makeNormalizedIntegral(hBackground, true);
  auto c2 = new TCanvas();
  hSignalInt->Draw();
  hBackgroundInt->Draw("same");
  c2->SetLogy(1);
  hSignalInt->SetMinimum(1e-9);
  TF1* fBackgroundFit = new TF1("fBackgroundFit", "exp([constant]+[slope]*x)", -2, 1);
  fBackgroundFit->SetNpx(1000);
  fBackgroundFit->SetLineStyle(3);
  fBackgroundFit->SetLineColor(kBlue);
  const double startFit = -1.6;
  const double endFit = -1.25;
  hBackgroundInt->Fit(fBackgroundFit,"N", "", startFit, endFit);
  fBackgroundFit->Draw("lsame");

  std::cerr << "Fitting range of function..." << std::endl;
  std::cerr << fBackgroundFit->Eval(startFit) << "\t" << fBackgroundFit->Eval(endFit) << std::endl;



  // so now we grab the params and invert this mother
  Double_t p0 = fBackgroundFit->GetParameter("constant");
  Double_t p1 = fBackgroundFit->GetParameter("slope");

  // y = exp(a + bx)
  // ln(y)  = a + bx
  // (ln(y) - a)/b = x
  // want x for y = 1e-7ish, maybe
  double cutValFisher = (TMath::Log(desiredBackgroundAcceptance) - p0)/p1;
  std::cerr << "cutValFisher = " << cutValFisher
	    << "\tdesiredBackgroundAcceptance = " << desiredBackgroundAcceptance << std::endl;
  TGraph* grCutLine = new TGraph();
  grCutLine->SetPoint(grCutLine->GetN(), cutValFisher, 1e-10);
  grCutLine->SetPoint(grCutLine->GetN(), cutValFisher, 10);
  grCutLine->SetLineColor(kMagenta);
  grCutLine->SetLineStyle(2);
  grCutLine->Draw("lsame");
  TGraph* grCutLine2 = new TGraph();
  grCutLine2->SetPoint(grCutLine2->GetN(), minFischer, desiredBackgroundAcceptance);
  grCutLine2->SetPoint(grCutLine2->GetN(), maxFischer, desiredBackgroundAcceptance);
  grCutLine2->SetLineColor(kMagenta);
  grCutLine2->SetLineStyle(2);
  grCutLine2->Draw("lsame");

  TLegend* l2 = new TLegend(0.8, 0.8, 1 , 1);
  l2->AddEntry(hSignal, "WAIS pulses", "l");
  l2->AddEntry(hBackground, "Min Bias (non-sun pointing events)", "l");
  l2->AddEntry(fBackgroundFit, "Exponential fit to background", "l");
  l2->AddEntry(grCutLine, TString::Format("w_{0} at %.2e background acceptance", desiredBackgroundAcceptance), "l");
  l2->Draw();


  // return;


  hSignal2->SetMarkerSize(1);
  hSignal2->Scale(1./hSignal2->Integral());
  hSignal2->SetMarkerStyle(1);
  hSignal2->GetXaxis()->SetRangeUser(0, 0.5);
  hSignal2->GetYaxis()->SetRangeUser(0, 200);

  hBackground2->Scale(1./hBackground2->Integral());
  hBackground2->GetXaxis()->SetRangeUser(0, 0.5);
  hBackground2->GetYaxis()->SetRangeUser(0, 200);

  const int numCont = 40;
  gStyle->SetNumberContours(numCont);
  TCanvas *c3  = new TCanvas("fancyCan", "fancyCan", 1500, 800);

  TPad *pad1 = new TPad("pad1", "",0, 0, 0.5, 1, 0);
  TPad *pad2 = new TPad("pad2", "",0.5, 0, 1, 1, 0);
  pad1->Draw();
  pad2->Draw();
  pad2->cd();

  pad2->SetLogz(1);
  hBackground2->Draw("colz");

  hBackground2->SetTitle(";;;Min Bias events per bin");
  hSignal2->SetTitle("Rotated Cross-Correlation; Image Peak (no units); Hilbert Peak (mV); WAIS events per bin");
  TExec *ex1 = new TExec("ex1","Pal1();");
  ex1->Draw();

  hBackground2->Draw("cont z list same");

  // hBackground2->Draw("cont z list same");
  // return;

  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hBackground2->GetListOfFunctions()->FindObject("palette");
  TObjArray *contours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");

  std::vector<std::vector<TGraph*> > grs;
  Int_t ncontours = contours->GetSize();

  for(int i=0; i < ncontours; i++){
    TList *list         = (TList*)contours->At(i);
    grs.push_back(std::vector<TGraph*>(0));


    TGraph* gr = (TGraph*)list->First();
    for(int j=0; j < list->GetSize(); j++){

      TGraph* gr = (TGraph*) list->At(j);
      TGraph* gr2 = (TGraph*) gr->Clone();


      Double_t val = hBackground2->GetContourLevel(i);
      Int_t col = palette->GetValueColor(val);
      gr2->SetLineColor(col);


      grs.back().push_back(gr2);
    }
  }
  gPad->Update();


  pad1->cd();
  hSignal2->Draw("colz");
  TExec *ex2 = new TExec("ex2","Pal2();");
  ex2->Draw();
  hSignal2->Draw("colzsame");

  UInt_t theN = grs.size();
  std::vector<Color_t> theCols(theN);

  for(Int_t i=0; i < theN; i++){
    // std::cout << i << "\t" << theN - i - 1 << "\t" << grs.at(theN - i - 1).size() <<  std::endl;
    theCols.at(i) = grs.at(theN - i - 1).at(0)->GetLineColor();
  }
  Int_t drawEvery = 1;
  for(UInt_t i=0; i < grs.size(); i++){
    if((i % drawEvery) == 0){
      for(UInt_t j=0; j < grs.at(i).size(); j++){
	// if(j==0){std::cout << i << std::endl;}

  	grs.at(i).at(j)->SetLineColor(theCols.at(i));
  	grs.at(i).at(j)->Draw("lsame");
      }
    }
  }

  gPad->Update();





  TF1* fLine = new TF1("fLine", "[0]*x + [1]", 0, 1);
  fLine->SetLineColor(kMagenta);
  fLine->SetLineStyle(2);
  fLine->SetLineWidth(2);
  fLine->SetNpx(10000);
  double grad = -weights[1]/weights[2];
  double intercept = (cutValFisher-weights[0])/weights[2];
  // double intercept = -(weights[0]-cutValFisher)/weights[2];
  // double intercept = -cutValFisher/weights[2];
  fLine->SetParameter(0, grad);
  fLine->SetParameter(1, intercept);
  gStyle->SetOptFit(1);
  TLegend* lFancyCan = new TLegend(0.8, 0.1, 1, 0.3);
  lFancyCan->AddEntry(fLine, TString::Format("Fisher discriminant extrapolated to accept %.2e Min Bias events", desiredBackgroundAcceptance), "l");

  fLine->Draw("lsame");
  lFancyCan->Draw();


  auto c1 = new TCanvas();
  hSignal->SetLineColor(kRed);
  hSignal->Draw();
  hBackground->Draw("same");

  if(hBackground->GetMaximum() > hSignal->GetMaximum()){
    hSignal->SetMaximum(1.1*hBackground->GetMaximum());
  }

  hSignal->SetTitle("Fisher Discriminant for Signal and Background; F (no units); Events per bin");
  auto l1 = new TLegend(0.8, 0.8, 1, 1);
  //c1->BuildLegend();
  l1->AddEntry(hSignal, "WAIS Pulses", "l");
  l1->AddEntry(hBackground, "Min Bias (after cuts)", "l");
  l1->Draw();

  c1->SetLogy(1);




  return;

  // this is shit
  TGraph* gr = new TGraph();
  for(int i=1; i <= hSignalInt->GetNbinsX(); i++){
    double signalRejection = 1-hSignalInt->GetBinContent(i);
    double backgroundAcceptance = hBackgroundInt->GetBinContent(i);
    if(signalRejection > 0 && backgroundAcceptance > 0){
      gr->SetPoint(gr->GetN(), signalRejection, backgroundAcceptance);
    }
  }
  auto c4 = new TCanvas();
  gr->Draw("al");
  gr->SetLineColor(kRed);
  // gr->SetMarkerStyle(8);
  // gr->SetMarkerColor(kBlack);

  gr->SetName("grFisherScan");
  gr->SetTitle("Scan through Fisher Discriminant; Signal Rejection; Backround Acceptance");
  c4->SetLogx(1);
  c4->SetLogy(1);
  gr->SetMinimum(1e-5);
  gr->GetXaxis()->SetRangeUser(1e-5, 1);
}


















TH1D* makeNormalizedIntegral(TH1D* h, bool invert){
  auto hInt = (TH1D*) h->Clone(TString::Format("%sInt", h->GetName()));
  hInt->Scale(1./hInt->Integral());
  Double_t* signalInt = hInt->GetIntegral();
  if(invert==false){
    for(int i=0; i <= hInt->GetNbinsX(); i++){
      hInt->SetBinContent(i, signalInt[i]);
    }
  }
  else{
    for(int i=10; i <= hInt->GetNbinsX(); i++){
      hInt->SetBinContent(i, 1-signalInt[i]);
      // std::cerr << i << signalInt[i] << std::endl;
    }
  }
  // delete [] signalInt;
  return hInt;
}












std::vector<Float_t> getWeightsFromXml(TString xmlWeightFileName){
  // First create engine
  TXMLEngine* xml = new TXMLEngine;

  // Now try to parse xml file
  // Only file with restricted xml syntax are supported
  XMLDocPointer_t xmldoc = xml->ParseFile(xmlWeightFileName.Data());
  if (xmldoc==0) {
    delete xml;
    std::cerr << "balls" << std::endl;
    std::vector<Float_t> a;
    return a;
  }

  // take access to main node
  XMLNodePointer_t node = xml->DocGetRootElement(xmldoc);
  node = xml->GetChild(node); // descend one level

  // loop until we get to weights
  while(strcmp(xml->GetNodeName(node), "Weights")!=0){
    // std::cout << xml->GetNodeName(node) << std::endl;
    node = xml->GetNext(node);
    // std::cout << node << std::endl;
  }

  XMLAttrPointer_t nCoefNode = xml->GetFirstAttr(node);
  // std::cout << "at node: " << xml->GetNodeName(node) << std::endl;
  const int nCoef = atoi(xml->GetAttrValue(nCoefNode));
  // std::cout << "I found " << nCoef << " weights" << std::endl;

  std::vector<Float_t> theWeights(nCoef, 0);

  // loop through children
  node = xml->GetChild(node);
  while(node!=0){
    // std::cout << "at node: " << xml->GetNodeName(node) << std::endl;
    XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    int currentIndex = -1; // to guarentee crash?
    while(attr!=0){
      // std::cout << "\t attr: " << xml->GetAttrName(attr) << " is " << xml->GetAttrValue(attr) << std::endl;
      if(strcmp(xml->GetAttrName(attr), "Index")==0){
	currentIndex = atoi(xml->GetAttrValue(attr));
      }
      else if(strcmp(xml->GetAttrName(attr), "Value")==0){
	theWeights.at(currentIndex) = atof(xml->GetAttrValue(attr));
      }
      attr = xml->GetNextAttr(attr);
    }
    node = xml->GetNext(node);
  }

  // display recursively all nodes and subnodes
  // DisplayNode(xml, node, 1);

  // Release memory before exit
  xml->FreeDoc(xmldoc);
  delete xml;
  return theWeights;
}




void Pal1()
{

  //gStyle->SetPalette(1);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor color;
  color.InitializeColors();
  color.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  // gStyle->SetPalette(50,colors
   // static Int_t  colors[50];
   // static Bool_t initialized = kFALSE;

   // Double_t Red[3]    = { 1.00, 0.00, 0.00};
   // Double_t Green[3]  = { 0.00, 1.00, 0.00};
   // Double_t Blue[3]   = { 1.00, 0.00, 1.00};
   // Double_t Length[3] = { 0.00, 0.50, 1.00 };

   // if(!initialized){
   //    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
   //    for (int i=0; i<50; i++) colors[i] = FI+i;
   //    initialized = kTRUE;
   //    return;
   // }
   // gStyle->SetPalette(50,colors);

}

void Pal2()
{
   // static Int_t  colors[50];
   // static Bool_t initialized = kFALSE;

   // Double_t Red[3]    = { 1.00, 0.50, 0.00};
   // Double_t Green[3]  = { 0.50, 0.00, 1.00};
   // Double_t Blue[3]   = { 1.00, 0.00, 0.50};
   // Double_t Length[3] = { 0.00, 0.50, 1.00 };

   // if(!initialized){
   //    Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
   //    for (int i=0; i<50; i++) colors[i] = FI+i;
   //    initialized = kTRUE;
   //    return;
   // }
   // gStyle->SetPalette(50,colors);
   gStyle->SetPalette(kDarkBodyRadiator);
}
