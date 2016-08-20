#include "OutputConvention.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
// #include "RootTools.h"
#include "TXMLEngine.h"
#include "TF1.h"
#include "TGraph.h"

TH1D* makeNormalizedIntegral(TH1D* h, bool invert=false);
std::vector<Float_t> getWeightsFromXml(TString xmlWeightFileName);





void plotDefineThermalCut(){  

  TFile* f = OutputConvention::getFile("defineThermalCutPlots*.root");
  TTree* t = (TTree*) f->Get("TestTree");
  t->Show(0);
  const double minFischer = -20;
  const double maxFischer = 20;
  const int nBins = 2048;

  TString isSignal = "classID==0";
  TString isBackground = "classID==1";  


  
  TH1D* hSignal = new TH1D("hSignal", "Signal", nBins, minFischer, maxFischer);

  hSignal->Sumw2();
  t->Draw("Fisher>>hSignal", isSignal, "goff");

  TH1D* hBackground = new TH1D("hBackground", "Background", nBins, minFischer, maxFischer);
  hBackground->Sumw2();
  t->Draw("Fisher>>hBackground", isBackground, "goff");


  TH2D* hSignal2 = new TH2D("hSignal2", "Signal", nBins, 0, 1, nBins, 0, 1024);
  t->Draw("hilbertPeak:imagePeak>>hSignal2", isSignal, "goff");
  TH2D* hBackground2 = new TH2D("hBackground2", "Background", nBins, 0, 1, nBins, 0, 1024);
  t->Draw("hilbertPeak:imagePeak>>hBackground2", isBackground, "goff");  


  hSignal->SetLineColor(kRed);
  hSignal2->SetLineColor(kRed);
  
  hBackground->SetLineColor(kBlue);
  hBackground2->SetLineColor(kBlue);

  // yFi(i) = W0 + W1*x1 + W2*x2;
  std::vector<Float_t> weights = getWeightsFromXml("weights/discriminationTraining_Fisher.weights.xml");
  std::cout << weights[0] << "\t" << weights[1] << "\t" << weights[2] << std::endl;  
  
  TF1* fLine = new TF1("fLine", "[0]*x + [1]", 0, 1);
  double grad = -weights[1]/weights[2];
  double intercept = -weights[0]/weights[2];
  fLine->SetParameter(0, grad);
  fLine->SetParameter(1, intercept);

  TF1* fLine2 = new TF1("fLine2", "[0]*x + [1]", 0, 1);
  double grad2 = -weights[1]/(weights[2]);
  double intercept2 = -(2*weights[0])/(weights[2]);
  fLine2->SetParameter(0, grad2);
  fLine2->SetParameter(1, intercept2);  
  

  auto c0 = new TCanvas();
  hSignal2->Draw("box");
  hBackground2->Draw("boxsame");  
  fLine->Draw("lsame");
  fLine2->Draw("lsame");  
  // return;
  
  auto c1 = new TCanvas();
  hSignal->SetLineColor(kRed);
  hSignal->Draw();
  hBackground->Draw("same");

  if(hBackground->GetMaximum() > hSignal->GetMaximum()){
    hSignal->SetMaximum(1.1*hBackground->GetMaximum());
  }
  
  auto l1 = c1->BuildLegend();

  c1->SetLogy(1);

  auto hSignalInt = makeNormalizedIntegral(hSignal, true);
  auto hBackgroundInt = makeNormalizedIntegral(hBackground, true);
  auto c2 = new TCanvas();
  hSignalInt->Draw();
  hBackgroundInt->Draw("same");
  c2->SetLogy(1);
  hSignalInt->SetMinimum(1e-8);
  hBackgroundInt->Fit("expo","", "", -1, 10);

  TGraph* gr = new TGraph();
  for(int i=1; i <= hSignalInt->GetNbinsX(); i++){
    double signalRejection = 1-hSignalInt->GetBinContent(i);
    double backgroundAcceptance = hBackgroundInt->GetBinContent(i);    
    if(signalRejection > 0 && backgroundAcceptance > 0){
      gr->SetPoint(gr->GetN(), signalRejection, backgroundAcceptance);
    }
  }
  auto c3 = new TCanvas();
  gr->Draw("al");
  gr->SetLineColor(kRed);
  // gr->SetMarkerStyle(8);
  // gr->SetMarkerColor(kBlack);
   
  gr->SetName("grFisherScan");
  gr->SetTitle("Scan through Fisher Discriminant; Signal Rejection; Backround Acceptance");
  c3->SetLogx(1);
  c3->SetLogy(1);
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

