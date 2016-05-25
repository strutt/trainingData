#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TString.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMVA/Factory.h>
#include <TMVA/Reader.h>

#include "OutputConvention.h"
#include "AnitaEventSummary.h"
#include "ProgressBar.h"
#include "TXMLEngine.h"

std::vector<Float_t> getWeightsFromXml(TString xmlWeightFileName);

int main(int argc, char *argv[]){

  // Create ouput file, factory object and open the input file

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  const double minFish = -5;
  const double maxFish = 5;
  const int numFishBins = 1024;
  TH1D* hFishSig = new TH1D("hFishWais", "Wais events", numFishBins, minFish, maxFish);
  TH1D* hFishBkg = new TH1D("hFishQuietMinBias", "Quiet min bias events", numFishBins, minFish, maxFish);

  TString  name = "hImagePeakHilbertPeakWais";
  TString title = "Image peak vs. Hilbert Peak ";
  const Int_t numImagePeakBins = 1024;
  const Int_t numHilbertPeakBins = 1024;
  // const Double_t maxHilbertPeak = 1000;
  const Double_t maxHilbertPeak = 400;

  // <Coefficient Index="0" Value="-2.2592648836783216e+00"/>
  // <Coefficient Index="1" Value="2.6308202243934204e+01"/>
  // <Coefficient Index="2" Value="-4.1458971509719778e-02"/>
  
  TH2D* hImagePeakHilbertPeakWais = new TH2D(name, title,
					     numImagePeakBins, 0, 1,
					     numHilbertPeakBins, 0, maxHilbertPeak);

  name = "hImagePeakHilbertPeakQuietMinBias";
  title = "Image peak vs. Hilbert Peak ";
  TH2D* hImagePeakHilbertPeakMinBias = new TH2D(name, title,
						numImagePeakBins, 0, 1,
						numHilbertPeakBins, 0, maxHilbertPeak);
  
  TChain* waisChain = new TChain("eventSummaryTree");
  waisChain->Add("../waisDistributionsPlots/*.root");
  TChain* minBiasChain = new TChain("eventSummaryTree");
  minBiasChain->Add("quietHPolEventFile.root");

// Set up the TMVA Reader object.
// The names in AddVariable must be same as in the input (weight) file.

  TMVA::Reader* reader = new TMVA::Reader();
  Float_t ip;
  Float_t ph;
  reader->AddVariable("peak[0][1].value", &ip);
  reader->AddVariable("coherent[0][1].peakHilbert", &ph);
  std::string dir    = "weights/";
  std::string prefix = "discriminationTraining";
  TString xmlWeightFileName = dir + prefix + "_Fisher.weights.xml";  
  
  reader->BookMVA("Fisher", xmlWeightFileName);

  std::vector<Float_t> theWeights = getWeightsFromXml(xmlWeightFileName);
  
  TF1* fCutLine = new TF1("fCutLine", "[0]*x + [1]");
  fCutLine->SetParameter(0, theWeights.at(1));
  fCutLine->SetParameter(1, theWeights.at(2));
  // fCutLine->SetParameter(1, -1*-2.2592648836783216e+00);  
  fCutLine->Write();
  
  // get the TTree objects from the input file

  // int nSig = waisChain->GetEntries();
  // int nBkg = minBiasChain->GetEntries();

  std::vector<TTree*> chainVec;
  chainVec.push_back(waisChain);
  chainVec.push_back(minBiasChain);

  Long64_t sigCount = 0;
  Long64_t bkgCount = 0;
  Long64_t numEntries[2]={waisChain->GetEntries(), minBiasChain->GetEntries()};
// Loop over TTrees

  for (UInt_t i=0; i<chainVec.size(); i++){

    // chainVec[i]->Print();

    // Double_t hilbertPeak;
    // Double_t imagePeak;
    AnitaEventSummary* eventSummary = new AnitaEventSummary();
    chainVec[i]->SetBranchAddress("eventSummary", &eventSummary);
    // chainVec[i]->SetBranchAddress("coherent[0][1].peakHilbert", &hilbertPeak);
    // chainVec[i]->SetBranchAddress("peak[0][1].value", &imagePeak);

    numEntries[i] = chainVec[i]->GetEntries();
    
    std::cout << "number of entries = " << numEntries[i] << std::endl;
    ProgressBar p(numEntries[i]);
    

    // Loop over events.  The test statistic is identified by the first 
    // argument used above in BookMVA (below, e.g., "Fisher").

    for (Long64_t j=0; j<numEntries[i]; j++){
      chainVec[i]->GetEntry(j);                // sets evt
      ph = static_cast<Float_t>(eventSummary->coherent[0][1].peakHilbert);
      ip = static_cast<Float_t>(eventSummary->peak[0][1].value);

      double tFisher = reader->EvaluateMVA("Fisher");
      // std::cout << ph << "  " << ip << std::endl;
      // std::cout << "tFisher = " << tFisher << std::endl;

      if ( i == 0 ){ 
        hFishSig->Fill(tFisher);
	hImagePeakHilbertPeakWais->Fill(ip, ph);	
	if ( tFisher > 0) sigCount++;
      }
      else if ( i == 1 ){
        hFishBkg->Fill(tFisher);
	hImagePeakHilbertPeakMinBias->Fill(ip, ph);
	if ( tFisher > 0) bkgCount++;
      }
      p.inc(j, numEntries[i]);
    }
    delete eventSummary;
    eventSummary = NULL;
  }

  std::cout << std::endl;
  std::cout << "Total signal events = " << numEntries[0] << std::endl;
  std::cout << "Total background events = " << numEntries[1] << std::endl;
  std::cout << "Selected signal events (tFisher>0) = " << sigCount << std::endl;
  std::cout << "Selected background events (tFisher>0)  = " << bkgCount << std::endl;

  double sigEfficiency = static_cast<double>(sigCount)/numEntries[0];
  double bkgEfficiency = static_cast<double>(bkgCount)/numEntries[1];

  std::cout << "Signal efficiency  = " << sigEfficiency  << std::endl;
  std::cout << "Background efficiency  = " << bkgEfficiency  << std::endl;

  double sigPrior = 0.5;
  double bkgPrior = 0.5;
  double purity = sigEfficiency*sigPrior/(sigEfficiency*sigPrior+bkgEfficiency*bkgPrior);

  std::cout << "Signal purity = " << purity << std::endl;

  outFile->Write();
  outFile->Close();
  delete reader;
  
  return 0;
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
