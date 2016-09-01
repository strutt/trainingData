// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
********************************************************************************************************* */

#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "TGraph2D.h"
#include "TRandom3.h"

#include "RawAnitaHeader.h"
#include "UsefulAdu5Pat.h"
#include "UsefulAnitaEvent.h"
#include "CalibratedAnitaEvent.h"
#include "AnitaEventCalibrator.h"

#include "ProgressBar.h"
#include "CrossCorrelator.h"
#include "OutputConvention.h"
#include "AnitaEventSummary.h"
#include "FFTtools.h"
#include "AntarcticaMapPlotter.h"

#include "AnitaClusterer.h"


int main(int argc, char *argv[])
{

    if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }

  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;

  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  // TChain* headChain = new TChain("headTree");
  // TChain* indexedHeadChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }
    if(run == 198 || run == 287){
      continue;
    }

    TString fileName = TString::Format("projectionPlots/projectionDecimatedPlots_%d_*", run);
    eventSummaryChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    // headChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);
    // indexedHeadChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
  }

  if(eventSummaryChain->GetEntries()==0){
    std::cerr << "Unable to find eventSummary files!" << std::endl;
    return 1;
  }

  if(gpsChain->GetEntries()==0){
    std::cerr << "Unable to find gps files!" << std::endl;
    return 1;
  }

  // if(indexedHeadChain->GetEntries()==0){
  //   std::cerr << "Unable to find header files!" << std::endl;
  //   return 1;
  // }
  // indexedHeadChain->BuildIndex("eventNumber");


  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  // UInt_t eventNumberPat = 0;
  // gpsChain->SetBranchAddress("eventNumber", &eventNumberPat);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);


  const int K = 50;
  const int numIterations = 10;
  AnitaClusterer clusterer(K, numIterations);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);

    // Long64_t entry2 = indexedHeadChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
    gpsChain->GetEntry(entry);

    // if(eventSummary->eventNumber!=eventNumberPat){
    //   std::cerr << "fuck" << std::endl;
    // }

    for(Int_t polInd=0; polInd < NUM_POL; polInd++){
      for(int peakInd=0; peakInd < 5; peakInd++){
	Double_t sourceLat = eventSummary->peak[polInd][peakInd].latitude;
	Double_t sourceLon = eventSummary->peak[polInd][peakInd].longitude;
	Double_t sourceAlt = eventSummary->peak[polInd][peakInd].altitude;
	if(sourceLat > -999 && sourceLon > -999){
	  // std::cout << (eventSummary->peak[polInd][peakInd].distanceToSource < 1e6 )<< std::endl;
	  if(eventSummary->peak[polInd][peakInd].theta < 0 && eventSummary->peak[polInd][peakInd].distanceToSource < 1e6){

	    clusterer.addPoint(pat, sourceLat,sourceLon,sourceAlt, eventSummary->run, eventSummary->eventNumber, 0.2, 0.5);
	    // Int_t n = clusterer.addPoint(sourceLat,sourceLon,sourceAlt);
	    // std::cout << n << std::endl;
	  }
	}
      }
    }
    p++;
  }

  // clusterer.kMeansCluster(1);
  clusterer.initializeCOMNAP();

  outFile->cd();
  for(int clusterInd=0; clusterInd < K; clusterInd++){
    TGraph* gr = clusterer.makeClusterSummaryTGraph(clusterInd);

    if(gr){
      gr->Write();
      delete gr;
    }
    else{
      std::cerr << "??????? " << clusterInd << std::endl;
    }
  }

  // no need to write?
  TTree* clusterTree = clusterer.makeClusterSummaryTree(outFile);

  std::cout << "Made clusterTree with " << clusterTree->GetEntries() << " entries." << std::endl;






  outFile->Write();
  outFile->Close();




  return 0;
}
