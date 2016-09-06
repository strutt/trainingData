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
#include "TVirtualIndex.h"
#include "TTreeIndex.h"
#include "TChainIndex.h"

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
#include "AnalysisCuts.h"


int main(int argc, char *argv[])
{

  if(!(argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;
    return 1;
  }

  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = firstRun;

  std::cout << firstRun << "\t" << lastRun << std::endl;

  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  TChain* dataQualityChain = new TChain("dataQualityTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    TString fileName = TString::Format("filter260-370-400-762-5peaks/reconstructDecimatedPlots_%d_*.root", run);
    eventSummaryChain->Add(fileName);

    fileName = TString::Format("finalDataQuality/makeDecimatedDataQualityTreesPlots_%d_*.root", run);
    // fileName = TString::Format("filter260-370-400-762/makeDecimatedDataQualityTreesPlots_%d_*.root", run);
    dataQualityChain->Add(fileName);

  }
  // return 1;
  if(eventSummaryChain->GetEntries()==0){
    std::cerr << "Unable to find eventSummary files!" << std::endl;
    return 1;
  }
  if(dataQualityChain->GetEntries()==0){
    std::cerr << "Unable to find dataQualityFiles files!" << std::endl;
    return 1;
  }

  Double_t peakToPeak[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("peakToPeak", peakToPeak);
  UInt_t eventNumberDQ;
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumberDQ);
  Double_t maxVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("maxVolts", maxVolts);
  Short_t numPoints[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("numPoints", numPoints);

  Double_t minVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("minVolts", minVolts);


  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  TTree* dataQualityTree = new TTree("dataQualityTree", "dataQualityTree");
  UInt_t eventNumber = 0;
  dataQualityTree->Branch("eventNumber", &eventNumber);

  Double_t maxPeakToPeakRatio[NUM_POL];
  dataQualityTree->Branch("maxPeakToPeakRatio",
			  maxPeakToPeakRatio,
			  TString::Format("maxPeakToPeakRatio[%d]/D",
					  NUM_POL));

  Short_t numPoints2[NUM_POL][NUM_SEAVEYS];
  dataQualityTree->Branch("numPoints",
			  numPoints2,
			  TString::Format("numPoints[%d][%d]/S",
					  NUM_POL, NUM_SEAVEYS));
  Double_t theMaxVolts[NUM_POL];
  dataQualityTree->Branch("theMaxVolts",
			  theMaxVolts,
			  TString::Format("theMaxVolts[%d]/D",
					  NUM_POL));
  Double_t theMinVolts[NUM_POL];
  dataQualityTree->Branch("theMinVolts",
			  theMinVolts,
			  TString::Format("theMinVolts[%d]/D",
					  NUM_POL));

  std::cerr << "building index" << std::endl;
  dataQualityChain->BuildIndex("eventNumber");
  std::cerr << "done" << std::endl;

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    Int_t dqEntry = dataQualityChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
    if(dqEntry < 0){
      eventNumber = eventSummary->eventNumber;
      for(int pol=0; pol<NUM_POL; pol++){
	maxPeakToPeakRatio[pol] = -9999;
	for(int ant=0; ant < NUM_SEAVEYS; ant++){
	  numPoints2[pol][ant] = -99;
	}
	theMaxVolts[pol] = -9999;
	theMinVolts[pol] = -9999;
      }
    }
    else{
      dataQualityChain->GetEntry(dqEntry);
      if(eventSummary->eventNumber!=eventNumberDQ){
	std::cerr << "?????????????\t" << eventSummary->eventNumber << "\t" << eventNumberDQ << std::endl;
      }
      eventNumber = eventSummary->eventNumber;
      for(int pol=0; pol<NUM_POL; pol++){
	theMaxVolts[pol] = 0;
	theMinVolts[pol] = 0;
	for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){

	  numPoints2[pol][ant] = numPoints[pol][ant];

	  if(maxVolts[pol][ant] > theMaxVolts[pol]){
	    theMaxVolts[pol] = maxVolts[pol][ant] ;
	  }
	  if(minVolts[pol][ant]  < theMinVolts[pol]){
	    theMinVolts[pol] = minVolts[pol][ant];
	  }
	}
	maxPeakToPeakRatio[pol] = -1;
	for(int phi=0; phi < NUM_PHI; phi++){
	  if(pol==AnitaPol::kVertical && phi==7){
	    continue;
	  }
	  Double_t ratio = peakToPeak[pol][phi+2*NUM_PHI]/peakToPeak[pol][phi];
	  if(ratio > maxPeakToPeakRatio[pol]){
	    maxPeakToPeakRatio[pol] = ratio;
	  }
	}
      }
    }
    dataQualityTree->Fill();
    p.inc(entry, maxEntry);
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
