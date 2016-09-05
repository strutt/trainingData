// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Reconstruct decimated data set.
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


int main(int argc, char *argv[]){

  if(!(argc==3 || argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [run]" << std::endl;
    std::cerr << "Usage 2: " << argv[0] << " [firstRun] [lastRun]" << std::endl;
    return 1;
  }

  std::cout << argv[0] << "\t" << argv[1];
  if(argc==3){std::cout << "\t" << argv[2];}
  std::cout << std::endl;
  const Int_t firstRun = atoi(argv[1]);
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;

  TChain* headChain = new TChain("headTree");
  // TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);
    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    // gpsChain->Add(fileName);
    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    // MagicDisplay *magicPtr = new MagicDisplay("https://anita:IceRadio@www.hep.ucl.ac.uk/uhen/anita/private/anita3/flight1415/root/",run,WaveCalType::kDefault);

    fileName = TString::Format("https://anita:IceRadio@www.hep.ucl.ac.uk/uhen/anita/private/anita3/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }
  // if(gpsChain->GetEntries()==0){
  //   std::cerr << "Unable to find gps files!" << std::endl;
  //   return 1;
  // }
  if(calEventChain->GetEntries()==0){
    std::cerr << "Unable to find calEvent files!" << std::endl;
    return 1;
  }

  calEventChain->BuildIndex("eventNumber");
  // gpsChain->BuildIndex("eventNumber");

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  // Adu5Pat* pat = NULL;
  // gpsChain->SetBranchAddress("pat", &pat);
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }


  TTree* dataQualityTree = new TTree("dataQualityTree", "dataQualityTree");
  // AnitaEventSummary* dataQuality = new AnitaEventSummary();

  UInt_t eventNumber = 0;
  dataQualityTree->Branch("eventNumber", &eventNumber);

  Double_t peakToPeak[NUM_POL][NUM_SEAVEYS];
  dataQualityTree->Branch("peakToPeak",
			  peakToPeak,
			  TString::Format("peakToPeak[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));
  Short_t numPoints[NUM_POL][NUM_SEAVEYS];
  dataQualityTree->Branch("numPoints",
			  numPoints,
			  TString::Format("numPoints[%d][%d]/S",
					  NUM_POL, NUM_SEAVEYS));

  Double_t power[NUM_POL][NUM_SEAVEYS];
  dataQualityTree->Branch("power",
			  power,
			  TString::Format("power[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));

  Double_t maxVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityTree->Branch("maxVolts",
			  maxVolts,
			  TString::Format("maxVolts[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));

  Double_t minVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityTree->Branch("minVolts",
			  minVolts,
			  TString::Format("minVolts[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));


  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    Int_t isMinBias = RootTools::isMinBiasSampleEvent(header);
    if(isMinBias == 0){ // don't want to do min bias events twice
      calEventChain->GetEntryWithIndex(header->eventNumber);

      UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
      // usefulEvent->setAlfaFilterFlag(false);

      eventNumber = header->eventNumber;

      for(int pol=0; pol<NUM_POL; pol++){
	for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){

	  TGraph* gr = usefulEvent->getGraph(ant, (AnitaPol::AnitaPol_t) pol);


	  numPoints[pol][ant] = gr->GetN();

	  Double_t maxY;
	  Double_t maxX;
	  Double_t minY;
	  Double_t minX;
	  RootTools::getLocalMaxToMin(gr, maxY, maxX, minY, minX);

	  peakToPeak[pol][ant] = maxY - minY;
	  maxVolts[pol][ant] = -9999;
	  minVolts[pol][ant] = 9999;
	  power[pol][ant] = 0;
	  for(int samp=0; samp < numPoints[pol][ant]; samp++){
	    Double_t v = gr->GetY()[samp];
	    power[pol][ant] += v*v*gr->GetX()[samp];
	    if(v > maxVolts[pol][ant]){
	      maxVolts[pol][ant] = v;
	    }
	    if(v < minVolts[pol][ant]){
	      minVolts[pol][ant] = v;
	    }
	  }

	  delete gr;
	}
      }

      dataQualityTree->Fill();
      delete usefulEvent;
    }
    p.inc(entry, maxEntry);
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
