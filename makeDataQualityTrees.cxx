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

  CrossCorrelator* cc = new CrossCorrelator();


  TChain* headChain = new TChain("headTree");
  // TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);
    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    // gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
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

  cc->kDoSimpleSatelliteFiltering = 1;
  if(cc->kDoSimpleSatelliteFiltering > 0){
    TNamed* comments = new TNamed("comments", "Applied simple, static notch at 260#pm26 MHz and 370#pm20");
    comments->Write();
    delete comments;
  }
    
  TTree* dataQualityTree = new TTree("dataQualityTree", "dataQualityTree");
  // AnitaEventSummary* dataQuality = new AnitaEventSummary();

  UInt_t eventNumber = 0;
  dataQualityTree->Branch("eventNumber", &eventNumber);

  Double_t maxVolts[NUM_POL][NUM_SEAVEYS];  
  dataQualityTree->Branch("maxVolts",
			  maxVolts,
			  TString::Format("maxVolts[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));

  Double_t maxNegSecondDeriv[NUM_POL][NUM_SEAVEYS];  
  dataQualityTree->Branch("maxNegSecondDeriv",
			  maxVolts,
			  TString::Format("maxNegSecondDeriv[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));

  Double_t sumVoltsSquared[NUM_POL][NUM_SEAVEYS];  
  dataQualityTree->Branch("sumVoltsSquared",
			  sumVoltsSquared,
			  TString::Format("sumVoltsSquared[%d][%d]/D",
					  NUM_POL, NUM_SEAVEYS));

  Int_t numPoints[NUM_POL][NUM_SEAVEYS];  
  dataQualityTree->Branch("numPoints",
			  numPoints,
			  TString::Format("numPoints[%d][%d]/I",
					  NUM_POL, NUM_SEAVEYS));  
  
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);


  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    Int_t isMinBias = header->getTriggerBitADU5() || header->getTriggerBitG12() || header->getTriggerBitSoftExt();

    if(isMinBias > 0){

      calEventChain->GetEntryWithIndex(header->eventNumber);

      UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

      eventNumber = header->eventNumber;

      const int numFreqs = FancyFFTs::getNumFreqs(cc->numSamples);

      for(int pol=0; pol<NUM_POL; pol++){
	cc->getNormalizedInterpolatedTGraphs(usefulEvent, (AnitaPol::AnitaPol_t) pol);

	for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	  FancyFFTs::doFFT(cc->numSamples, cc->grsResampled[pol][ant]->GetY(), cc->ffts[pol][ant]);

	  cc->simple260MHzSatelliteNotch((AnitaPol::AnitaPol_t) pol, ant);
	  cc->simple370MHzSatelliteNotch((AnitaPol::AnitaPol_t) pol, ant);

	  Double_t vSq = 0;
	  for(int freqInd=0; freqInd < numFreqs; freqInd++){
	    Double_t theNorm = std::norm(cc->ffts[pol][ant][freqInd]);
	    theNorm = freqInd == numFreqs - 1 ? 0.5*theNorm : theNorm;
	    vSq += theNorm;
	  }
	  vSq /= (cc->numSamples*cc->numSamples*0.5);
	  vSq *= cc->interpRMS[pol][ant]*cc->interpRMS[pol][ant];
	  vSq *= cc->numSamples;
	  // std::cout << vSq << std::endl;

	  TGraph* gr = usefulEvent->getGraph(ant, (AnitaPol::AnitaPol_t) pol);

	  sumVoltsSquared[pol][ant] = vSq;
	  maxNegSecondDeriv[pol][ant] = -9999;
	  numPoints[pol][ant] = gr->GetN();
	  maxVolts[pol][ant] = -9999;

	  for(int samp = 0; samp < numPoints[pol][ant]; samp++){

	    Double_t V = gr->GetY()[samp];

	    if(samp < numPoints[pol][ant] - 2){
	      Double_t V2 = gr->GetY()[samp+1];
	      Double_t V3 = gr->GetY()[samp+2];
	      double thisNegSecondDeriv = 2*V2 - V3 - V;
	      if(thisNegSecondDeriv > maxNegSecondDeriv[pol][ant]){
		maxNegSecondDeriv[pol][ant] = thisNegSecondDeriv;
	      }
	    }

	    if(V > maxVolts[pol][ant]){
	      maxVolts[pol][ant] = V;
	    }
	  }

	  if(maxVolts[pol][ant] < 0){
	    std::cerr << eventNumber << "\t" << maxVolts[pol][ant] << std::endl;
	    for(int samp=0; samp < numPoints[pol][ant]; samp++){
	      std::cerr << gr->GetY()[samp] << ", ";
	    }
	    std::cerr << std::endl;
	  }

	  delete gr;

	}
      }

      dataQualityTree->Fill();
      delete usefulEvent;      
    }
    p++;
  }

  outFile->Write();
  outFile->Close();  
  
  return 0;
}
