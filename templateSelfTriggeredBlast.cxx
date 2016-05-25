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
  const Int_t lastRun = firstRun; //argc==3 ? atoi(argv[2]) : firstRun;

  TChain* headChain = new TChain("headTree");
  TChain* calEventChain = new TChain("eventTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }
  if(calEventChain->GetEntries()==0){
    std::cerr << "Unable to find calEvent files!" << std::endl;
    return 1;
  }
  
  calEventChain->BuildIndex("eventNumber");

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);

  // const UInt_t templateEventNumber = 61072528;
  // calEventChain->GetEntryWithIndex(templateEventNumber);
  // UsefulAnitaEvent* templateEvent = new UsefulAnitaEvent(calEvent);

  // CrossCorrelator* ccTemplate = new CrossCorrelator();
  // ccTemplate->kDoSimpleSatelliteFiltering = 0;
  // for(int polInd=0; polInd < NUM_POL; polInd++){
  //   AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
  //   ccTemplate->getNormalizedInterpolatedTGraphs(templateEvent, pol);
  //   // ccTemplate->doFFTs(pol);
  // }
  
  // delete templateEvent;  

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  
    
  TTree* selfTriggeredBlastTree = new TTree("selfTriggeredBlastTree", "selfTriggeredBlastTree");
  // AnitaEventSummary* eventSummary = new AnitaEventSummary();

  Double_t corrFactors[NUM_POL][NUM_SEAVEYS];
  selfTriggeredBlastTree->Branch("corrFactors",
				 corrFactors,
				 TString::Format("corrFactors[%d][%d]/D",
						 NUM_POL, NUM_SEAVEYS));
  UInt_t eventNumber = 0;
  selfTriggeredBlastTree->Branch("eventNumber", &eventNumber);

  Int_t isCandidate = 0;
  selfTriggeredBlastTree->Branch("isCandidate", &isCandidate);

  Int_t maxAnt = 0;
  selfTriggeredBlastTree->Branch("maxAnt", &maxAnt);
  Int_t maxPol = 0;
  selfTriggeredBlastTree->Branch("maxPol", &maxPol);  
  Double_t maxRms = 0;
  selfTriggeredBlastTree->Branch("maxRms", &maxRms);

  Int_t topRingAnt = 0;
  selfTriggeredBlastTree->Branch("topRingAnt", &topRingAnt);
  Double_t topRingRms = 0;
  selfTriggeredBlastTree->Branch("topRingRms", &topRingRms);


  Double_t xMinAddPeds;
  selfTriggeredBlastTree->Branch("xMinAddPeds", &xMinAddPeds);
  Double_t xMaxAddPeds;  
  selfTriggeredBlastTree->Branch("xMaxAddPeds", &xMaxAddPeds);

  Double_t chanXMax[NUM_POL][NUM_SEAVEYS];
  selfTriggeredBlastTree->Branch("chanXMax",
				 chanXMax,
				 TString::Format("chanXMax[%d][%d]/D",
						 NUM_POL, NUM_SEAVEYS));
  Double_t chanXMin[NUM_POL][NUM_SEAVEYS];
  selfTriggeredBlastTree->Branch("chanXMin",
				 chanXMin,
				 TString::Format("chanXMin[%d][%d]/D",
						 NUM_POL, NUM_SEAVEYS));  
  

  Double_t maxSumBottomMidSix[NUM_POL] = {0};
  selfTriggeredBlastTree->Branch("maxSumBottomMidSix",
				 maxSumBottomMidSix,
				 TString::Format("maxSumBottomMidSix[%d]/D",
						 NUM_POL));
  Double_t maxSumTopSix[NUM_POL] = {0};
  selfTriggeredBlastTree->Branch("maxSumTopSix",
				 maxSumTopSix,
				 TString::Format("maxSumTopSix[%d]/D",
						 NUM_POL));

  Int_t maxPhi[NUM_POL] = {0};
  selfTriggeredBlastTree->Branch("maxPhi",
				 maxPhi,
				 TString::Format("maxPhi[%d]/I",
						 NUM_POL));    
    
  CrossCorrelator* cc = new CrossCorrelator();
  cc->kDoSimpleSatelliteFiltering = 0;

  // Double_t* theCrossCorr = new Double_t[cc->numSamples];
  // const double threshVal = 0.5;

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    Int_t isMinBias = RootTools::isMinBiasSampleEvent(header);
    
    if(!(isMinBias > 0)){
      p++;
      continue;
    }

    eventNumber = header->eventNumber;
    isCandidate = 0;
    calEventChain->GetEntryWithIndex(header->eventNumber);

    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
    UsefulAnitaEvent* pedEvent = new UsefulAnitaEvent(calEvent, WaveCalType::kAddPeds);

    xMaxAddPeds = -DBL_MAX;
    xMinAddPeds = DBL_MAX;
    
    for(int polInd=0; polInd < NUM_POL; polInd++){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	TGraph* gr = pedEvent->getGraph(ant, (AnitaPol::AnitaPol_t) polInd);

	Double_t xMax = -DBL_MAX;
	Double_t xMin = DBL_MAX;
	for(int samp=0; samp < gr->GetN(); samp++){
	  Double_t y = gr->GetY()[samp];
	  if(y > xMax){
	    xMax = y;
	  }
	  if(y < xMin){
	    xMin = y;
	  }	  
	}

	chanXMax[polInd][ant] = xMax;
	chanXMin[polInd][ant] = xMin;	

	if(xMax > xMaxAddPeds){
	  xMaxAddPeds = xMax;
	}
	if(xMin < xMinAddPeds){
	  xMinAddPeds = xMax;
	}
	
	delete gr;
      }
    }
    

    for(int polInd=0; polInd < NUM_POL; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      cc->getNormalizedInterpolatedTGraphs(usefulEvent, pol);
      // cc->doFFTs(pol);
    }

    maxRms = 0;
    for(int polInd=0; polInd < NUM_POL; polInd++){
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	if(cc->interpRMS[polInd][ant] > maxRms){
	  maxRms = cc->interpRMS[polInd][ant];
	  maxAnt = ant;
	  maxPol = polInd;
	}
      }
    }
    topRingAnt = maxAnt % NUM_PHI;
    topRingRms = cc->interpRMS[maxPol][topRingAnt];

    for(int polInd=0; polInd < NUM_POL; polInd++){
      maxSumBottomMidSix[polInd] = 0;
      maxSumTopSix[polInd] = 0;      
      for(int phi=0; phi<NUM_PHI; phi++){
	double sumBottomMidSix = 0;
	for(int dPhi=0; dPhi < 6; dPhi++){
	  int thisPhi = phi + dPhi;
	  thisPhi = thisPhi >= NUM_PHI ? thisPhi - NUM_PHI : thisPhi;	  

	  sumBottomMidSix += cc->interpRMS[polInd][thisPhi+NUM_PHI];
	  sumBottomMidSix += cc->interpRMS[polInd][thisPhi+2*NUM_PHI];
	}

	if(sumBottomMidSix > maxSumBottomMidSix[polInd]){
	  maxSumBottomMidSix[polInd] = sumBottomMidSix;
	  maxPhi[polInd] = phi;
	}
      }
      for(int dPhi=0; dPhi < 6; dPhi++){
	int thisPhi = maxPhi[polInd] + dPhi;
	thisPhi = thisPhi >= NUM_PHI ? thisPhi - NUM_PHI : thisPhi;

	maxSumTopSix[polInd] += cc->interpRMS[polInd][thisPhi];
      }
    }
    // for(int polInd=0; polInd < NUM_POL; polInd++){
    //   for(int ant=0; ant < NUM_SEAVEYS; ant++){
    // 	FancyFFTs::crossCorrelate(cc->numSamples, cc->ffts[polInd][ant],
    // 				  ccTemplate->ffts[polInd][ant],
    // 				  theCrossCorr);

    // 	corrFactors[polInd][ant] = TMath::MaxElement(cc->numSamples, theCrossCorr);
    // 	if(ant >= NUM_PHI){
    // 	  if(corrFactors[polInd][ant] >= threshVal){
    // 	    isCandidate++;
    // 	  }
    // 	}
    //   }
    // }



    delete usefulEvent;
    delete pedEvent;
    selfTriggeredBlastTree->Fill();

    p++;    
  }

  delete cc;
  // delete ccTemplate;
  // delete [] theCrossCorr;
  
  outFile->Write();
  outFile->Close();

  return 0;
}
