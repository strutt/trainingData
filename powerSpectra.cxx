// -*- C++ -*-.
/*****************************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to find peak cross correlation offsets between antenna pairs in pulses from Wais Divide.
*************************************************************************************************************** */

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
#include "AveragePowerSpectrum.h"

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

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);
  calEventChain->BuildIndex("eventNumber");
  
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const Double_t deltaT = 1./2.6;
  const Int_t numSamps = 256;  
  AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS];
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      avePowSpecs[polInd][ant] = new AveragePowerSpectrum(deltaT, numSamps);
    }
  }
  const Int_t numEvents = 100;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->realTime);
    calEventChain->GetEntryWithIndex(header->eventNumber);        

    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

    for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	TGraph* gr = usefulEvent->getGraph(ant, pol);
	TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], deltaT, numSamps);
	avePowSpecs[polInd][ant]->add(grInterp);
	if((entry % numEvents) == 0){
	  TString name = TString::Format("grAvePs_%d_%d_%u", pol, ant, header->eventNumber);
	  TGraph* grAvePs = avePowSpecs[polInd][ant]->get(name, name);
	  grAvePs->Write();
	  delete grAvePs;

	  avePowSpecs[polInd][ant]->empty();
	}
	delete grInterp;
	delete gr;
      }
    }

    
    delete usefulEvent;

    p++;
  }
  
  outFile->Write();
  outFile->Close();

  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      delete avePowSpecs[polInd][ant];
    }
  }

  return 0;
}
