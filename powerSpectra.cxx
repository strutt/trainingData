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
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run); 
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

  TH2D* hHeadingTest = new TH2D("hHeadingTest", "Heading test; Heading (Degrees); Northest phi-sector",
				360, 0, 360,
				16, 0, 16);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //5000; //20000; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const Double_t deltaT = 1./2.6;
  const Int_t upsampleFactor = 16;
  const Int_t numSamps = 256*upsampleFactor;

  const char* dirNames[NUM_PHI] = {"N", "NEN", "NE", "NEE",
				   "E", "SEE", "SE", "SES",
				   "S", "SWS", "SW", "SWW",
				   "W", "NWW", "NW", "NWN"};

  AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS];
  AveragePowerSpectrum* directionalPowSpecs[AnitaPol::kNotAPol][NUM_PHI];
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      TString histBaseName = TString::Format("avePowSpec_%d_%d", polInd, ant);
      avePowSpecs[polInd][ant] = new AveragePowerSpectrum(histBaseName, "test",
							  deltaT, numSamps,
							  AveragePowerSpectrum::kSummed);

    }
    for(int phi=0; phi<NUM_PHI; phi++){
      TString histBaseName = TString::Format("aveDirPowSpec_%d_%s", polInd, dirNames[phi]);
      directionalPowSpecs[polInd][phi] = new AveragePowerSpectrum(histBaseName, dirNames[phi],
								  deltaT, numSamps,
								  AveragePowerSpectrum::kSummed);
    }
    
  }

  // const Int_t numEvents = maxEntry;
  // headChain->GetEntry(0);
  // const UInt_t startTime = header->realTime;
  // headChain->GetEntry(maxEntry-1);
  // const UInt_t endTime = header->realTime;
  // const Int_t numTimeBins = 100;

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    // Get header.
    headChain->GetEntry(entry);

    // If event is RF trigger then skip.
    if((header->trigType & 1)==1) {
      p++;
      continue;
    }

    // If have bad heading info then skip.
    gpsChain->GetEntryWithIndex(header->realTime);
    if(pat->heading < 0){
      p++;
      continue;
    }

    calEventChain->GetEntryWithIndex(header->eventNumber);
    
    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
    Double_t heading = pat->heading;

    Int_t northestPhiSector = -1;
    Double_t northestPhiDeg = 180;

    Int_t directions[NUM_PHI];
    // What do I mean by directions[NUM_PHI]?
    // In this array 0 means the "northest" phi-sector,
    // 1 means NW by N
    // 2 means NW
    // 3 mean NW by W
    // 4 means W
    // 5 means SW by W
    // 6 means SW     
    // 7 means SW by S
    // 8 means S 
    // 9 means SE by S     
    // 10 mean SE  
    // 11 means SE by E
    // 12 means E
    // 13 means NE by E
    // 14 means NE
    // 15 means NE by N    

    // heading = (entry % 180)*2;
    // std::cout << std::endl << std::endl;    
    // std::cout << heading << ": ";
    

    // Here I figure out which phi-sector points closest to north.
    for(Int_t phiSector=0; phiSector < NUM_PHI; phiSector++){
      Double_t antPhiDeg = geom->getAntPhiPositionRelToAftFore(phiSector)*TMath::RadToDeg();
      while(antPhiDeg >= 180){
	antPhiDeg -= 360;
      }
      while(antPhiDeg < -180){
	antPhiDeg += 360;
      }

      Double_t thisAntHeading = RootTools::getDeltaAngleDeg(heading, antPhiDeg);

      // std::cout << heading << "\t" << phiSector << "\t" << thisAntHeading << "\t" << std::endl;
      
      if(TMath::Abs(thisAntHeading) < northestPhiDeg){
	northestPhiSector = phiSector;
	northestPhiDeg = TMath::Abs(thisAntHeading);
      }

      if(TMath::Abs(thisAntHeading) > 180){
	std::cerr << "Think harder! " << thisAntHeading << std::endl;
      }
    }

    hHeadingTest->Fill(heading, northestPhiSector);
    
    // Then loop over phi-sectors and fill in the direction array.
    for(Int_t phi=0; phi<NUM_PHI; phi++){
      directions[phi] = (northestPhiSector + phi)%NUM_PHI;
      // std::cout << directions[phi] << ", ";      
    }
    // std::cout << std::endl << std::endl;

    // Now I do the heavy lifting of making the average power spectra.
    for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	Int_t phiSector = ant%NUM_PHI;

	TGraph* gr = usefulEvent->getGraph(ant, pol);
	TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], deltaT, numSamps);
	delete gr;

	avePowSpecs[polInd][ant]->add(grInterp);
	// std::cout << pol << "\t" << ant
	// 	  << avePowSpecs[polInd][ant]->getRayleighHistogramFromFrequencyMHz(250)->Integral()
	// 	  << std::endl;

	Int_t direction = -1000000000;
	for(Int_t dirInd=0; dirInd<NUM_PHI; dirInd++){
	  if(directions[dirInd]==phiSector){
	    direction = dirInd;
	    break;
	  }
	}
	if(!(ant==4 && pol==AnitaPol::kHorizontal)){
	  directionalPowSpecs[polInd][direction]->add(grInterp);
	}
	delete grInterp;
      }
    }
    
    delete usefulEvent;

    p++;
  }

  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      avePowSpecs[polInd][ant]->Write();
    }
    for(int phi=0; phi<NUM_PHI; phi++){
      directionalPowSpecs[polInd][phi]->Write();
    }
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
