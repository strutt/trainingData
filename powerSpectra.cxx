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

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //20000; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const Double_t deltaT = 1./2.6;
  const Int_t numSamps = 256;
  // AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS];
  AveragePowerSpectrum* avePowSpecs[AnitaPol::kNotAPol][NUM_SEAVEYS];
  // AveragePowerSpectrum* northPs[AnitaPol::kNotAPol][NUM_SEAVEYS];
  AveragePowerSpectrum* northPs[AnitaPol::kNotAPol];  
  // AveragePowerSpectrum* southPs[AnitaPol::kNotAPol][NUM_SEAVEYS];
  AveragePowerSpectrum* southPs[AnitaPol::kNotAPol];
  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    TString histBaseName = TString::Format("hNorthPs_%d", polInd);
    northPs[polInd] = new AveragePowerSpectrum(histBaseName, deltaT, numSamps,
						    AveragePowerSpectrum::kSummed);
    histBaseName = TString::Format("hSouthPs_%d", polInd);
    southPs[polInd] = new AveragePowerSpectrum(histBaseName, deltaT, numSamps,
					       AveragePowerSpectrum::kSummed);
    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      histBaseName = TString::Format("hAvePowSpec_%d_%d", polInd, ant);
      avePowSpecs[polInd][ant] = new AveragePowerSpectrum(histBaseName, deltaT, numSamps,
							  AveragePowerSpectrum::kSummed);
      // histBaseName = TString::Format("hRollingAvePowSpec_%d_%d", polInd, ant);
      // avePowSpecs[polInd][ant] = new AveragePowerSpectrum(histBaseName, deltaT, numSamps,
      // // 							  AveragePowerSpectrum::kRolling);
      // histBaseName = TString::Format("hNorthPs_%d_%d", polInd, ant);
      // northPs[polInd][ant] = new AveragePowerSpectrum(histBaseName, deltaT, numSamps,
      // 						      AveragePowerSpectrum::kSummed);
      // histBaseName = TString::Format("hSouthPs_%d_%d", polInd, ant);
      // southPs[polInd][ant] = new AveragePowerSpectrum(histBaseName, deltaT, numSamps,
      // 						      AveragePowerSpectrum::kSummed);
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
    headChain->GetEntry(entry);

    if((header->trigType & 1)==1) {
      p++;
      continue;
    }

    gpsChain->GetEntryWithIndex(header->realTime);
    if(pat->heading < 0){
      p++;
      continue;
    }

    calEventChain->GetEntryWithIndex(header->eventNumber);

    
    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
    Double_t heading = pat->heading;

    Int_t northestPhiSector = -1;
    Int_t southestPhiSector = -1;

    Double_t northestPhiDeg = 180;
    Double_t southestPhiDeg = 0;
    for(Int_t phiSector=0; phiSector < NUM_PHI; phiSector++){
      Double_t antPhiDeg = geom->getAntPhiPositionRelToAftFore(phiSector)*TMath::RadToDeg();
      Double_t thisAntHeading = RootTools::getDeltaAngleDeg(heading, antPhiDeg);

      if(TMath::Abs(thisAntHeading) > southestPhiDeg){
	southestPhiSector = phiSector;
	southestPhiDeg = thisAntHeading;
      }
      if(TMath::Abs(thisAntHeading) < northestPhiDeg){
	northestPhiSector = phiSector;
	northestPhiDeg = thisAntHeading;
      }

      if(TMath::Abs(thisAntHeading) > 180){
	std::cerr << "think harder! " << thisAntHeading << std::endl;
      }
    }
    // std::cout << heading << "\t" << northestPhiSector << "\t" << southestPhiSector << std::endl;

    for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;
      for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
	Int_t phiSector = ant%NUM_PHI;

	TGraph* gr = usefulEvent->getGraph(ant, pol);
	TGraph* grInterp = RootTools::interpolateWithStartTime(gr, gr->GetX()[0], deltaT, numSamps);
	delete gr;

	avePowSpecs[pol][ant]->add(grInterp);
	if(!(ant==4 && pol==AnitaPol::kHorizontal)){
	  if(phiSector==northestPhiSector){
	    // for(int i=0; i<129; i++){
	    //   std::cerr << grInterp->GetY()[i] << ", ";
	    // }
	    // std::cerr << std::endl;
	    northPs[pol]->add(grInterp);
	  }
	  else if(phiSector==southestPhiSector){
	    southPs[pol]->add(grInterp);
	  }
	}
	delete grInterp;
      }
    }
    
    delete usefulEvent;

    p++;
  }

  for(Int_t polInd = 0; polInd < AnitaPol::kNotAPol; polInd++){
    TString name = TString::Format("grNorthAvePs_%d_%d_%d", polInd, firstRun, lastRun);
    TGraph* grAvePs = northPs[polInd]->getScaled(name, name);
    // grAvePs = northPs[polInd]->get(name, name);      
    grAvePs->Write();
    delete grAvePs;
    delete northPs[polInd];

    name = TString::Format("grSouthAvePs_%d_%d_%d", polInd, firstRun, lastRun);
    grAvePs = southPs[polInd]->getScaled(name, name);
    grAvePs->Write();
    delete grAvePs;
    delete southPs[polInd];      

    for(Int_t ant=0; ant<NUM_SEAVEYS; ant++){
      name = TString::Format("grAvePs_%d_%d_%d_%d", polInd, ant, firstRun, lastRun);
      grAvePs = avePowSpecs[polInd][ant]->getScaled(name, name);
      grAvePs->Write();
      delete grAvePs;
      delete avePowSpecs[polInd][ant];
     
      // name = TString::Format("grNorthAvePs_%d_%d_%d_%d", polInd, ant, firstRun, lastRun);
      // // grAvePs = northPs[polInd][ant]->getScaled(name, name);
      // grAvePs = northPs[polInd][ant]->get(name, name);      
      // grAvePs->Write();
      // delete grAvePs;
      // delete northPs[polInd][ant];

      // name = TString::Format("grSouthAvePs_%d_%d_%d_%d", polInd, ant, firstRun, lastRun);
      // grAvePs = southPs[polInd][ant]->getScaled(name, name);
      // grAvePs->Write();
      // delete grAvePs;
      // delete southPs[polInd][ant];      
    }
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
