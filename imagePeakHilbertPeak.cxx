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
  const Double_t maxDeltaTriggerTimeNs = 1200;

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  // geom->useKurtAnitaIIINumbers(1);
  // AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
  // for(int surf=0; surf<NUM_SURF; surf++){
  //   for(int chan=0; chan<NUM_CHAN; chan++){
  //     cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
  //   }
  // }

  for(Int_t run=firstRun; run<=lastRun; run++){
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
  TTree* snrTree = new TTree("snrTree", "snrTree");

  UInt_t eventNumber = 0;
  Double_t zoomPeak = 0;
  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t hilbertPeak = 0;
  Double_t hilbertPeakTime = 0;
  Double_t coherentPeakToPeak = 0;
  Double_t rms = 0;
  Double_t snr = 0;
  // std::vector<Double_t>* deltaPhiDeg = NULL;

  snrTree->Branch("eventNumber", &eventNumber);
  snrTree->Branch("zoomPeak", &zoomPeak);
  snrTree->Branch("zoomPhiDeg", &zoomPhiDeg);
  snrTree->Branch("zoomThetaDeg", &zoomThetaDeg);
  snrTree->Branch("hilbertPeak", &hilbertPeak);
  snrTree->Branch("hilbertPeakTime", &hilbertPeakTime);
  snrTree->Branch("coherentPeakToPeak", &coherentPeakToPeak);
  snrTree->Branch("rms", &rms);
  snrTree->Branch("snr", &snr);
  
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  CrossCorrelator* cc = new CrossCorrelator();

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const Double_t amountOfBackTime = 10;
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    calEventChain->GetEntryWithIndex(header->eventNumber);
    gpsChain->GetEntryWithIndex(header->realTime);    
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      UInt_t triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      UInt_t triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
	eventNumber = header->eventNumber;

	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	cc->correlateEvent(usefulEvent, pol);
	UShort_t l3TrigPattern = header->getL3TrigPattern(pol);

	Double_t triggeredPeak, triggeredPhiDeg,triggeredThetaDeg;
	TH2D* hTriggeredImageH = cc->makeTriggeredImage(pol, triggeredPeak, triggeredPhiDeg,
							triggeredThetaDeg, l3TrigPattern);

	TH2D* hZoomedImageH = cc->makeZoomedImage(pol, zoomPeak, zoomPhiDeg,
						 zoomThetaDeg, l3TrigPattern,
						 triggeredPhiDeg, triggeredThetaDeg);
	
	zoomPhiDeg = zoomPhiDeg < 0 ? zoomPhiDeg + 360 : zoomPhiDeg;
	zoomPhiDeg = zoomPhiDeg >= 360 ? zoomPhiDeg - 360 : zoomPhiDeg;

	TGraph* grCoherent = cc->makeCoherentlySummedWaveform(pol, zoomPhiDeg, zoomThetaDeg, l3TrigPattern);

	Double_t maxY = 0, maxX = 0, minY = 0, minX = 0;
	rms = -1000;
	snr = -1000;
	hilbertPeak = -1000;
	hilbertPeakTime = -1000;
	coherentPeakToPeak = -1000;
	if(grCoherent!=NULL){
	  RootTools::getLocalMaxToMin(grCoherent, maxY, maxX, minY, minX);
	  coherentPeakToPeak = maxY - minY;

	  TGraph* grHilbert = FFTtools::getHilbertEnvelope(grCoherent);

	  Int_t maxSamp = TMath::LocMax(grHilbert->GetN(),grHilbert->GetY());
	  hilbertPeak = grHilbert->GetY()[maxSamp];
	  hilbertPeakTime = grHilbert->GetX()[maxSamp];

	  Double_t sum=0;
	  Double_t square=0;
	  Int_t numSampsCounted=0;
	  for(Int_t phiSect=0; phiSect < NUM_PHI; phiSect++){
	    if(RootTools::getBit(phiSect, l3TrigPattern)){
	      for(Int_t ring=0; ring<NUM_RING; ring++){
		Int_t ant = ring*NUM_PHI + phiSect;

		// Don't delete
		TGraph* gr = cc->grs[pol][ant];

		Int_t numSamps = gr->GetN();
		for(Int_t samp=numSamps-1; samp >= 0; --samp){
		  sum += gr->GetY()[samp];
		  square += gr->GetY()[samp]*gr->GetY()[samp];
		  numSampsCounted++;
		  if(gr->GetX()[numSamps-1] - gr->GetX()[samp] > amountOfBackTime){
		    break;
		  }
		}
	      }
	    }
	  }
	  Double_t mean = sum/numSampsCounted;
	  rms = TMath::Sqrt(square/numSampsCounted - mean*mean);
	  snr = coherentPeakToPeak/(2*rms);

	  delete grCoherent;
	  delete grHilbert;
	}      
	delete hTriggeredImageH;
	delete hZoomedImageH;
	
	snrTree->Fill();

	delete usefulEvent;
	
      }
    }
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
