// -*- C++ -*-.
/***********************************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             Program to reconstruct HPol pulses from Wais Divide.
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

  CrossCorrelator* cc = new CrossCorrelator();
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");
  
  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);    
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
  calEventChain->BuildIndex("eventNumber");
  gpsChain->BuildIndex("eventNumber");  

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  CalibratedAnitaEvent* calEvent = NULL;
  calEventChain->SetBranchAddress("event", &calEvent);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  TTree* eventSummaryTree = new TTree("eventSummaryTree", "eventSummaryTree");
  AnitaEventSummary* eventSummary;
  eventSummaryTree->Branch("eventSummary", &eventSummary);
  
  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //100;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    if(header->eventNumber != 60832108){
      p++;
      continue;
    }
				       
    // gpsChain->GetEntry(entry);
    // calEventChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->eventNumber);
    calEventChain->GetEntryWithIndex(header->eventNumber);

    UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);
    UsefulAdu5Pat usefulPat(pat);
    cc->correlateEvent(usefulEvent);

    // cc->reconstructEvent(usefulEvent, header);
    
    eventSummary = new AnitaEventSummary(header, pat);
    // // eventSummary = new AnitaEventSummary(); //header, pat);
    // eventSummary->eventNumber = header->eventNumber;

    // Recontruction...
    for(Int_t polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      Int_t hypInd = 0;      
      AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) polInd;      
      const Int_t coherentDeltaPhi = 2;
      Double_t minY = 0;

      // cc->getPeakInfoTriggered(pol,
      // 			       eventSummary->peak[pol][hypInd].value,
      // 			       eventSummary->peak[pol][hypInd].phi,
      // 			       eventSummary->peak[pol][hypInd].theta);
      TH2D* hTriggeredImage = cc->makeTriggeredImage(pol,
						     eventSummary->peak[pol][hypInd].value,
						     eventSummary->peak[pol][hypInd].phi,
						     eventSummary->peak[pol][hypInd].theta,
						     header->getL3TrigPattern(pol));
			     
      hTriggeredImage->Write();
			       
      
      usefulPat.getSourceLonAndLatAtAlt(eventSummary->peak[pol][hypInd].phi*TMath::DegToRad(),
					eventSummary->peak[pol][hypInd].theta*TMath::DegToRad(),
					eventSummary->peak[pol][hypInd].longitude,
					eventSummary->peak[pol][hypInd].latitude,
					eventSummary->peak[pol][hypInd].altitude);

      TGraph* grTriggered = cc->makeCoherentlySummedWaveform(pol,
							     eventSummary->peak[pol][hypInd].phi,
							     eventSummary->peak[pol][hypInd].theta,
							     coherentDeltaPhi,
							     eventSummary->peak[pol][hypInd].snr);

      grTriggered->Write();
      
      TGraph* grTriggeredHilbert = FFTtools::getHilbertEnvelope(grTriggered);

      
      RootTools::getMaxMin(grTriggeredHilbert, eventSummary->coherent[pol][hypInd].peakHilbert, minY);

      delete hTriggeredImage;      
      delete grTriggered;
      delete grTriggeredHilbert;

      hypInd++;

      // cc->getPeakInfoZoom(pol,
      // 			  eventSummary->peak[pol][hypInd].value,
      // 			  eventSummary->peak[pol][hypInd].phi,
      // 			  eventSummary->peak[pol][hypInd].theta);
      TH2D* hZoomedImage = cc->makeZoomedImage(pol,
					       eventSummary->peak[pol][hypInd].value,
					       eventSummary->peak[pol][hypInd].phi,
					       eventSummary->peak[pol][hypInd].theta,
					       eventSummary->peak[pol][hypInd-1].phi,
					       eventSummary->peak[pol][hypInd-1].theta);	

      hZoomedImage->Write();
      
      usefulPat.getSourceLonAndLatAtAlt(eventSummary->peak[pol][hypInd].phi*TMath::DegToRad(),
					eventSummary->peak[pol][hypInd].theta*TMath::DegToRad(),
					eventSummary->peak[pol][hypInd].longitude,
					eventSummary->peak[pol][hypInd].latitude,
					eventSummary->peak[pol][hypInd].altitude);
      
      TGraph* grZoomed = cc->makeUpsampledCoherentlySummedWaveform(pol,
								   eventSummary->peak[pol][hypInd].phi,
								   eventSummary->peak[pol][hypInd].theta,
								   coherentDeltaPhi,
								   eventSummary->peak[pol][hypInd].snr);

      TGraph* grZoomedHilbert = FFTtools::getHilbertEnvelope(grZoomed);

      // Double_t minY = 0;
      RootTools::getMaxMin(grZoomedHilbert, eventSummary->coherent[pol][hypInd].peakHilbert, minY);


      delete hZoomedImage;
      delete grZoomed;
      delete grZoomedHilbert;

      hypInd++;
    }

    // Flags
    eventSummary->flags.isGood = 1;
    eventSummary->flags.isRF = header->getTriggerBitRF();
    eventSummary->flags.isSoftwareTrigger = header->getTriggerBitADU5() | header->getTriggerBitG12() | header->getTriggerBitSoftExt();
    eventSummary->flags.isPayloadBlast = 0; //!< To be determined
    eventSummary->flags.nadirFlag = 0; //!< Not sure I will use this.
    eventSummary->flags.strongCWFlag = 0; //!< Not sure I will use this.
    eventSummary->flags.isVarner = 0; //!< Not sure I will use this.
    eventSummary->flags.isVarner2 = 0; //!< Not sure I will use this.
    eventSummary->flags.pulser = AnitaEventSummary::EventFlags::NONE; //!< Not yet
    
    delete usefulEvent;

    eventSummaryTree->Fill();
    delete eventSummary;
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
