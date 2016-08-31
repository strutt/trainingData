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
#include "AnalysisCuts.h"

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
  // const Int_t lastRun = argc==3 ? atoi(argv[2]) : firstRun;
  const Int_t lastRun = firstRun;

  const int cutStep = 4;  

  TChain* headChain = new TChain("headTree");
  TChain* indexedHeadChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  TChain* dataQualityChain = new TChain("dataQualityTree");  

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }
    if(run == 198 || run == 287){
      continue;
    }

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);
    indexedHeadChain->Add(fileName);
    
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    
    // fileName = TString::Format("~/UCL/ANITA/anita3Analysis/trainingData/filter260-370-400-762-speed/reconstructWaisPlots_%d_*.root", run);
    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/trainingData/filter260-370-400-762-speed/reconstructDecimatedPlots_%d_*.root", run);    
    eventSummaryChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/trainingData/filter260-370-400-762/makeDecimatedDataQualityTreesPlots_%d_*.root", run);
    dataQualityChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }
  if(gpsChain->GetEntries()==0){
    std::cerr << "Unable to find gps files!" << std::endl;
    return 1;
  }
  if(eventSummaryChain->GetEntries()==0){
    std::cerr << "Unable to find eventSummary files!" << std::endl;
    return 1;
  }
  if(dataQualityChain->GetEntries()==0){
    std::cerr << "Unable to find data quality files!" << std::endl;
    return 1;
  }

  indexedHeadChain->BuildIndex("eventNumber");    
  // headChain->BuildIndex("eventNumber");
  // gpsChain->BuildIndex("eventNumber");

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  Double_t peakToPeak[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("peakToPeak", peakToPeak);
  UInt_t eventNumberDQ;
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumberDQ);
  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);  

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  TTree* outTree = new TTree("eventSummaryTree", "eventSummaryTree");
  AnitaEventSummary* outEventSummary = NULL;
  outTree->Branch("eventSummary", &outEventSummary);

  TTree* outTree2 = new TTree("adu5PatTree", "adu5PatTree");
  Adu5Pat* pat2 = NULL;
  outTree2->Branch("pat", &pat2);
  
  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    Int_t entry2 = indexedHeadChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
    headChain->GetEntry(entry2);
    // headChain->GetEntryWithIndex(eventSummary->eventNumber);
    // gpsChain->GetEntryWithIndex(eventSummary->eventNumber);


      
    AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
    if(eventSummary->peak[AnitaPol::kHorizontal][0].value > eventSummary->peak[AnitaPol::kVertical][0].value){
      pol = AnitaPol::kHorizontal;
    }

    dataQualityChain->GetEntry(entry);    

    if(eventSummary->eventNumber != eventNumberDQ){
      std::cerr << "???" << eventSummary->eventNumber << "\t" << eventNumberDQ << std::endl;
    }

      
      
    
    Double_t maxRatio;
    AnalysisCuts::Status_t selfTriggeredBlastCut;
    selfTriggeredBlastCut = AnalysisCuts::applyBottomToTopRingPeakToPeakRatioCut(pol, peakToPeak[pol], maxRatio);
    if(cutStep >= 1 && selfTriggeredBlastCut==AnalysisCuts::kFail){
      p.inc(entry, maxEntry);	
      continue;
    }
      
      
    // Get event info
    const int peakInd = 0;
      
    Double_t recoPhiDeg = eventSummary->peak[pol][peakInd].phi;
    recoPhiDeg += recoPhiDeg < 0 ? DEGREES_IN_CIRCLE : 0;      
    Double_t recoThetaDeg = eventSummary->peak[pol][peakInd].theta;
    Double_t imagePeak = eventSummary->peak[pol][peakInd].value;      
    Double_t hilbertPeak = eventSummary->coherent[pol][peakInd].peakHilbert;

      

      
    // CUT FLOW
    // Step 2: cut phi-sector angle triggers
    Int_t deltaPhiSect = NUM_PHI/2;

    AnalysisCuts::Status_t l3TriggerCut;
    l3TriggerCut = AnalysisCuts::L3TriggerDirectionCut(pol, header, recoPhiDeg, deltaPhiSect);
    if(cutStep >= 2 && l3TriggerCut==AnalysisCuts::kFail){
      p.inc(entry, maxEntry);
      continue;
    }

    gpsChain->GetEntry(entry2);
    

    // CUT FLOW
    // Step 3: cut phi-direction relative to sun      
    Double_t solarPhiDeg = eventSummary->sun.phi;
    // Double_t solarThetaDeg = -1*eventSummary->sun.theta;

    Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(recoPhiDeg, solarPhiDeg);
    solarPhiDeg = solarPhiDeg < 0 ? solarPhiDeg + 360 : solarPhiDeg;
    // Double_t deltaSolarThetaDeg = recoThetaDeg - solarThetaDeg;
      
    AnalysisCuts::Status_t sunCut;
    sunCut = AnalysisCuts::applySunPointingCut(deltaSolarPhiDeg);
    if(cutStep>=3 && sunCut==AnalysisCuts::kFail){
      p.inc(entry, maxEntry);	
      continue;
    }

    AnalysisCuts::Status_t thermalCut;
    Double_t fisher;
    thermalCut = AnalysisCuts::applyThermalBackgroundCut(imagePeak, hilbertPeak, fisher);
    if(cutStep>=4 && thermalCut==AnalysisCuts::kFail){
      p.inc(entry, maxEntry);	
      continue;
    }
    
    UsefulAdu5Pat usefulPat(pat);

    Int_t goodPeak = 0; // for now

    // AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
    // if(eventSummary->peak[AnitaPol::kHorizontal][0].value > eventSummary->peak[AnitaPol::kVertical][0].value){
    //   pol = AnitaPol::kHorizontal;
    // }

    // assignment constructor works?
    outEventSummary = eventSummary;
    outEventSummary->realTime = header->realTime;

    for(int polInd=0; polInd < AnitaPol::kNotAPol; polInd++){
      for(int peak=0; peak < AnitaEventSummary::maxDirectionsPerPol; peak++){
	outEventSummary->peak[polInd][peak].latitude = -9999;
	outEventSummary->peak[polInd][peak].longitude = -9999;
      }
    }

    if(recoPhiDeg < 0) recoPhiDeg += 360;
    if(recoPhiDeg >= 360) recoPhiDeg -= 360;
    Double_t phiWave = recoPhiDeg*TMath::DegToRad();
    Double_t thetaWave = -1*recoThetaDeg*TMath::DegToRad();
    Double_t sourceLon, sourceLat, sourceAlt = 0; // altitude zero for now
    int success = usefulPat.getSourceLonAndLatAtAlt(phiWave, thetaWave,
						    sourceLon, sourceLat, sourceAlt);

    if(success){
      outEventSummary->peak[pol][goodPeak].latitude = sourceLat;
      outEventSummary->peak[pol][goodPeak].longitude = sourceLon;
      outEventSummary->peak[pol][goodPeak].altitude = sourceAlt;      
      outEventSummary->peak[pol][goodPeak].distanceToSource = SPEED_OF_LIGHT_NS*usefulPat.getTriggerTimeNsFromSource(sourceLat, sourceLon, sourceAlt);
    }

    outTree->Fill();
    
    pat2 = pat;
    
    outTree2->Fill();

    // delete pat2;
    
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
