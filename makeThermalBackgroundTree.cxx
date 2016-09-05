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


  const int numGoodTimes = 1;
  UInt_t goodTimesStart[numGoodTimes] = {1419220000};
  UInt_t goodTimesEnd[numGoodTimes] = {1419500000};

  // UInt_t goodTimesStart[numGoodTimes] = {1419220000};
  // UInt_t goodTimesEnd[numGoodTimes] = {1419690000};
  // UInt_t goodTimesStart[numGoodTimes] = {1419100000};
  // UInt_t goodTimesEnd[numGoodTimes] = {1419500000};


  TChain* headChain = new TChain("headTree");
  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  TChain* dataQualityChain = new TChain("dataQualityTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);


    fileName = TString::Format("filter260-370-400-762-5peaks/reconstructMinBiasPlots_%d_*.root", run);
    // fileName = TString::Format("filter260-370-400-762-5peaks/reconstructDecimatedPlots_%d_*.root", run);
    // TString fileName = TString::Format("test400MHzExt/reconstructDecimatedPlots_%d_*.root", run);
    eventSummaryChain->Add(fileName);


    fileName = TString::Format("filter260-370-400-762/makeMinBiasDataQualityTreesPlots_%d*.root", run);
    // fileName = TString::Format("testNewDataQuality/makeDecimatedDataQualityTreesPlots_%d*.root", run);
    dataQualityChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }
  headChain->BuildIndex("eventNumber");
  if(eventSummaryChain->GetEntries()==0){
    std::cerr << "Unable to find eventSummary files!" << std::endl;
    return 1;
  }
  if(dataQualityChain->GetEntries()==0){
    std::cerr << "Unable to find dataQualityFiles files!" << std::endl;
    return 1;
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
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

  TTree* outTree = new TTree("thermalTree", "HQ Thermal Background Events");
  UInt_t eventNumber;
  outTree->Branch("eventNumber", &eventNumber);
  Double_t imagePeak;
  outTree->Branch("imagePeak", &imagePeak);
  Double_t hilbertPeak;
  outTree->Branch("hilbertPeak", &hilbertPeak);

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);

    Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber, 0);
    if(headEntry < 0){
      std::cerr << "Now what!?\t" << headEntry << "\t" << eventSummary->eventNumber << std::endl;
    }
    else{
      headChain->GetEntry(headEntry);

      bool isGoodTime = false;
      for(int i=0; i < numGoodTimes; i++){
	if(header->realTime >= goodTimesStart[i] && header->realTime < goodTimesEnd[i]){
	  isGoodTime = true;
	}
      }
      if(isGoodTime==false){
	p.inc(entry, maxEntry);
	continue;
      }



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
      if(selfTriggeredBlastCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }

      // Get event info
      const int peakInd = 0;

      Double_t recoPhiDeg = eventSummary->peak[pol][peakInd].phi;
      recoPhiDeg += recoPhiDeg < 0 ? DEGREES_IN_CIRCLE : 0;
      imagePeak = eventSummary->peak[pol][peakInd].value;
      hilbertPeak = eventSummary->coherent[pol][peakInd].peakHilbert;




      // CUT FLOW
      // Step 2: cut phi-sector angle triggers
      Int_t deltaPhiSect = NUM_PHI/2;

      AnalysisCuts::Status_t l3TriggerCut;
      l3TriggerCut = AnalysisCuts::L3TriggerDirectionCut(pol, header, recoPhiDeg, deltaPhiSect);
      if(l3TriggerCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }



      // CUT FLOW
      // Step 3: cut phi-direction relative to sun
      Double_t solarPhiDeg = eventSummary->sun.phi;
      Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(recoPhiDeg, solarPhiDeg);
      AnalysisCuts::Status_t sunCut;
      sunCut = AnalysisCuts::applySunPointingCut(deltaSolarPhiDeg);
      if(sunCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }

      if(imagePeak > 0.075){
	std::cerr << "\nIP: run " << header->run << "\teventNumber: " << header->eventNumber << std::endl;
      }
      else if(hilbertPeak > 75){
	std::cerr << "\nHP: run " << header->run << "\teventNumber: " << header->eventNumber << std::endl;
      }

      imagePeak = eventSummary->peak[pol][peakInd].value;
      hilbertPeak = eventSummary->coherent[pol][peakInd].peakHilbert;

      outTree->Fill();
    }
    p.inc(entry, maxEntry);
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
