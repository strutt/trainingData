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

  
  const double ratioCut = 2.8;

  const int polInd = AnitaPol::kVertical;

  const int numGoodTimes = 1;
  UInt_t goodTimesStart[numGoodTimes] = {1419100000};
  UInt_t goodTimesEnd[numGoodTimes] = {1419500000};

  
  TChain* headChain = new TChain("headTree");
  TChain* eventSummaryChain = new TChain("eventSummaryTree");  
  TChain* dataQualityChain = new TChain("dataQualityTree");    

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);


    fileName = TString::Format("filter260-370-400-762/reconstructDecimatedPlots_%d_*.root", run);
    // TString fileName = TString::Format("test400MHzExt/reconstructDecimatedPlots_%d_*.root", run); 
    eventSummaryChain->Add(fileName);

    
    fileName = TString::Format("filter260-370-400-762/makeDecimatedDataQualityTreesPlots_%d*.root", run);
    // fileName = TString::Format("testNewDataQuality/makeDecimatedDataQualityTreesPlots_%d*.root", run);    
    dataQualityChain->Add(fileName);
  }

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
    headChain->GetEntry(entry);

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
    
    
    dataQualityChain->GetEntry(entry);
    
    Double_t maxRatio = 0;
    for(Int_t polInd=0; polInd < NUM_POL; polInd++){      
      for(int phi=0; phi < NUM_PHI; phi++){
	if(polInd==1 && phi==7){
	  continue;
	}
	Double_t ratio = peakToPeak[polInd][phi+32]/peakToPeak[polInd][phi];
	if(ratio > maxRatio){
	  maxRatio = ratio;
	}
      }
    }
    if(maxRatio > ratioCut){
      p.inc(entry, maxEntry);
      continue;
    }

    eventSummaryChain->GetEntry(entry);

    if(eventSummary->eventNumber != eventNumberDQ){
      std::cerr << "???" << eventSummary->eventNumber << "\t" << eventNumberDQ << std::endl;
    }
    
    const double deltaSolarPhiDegCut = 10; // degrees    
    Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;
    Double_t solarPhiDeg  = eventSummary->sun.phi;
    Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(solarPhiDeg, recoPhiDeg);    
    Int_t sun = TMath::Abs(deltaSolarPhiDeg) > deltaSolarPhiDegCut ? 0 : 1;

    if(sun==1){

      p.inc(entry, maxEntry);
      continue;
    }

    // Double_t deltaSolarThetaDeg = solarThetaDeg - recoThetaDeg;

    
    imagePeak = eventSummary->peak[polInd][1].value;
    hilbertPeak = eventSummary->coherent[polInd][1].peakHilbert;

    outTree->Fill();
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}

