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


  if(argc!=2){
    std::cerr << "Usage: " << argv[0] << " [useTimeCut] " << std::endl;
    return 1;
  }

  const bool useTimeCut = atoi(argv[1]) > 0 ? true : false;
  const Int_t firstRun = useTimeCut ? 202 : 130;
  const Int_t lastRun = useTimeCut ? 251 : 434;

  std::cout << firstRun << "\t" << lastRun << std::endl;

  const int numGoodTimes = 1;
  UInt_t goodTimesStart[numGoodTimes] = {1419220000};
  UInt_t goodTimesEnd[numGoodTimes] = {1419500000};



  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  TChain* dataQualityChain = new TChain("dataQualityTree");


  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    // if(filterInd > 0 && filterInd < 4){
    //   fileName = TString::Format("rayleigh/reconstructFilteredMinBiasPlots_%d_%d_*.root", run, filterInd);
    // }
    // else{
    fileName = TString::Format("filter260-370-400-762-5peaks/reconstructMinBiasPlots_%d_*.root", run);
    // }
    eventSummaryChain->Add(fileName);

    fileName = TString::Format("finalDataQuality/makeSlimMinBiasDataQualityTreesPlots_%d_*.root", run);
    dataQualityChain->Add(fileName);
  }

  if(headChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
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
    std::cerr << "Unable to find dataQualityFiles files!" << std::endl;
    return 1;
  }

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);


  Double_t maxPeakToPeakRatio[NUM_POL];
  dataQualityChain->SetBranchAddress("maxPeakToPeakRatio", maxPeakToPeakRatio);
  Double_t theMaxVolts[NUM_POL];
  dataQualityChain->SetBranchAddress("theMaxVolts", theMaxVolts);
  Double_t theMinVolts[NUM_POL];
  dataQualityChain->SetBranchAddress("theMinVolts",theMinVolts);


  UInt_t eventNumberDQ;
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumberDQ);
  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){;
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  TTree* cutTree = new TTree("cutTree", "cutTree");
  UInt_t eventNumber = 0;
  cutTree->Branch("eventNumber", &eventNumber);


  AnitaPol::AnitaPol_t pol;
  cutTree->Branch("pol", (Int_t*)&pol);

  Double_t maxV, minV, absSumMaxMin;
  cutTree->Branch("maxV", &maxV);
  cutTree->Branch("minV",&minV);
  cutTree->Branch("absSumMaxMin", &absSumMaxMin);
  AnalysisCuts::Status_t surfSaturation;
  cutTree->Branch("surfSaturation",(Int_t*) &surfSaturation);

  Double_t theMaxPeakToPeakRatio;
  cutTree->Branch("theMaxPeakToPeakRatio", &theMaxPeakToPeakRatio);
  AnalysisCuts::Status_t selfTriggeredBlastCut;
  cutTree->Branch("selfTriggeredBlastCut",(Int_t*) &selfTriggeredBlastCut);

  Int_t deltaPhiSect;
  cutTree->Branch("deltaPhiSect", &deltaPhiSect);
  AnalysisCuts::Status_t l3TriggerCut;
  cutTree->Branch("l3TriggerCut", (Int_t*)&l3TriggerCut);


  Double_t deltaSolarPhiDeg, deltaSolarThetaDeg;
  cutTree->Branch("deltaSolarPhiDeg", &deltaSolarPhiDeg);
  cutTree->Branch("deltaSolarThetaDeg", &deltaSolarThetaDeg);
  AnalysisCuts::Status_t sunCut;
  cutTree->Branch("sunCut", (Int_t*)&sunCut);


  Double_t imagePeak, hilbertPeak, fisher;
  cutTree->Branch("imagePeak", &imagePeak);
  cutTree->Branch("hilbertPeak", &hilbertPeak);
  cutTree->Branch("fisher", &fisher);
  AnalysisCuts::Status_t thermalCut;
  cutTree->Branch("thermalCut", (Int_t*)&thermalCut);


  Double_t peakRatio ,imagePeak2;
  cutTree->Branch("peakRatio", &peakRatio);
  cutTree->Branch("imagePeak2", &imagePeak2);
  AnalysisCuts::Status_t peakRatioCut;
  cutTree->Branch("peakRatioCut", (int*)&peakRatioCut);

  Double_t recoThetaDeg;
  cutTree->Branch("recoThetaDeg", &recoThetaDeg);
  AnalysisCuts::Status_t thetaAngleCut;
  cutTree->Branch("thetaAngleCut", (int*)&thetaAngleCut);

  Double_t recoPhiDeg;
  cutTree->Branch("recoPhiDeg", &recoPhiDeg);
  Double_t directionWrtNorth;
  cutTree->Branch("directionWrtNorth", &directionWrtNorth);
  Adu5Pat* pat2;
  cutTree->Branch("pat", &pat2);
  RawAnitaHeader* header2;
  cutTree->Branch("header", &header2);


  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);


  std::cerr << "building index" << std::endl;
  headChain->BuildIndex("eventNumber");
  dataQualityChain->BuildIndex("eventNumber");
  std::cerr << "done" << std::endl;

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);

    Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber, 0);
    if(headEntry < 0){
      std::cerr << "Now what!?\t" << headEntry << "\t" << eventSummary->eventNumber << std::endl;
    }
    else{
      headChain->GetEntry(headEntry);
      eventNumber = eventSummary->eventNumber;

      pol = AnitaPol::kVertical;
      if(eventSummary->peak[AnitaPol::kHorizontal][0].value > eventSummary->peak[AnitaPol::kVertical][0].value){
	pol = AnitaPol::kHorizontal;
      }

      // timing selection...
      bool isGoodTime = false;
      for(int i=0; i < numGoodTimes; i++){
	if(header->realTime >= goodTimesStart[i] && header->realTime < goodTimesEnd[i]){
	  isGoodTime = true;
	}
      }
      if(useTimeCut==true && isGoodTime==false){
	p.inc(entry, maxEntry);
	continue;
      }


      // dataQualityChain->GetEntry(entry);
      dataQualityChain->GetEntryWithIndex(header->eventNumber);

      if(eventSummary->eventNumber != eventNumberDQ){
	std::cerr << "??? " << eventSummary->eventNumber << "\t" << eventNumberDQ << std::endl;
      }




      maxV = TMath::Max(theMaxVolts[0],  theMaxVolts[1]);
      minV = TMath::Max(theMinVolts[0],  theMinVolts[1]);
      absSumMaxMin = maxV + minV;
      surfSaturation = AnalysisCuts::applySurfSaturationCutBetter(maxV, minV, absSumMaxMin);

      theMaxPeakToPeakRatio = TMath::Max(maxPeakToPeakRatio[0], maxPeakToPeakRatio[1]);
      selfTriggeredBlastCut = AnalysisCuts::applyBottomToTopRingPeakToPeakRatioCut(theMaxPeakToPeakRatio);


      // Get event info
      const int peakInd = 0;

      recoPhiDeg = eventSummary->peak[pol][peakInd].phi;
      recoPhiDeg += recoPhiDeg < 0 ? DEGREES_IN_CIRCLE : 0;
      recoThetaDeg = eventSummary->peak[pol][peakInd].theta;
      imagePeak = eventSummary->peak[pol][peakInd].value;
      hilbertPeak = eventSummary->coherent[pol][peakInd].peakHilbert;




      // CUT FLOW
      // Step 2: cut phi-sector angle triggers
      deltaPhiSect = NUM_PHI/2;

      l3TriggerCut = AnalysisCuts::L3TriggerDirectionCut(pol, header, recoPhiDeg, deltaPhiSect);


      gpsChain->GetEntry(headEntry);




      // CUT FLOW
      // Step 3: cut phi-direction relative to sun
      Double_t solarPhiDeg = eventSummary->sun.phi;
      Double_t solarThetaDeg = -1*eventSummary->sun.theta;

      deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(recoPhiDeg, solarPhiDeg);
      solarPhiDeg = solarPhiDeg < 0 ? solarPhiDeg + 360 : solarPhiDeg;
      deltaSolarThetaDeg = recoThetaDeg - solarThetaDeg;


      sunCut = AnalysisCuts::applySunPointingCut(deltaSolarPhiDeg);



      thermalCut = AnalysisCuts::applyThermalBackgroundCut(imagePeak, hilbertPeak, fisher);


      if(recoPhiDeg < 0) recoPhiDeg += 360;
      else if(recoPhiDeg >= 360) recoPhiDeg -= 360;


      if(recoPhiDeg < 0){
	std::cerr << recoPhiDeg <<"\t" << std::endl;
      }


      if(recoPhiDeg < 0) recoPhiDeg += 360;
      else if(recoPhiDeg >= 360) recoPhiDeg -= 360;


      if(recoPhiDeg < 0){
	std::cerr << recoPhiDeg <<"\t" << std::endl;
      }


      directionWrtNorth = RootTools::getDeltaAngleDeg(pat->heading, recoPhiDeg);

      directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
      directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;
      imagePeak2 = eventSummary->peak[pol][1].value;



      peakRatioCut = AnalysisCuts::applyImagePeakRatioCut(imagePeak, imagePeak2, peakRatio);


      thetaAngleCut = AnalysisCuts::applyThetaAngleCut(recoThetaDeg);


      header2 = (RawAnitaHeader*) header->Clone();
      pat2 = (Adu5Pat*) pat->Clone();

      cutTree->Fill();

    }

    p.inc(entry, maxEntry);

  }

  outFile->Write();
  outFile->Close();

  return 0;
}
