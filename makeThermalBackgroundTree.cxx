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

  TFile* inFile = OutputConvention::getFile("plotReconstructedMinBiasPlots_*.root");
  TTree* cutTree = (TTree*) inFile->Get("cutTree");

  UInt_t eventNumber = 0;
  cutTree->SetBranchAddress("eventNumber", &eventNumber);

  AnitaPol::AnitaPol_t pol;
  cutTree->SetBranchAddress("pol", (Int_t*)&pol);

  Double_t maxV, minV, absSumMaxMin;
  cutTree->SetBranchAddress("maxV", &maxV);
  cutTree->SetBranchAddress("minV",&minV);
  cutTree->SetBranchAddress("absSumMaxMin", &absSumMaxMin);
  AnalysisCuts::Status_t surfSaturation;
  cutTree->SetBranchAddress("surfSaturation",(Int_t*) &surfSaturation);

  Double_t theMaxPeakToPeakRatio;
  cutTree->SetBranchAddress("theMaxPeakToPeakRatio", &theMaxPeakToPeakRatio);
  AnalysisCuts::Status_t selfTriggeredBlastCut;
  cutTree->SetBranchAddress("selfTriggeredBlastCut",(Int_t*) &selfTriggeredBlastCut);

  Int_t deltaPhiSect;
  cutTree->SetBranchAddress("deltaPhiSect", &deltaPhiSect);
  AnalysisCuts::Status_t l3TriggerCut;
  cutTree->SetBranchAddress("l3TriggerCut", (Int_t*)&l3TriggerCut);


  Double_t deltaSolarPhiDeg, deltaSolarThetaDeg;
  cutTree->SetBranchAddress("deltaSolarPhiDeg", &deltaSolarPhiDeg);
  cutTree->SetBranchAddress("deltaSolarThetaDeg", &deltaSolarThetaDeg);
  AnalysisCuts::Status_t sunCut;
  cutTree->SetBranchAddress("sunCut", (Int_t*)&sunCut);


  Double_t imagePeak, hilbertPeak, fisher;
  cutTree->SetBranchAddress("imagePeak", &imagePeak);
  cutTree->SetBranchAddress("hilbertPeak", &hilbertPeak);
  cutTree->SetBranchAddress("fisher", &fisher);
  AnalysisCuts::Status_t thermalCut;
  cutTree->SetBranchAddress("thermalCut", (Int_t*)&thermalCut);


  Double_t peakRatio ,imagePeak2;
  cutTree->SetBranchAddress("peakRatio", &peakRatio);
  cutTree->SetBranchAddress("imagePeak2", &imagePeak2);
  AnalysisCuts::Status_t peakRatioCut;
  cutTree->SetBranchAddress("peakRatioCut", (int*)&peakRatioCut);

  Double_t recoThetaDeg;
  cutTree->SetBranchAddress("recoThetaDeg", &recoThetaDeg);
  AnalysisCuts::Status_t thetaAngleCut;
  cutTree->SetBranchAddress("thetaAngleCut", (int*)&thetaAngleCut);

  Double_t recoPhiDeg;
  cutTree->SetBranchAddress("recoPhiDeg", &recoPhiDeg);
  Double_t directionWrtNorth;
  cutTree->SetBranchAddress("directionWrtNorth", &directionWrtNorth);
  Adu5Pat* pat2;
  cutTree->SetBranchAddress("pat", &pat2);
  RawAnitaHeader* header2;
  cutTree->SetBranchAddress("header", &header2);




  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  TTree* outTree = new TTree("thermalTree", "Quiet Min Bias events");
  // UInt_t eventNumber;
  outTree->Branch("eventNumber", &eventNumber);
  // Double_t imagePeak;
  outTree->Branch("imagePeak", &imagePeak);
  // Double_t hilbertPeak;
  outTree->Branch("hilbertPeak", &hilbertPeak);

  Long64_t nEntries = cutTree->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    cutTree->GetEntry(entry);

    if(l3TriggerCut==0 && surfSaturation==0 && selfTriggeredBlastCut==0 && sunCut==0){
      outTree->Fill();
    }
    p.inc(entry, maxEntry);
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
