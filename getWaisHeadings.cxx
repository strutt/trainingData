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

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsFile%d.root", run, run);
    gpsChain->Add(fileName);
  }
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  gpsChain->BuildIndex("realTime");
  
  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  TTree* snrFriendTree = new TTree("snrFriendTree", "snrFriendTree");

  UInt_t eventNumber = 0;
  UInt_t triggerTimeNs = 0;
  UInt_t triggerTimeNsExpected = 0;
  Double_t heading = 0;
  Double_t realTime = 0;
  Double_t distanceFromWais = 0;
  
  snrFriendTree->Branch("triggerTimeNs", &triggerTimeNs);
  snrFriendTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);  
  snrFriendTree->Branch("eventNumber", &eventNumber);
  snrFriendTree->Branch("heading", &heading);
  snrFriendTree->Branch("realTime", &realTime);
  snrFriendTree->Branch("distanceFromWais", &distanceFromWais);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //5000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    headChain->GetEntry(entry);
    gpsChain->GetEntryWithIndex(header->realTime);    
    if((header->trigType & 1)==1){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
	eventNumber = header->eventNumber;
	heading = usefulPat.heading;
	realTime = header->realTime;
	distanceFromWais = triggerTimeNsExpected*1e-9*C_LIGHT/1e3;
	snrFriendTree->Fill();

      }
    }
    p++;
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
