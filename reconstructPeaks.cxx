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
  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  // TChain* dataQualityChain = new TChain("dataQualityTree");    

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    // if(run >= 211 && run <= 263){
    //   continue;
    // }
    // if(run == 241){
    //   continue;
    // }
      
    
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/decimatedDistributions/initialDistributions260MHzAnd370MHzFiltered/initialDistributionsPlots_%d_*.root", run);
    eventSummaryChain->Add(fileName);

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

  ProgressBar p2(1);
  headChain->BuildIndex("eventNumber");
  p2++;
  
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  UInt_t gpsEventNumber = 0;
  gpsChain->SetBranchAddress("eventNumber", &gpsEventNumber);
  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  AntarcticaMapPlotter* MapVPol = new AntarcticaMapPlotter("ampHPol", "HPol", 1024, 1024);
  AntarcticaMapPlotter* MapHPol = new AntarcticaMapPlotter("ampVPol", "VPol", 1024, 1024);

  TTree* t = new TTree("recoTree", "recoTree");
  Double_t sourceLon[NUM_POL];
  t->Branch("sourceLon", &sourceLon[0], "sourceLon[2]/D");
  Double_t sourceLat[NUM_POL];
  t->Branch("sourceLat", &sourceLat[0], "sourceLat[2]/D");
  Double_t sourceAlt[NUM_POL];
  t->Branch("sourceAlt", &sourceAlt[0], "sourceAlt[2]/D");
  Double_t imagePeak[NUM_POL];
  t->Branch("imagePeak", &imagePeak[0], "imagePeak[2]/D");
  Double_t thetaDeg[NUM_POL];
  t->Branch("thetaDeg", &thetaDeg[0], "thetaDeg[2]/D");
  Double_t phiDeg[NUM_POL];
  t->Branch("phiDeg", &phiDeg[0], "phiDeg[2]/D");  
  
  Double_t distMeters[NUM_POL];
  t->Branch("distMeters", &distMeters[0], "distMeters[2]/D");
  
  UInt_t eventNumber;
  t->Branch("eventNumber", &eventNumber);
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    // Int_t headBytes = headChain->GetEntryWithIndex(eventSummary->eventNumber);
    // Int_t patBytes = gpsChain->GetEntryWithIndex(eventSummary->eventNumber);

    Int_t entry2 = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber);

    if(entry2 < 0){
      std::cerr << entry2 << "\t" << eventSummary->eventNumber << "\t"
    		<< header->eventNumber << "\t" << header->realTime << "\t" << pat->realTime << "\t"
    		<< std::endl;
      p++;
      continue;
    }

    headChain->GetEntry(entry2);
    gpsChain->GetEntry(entry2);
    
    for(Int_t polInd=0; polInd < NUM_POL; polInd++){

      Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;
      phiDeg[polInd] = recoPhiDeg;
      Double_t recoThetaDeg = eventSummary->peak[polInd][1].theta;
      thetaDeg[polInd] = recoThetaDeg;
      UsefulAdu5Pat usefulPat(pat);

      usefulPat.getSourceLonAndLatAtAlt(recoPhiDeg*TMath::DegToRad(), recoThetaDeg*TMath::DegToRad(),
					sourceLon[polInd], sourceLat[polInd], sourceAlt[polInd]);

      distMeters[polInd] = usefulPat.getDistanceFromSource(sourceLat[polInd], sourceLon[polInd], sourceAlt[polInd]);

      // if(sourceAlt[polInd] > -9000){
      // if(thetaDeg[polInd] >= 6){
      if(thetaDeg[polInd] >= 6 && sourceAlt[polInd] > -9000){
	if(polInd==0){
	  MapHPol->Fill(sourceLat[polInd], sourceLon[polInd], imagePeak[polInd]);
	}
	else{
	  MapVPol->Fill(sourceLat[polInd], sourceLon[polInd], imagePeak[polInd]);
	}
      }
      imagePeak[polInd] = eventSummary->peak[polInd][1].value;
    }
    eventNumber = eventSummary->eventNumber;    
    t->Fill();
    p.inc(entry, maxEntry);
  }
  // grSunPhiDeg->SetName("grSunPhiDeg");
  // grSunPhiDeg->SetTitle("Predicted sun #phi (Degrees)");
  // grSunThetaDeg->SetName("grSunThetaDeg");
  // grSunThetaDeg->SetTitle("Predicted sun #theta (Degrees)");
  // grSunPhiDeg->Write();  
  // grSunThetaDeg->Write();
  MapVPol->Write();
  delete MapVPol;
  MapHPol->Write();
  delete MapHPol;
  

  outFile->Write();
  outFile->Close();

  return 0;
}




