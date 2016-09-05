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

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }
    if(run == 198 || run == 287){
      continue;
    }
    
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/trainingData/filter260-370-400-762-speed/reconstructedWaisPlots_%d_*.root", run);
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
  
  headChain->BuildIndex("eventNumber");
  gpsChain->BuildIndex("eventNumber");

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
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


  AntarcticaMapPlotter* hs[NUM_POL];
  
  for(int pol=0; pol<NUM_POL; pol++){

    TString name = "hLatLon_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    TString title = "Lat/Long";
    title += pol == AnitaPol::kHorizontal ? "HPol" : "VPol"; 

    hs[pol] = new AntarcticaMapPlotter(name, title,1024, 1024);
  }  
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    headChain->GetEntryWithIndex(eventSummary->eventNumber);
    gpsChain->GetEntryWithIndex(eventSummary->eventNumber);

    UsefulAdu5Pat usefulPat(pat);
    
    for(Int_t polInd=0; polInd < NUM_POL; polInd++){

      Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;

      if(recoPhiDeg < 0) recoPhiDeg += 360;
      if(recoPhiDeg >= 360) recoPhiDeg -= 360;      


      Double_t recoThetaDeg = eventSummary->peak[polInd][1].theta;

      
      Double_t phiWave = recoPhiDeg*TMath::DegToRad();
      Double_t thetaWave = recoThetaDeg*TMath::DegToRad();      
      Double_t sourceLon, sourceLat;
      int success = usefulPat.getSourceLonAndLatAltZero(phiWave, thetaWave,
							sourceLon, sourceLat);

      if(success){
	hs[polInd]->Fill(sourceLat, sourceLon);
      }
    }
    
    p++;
  }

  hs[0]->Write();
  hs[1]->Write();  
  
  outFile->Write();
  outFile->Close();

  return 0;
}
