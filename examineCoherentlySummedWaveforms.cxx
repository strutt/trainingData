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
#include "TGlobal.h"
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

  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* eventChain = new TChain("eventTree");  
  TChain* eventSummaryChain = new TChain("eventSummaryTree");

  // gDirectory->pwd();
  // // // std::cout << gDirectory << std::endl;  
  // auto thecal = AnitaEventCalibrator::Instance();
  // std::cout << thecal->relativeCableDelays[0][0][0] << std::endl;
  // gDirectory->pwd();

  
  // std::cout << gDirectory << std::endl;  
  // return 0;
  const double maxDirWrtNorth = 180;
  const double deltaSolarPhiDegCut = 10; // degrees
  
  
  // const double minDirWrtNorth = -180;

  int previousRun = 0;
  
  const int numRuns = lastRun - firstRun + 1;
  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);    
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    eventChain->Add(fileName);    

    fileName = TString::Format("filter260and370/reconstructMinBiasPlots_%d_*.root", run);
    eventSummaryChain->Add(fileName);
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

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  CalibratedAnitaEvent* calEvent = NULL;
  eventChain->SetBranchAddress("event", &calEvent);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);

  CrossCorrelator* cc = new CrossCorrelator();  
  CrossCorrelator::SimpleNotch notch260("n260Notch", "260MHz Satellite And 200MHz Notch Notch",
					260-26, 260+26);
  CrossCorrelator::SimpleNotch notch370("n370Notch", "370MHz Satellite Notch",
					370-26, 370+26);
  cc->addNotch(notch260);
  cc->addNotch(notch370);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }


  headChain->GetEntry(0);
  UInt_t firstRealTime = header->realTime;
  headChain->GetEntry(headChain->GetEntries()-2);
  UInt_t lastRealTime = header->realTime;

  // std::cerr << firstRealTime << "\t" << lastRealTime << std::endl;


  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  // const int numTrigTypes = 2; // rf=0 and min bias=1

  headChain->BuildIndex("eventNumber");

  Int_t numGr[NUM_POL] = {0};
  TGraph* grAve[NUM_POL] = {NULL};

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    
    Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber, 0);    
    if(headEntry < 0){
      std::cerr << "Now what!?\t" << headEntry << "\t" << eventSummary->eventNumber << std::endl;
    }
    else{
      headChain->GetEntry(headEntry);
      gpsChain->GetEntry(headEntry);

      Double_t solarPhiDeg, solarThetaDeg;

      solarPhiDeg = eventSummary->sun.phi;
      solarThetaDeg = eventSummary->sun.theta;


      UsefulAnitaEvent* usefulEvent = NULL;
      
      for(Int_t polInd=0; polInd < NUM_POL; polInd++){

	Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;

	if(recoPhiDeg < 0) recoPhiDeg += 360;
	else if(recoPhiDeg >= 360) recoPhiDeg -= 360;

	Double_t recoThetaDeg = eventSummary->peak[polInd][1].theta;	
	//std::cerr << "\t" << pat->heading << "\t" << recoPhiDeg << "\t" << directionWrtNorth << std::endl;

	
	Double_t directionWrtNorth = RootTools::getDeltaAngleDeg(pat->heading, recoPhiDeg);
	directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
	directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;

	// std::cerr << "\t" << header->realTime - firstRealTime << "\t" << lastRealTime - header->realTime << "\t" << (header->realTime >= firstRealTime) << "\t" << (header->realTime < lastRealTime) << std::endl;	
	

	Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(solarPhiDeg, recoPhiDeg);
	Double_t deltaSolarThetaDeg = solarThetaDeg - recoThetaDeg;


	if(TMath::Abs(directionWrtNorth) > maxDirWrtNorth){
	  std::cerr << recoPhiDeg << "\t" << pat->heading << std::endl;
	}

	Int_t sun = TMath::Abs(deltaSolarPhiDeg) > deltaSolarPhiDegCut ? 0 : 1;

	// hDeltaSolarPhiDeg[polInd]->Fill(deltaSolarPhiDeg);
	// hDeltaSolarThetaDeg[polInd]->Fill(deltaSolarThetaDeg);
	// hDeltaSolarPhiDegZoom[polInd]->Fill(solarPhiDeg, deltaSolarPhiDeg);
	// hDeltaSolarThetaDegZoom[polInd]->Fill(solarPhiDeg, deltaSolarThetaDeg);

	// hThetaDeg[polInd]->Fill(recoThetaDeg);

	if(sun){
	  // hDeltaSolarThetaDegZoomVsTheta[polInd]->Fill(solarThetaDeg, deltaSolarThetaDeg);	
	  // hDeltaSolarThetaDegZoomVsPhi[polInd]->Fill(solarPhiDeg, deltaSolarThetaDeg);
	  // hDeltaSolarThetaDegZoomVsTimeOfDay[polInd]->Fill(pat->timeOfDay/1000, deltaSolarThetaDeg);
	  // hDeltaSolarPhiDegZoomVsTheta[polInd]->Fill(solarThetaDeg, deltaSolarPhiDeg);
	  // hDeltaSolarPhiDegZoomVsPhi[polInd]->Fill(solarPhiDeg, deltaSolarPhiDeg);
	  // hDeltaSolarPhiDegZoomVsTimeOfDay[polInd]->Fill(pat->timeOfDay/1000, deltaSolarPhiDeg);
	}
	else{

	  /// pick event numbers
	  const double deltaX = 1420050000-1419789000;
	  const double deltaY = 53;
	  const Double_t timeAtZeroHeading = 1420050000; //ish
	  const Double_t timeAtHeadingEqMinus53 = 1419789000; //ish
	  const Double_t theY = (header->realTime - timeAtZeroHeading)*deltaY/deltaX;	  
	  
	  if(header->realTime >= timeAtHeadingEqMinus53 && header->realTime < timeAtHeadingEqMinus53 + 2*deltaX){
	    if(TMath::Abs(RootTools::getDeltaAngleDeg(directionWrtNorth, theY)) < 10 && recoThetaDeg > -10){
	      // yes

	      
	      if(usefulEvent == NULL){
		eventChain->GetEntry(headEntry);
		usefulEvent = new UsefulAnitaEvent(calEvent);      	      
	      }
	      
	      cc->getNormalizedInterpolatedTGraphs(usefulEvent, AnitaPol::AnitaPol_t (polInd));
	      cc->doFFTs(AnitaPol::AnitaPol_t (polInd));

	      
	      if(header->run != previousRun){
		std::cout << "At run " << header->run << std::endl;
		previousRun = header->run;
	      }

	      Double_t snr = 0;
	      TGraph* gr = cc->makeCoherentlySummedWaveform((AnitaPol::AnitaPol_t) polInd,
							    recoPhiDeg,
							    recoThetaDeg,
							    0, snr);

	      TGraph* grPs = FFTtools::makePowerSpectrumMilliVoltsNanoSeconds(gr);
	      delete gr;	
	      if(grAve[polInd] == NULL){
		grAve[polInd] = grPs;
		numGr[polInd]++;
	      }
	      else{
		for(int freqInd=0; freqInd < grAve[polInd]->GetN(); freqInd++){
		  grAve[polInd]->GetY()[freqInd] += grPs->GetY()[freqInd];
		}
		delete grPs;
		numGr[polInd]++;		
	      }
	    }
	  }
	}
      }
      if(usefulEvent != NULL){
	delete usefulEvent;
	usefulEvent = NULL;
      }
      // if(pat->heading < 0 || pat->heading >= 360){
      // 	std::cerr << "bad heading ? " << pat->heading << "\t" << header->run << "\t" << header->eventNumber << std::endl;
      // }
    }
    // p++;
    p.inc(entry, maxEntry);
  }

  // grSunPhiDeg->SetName("grSunPhiDeg");
  // grSunPhiDeg->SetTitle("Predicted sun #phi (Degrees)");
  // grSunThetaDeg->SetName("grSunThetaDeg");
  // grSunThetaDeg->SetTitle("Predicted sun Elevation (Degrees)");
  // grSunPhiDeg->Write();
  // grSunThetaDeg->Write();
  outFile->cd();


  for(Int_t polInd=0; polInd < NUM_POL; polInd++){
    TString title = "Average power spectrum of coherently summed wave";
    title += polInd == 0 ? "HPol" : "VPol";
    grAve[polInd]->SetTitle(title);

    for(int freqInd=0; freqInd < grAve[polInd]->GetN(); freqInd++){
      grAve[polInd]->GetY()[freqInd] /= numGr[polInd];
      grAve[polInd]->GetY()[freqInd] /= numGr[polInd];

      grAve[polInd]->GetY()[freqInd] = 10*TMath::Log10(grAve[polInd]->GetY()[freqInd]);
    }

    TString name = "grAvePowSpecCohSumWave";
    name += polInd == 0 ? "HPol" : "VPol";
    grAve[polInd]->SetName(name);

    grAve[polInd]->Write();
    delete grAve[polInd];
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}

