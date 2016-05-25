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
  const Double_t maxDeltaTriggerTimeNs = 1200;

  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;


  TString lindaFileName = "../antennaPositionCalib/newLindaNumbers_4steps_WAISHPOL_NEW11_cosminV3_nfixedBug_2016_02_05_time_15_42_15.txt";

  Int_t insertion = CrossCorrelator::directlyInsertGeometry(lindaFileName, pol);  
  if(insertion > 0){
    std::cerr << "Couldn't find file " << lindaFileName.Data() << std::endl;
    return 1;
  }
  
  CrossCorrelator* cc = new CrossCorrelator();
  // cc->kDoSimpleSatelliteFiltering = 1;
  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* calEventChain = new TChain("eventTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);
    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/calEventFile%d.root", run, run);
    calEventChain->Add(fileName);
  }
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

  TTree* angResTree = new TTree("angResTree", "angResTree");

  Double_t globalPeak = 0;
  Double_t globalPhiDeg = 0;
  Double_t globalThetaDeg = 0;
  Double_t zoomPeak = 0;
  Double_t zoomPhiDeg = 0;
  Double_t zoomThetaDeg = 0;
  Double_t deltaThetaDeg = 0;
  Double_t deltaPhiDeg = 0;

  Double_t thetaExpected = 0;
  Double_t phiExpected = 0;
  UInt_t eventNumber = 0;
  UInt_t triggerTimeNs = 0;
  UInt_t triggerTimeNsExpected = 0;
  Double_t heading = 0;
  UInt_t l3TrigPattern = 0;
  UInt_t l3TrigPatternH = 0;
  Int_t run = 0;
  UInt_t realTime = 0;
  Double_t snr = 0;
  
  angResTree->Branch("globalPeak", &globalPeak);
  angResTree->Branch("globalPhiDeg", &globalPhiDeg);
  angResTree->Branch("globalThetaDeg", &globalThetaDeg);

  angResTree->Branch("zoomPeak", &zoomPeak);
  angResTree->Branch("zoomPhiDeg", &zoomPhiDeg);
  angResTree->Branch("zoomThetaDeg", &zoomThetaDeg);
  angResTree->Branch("deltaPhiDeg", &deltaPhiDeg);
  angResTree->Branch("deltaThetaDeg", &deltaThetaDeg);
  
  angResTree->Branch("thetaExpected", &thetaExpected);
  angResTree->Branch("phiExpected", &phiExpected);
  angResTree->Branch("triggerTimeNs", &triggerTimeNs);
  angResTree->Branch("triggerTimeNsExpected", &triggerTimeNsExpected);
  angResTree->Branch("heading", &heading);
  angResTree->Branch("l3TrigPattern", &l3TrigPattern);
  angResTree->Branch("l3TrigPatternH", &l3TrigPatternH);
  angResTree->Branch("eventNumber", &eventNumber);
  angResTree->Branch("run", &run);

  Double_t lat, lon, alt;
  angResTree->Branch("lat", &lat);
  angResTree->Branch("lon", &lon);
  angResTree->Branch("alt", &alt);
  angResTree->Branch("snr", &snr);
  angResTree->Branch("realTime", &realTime);

  Long64_t nEntries = headChain->GetEntries();
  Long64_t maxEntry = 0; //1000; //10000;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  Int_t numSaved = 0;
  Int_t maxToSave = 5;
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    headChain->GetEntry(entry);

    gpsChain->GetEntry(entry);        
    if((header->trigType & 1)==1 && header->l3TrigPatternH > 0){
      UsefulAdu5Pat usefulPat(pat);
      triggerTimeNsExpected = usefulPat.getWaisDivideTriggerTimeNs();
      triggerTimeNs = header->triggerTimeNs;
      Int_t deltaTriggerTimeNs = Int_t(triggerTimeNs) - Int_t(triggerTimeNsExpected);
      if(TMath::Abs(deltaTriggerTimeNs) < maxDeltaTriggerTimeNs){
      // if(true){
	eventNumber = header->eventNumber;
	run = header->run;

	calEventChain->GetEntry(entry);
	heading = usefulPat.heading;
	realTime = header->realTime;
	l3TrigPattern = header->l3TrigPattern;
	l3TrigPatternH = header->l3TrigPatternH;
	UsefulAnitaEvent* usefulEvent = new UsefulAnitaEvent(calEvent);

	usefulPat.getThetaAndPhiWaveWaisDivide(thetaExpected, phiExpected);
	phiExpected *= TMath::RadToDeg();
	thetaExpected *= TMath::RadToDeg();

	phiExpected = phiExpected < 0 ? phiExpected + 360 : phiExpected;
	phiExpected = phiExpected >= 360 ? phiExpected - 360 : phiExpected;

	cc->reconstructEvent(usefulEvent, 1, 1);
	cc->getCoarsePeakInfo(pol, 0, globalPeak, globalPhiDeg, globalThetaDeg);
	cc->getFinePeakInfo(pol, 0, zoomPeak, zoomPhiDeg, zoomThetaDeg);

	deltaPhiDeg = RootTools::getDeltaAngleDeg(zoomPhiDeg, phiExpected);
	deltaThetaDeg = zoomThetaDeg - thetaExpected;

	TGraph * gr = cc->makeUpsampledCoherentlySummedWaveform(pol, zoomPhiDeg, zoomThetaDeg, 1, snr);

	usefulPat.getSourceLonAndLatAtAlt(zoomPhiDeg*TMath::DegToRad(), zoomThetaDeg*TMath::DegToRad(),
					  lon, lat, alt);

	if(numSaved < maxToSave){
	  Double_t peakValue;
	  Double_t peakPhiDeg;
	  Double_t peakThetaDeg;	  
	  
	  TH2D* hGlobal = cc->getMap(pol, peakValue, peakPhiDeg, peakThetaDeg);
	  hGlobal->Write();
	  delete hGlobal;

	  TH2D* hZoomed = cc->getZoomMap(pol);
	  hZoomed->Write();	  
	  delete hZoomed;

	  gr->Write();
	  
	  numSaved++;
	}
	else{
	  delete gr;	
	}

	angResTree->Fill();

	delete usefulEvent;
      }
    }
    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}
