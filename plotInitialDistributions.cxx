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
    // fileName = TString::Format("~/UCL/ANITA/anita3Analysis/decimatedDistributions/initialDistributions260MHzAnd370MHzFiltered/initialDistributionsPlots_%d_*.root", run);
    fileName = TString::Format("~/UCL/ANITA/anita3Analysis/decimatedDistributions/initialDistributionPlots/initialDistributionsPlots_%d_*.root", run);    
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

  headChain->GetEntry(0);
  UInt_t firstRealTime = header->realTime;
  headChain->GetEntry(headChain->GetEntries()-1);
  UInt_t lastRealTime = header->realTime;  

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const int numTrigTypes = 2; // rf=0 and min bias=1
  const int numSunDirs = 2; // away from the sun=0, towards the sun=1
  // const int numUpDownDirs = 2; // downwards=0 and upwards=1;
  // const int numNoiseStates= 2; // 0 is quiet, 1 is noisy
  
  // TH2D* hPeakDirWrtNorth[NUM_POL][numTrigTypes][numSunDirs][numUpDownDirs][numNoiseStates];
  // TH2D* hPeakTheta[NUM_POL][numTrigTypes][numSunDirs][numUpDownDirs][numNoiseStates];
  // TH2D* hDeltaSolarPhiDeg[NUM_POL][numTrigTypes][numSunDirs][numUpDownDirs][numNoiseStates];
  // TH2D* hImagePeakHilbertPeak[NUM_POL][numTrigTypes][numSunDirs][numUpDownDirs][numNoiseStates];
  // TProfile2D* hpDeltaSolarPhiDeg[NUM_POL][numTrigTypes][numSunDirs][numUpDownDirs][numNoiseStates];

  TH2D* hPeakDirWrtNorth[NUM_POL][numTrigTypes][numSunDirs];
  TH2D* hPeakTheta[NUM_POL][numTrigTypes][numSunDirs];
  TH2D* hImagePeakHilbertPeak[NUM_POL][numTrigTypes][numSunDirs];
  
  const int numBinsPhi = 420;
  const int numBinsTheta = numBinsPhi/2;  
  const int numTimeBins = 1024*2;
  const double maxDirWrtNorth = 180;
  const double minDirWrtNorth = -180;
  const double maxTheta = 90;
  const double minTheta = -90;

  // const double deltaSolarPhiDegCut = 20; // degrees
  const double deltaSolarPhiDegCut = 10; // degrees  

  const Int_t numImagePeakBins = 1024;
  const Int_t numHilbertPeakBins = 1024;
  const Double_t maxHilbertPeak = 1000;

  
  for(int pol=0; pol<NUM_POL; pol++){
    for(int trig=0; trig < numTrigTypes; trig++){
      for(int sun=0; sun < numSunDirs; sun++){

	TString name = "hPeakDirWrtNorth_";
	name += trig == 0 ? "" : "minBias_";
	name += sun == 0 ? "awayFromSun_" : "towardsSun_";
	name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

	TString title = trig == 0 ? "RF triggered" : "Min bias";
	title += " peak direction w.r.t North vs. realTime ";
	title += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
	title += "; realTime; peak direction wrt. north (Degrees)";

	hPeakDirWrtNorth[pol][trig][sun] = new TH2D(name, title,
						    numTimeBins,
						    firstRealTime, lastRealTime+1,
						    numBinsPhi,
						    minDirWrtNorth, maxDirWrtNorth);

	name = "hPeakTheta_";
	name += trig == 0 ? "" : "minBias_";
	name += sun == 0 ? "awayFromSun_" : "towardsSun_";
	name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

	title = trig == 0 ? "RF triggered" : "Min bias";
	title += " Peak direction #theta vs. realTime ";
	title += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
	title += "; realTime; #theta (Degrees)";
	
	hPeakTheta[pol][trig][sun] = new TH2D(name, title,
					      numTimeBins, firstRealTime, lastRealTime+1,
					      numBinsTheta, minTheta, maxTheta);

	name = "hImagePeakHilbertPeak_";
	name += trig == 0 ? "" : "minBias_";
	name += sun == 0 ? "awayFromSun_" : "towardsSun_";
	name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

	title = trig == 0 ? " RF triggered" : " Min bias";
	title += " image peak vs. Hilbert Peak ";
	title += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";
	title += "; Image peak (no units); Hilbert peak (mV); Number of events";

	hImagePeakHilbertPeak[pol][trig][sun] = new TH2D(name, title,
							 numImagePeakBins,
							 0, 1,
							 numHilbertPeakBins,
							 0, maxHilbertPeak);
      }
    }
  }
  
  TH2D* hHeading = new TH2D("hHeading", "ANITA heading vs. realTime; realTime; heading (Degrees)",
			    1024, firstRealTime, lastRealTime+1,
			    numBinsPhi, 0, 360);

  TH1D* hDeltaSolarPhiDeg[NUM_POL];
  TH1D* hDeltaSolarThetaDeg[NUM_POL];
  TH1D* hDeltaSolarPhiDegZoom[NUM_POL];
  TH1D* hDeltaSolarThetaDegZoom[NUM_POL];

  for(int pol=0; pol<NUM_POL; pol++){
    TString name = pol==AnitaPol::kHorizontal ? "hDeltaSolarPhiDegHPol" : "hDeltaSolarPhiDegVPol";
    TString title = "#delta#phi_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #delta#phi_{sun} (Degrees); Number of events / bin";

    hDeltaSolarPhiDeg[pol] = new TH1D(name, title,
				      numBinsPhi, -180, 180);


    name = pol==AnitaPol::kHorizontal ? "hDeltaSolarThetaDegHPol" : "hDeltaSolarThetaDegVPol";
    title = "#delta#theta_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #delta#theta_{sun} (Degrees); Number of events / bin";

    hDeltaSolarThetaDeg[pol] = new TH1D(name, title,
					numBinsPhi, -180, 180);

    name = pol==AnitaPol::kHorizontal ? "hDeltaSolarPhiDegZoomHPol" : "hDeltaSolarPhiDegZoomVPol";
    title = "#delta#phi_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #delta#phi_{sun} (Degrees); Number of events / bin";

    hDeltaSolarPhiDegZoom[pol] = new TH1D(name, title,
				      numBinsPhi, -20, 20);


    name = pol==AnitaPol::kHorizontal ? "hDeltaSolarThetaDegZoomHPol" : "hDeltaSolarThetaDegZoomVPol";
    title = "#delta#theta_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #delta#theta_{sun} (Degrees); Number of events / bin";

    hDeltaSolarThetaDegZoom[pol] = new TH1D(name, title,
					numBinsPhi, -20, 20);

  }



  TGraph* grSunPhiDeg = new TGraph();
  TGraph* grSunThetaDeg = new TGraph();
  
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    headChain->GetEntryWithIndex(eventSummary->eventNumber);
    gpsChain->GetEntryWithIndex(eventSummary->eventNumber);

    UsefulAdu5Pat usefulPat(pat);
    // Double_t sunHeading = -1*usefulPat.getAzimuthOfSunRelativeToNorth();
    // Double_t solarPhiDeg = usefulPat.getAzimuthOfSun();
    Double_t solarPhiDeg, solarThetaDeg;
    usefulPat.getSunPosition(solarPhiDeg, solarThetaDeg);

    // std::cerr << pat->heading - solarPhiDeg << "\t" << usefulPat.getAzimuthOfSunRelativeToNorth() << std::endl;
    if((entry % 16)==0){
      grSunPhiDeg->SetPoint(grSunPhiDeg->GetN(), header->realTime, pat->heading - solarPhiDeg);
      grSunThetaDeg->SetPoint(grSunThetaDeg->GetN(), header->realTime, solarThetaDeg);    
    }
    
    for(Int_t polInd=0; polInd < NUM_POL; polInd++){

      Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;
      Double_t recoThetaDeg = eventSummary->peak[polInd][1].theta;

      Double_t directionWrtNorth = pat->heading - recoPhiDeg;
      directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
      directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;

      Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(solarPhiDeg, recoPhiDeg);
      Double_t deltaSolarThetaDeg = solarThetaDeg - recoThetaDeg;

      Double_t deltaSolarAngleDeg = TMath::Sqrt(deltaSolarThetaDeg*deltaSolarThetaDeg + deltaSolarPhiDeg*deltaSolarPhiDeg);
      
      if(TMath::Abs(directionWrtNorth) > maxDirWrtNorth){
	std::cerr << recoPhiDeg << "\t" << pat->heading << std::endl;
      }

      Int_t trig = eventSummary->flags.isMinBiasTrigger > 0 ? 1 : 0;
      // Int_t sun = TMath::Abs(deltaSolarPhiDeg) > deltaSolarPhiDegCut ? 0 : 1;
      // Int_t sun = TMath::Abs(deltaSolarPhiDeg) > deltaSolarPhiDegCut ? 0 : 1;
      Int_t sun = TMath::Abs(deltaSolarAngleDeg) > deltaSolarPhiDegCut ? 0 : 1;
      // Int_t upDir = recoThetaDeg > 0 ? 0 : 1

      // UInt_t realTime = header->realTime;
      // Int_t noisy = 1;

      // if((realTime >= timeCut1Low && realTime < timeCut1High) || (realTime >= timeCut2Low && realTime < timeCut2High)){
      // 	noisy = 0;
      // }
      
      // hPeakDirWrtNorth[polInd][trig][sun][upDir][noisy]->Fill(header->realTime,
      // 							      directionWrtNorth);
      // hPeakTheta[polInd][trig][sun][upDir][noisy]->Fill(header->realTime,
      // 							recoThetaDeg);
      // hDeltaSolarPhiDeg[polInd][trig][sun][upDir][noisy]->Fill(deltaSolarPhiDeg, deltaSolarThetaDeg);
      // hpDeltaSolarPhiDeg[polInd][trig][sun][upDir][noisy]->Fill(deltaSolarPhiDeg,
      // 								deltaSolarThetaDeg,
      // 								eventSummary->peak[polInd][1].value);
      // hImagePeakHilbertPeak[polInd][trig][sun][upDir][noisy]->Fill(eventSummary->peak[polInd][1].value,
      // 								   eventSummary->coherent[polInd][1].peakHilbert);

      hPeakDirWrtNorth[polInd][trig][sun]->Fill(header->realTime,
						directionWrtNorth);
      hPeakTheta[polInd][trig][sun]->Fill(header->realTime,
					  recoThetaDeg);

      // const UInt_t timeCut1Low  = 1419000000;
      // const UInt_t timeCut1High = 1419600000;
      // const UInt_t timeCut2Low  = 1419800000;
      // const UInt_t timeCut2High = 1420280000;      
      // UInt_t realTime = header->realTime;
      // if((realTime >= timeCut1Low && realTime < timeCut1High) ||
      // 	 (realTime >= timeCut2Low && realTime < timeCut2High)){
      {
	hImagePeakHilbertPeak[polInd][trig][sun]->Fill(eventSummary->peak[polInd][1].value,
						       eventSummary->coherent[polInd][1].peakHilbert);
      }

      if(eventSummary->flags.isMinBiasTrigger > 0){
	// std::cerr << deltaSolarPhiDeg << "\t" << deltaSolarThetaDeg << std::endl;
	hDeltaSolarPhiDeg[polInd]->Fill(deltaSolarPhiDeg);
	hDeltaSolarThetaDeg[polInd]->Fill(deltaSolarThetaDeg);

	hDeltaSolarPhiDegZoom[polInd]->Fill(deltaSolarPhiDeg);
	hDeltaSolarThetaDegZoom[polInd]->Fill(deltaSolarThetaDeg);	
      }
    }
    hHeading->Fill(header->realTime, pat->heading);
    
    p++;
  }


  grSunPhiDeg->SetName("grSunPhiDeg");
  grSunPhiDeg->SetTitle("Predicted sun #phi (Degrees)");
  grSunThetaDeg->SetName("grSunThetaDeg");
  grSunThetaDeg->SetTitle("Predicted sun #theta (Degrees)");
  grSunPhiDeg->Write();
  grSunThetaDeg->Write();
  
  outFile->Write();
  outFile->Close();

  return 0;
}
