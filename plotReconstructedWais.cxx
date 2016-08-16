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
  // TChain* dataQualityChain = new TChain("dataQualityTree");    

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("filter260and370/reconstructWaisPlots_%d_*.root", run);
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
  // if(dataQualityChain->GetEntries()==0){
  //   std::cerr << "Unable to find data quality files!" << std::endl;
  //   return 1;
  // }

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

  headChain->BuildIndex("eventNumber");
  
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

  TH2D* hPeakDirWrtNorth[NUM_POL];
  TH2D* hPeakTheta[NUM_POL];
  TH2D* hImagePeakHilbertPeak[NUM_POL];
  TH2D* hImagePeakTime[NUM_POL];
  TH2D* hImagePeakPhi[NUM_POL];
  TH2D* hImagePeakPhi0[NUM_POL];
  TH2D* hImagePeakDeltaPhi[NUM_POL];
  TH2D* hPhi0DeltaPhi[NUM_POL];
  TH2D* hPhi0DeltaPhi2[NUM_POL];
  TH2D* hPhi0Theta0[NUM_POL];
  TH2D* hHilbertPeakTime[NUM_POL];

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
  // const Double_t maxHilbertPeak = 1000;
  const Double_t maxHilbertPeak = 2048;  

  // const Double_t ipCenter = 0.04;
  // const Double_t hpCenter = 20;
  // const Double_t ipScale = 0.55 - ipCenter;
  // const Double_t hpScale = 30 - hpCenter;

  for(int pol=0; pol<NUM_POL; pol++){

    TString name = "hPeakDirWrtNorth_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    TString title = " Peak direction w.r.t North vs. time for all min bias events ";
    title += "; Time; Direction wrt. north (Degrees)";

    hPeakDirWrtNorth[pol] = new TH2D(name, title,
				     numTimeBins,
				     firstRealTime, lastRealTime+1,
				     numBinsPhi,
				     minDirWrtNorth, maxDirWrtNorth);

    name = "hPeakTheta_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title += "Peak direction #theta vs. time for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; Time; #theta (Degrees)";
	
    hPeakTheta[pol] = new TH2D(name, title,
			       numTimeBins, firstRealTime, lastRealTime+1,
			       numBinsTheta, minTheta, maxTheta);

    name = "hImagePeakHilbertPeak_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "Image peak vs. Hilbert Peak for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; Image peak (no units); Hilbert peak (mV); Number of events";

    hImagePeakHilbertPeak[pol] = new TH2D(name, title,
					  numImagePeakBins, 0, 1,
					  numHilbertPeakBins, 0, maxHilbertPeak);


    name = "hImagePeakPhi_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "Image peak vs. #Phi for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #Phi (Degrees); Image peak (no units); Number of events";

    hImagePeakPhi[pol] = new TH2D(name, title,
				  numBinsPhi, 0, 360,
				  numImagePeakBins, 0, 1);

    name = "hImagePeakPhi0_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "Image peak vs. #Phi for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #Phi (Degrees); Image peak (no units); Number of events";

    hImagePeakPhi0[pol] = new TH2D(name, title,
				   NUM_BINS_PHI*NUM_PHI, 0, 360,
				   numImagePeakBins, 0, 1);
	


    name = "hImagePeakDeltaPhi_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "Image peak vs. #Phi for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #delta#Phi (Degrees); Image peak (no units); Number of events";

    hImagePeakDeltaPhi[pol] = new TH2D(name, title,
				       1024, -6, 6,
				       numImagePeakBins, 0, 1);
	

    name = "hPhi0DeltaPhi_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "#delta#Phi vs. #Phi_{0} for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #Phi_{0} (Degrees); #delta#Phi (Degrees); Number of events";

    hPhi0DeltaPhi[pol] = new TH2D(name, title,
				  NUM_BINS_PHI*NUM_PHI, 0, 360,
				  1024, -6, 6);


    name = "hPhi0DeltaPhi2_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "#delta#Phi vs. #Phi_{0} for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #Phi_{0} (Degrees); #delta#Phi (Degrees); Number of events";

    hPhi0DeltaPhi2[pol] = new TH2D(name, title,
				   NUM_BINS_PHI*NUM_PHI, 0, 360,
				   1024, -6, 6);
	
    name = "hPhi0Theta0_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "#Phi_{0} vs #Theta_{0} for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #Phi_{0} (Degrees); #theta_{0} (Degrees); Number of events";

    hPhi0Theta0[pol] = new TH2D(name, title,
				NUM_BINS_PHI*NUM_PHI, 0, 360,
				NUM_BINS_THETA, -75, 75);

    name = "hImagePeakTime_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "Image peak vs. time for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; Time; Image peak (no units); Number of events";

    hImagePeakTime[pol] = new TH2D(name, title,
				   numTimeBins, firstRealTime, lastRealTime+1,
				   numImagePeakBins, 0, 1);
    name = "hHilbertPeakTime_";
    // name += trig == 0 ? "" : "goodTime_";
    // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
    name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

    title = "Hilbert peak vs. time for all min bias events ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; Time; Hilbert peak (mV); Number of events";

    hHilbertPeakTime[pol] = new TH2D(name, title,
				     numTimeBins, firstRealTime, lastRealTime+1,
				     numHilbertPeakBins, 0, maxHilbertPeak);
  }
  
  TH2D* hHeading = new TH2D("hHeading", "ANITA heading vs. realTime; realTime; heading (Degrees)",
			    1024, firstRealTime, lastRealTime+1,
			    numBinsPhi, 0, 360);

  TH1D* hDeltaSolarPhiDeg[NUM_POL];
  TH1D* hDeltaSolarThetaDeg[NUM_POL];
  TH2D* hDeltaSolarPhiDegZoom[NUM_POL];
  TH2D* hDeltaSolarThetaDegZoom[NUM_POL];
  TH1D* hThetaDeg[NUM_POL];

  TH2D* hDeltaSolarThetaDegZoomVsTheta[NUM_POL];
  TH2D* hDeltaSolarThetaDegZoomVsPhi[NUM_POL];
  TH2D* hDeltaSolarThetaDegZoomVsTimeOfDay[NUM_POL];
  TH2D* hDeltaSolarPhiDegZoomVsTheta[NUM_POL];
  TH2D* hDeltaSolarPhiDegZoomVsPhi[NUM_POL];
  TH2D* hDeltaSolarPhiDegZoomVsTimeOfDay[NUM_POL];
  TH2D* hImagePeakSolarTheta[NUM_POL];
  
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
    title += "; Measured #phi (Degrees); #delta#phi_{sun} (Degrees); Number of events / bin";

    hDeltaSolarPhiDegZoom[pol] = new TH2D(name, title,
					  numBinsPhi, 0, 360,
					  numBinsPhi, -10, 10);

    name = pol==AnitaPol::kHorizontal ? "hDeltaSolarThetaDegZoomHPol" : "hDeltaSolarThetaDegZoomVPol";
    title = "#delta#theta_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; Measured #phi (Degrees); #delta#theta_{sun} (Degrees); Number of events / bin";

    hDeltaSolarThetaDegZoom[pol] = new TH2D(name, title,
					    numBinsPhi, -0, 360,
					    numBinsPhi, -10, 10);

    name = pol==AnitaPol::kHorizontal ? "hThetaDegHPol" : "hThetaDegVPol";
    title = "Min bias measured #theta ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "#theta (Degrees); Number of events / bin";

    hThetaDeg[pol] = new TH1D(name, title,
			      numBinsPhi, -90, 90);


    name = "hDeltaSolarThetaDegZoomVsTheta";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "#delta#theta_{sun} vs. #theta_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #theta (Degrees); #delta#theta_{sun}; Number of events / bin";

    hDeltaSolarThetaDegZoomVsTheta[pol] = new TH2D(name, title,
						   numBinsPhi, -90, 90,
						   numBinsPhi, -10, 10);


    name = "hDeltaSolarThetaDegZoomVsPhi";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "#delta#theta_{sun} vs. #phi_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #phi (Degrees); #delta#theta_{sun}; Number of events / bin";

    hDeltaSolarThetaDegZoomVsPhi[pol] = new TH2D(name, title,
						 numBinsPhi, 0, 360,
						 numBinsPhi, -10, 10);

    name = "hDeltaSolarThetaDegZoomVsTimeOfDay";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "#delta#theta_{sun} vs. timeOfDay ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; Time of day; #delta#theta_{sun}; Number of events / bin";

    hDeltaSolarThetaDegZoomVsTimeOfDay[pol] = new TH2D(name, title,
						       1024, 0, 24*60*60,
						       numBinsPhi, -10, 10);

    name = "hDeltaSolarPhiDegZoomVsTheta";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "#delta#phi_{sun} vs. #theta_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #theta (Degrees); #delta#phi_{sun}; Number of events / bin";

    hDeltaSolarPhiDegZoomVsTheta[pol] = new TH2D(name, title,
						   numBinsPhi, -90, 90,
						   numBinsPhi, -10, 10);


    name = "hDeltaSolarPhiDegZoomVsPhi";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "#delta#phi_{sun} vs. #phi_{sun} ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; #phi (Degrees); #delta#phi_{sun}; Number of events / bin";

    hDeltaSolarPhiDegZoomVsPhi[pol] = new TH2D(name, title,
						 numBinsPhi, 0, 360,
						 numBinsPhi, -10, 10);

    name = "hDeltaSolarPhiDegZoomVsTimeOfDay";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "#delta#phi_{sun} vs. timeOfDay ";
    title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title += "; Time of day; #delta#phi_{sun}; Number of events / bin";

    hDeltaSolarPhiDegZoomVsTimeOfDay[pol] = new TH2D(name, title,
						       1024, 0, 24*60*60,
						       numBinsPhi, -10, 10);    

    name = "hImagePeakSolarTheta";
    name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
    title = "Image peak vs. time of day ";
    title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
    title += "; #theta_{sun} (Degrees); Image peak (no units)";
    
    hImagePeakSolarTheta[pol] = new TH2D(name, title,
					 numBinsTheta, minTheta, maxTheta,
					 numImagePeakBins, 0, 1);

  }

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
    // std::cerr << headEntry << "\t" << std::endl;
    if(headEntry < 0){
      std::cerr << "Now what!?" << std::endl;      
    }
    else{
      headChain->GetEntry(headEntry);
      gpsChain->GetEntry(headEntry);

      Double_t solarPhiDeg, solarThetaDeg;

      // usefulPat.getSunPosition(solarPhiDeg, solarThetaDeg);
      solarPhiDeg = eventSummary->sun.phi;
      solarThetaDeg = eventSummary->sun.theta;
    
      for(Int_t polInd=0; polInd < NUM_POL; polInd++){

	Double_t imagePeak = eventSummary->peak[polInd][1].value;
	Double_t hilbertPeak = eventSummary->coherent[polInd][1].peakHilbert;

	Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;
	Double_t recoPhiDeg0 = eventSummary->peak[polInd][0].phi;

	if(recoPhiDeg < 0) recoPhiDeg += 360;
	else if(recoPhiDeg >= 360) recoPhiDeg -= 360;

	if(recoPhiDeg0 < 0) recoPhiDeg0 += 360;
	else if(recoPhiDeg0 >= 360) recoPhiDeg0 -= 360;


	if(recoPhiDeg0 < 0 || recoPhiDeg < 0){
	  std::cerr << recoPhiDeg0 << "\t" << recoPhiDeg <<"\t" << std::endl;
	}

	Double_t recoThetaDeg = eventSummary->peak[polInd][1].theta;
	Double_t recoThetaDeg0 = eventSummary->peak[polInd][0].theta;

	// std::cerr << recoPhiDeg << "\t" << recoThetaDeg << "\t"
	// 	  << recoThetaDeg0 << "\t" << recoPhiDeg << std::endl;

	Double_t directionWrtNorth = RootTools::getDeltaAngleDeg(pat->heading, recoPhiDeg);

	// std::cerr << "\t" << pat->heading << "\t" << recoPhiDeg << "\t" << directionWrtNorth << std::endl;

	directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
	directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;

	// std::cerr << "\t" << header->realTime - firstRealTime << "\t" << lastRealTime - header->realTime << "\t" << (header->realTime >= firstRealTime) << "\t" << (header->realTime < lastRealTime) << std::endl;	
	

	Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(solarPhiDeg, recoPhiDeg);
	Double_t deltaSolarThetaDeg = solarThetaDeg - recoThetaDeg;


	if(TMath::Abs(directionWrtNorth) > maxDirWrtNorth){
	  std::cerr << recoPhiDeg << "\t" << pat->heading << std::endl;
	}

	Int_t sun = TMath::Abs(deltaSolarPhiDeg) > deltaSolarPhiDegCut ? 0 : 1;

	hDeltaSolarPhiDeg[polInd]->Fill(deltaSolarPhiDeg);
	hDeltaSolarThetaDeg[polInd]->Fill(deltaSolarThetaDeg);
	hDeltaSolarPhiDegZoom[polInd]->Fill(solarPhiDeg, deltaSolarPhiDeg);
	hDeltaSolarThetaDegZoom[polInd]->Fill(solarPhiDeg, deltaSolarThetaDeg);

	hThetaDeg[polInd]->Fill(recoThetaDeg);

	if(sun){
	  hDeltaSolarThetaDegZoomVsTheta[polInd]->Fill(solarThetaDeg, deltaSolarThetaDeg);	
	  hDeltaSolarThetaDegZoomVsPhi[polInd]->Fill(solarPhiDeg, deltaSolarThetaDeg);
	  hDeltaSolarThetaDegZoomVsTimeOfDay[polInd]->Fill(pat->timeOfDay/1000, deltaSolarThetaDeg);
	  hDeltaSolarPhiDegZoomVsTheta[polInd]->Fill(solarThetaDeg, deltaSolarPhiDeg);
	  hDeltaSolarPhiDegZoomVsPhi[polInd]->Fill(solarPhiDeg, deltaSolarPhiDeg);
	  hDeltaSolarPhiDegZoomVsTimeOfDay[polInd]->Fill(pat->timeOfDay/1000, deltaSolarPhiDeg);
	}

	double deltaPhiCoarseZoom = RootTools::getDeltaAngleDeg(recoPhiDeg0, recoPhiDeg);

	hPeakDirWrtNorth[polInd]->Fill(header->realTime,
				       directionWrtNorth);

	hPeakTheta[polInd]->Fill(header->realTime,
				 recoThetaDeg);

	hImagePeakHilbertPeak[polInd]->Fill(imagePeak,
					    hilbertPeak);
	hImagePeakPhi[polInd]->Fill(recoPhiDeg, imagePeak);
	hImagePeakPhi0[polInd]->Fill(recoPhiDeg0,
				     eventSummary->peak[polInd][0].value);
	hImagePeakDeltaPhi[polInd]->Fill(deltaPhiCoarseZoom,
					 imagePeak);
	hPhi0DeltaPhi[polInd]->Fill(recoPhiDeg0, deltaPhiCoarseZoom);

	if(deltaPhiCoarseZoom > -4.94 && deltaPhiCoarseZoom < 4.99){
	  hPhi0DeltaPhi2[polInd]->Fill(recoPhiDeg0, deltaPhiCoarseZoom);
	}
	
	hPhi0DeltaPhi[polInd]->Fill(recoPhiDeg0, deltaPhiCoarseZoom);
      
	hPhi0Theta0[polInd]->Fill(recoPhiDeg0,
				  recoThetaDeg0+0.0001);

	hImagePeakTime[polInd]->Fill(header->realTime,
				     imagePeak);

	hHilbertPeakTime[polInd]->Fill(header->realTime,
				       hilbertPeak);

	hImagePeakSolarTheta[polInd]->Fill(solarThetaDeg, imagePeak);
      }

      hHeading->Fill(header->realTime, pat->heading);
      if(pat->heading < 0 || pat->heading >= 360){
	std::cerr << "bad heading ? " << pat->heading << "\t" << header->run << "\t" << header->eventNumber << std::endl;
      }
    }
    // p++;
    p.inc(entry, maxEntry);
  }

  // grSunPhiDeg->SetName("grSunPhiDeg");
  // grSunPhiDeg->SetTitle("Predicted sun #phi (Degrees)");
  // grSunThetaDeg->SetName("grSunThetaDeg");
  // grSunThetaDeg->SetTitle("Predicted sun #theta (Degrees)");
  // grSunPhiDeg->Write();
  // grSunThetaDeg->Write();
  
  outFile->Write();
  outFile->Close();

  return 0;
}

