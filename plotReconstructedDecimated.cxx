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

  if(!(argc==2)){
    std::cerr << "Usage 1: " << argv[0] << " [cutStep]" << std::endl;
    return 1;
  }

  Int_t firstRun = 130; //atoi(argv[1]);
  Int_t lastRun = 434;
  const bool useTimeCut = false; //true;

  std::cout << firstRun << "\t" << lastRun << std::endl;
  const int cutStep = atoi(argv[1]);
  std::cout << "cutStep = " << cutStep << std::endl;


  // quite complicated, let's try and assign the cut values based on a passed variable


  // for image peak / hilbert peak histograms
  const double deltaSolarPhiDegClose = 20;
  const double deltaSolarThetaClose = 5;



  if(useTimeCut==true){
    firstRun = 202; //atoi(argv[1]);
    lastRun = 251;
  }
  const int numGoodTimes = 1;
  UInt_t goodTimesStart[numGoodTimes] = {1419220000};
  UInt_t goodTimesEnd[numGoodTimes] = {1419500000};

  const double findHilbert = 0; //120; //50;
  const double findImage = 0; //0.074; //0.1; //0.06;

  // const double cutHilbert = 100; //50;
  // const double cutImage = 0.1; //0.06;
  const double cutHilbert = 0; //100; //50;
  const double cutImage = 0; //0.088; //0; //0.1; //0.06;


  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  TChain* dataQualityChain = new TChain("dataQualityTree");


  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);

    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("filter260-370-400-762-5peaks/reconstructDecimatedPlots_%d_*.root", run);
    eventSummaryChain->Add(fileName);

    fileName = TString::Format("finalDataQuality/slimTreePlots_%d_*.root", run);
    // fileName = TString::Format("finalDataQuality/makeDecimatedDataQualityTreesPlots_%d_*.root", run);
    // fileName = TString::Format("filter260-370-400-762/makeDecimatedDataQualityTreesPlots_%d_*.root", run);
    dataQualityChain->Add(fileName);

  }
  // return 1;
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

  Int_t goodTimeRun = -1;

  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  // Double_t peakToPeak[NUM_POL][NUM_SEAVEYS];
  // dataQualityChain->SetBranchAddress("peakToPeak", peakToPeak);
  UInt_t eventNumberDQ;
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumberDQ);
  // Double_t maxVolts[NUM_POL][NUM_SEAVEYS];
  // dataQualityChain->SetBranchAddress("maxVolts", maxVolts);
  // Double_t minVolts[NUM_POL][NUM_SEAVEYS];
  // dataQualityChain->SetBranchAddress("minVolts", minVolts);

  Double_t maxPeakToPeakRatio[NUM_POL];
  dataQualityChain->SetBranchAddress("maxPeakToPeakRatio", maxPeakToPeakRatio);
  Double_t theMaxVolts[NUM_POL];
  dataQualityChain->SetBranchAddress("theMaxVolts", theMaxVolts);
  Double_t theMinVolts[NUM_POL];
  dataQualityChain->SetBranchAddress("theMinVolts",theMinVolts);


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
  headChain->GetEntry(headChain->GetEntries()-2);
  UInt_t lastRealTime = header->realTime;

  // std::cerr << firstRealTime << "\t" << lastRealTime << std::endl;

  Long64_t nEntries = eventSummaryChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  const int numBinsPhi = 420;
  const int numBinsTheta = numBinsPhi/2;
  const int numTimeBins = 1024*2;
  const double maxDirWrtNorth = 180;
  const double minDirWrtNorth = -180;
  const double maxTheta = 90;
  const double minTheta = -90;

  const Int_t numImagePeakBins = 1024;
  const Int_t numHilbertPeakBins = 1024;
  const Double_t maxHilbertPeak = 2048;


  // const int numPeaks = 1;
  const int numPeaks = 1;

  TH1D* hMaxVolts = new TH1D("hMaxVolts", "Maximum amplitude in event", 1024, 0, 5000);
  TH1D* hMinVolts = new TH1D("hMinVolts", "Minimum amplitude in event", 1024, -5000, 0);
  TH1D* hAbsSumMaxMin = new TH1D("hAbsSumMaxMin", "hAbsSumMaxMin", 1024, 0, 5000);


  TH2D* hPeakHeading = new TH2D("hPeakHeading",
				"Peak Heading; Time; Peak Heading (Degrees)",
				numTimeBins,
				firstRealTime, lastRealTime+1,
				numBinsPhi,
				minDirWrtNorth, maxDirWrtNorth);

  TH1D* hGoodPeak = new TH1D("hGoodPeak",
			     "Which peak is good; Peak Number; Number of events",
			     numPeaks+1, -1, numPeaks);
  hGoodPeak->GetXaxis()->SetBinLabel(1, "None");
  for(int peakInd=0; peakInd < numPeaks; peakInd++){
    hGoodPeak->GetXaxis()->SetBinLabel(peakInd+2, TString::Format("%d", peakInd+1));
  }

  TH2D* hPeakElevation = new TH2D("hPeakElevation",
				  "Peak Elevation; Time; Peak Elevation (Degrees)",
				  numTimeBins, firstRealTime, lastRealTime+1,
				  numBinsTheta, minTheta, maxTheta);

  TProfile2D* pPeakHeading = new TProfile2D("pPeakHeading",
					    "Profile of Image Peak vs. Peak Heading vs. Time; Time; Peak Heading (Degrees); Image peak (no units)",
					    numTimeBins,
					    firstRealTime, lastRealTime+1,
					    numBinsPhi,
					    minDirWrtNorth, maxDirWrtNorth);

  TProfile2D* pPeakElevation = new TProfile2D("pPeakElevation",
					      "Profile of Image Peak vs. Peak Elevation vs. Time; Time; Peak Elevation (Degrees); Image peak (no units)",
					      numTimeBins, firstRealTime, lastRealTime+1,
					      numBinsTheta, minTheta, maxTheta);


  TH1D* hThetaDeg = new TH1D("hThetaDeg",
			     "Peak Elevation; #theta_{peak} (Degrees); Events per bin",
			     numBinsPhi, -90, 90);

  TH1D* hDeltaPhiSect = new TH1D("hDeltaPhiSect",
				 "Number of #Phi-sectors between image peak and nearest L3 trigger; #delta#Phi-sector; Events per bin",
				 NUM_PHI/2+1, 0, NUM_PHI/2 + 1); // last bin == no L3 triggers

  TH1D* hMaxBottomToTopPeakToPeakRatio = new TH1D("hBottomToTopPeakToPeakRatio",
						  "Maximum peak-to-peak ratio between top and bottom rings; Ratio of peak-to-peak value (no units); Number of events",
						  1024, 0, 10);

  TH2D* hImagePeak1Peak2 = new TH2D("hImagePeak1Peak2",
				    "Image peak 1 vs. Image Peak 2 ",
				    numImagePeakBins, 0, 1,
				    numImagePeakBins, 0, 1);

  TH2D* hImagePeakHilbertPeak = new TH2D("hImagePeakHilbertPeak",
					 "Image peak vs. Hilbert Peak ",
					 numImagePeakBins, 0, 1,
					 numHilbertPeakBins, 0, maxHilbertPeak);

  TH2D* hImagePeakHilbertPeakSunPhi = new TH2D("hImagePeakHilbertPeakSunPhi",
					       "Image peak vs. Hilbert Peak ",
					       numImagePeakBins, 0, 1,
					       numHilbertPeakBins, 0, maxHilbertPeak);

  TH2D* hImagePeakHilbertPeakSunPhiTheta = new TH2D("hImagePeakHilbertPeakSunPhiTheta",
						    "Image peak vs. Hilbert Peak ",
						    numImagePeakBins, 0, 1,
						    numHilbertPeakBins, 0, maxHilbertPeak);


  TH2D* hImagePeakPhi = new TH2D("hImagePeakPhi",
				 "Image peak vs. #Phi ",
				 numBinsPhi, 0, 360,
				 numImagePeakBins, 0, 1);

  TH2D* hImagePeakTime = new TH2D("hImagePeakTime",
				  "Image peak vs. Time ",
				  numTimeBins, firstRealTime, lastRealTime+1,
				  numImagePeakBins, 0, 1);

  TH2D* hHilbertPeakTime = new TH2D("hHilbertPeakTime",
				    "Hilbert peak vs. Time ",
				    numTimeBins, firstRealTime, lastRealTime+1,
				    numHilbertPeakBins, 0, maxHilbertPeak);

  TH2D* hHeading = new TH2D("hHeading",
			    "ANITA heading vs. realTime; realTime; heading (Degrees)",
			    1024, firstRealTime, lastRealTime+1,
			    numBinsPhi, 0, 360);

  TH1D* hDeltaSolarPhiDeg[numPeaks];
  TH1D* hDeltaSolarThetaDeg[numPeaks];

  for(int peakInd=0; peakInd < numPeaks; peakInd++){
    TString name = TString::Format("hDeltaSolarPhiDeg%d", peakInd);
    hDeltaSolarPhiDeg[peakInd] = new TH1D(name,
					  "#delta#phi_{sun}; #delta#phi_{sun} (Degrees); Events per bin",
					  numBinsPhi, -180, 180);

    name = TString::Format("hDeltaSolarThetaDeg%d", peakInd);
    hDeltaSolarThetaDeg[peakInd] = new TH1D(name,
					    "#deltaElevation_{sun}; #deltaElevation_{sun} (Degrees); Events per bin",
					    numBinsPhi, -180, 180);


  }



  TH2D* hDeltaSolarThetaDegVsTheta = new TH2D("hDeltaSolarThetaDegVsTheta",
					      "#delta#theta_{sun} vs. Elevation_{sun}; Elevation (Degrees); #delta#theta (Degrees)",
					      numBinsPhi, -90, 90,
					      numBinsPhi, -10, 10);


  TH2D* hDeltaSolarThetaDegVsPhi = new TH2D("hDeltaSolarThetaDegVsPhi",
					    "#deltaElevation_{sun} vs. #phi_{sun}; ",
					    numBinsPhi, 0, 360,
					    numBinsPhi, -10, 10);


  TH2D* hDeltaSolarThetaDegVsTimeOfDay = new TH2D("hDeltaSolarThetaDegVsTimeOfDay",
						  "#deltaElevation_{sun} vs. TimeOfDay; Time of day; #deltaElevation_{sun}; Events per bin",
						  1024, 0, 24*60*60,
						  numBinsPhi, -10, 10);


  TH2D* hDeltaSolarPhiDegVsTheta = new TH2D("hDeltaSolarPhiDegVsTheta",
					    "#delta#phi_{sun} vs. Elevation_{sun} ",
					    numBinsPhi, -90, 90,
					    numBinsPhi, -10, 10);



  TH2D* hDeltaSolarPhiDegVsPhi = new TH2D("hDeltaSolarPhiDegVsPhi",
					  "#delta#phi_{sun} vs. #phi_{sun} ",
					  numBinsPhi, 0, 360,
					  numBinsPhi, -10, 10);


  TH2D* hDeltaSolarPhiDegVsTimeOfDay = new TH2D("hDeltaSolarPhiDegVsTimeOfDay",
						"#delta#phi_{sun} vs. TimeOfDay; Time of day; #delta#phi_{sun}; Events per bin",
						1024, 0, 24*60*60,
						numBinsPhi, -10, 10);


  TH2D* hDeltaSolarPhiDegVsDeltaSolarThetaDeg = new TH2D("hDeltaSolarPhiDegVsDeltaSolarThetaDeg",
						       "#delta#theta_{sun} vs. #delta#phi_{sun}; #delta#phi_{sun} (Degrees); #delta#theta_{sun} (Degrees); Events per bin",
						       numBinsPhi, -180, 180,
						       2*numBinsTheta, -180, 180);
  TProfile2D* pImagePeakVsDeltaSolarPhiDegVsDeltaSolarThetaDeg = new TProfile2D("pImagePeakVsDeltaSolarPhiDegVsDeltaSolarThetaDeg",
									      "Profile of Image Peak vs. #delta#theta_{sun} vs. #delta#phi_{sun}; #delta#phi_{sun} (Degrees); #delta#theta_{sun} (Degrees); Image peak (no units)",
									      numBinsPhi, -180, 180,
									      2*numBinsTheta, -180, 180);


  TProfile2D* pImagePeakVsDeltaSolarPhiDegVsSunTheta = new TProfile2D("pImagePeakVsDeltaSolarPhiDegVsSunTheta",
								      "Profile of Image Peak vs. #delta#theta_{sun} vs. #theta_{sun}; #theta_{sun} (Degrees); #delta#theta_{sun} (Degrees); Image peak (no units)",
								      numBinsPhi, 10, 40,
								      2*numBinsTheta, -180, 180);

  TH2D* hPeakRatio = new TH2D("hPeakRatio", "Map peaks; P1; P2/P1", 1024, 0, 1, 1024, 0, 1);

  std::cerr << "building index" << std::endl;
  headChain->BuildIndex("eventNumber");
  // dataQualityChain->BuildIndex("eventNumber");
  std::cerr << "done" << std::endl;


  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);

    Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber, 0);
    if(headEntry < 0){
      std::cerr << "Now what!?\t" << headEntry << "\t" << eventSummary->eventNumber << std::endl;
    }
    else{
      headChain->GetEntry(headEntry);

      AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
      if(eventSummary->peak[AnitaPol::kHorizontal][0].value > eventSummary->peak[AnitaPol::kVertical][0].value){
	pol = AnitaPol::kHorizontal;
      }


      // timing selection...
      bool isGoodTime = false;
      for(int i=0; i < numGoodTimes; i++){
	if(header->realTime >= goodTimesStart[i] && header->realTime < goodTimesEnd[i]){

	  if(goodTimeRun==-1){
	    std::cout << "First good time is in run " << header->run << std::endl;
	  }
	  goodTimeRun = header->run;
	  isGoodTime = true;
	}
      }
      if(useTimeCut==true && isGoodTime==false){
	p.inc(entry, maxEntry);
	continue;
      }



      dataQualityChain->GetEntry(entry);
      if(eventSummary->eventNumber != eventNumberDQ){
	std::cerr << "???" << eventSummary->eventNumber << "\t" << eventNumberDQ << std::endl;
      }


      AnalysisCuts::Status_t surfSaturation;
      Double_t absSumMaxMin = theMaxVolts[pol] + theMinVolts[pol];
      surfSaturation = AnalysisCuts::applySurfSaturationCut(theMaxVolts[pol], theMinVolts[pol], absSumMaxMin);

      if(cutStep >= 1 && surfSaturation==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }
      hMaxVolts->Fill(theMaxVolts[pol]);
      hMinVolts->Fill(theMinVolts[pol]);
      hAbsSumMaxMin->Fill(absSumMaxMin);

      AnalysisCuts::Status_t selfTriggeredBlastCut;
      selfTriggeredBlastCut = AnalysisCuts::applyBottomToTopRingPeakToPeakRatioCut(maxPeakToPeakRatio[pol]);
      if(cutStep >= 1 && selfTriggeredBlastCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }
      hMaxBottomToTopPeakToPeakRatio->Fill(maxPeakToPeakRatio[pol]);


      if(maxPeakToPeakRatio[pol] < 0.9){
	std::cerr << "low ratio: " << header->run << "\t" << header->eventNumber << std::endl;
      }



      // Get event info
      const int peakInd = 0;

      Double_t recoPhiDeg = eventSummary->peak[pol][peakInd].phi;
      recoPhiDeg += recoPhiDeg < 0 ? DEGREES_IN_CIRCLE : 0;
      Double_t recoThetaDeg = eventSummary->peak[pol][peakInd].theta;
      Double_t imagePeak = eventSummary->peak[pol][peakInd].value;
      Double_t hilbertPeak = eventSummary->coherent[pol][peakInd].peakHilbert;


      // CUT FLOW
      // Step 2: cut phi-sector angle triggers
      Int_t deltaPhiSect = NUM_PHI/2;

      AnalysisCuts::Status_t l3TriggerCut;
      l3TriggerCut = AnalysisCuts::L3TriggerDirectionCut(pol, header, recoPhiDeg, deltaPhiSect);
      if(cutStep >= 2 && l3TriggerCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }

      hDeltaPhiSect->Fill(deltaPhiSect);


      gpsChain->GetEntry(headEntry);




      // CUT FLOW
      // Step 3: cut phi-direction relative to sun
      Double_t solarPhiDeg = eventSummary->sun.phi;
      Double_t solarThetaDeg = -1*eventSummary->sun.theta;

      Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(recoPhiDeg, solarPhiDeg);
      solarPhiDeg = solarPhiDeg < 0 ? solarPhiDeg + 360 : solarPhiDeg;
      Double_t deltaSolarThetaDeg = recoThetaDeg - solarThetaDeg;

      // std::cerr << recoThetaDeg << "\t" << solarThetaDeg << "\t" << deltaSolarThetaDeg << "\t";
      // std::cerr << recoPhiDeg << "\t" << solarPhiDeg << "\t" << deltaSolarPhiDeg << std::endl;

      // abs sun position relative to other variables
      pImagePeakVsDeltaSolarPhiDegVsSunTheta->Fill(solarThetaDeg,
						   deltaSolarThetaDeg,
						   imagePeak);

      if(TMath::Abs(deltaSolarPhiDeg) < deltaSolarPhiDegClose &&
	 TMath::Abs(deltaSolarThetaDeg) < deltaSolarThetaClose){
	hImagePeakHilbertPeakSunPhiTheta->Fill(imagePeak, hilbertPeak);
      }
      else if(TMath::Abs(deltaSolarPhiDeg) < deltaSolarPhiDegClose){
	hImagePeakHilbertPeakSunPhi->Fill(imagePeak, hilbertPeak);
      }


      AnalysisCuts::Status_t sunCut;
      sunCut = AnalysisCuts::applySunPointingCut(deltaSolarPhiDeg);
      if(cutStep>=3 && sunCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }

      hDeltaSolarThetaDegVsTheta->Fill(solarThetaDeg, deltaSolarThetaDeg);
      hDeltaSolarThetaDegVsPhi->Fill(solarPhiDeg, deltaSolarThetaDeg);
      hDeltaSolarThetaDegVsTimeOfDay->Fill(pat->timeOfDay/1000, deltaSolarThetaDeg);
      hDeltaSolarPhiDegVsTheta->Fill(solarThetaDeg, deltaSolarPhiDeg);
      hDeltaSolarPhiDegVsPhi->Fill(solarPhiDeg, deltaSolarPhiDeg);
      hDeltaSolarPhiDegVsTimeOfDay->Fill(pat->timeOfDay/1000, deltaSolarPhiDeg);
      // relative sun position
      hDeltaSolarPhiDeg[peakInd]->Fill(deltaSolarPhiDeg);
      hDeltaSolarThetaDeg[peakInd]->Fill(deltaSolarThetaDeg);
      hDeltaSolarPhiDegVsDeltaSolarThetaDeg->Fill(deltaSolarPhiDeg, deltaSolarThetaDeg);
      pImagePeakVsDeltaSolarPhiDegVsDeltaSolarThetaDeg->Fill(deltaSolarPhiDeg,
							     deltaSolarThetaDeg,
							     imagePeak);





      AnalysisCuts::Status_t thermalCut;
      Double_t fisher;
      thermalCut = AnalysisCuts::applyThermalBackgroundCut(imagePeak, hilbertPeak, fisher);
      if(cutStep>=4 && thermalCut==AnalysisCuts::kFail){
	p.inc(entry, maxEntry);
	continue;
      }

      // if(imagePeak < 0.075 || (imagePeak < 0.3 && hilbertPeak > 200)){
      // 	std::cout << std::endl << header->run << "\t" << header->eventNumber << "\t" << imagePeak << "\t" << hilbertPeak << "\t" << fisher << std::endl;
      // }










      if(imagePeak < cutImage || hilbertPeak < cutHilbert){
	p.inc(entry, maxEntry);
	continue;
      }

      if((findImage > 0 && imagePeak >= findImage) || (findHilbert > 0 && hilbertPeak >= findHilbert)){
      	std::cerr << header->run << "\t" << header->eventNumber << "\t" << imagePeak << "\t" << hilbertPeak << std::endl;;
      }


      if(recoPhiDeg < 0) recoPhiDeg += 360;
      else if(recoPhiDeg >= 360) recoPhiDeg -= 360;


      if(recoPhiDeg < 0){
	std::cerr << recoPhiDeg <<"\t" << std::endl;
      }



      Double_t directionWrtNorth = RootTools::getDeltaAngleDeg(pat->heading, recoPhiDeg);

      directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
      directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;



      if(TMath::Abs(directionWrtNorth) > maxDirWrtNorth){
	std::cerr << recoPhiDeg << "\t" << pat->heading << std::endl;
      }

      hImagePeak1Peak2->Fill(imagePeak, eventSummary->peak[pol][peakInd+1].value);


      hThetaDeg->Fill(recoThetaDeg);

      // std::cerr << solarThetaDeg << "\t" << solarPhiDeg << "\t"
      // 	  << deltaSolarThetaDeg << "\t" << deltaSolarPhiDeg << "\t"
      // 	  << std::endl;



      // event info (peak direction, and value)
      hPeakHeading->Fill(header->realTime,
			 directionWrtNorth);
      hPeakElevation->Fill(header->realTime,
			   recoThetaDeg);
      pPeakHeading->Fill(header->realTime,
			 directionWrtNorth,
			 imagePeak);
      pPeakElevation->Fill(header->realTime,
			   recoThetaDeg,
			   imagePeak);

      hImagePeakHilbertPeak->Fill(imagePeak,
				  hilbertPeak);

      hImagePeakTime->Fill(header->realTime,
			   imagePeak);
      hHilbertPeakTime->Fill(header->realTime,
			     hilbertPeak);
      Double_t p1 = eventSummary->peak[pol][0].value;
      Double_t p2 = eventSummary->peak[pol][1].value;

      hPeakRatio->Fill(p1, p2/p1);

      // reconstruction data quality
      hImagePeakPhi->Fill(recoPhiDeg, imagePeak);
      hHeading->Fill(header->realTime, pat->heading);
    }

    p.inc(entry, maxEntry);

  }
  std::cout << "Last good time run is " << goodTimeRun << std::endl;

  outFile->Write();
  outFile->Close();

  return 0;
}
