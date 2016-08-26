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

  
  // const double cutHilbert = 100; //50;
  // const double cutImage = 0.1; //0.06;
  const double cutHilbert = 0; //100; //50;
  const double cutImage = 0; //0.088; //0; //0.1; //0.06;  
  
  const double ratioCut = 2.8;

  const bool useTimeCut = false;
  const int numGoodTimes = 1;
  UInt_t goodTimesStart[numGoodTimes] = {1419100000};
  UInt_t goodTimesEnd[numGoodTimes] = {1419500000};

  
  TChain* headChain = new TChain("headTree");
  TChain* gpsChain = new TChain("adu5PatTree");
  TChain* eventSummaryChain = new TChain("eventSummaryTree");
  TChain* dataQualityChain = new TChain("dataQualityTree");    

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    headChain->Add(fileName);

    fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    gpsChain->Add(fileName);

    fileName = TString::Format("filter260-370-400-762/reconstructWaisPlots_%d_*.root", run);
    eventSummaryChain->Add(fileName);

    fileName = TString::Format("filter260-370-400-762/makeWaisDataQualityTreesPlots_%d_*.root", run);
    dataQualityChain->Add(fileName);    

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
  if(dataQualityChain->GetEntries()==0){
    std::cerr << "Unable to find data quality files!" << std::endl;
    return 1;
  }

  
  RawAnitaHeader* header = NULL;
  headChain->SetBranchAddress("header", &header);
  Adu5Pat* pat = NULL;
  gpsChain->SetBranchAddress("pat", &pat);
  Double_t peakToPeak[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("peakToPeak", peakToPeak);
  UInt_t eventNumberDQ;
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumberDQ);
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

  // const double deltaSolarPhiDegCut = 20; // degrees
  const double deltaSolarPhiDegCut = 10; // degrees

  const Int_t numImagePeakBins = 1024;
  const Int_t numHilbertPeakBins = 1024;
  const Double_t maxHilbertPeak = 2048;  

  TH2D* hPeakHeading = new TH2D("hPeakHeading",
				"Peak Heading; Time; Peak Heading (Degrees)",
				numTimeBins,
				firstRealTime, lastRealTime+1,
				numBinsPhi,
				minDirWrtNorth, maxDirWrtNorth);
	
  TH2D* hPeakElevation = new TH2D("hPeakElevation",
				  "Peak Heading; Time; Peak Heading (Degrees)",			      
				  numTimeBins, firstRealTime, lastRealTime+1,
				  numBinsTheta, minTheta, maxTheta);

  TProfile2D* pPeakHeading = new TProfile2D("pPeakHeading",
					    "Profile of Image Peak vs. Peak Heading vs. Time; Time; Peak Heading (Degrees); Image peak (no units)",
					    numTimeBins,
					    firstRealTime, lastRealTime+1,
					    numBinsPhi,
					    minDirWrtNorth, maxDirWrtNorth);

  TProfile2D* pPeakElevation = new TProfile2D("pPeakElevation",
					      "Profile of Image Peak vs. Peak Elevtaion vs. Time; Time; Peak Elevation (Degrees); Image peak (no units)",
					      numTimeBins, firstRealTime, lastRealTime+1,
					      numBinsTheta, minTheta, maxTheta);

  TH2D* hImagePeakHilbertPeak = new TH2D("hImagePeakHilbertPeak",
					 "Image peak vs. Hilbert Peak ",
					 numImagePeakBins, 0, 1,
					 numHilbertPeakBins, 0, maxHilbertPeak);
  

  TH2D* hImagePeakPhi = new TH2D("hImagePeakPhi",
				 "Image peak vs. #Phi ",
				 numBinsPhi, 0, 360,
				 numImagePeakBins, 0, 1);



  TH2D* hImagePeakPhi0 = new TH2D("hImagePeakPhi0",
				  "Image peak vs. #Phi ",
				  NUM_BINS_PHI*NUM_PHI, 0, 360,
				  numImagePeakBins, 0, 1);

  TH2D* hImagePeakDeltaPhi = new TH2D("hImagePeakDeltaPhi",
				      "Image peak vs. #Phi ",
				      1024, -6, 6,
				      numImagePeakBins, 0, 1);

  TH2D* hPhi0DeltaPhi = new TH2D("hPhi0DeltaPhi",
				 "#delta#Phi vs. #Phi_{0} ",
				 NUM_BINS_PHI*NUM_PHI, 0, 360,
				 1024, -6, 6);

  TH2D* hPhi0DeltaPhi2 = new TH2D("hPhi0DeltaPhi2",
				  "#delta#Phi vs. #Phi_{0} ",
				  NUM_BINS_PHI*NUM_PHI, 0, 360,
				  1024, -6, 6);

  TH2D* hPhi0Theta0 = new TH2D("hPhi0Theta0",
			       "#Phi_{0} vs Elevation_{0} ",
			       NUM_BINS_PHI*NUM_PHI, 0, 360,
			       NUM_BINS_THETA, -75, 75);

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


  TH1D* hDeltaSolarPhiDeg = new TH1D("hDeltaSolarPhiDeg",
				     "#delta#phi_{sun}; #delta#phi_{sun} (Degrees); Events per bin",
				     numBinsPhi, -180, 180);

  

  TH1D* hDeltaSolarThetaDeg = new TH1D("hDeltaSolarThetaDeg",
				       "#deltaElevation_{sun}; #deltaElevation_{sun} (Degrees); Events per bin",
				       numBinsPhi, -180, 180);

  TH2D* hDeltaSolarPhiDegZoom = new TH2D("hDeltaSolarPhiDegZoom",
					 "#delta#phi_{sun}; #delta#phi_{sun} (Degrees); Events per bin",
					 numBinsPhi, 0, 360,
					 numBinsPhi, -10, 10);


  TH2D* hDeltaSolarThetaDegZoom = new TH2D("hDeltaSolarThetaDegZoom",
					   "#deltaElevation_{sun}; #deltaElevation_{sun} (Degrees); Events per bin",
					   numBinsPhi, -0, 360,
					   numBinsPhi, -10, 10);

  TH1D* hThetaDeg = new TH1D("hThetaDeg",
			     "Peak Elevation; #theta_{peak} (Degrees); Events per bin",
			     numBinsPhi, -90, 90);



  TH2D* hDeltaSolarThetaDegZoomVsTheta = new TH2D("hDeltaSolarThetaDegZoomVsTheta",
						  "#delta#theta_{sun} vs. Elevation_{sun}; Elevation (Degrees); #delta#theta (Degrees)",
						  numBinsPhi, -90, 90,
						  numBinsPhi, -10, 10);


  TH2D* hDeltaSolarThetaDegZoomVsPhi = new TH2D("hDeltaSolarThetaDegZoomVsPhi",
						"#deltaElevation_{sun} vs. #phi_{sun}; ",
						numBinsPhi, 0, 360,
						numBinsPhi, -10, 10);


  TH2D* hDeltaSolarThetaDegZoomVsTimeOfDay = new TH2D("hDeltaSolarThetaDegZoomVsTimeOfDay",
						      "#deltaElevation_{sun} vs. TimeOfDay; Time of day; #deltaElevation_{sun}; Events per bin",
						      1024, 0, 24*60*60,
						      numBinsPhi, -10, 10);


  TH2D* hDeltaSolarPhiDegZoomVsTheta = new TH2D("hDeltaSolarPhiDegZoomVsTheta",
						"#delta#phi_{sun} vs. Elevation_{sun} ",
						numBinsPhi, -90, 90,
						numBinsPhi, -10, 10);



  TH2D* hDeltaSolarPhiDegZoomVsPhi = new TH2D("hDeltaSolarPhiDegZoomVsPhi",
					      "#delta#phi_{sun} vs. #phi_{sun} ",
					      numBinsPhi, 0, 360,
					      numBinsPhi, -10, 10);


  TH2D* hDeltaSolarPhiDegZoomVsTimeOfDay = new TH2D("hDeltaSolarPhiDegZoomVsTimeOfDay",
						    "#delta#phi_{sun} vs. TimeOfDay; Time of day; #delta#phi_{sun}; Events per bin",
						    1024, 0, 24*60*60,
						    numBinsPhi, -10, 10);


  TH1D* hDeltaPhiSect = new TH1D("hDeltaPhiSect",
				 "Peak distance to nearest #Phi-sector; #delta#Phi-sector; Events per bin",
				 NUM_PHI, -NUM_PHI/2, NUM_PHI/2);

  TH1D* hDeltaPhiSect0 = new TH1D("hDeltaPhiSect0",
				 "Peak distance to nearest #Phi-sector; #delta#Phi-sector; Events per bin",
				 NUM_PHI, -NUM_PHI/2, NUM_PHI/2);  


  TH2D* hDeltaSolarPhiDegVsDeltaSolarPhiDeg = new TH2D("hDeltaSolarPhiDegVsDeltaSolarPhiDeg",
						       "#delta#theta_{sun} vs. #delta#phi_{sun}; #delta#phi_{sun} (Degrees); #delta#theta_{sun} (Degrees); Events per bin",
						       numBinsPhi, -180, 180,
						       2*numBinsTheta, -180, 180);
  TProfile2D* hImagePeakVsDeltaSolarPhiDegVsDeltaSolarPhiDeg = new TProfile2D("hImagePeakVsDeltaSolarPhiDegVsDeltaSolarPhiDeg",
						       "Profile of Image Peak vs. #delta#theta_{sun} vs. #delta#phi_{sun}; #delta#phi_{sun} (Degrees); #delta#theta_{sun} (Degrees); Image peak (no units)",
						       numBinsPhi, -180, 180,
						       2*numBinsTheta, -180, 180);
  
  
  std::cerr << "building index" << std::endl;
  headChain->BuildIndex("eventNumber");
  std::cerr << "done" << std::endl;

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);
    dataQualityChain->GetEntry(entry);

    if(eventSummary->eventNumber != eventNumberDQ){
      std::cerr << "???" << eventSummary->eventNumber << "\t" << eventNumberDQ << std::endl;
    }

    AnitaPol::AnitaPol_t pol = AnitaPol::kVertical;
    if(eventSummary->peak[0][1].value > eventSummary->peak[1][1].value){
      pol = AnitaPol::kHorizontal;
    }
    
    Double_t maxRatio = 0;
    for(int phi=0; phi < NUM_PHI; phi++){
      if(pol==AnitaPol::kVertical && phi==7){
	continue;
      }
      Double_t ratio = peakToPeak[pol][phi+2*NUM_PHI]/peakToPeak[pol][phi];
      if(ratio > maxRatio){
	maxRatio = ratio;
      }
    }
    if(maxRatio > ratioCut){
      // std::cerr << eventNumberDQ << "\t" << maxRatio << std::endl;
      p.inc(entry, maxEntry);
      continue;
    }
    
    Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber, 0);
    if(headEntry < 0){
      std::cerr << "Now what!?\t" << headEntry << "\t" << eventSummary->eventNumber << std::endl;
    }
    else{
      headChain->GetEntry(headEntry);

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
      
      gpsChain->GetEntry(headEntry);

      Double_t solarPhiDeg, solarThetaDeg;

      solarPhiDeg = eventSummary->sun.phi;
      solarThetaDeg = -1*eventSummary->sun.theta;
    
      Double_t imagePeak = eventSummary->peak[pol][1].value;
      Double_t hilbertPeak = eventSummary->coherent[pol][1].peakHilbert;

      Double_t recoPhiDeg = eventSummary->peak[pol][1].phi;
      Double_t recoPhiDeg0 = eventSummary->peak[pol][0].phi;

      const double aftForeOffset = 45;
      const double phiRelativeToPhiSector0 = RootTools::getDeltaAngleDeg(recoPhiDeg, -aftForeOffset);
      Int_t phiSectorOfPeak = TMath::Nint(phiRelativeToPhiSector0/PHI_RANGE);


      const double phi0RelativeToPhiSector0 = RootTools::getDeltaAngleDeg(recoPhiDeg0, -aftForeOffset);
      Int_t phiSectorOfPeak0 = TMath::Nint(phi0RelativeToPhiSector0/PHI_RANGE);
      
      Int_t deltaPhiSect = NUM_PHI;
      Int_t deltaPhiSect0 = NUM_PHI;            
      for(int phi=0; phi<NUM_PHI; phi++){
	UInt_t phiMask = RootTools::getBit(phi, header->getL3TrigPattern(pol));

	if(phiMask > 0){
	  Int_t dPhiSect = phiSectorOfPeak - phi;
	  if(dPhiSect < -NUM_PHI/2){
	    dPhiSect += NUM_PHI;
	  }
	  else if(dPhiSect >= NUM_PHI/2){
	    dPhiSect -= NUM_PHI;
	  }
	  if(TMath::Abs(dPhiSect) < deltaPhiSect){
	    deltaPhiSect = dPhiSect;
	  }


	  Int_t dPhiSect0 = phiSectorOfPeak0 - phi;
	  if(dPhiSect0 < -NUM_PHI/2){
	    dPhiSect0 += NUM_PHI;
	  }
	  else if(dPhiSect0 >= NUM_PHI/2){
	    dPhiSect0 -= NUM_PHI;
	  }
	  if(TMath::Abs(dPhiSect) < deltaPhiSect0){
	    deltaPhiSect0 = dPhiSect0;
	  }
	  
	}
      }

      hDeltaPhiSect->Fill(deltaPhiSect);
      hDeltaPhiSect0->Fill(deltaPhiSect0);

      if(TMath::Abs(deltaPhiSect) > 1){
	p.inc(entry, maxEntry);	
	continue;
      }
      
      

      // if(imagePeak < cutImage && hilbertPeak < cutHilbert){
      if(imagePeak < cutImage || hilbertPeak < cutHilbert){
	p.inc(entry, maxEntry);	
	continue;
      }

      if(recoPhiDeg < 0) recoPhiDeg += 360;
      else if(recoPhiDeg >= 360) recoPhiDeg -= 360;

      if(recoPhiDeg0 < 0) recoPhiDeg0 += 360;
      else if(recoPhiDeg0 >= 360) recoPhiDeg0 -= 360;

	
      if(recoPhiDeg0 < 0 || recoPhiDeg < 0){
	std::cerr << recoPhiDeg0 << "\t" << recoPhiDeg <<"\t" << std::endl;
      }

      Double_t recoThetaDeg = eventSummary->peak[pol][1].theta;
      Double_t recoThetaDeg0 = eventSummary->peak[pol][0].theta;

      // if(header->run==352){
      // if(recoThetaDeg < -25 && recoThetaDeg > -33){

      // if(recoThetaDeg > 20 && recoThetaDeg < 30){
      //   std::cout << std::endl << header->run << "\t" << header->eventNumber << std::endl;
      // }
      // if(imagePeak > 0.2 && hilbertPeak < 100){
      //   std::cout << "coherent-low-power: " << header->run << "\t" << header->eventNumber << std::endl;
      // }
      // if(hilbertPeak > 80){
      //   std::cout << "high-power-incoherent: " << header->run << "\t" << header->eventNumber << std::endl;
      // }
      // if(imagePeak > 0.088){
      //   std::cout << "What are these: " << header->run << "\t" << header->eventNumber << std::endl;
      // }
      // }

      // std::cerr << recoPhiDeg << "\t" << recoThetaDeg << "\t"
      // 	  << recoThetaDeg0 << "\t" << recoPhiDeg << std::endl;

      Double_t directionWrtNorth = RootTools::getDeltaAngleDeg(pat->heading, recoPhiDeg);

      // std::cerr << "\t" << pat->heading << "\t" << recoPhiDeg << "\t" << directionWrtNorth << std::endl;
      directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
      directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;

      // std::cerr << "\t" << header->realTime - firstRealTime << "\t" << lastRealTime - header->realTime << "\t" << (header->realTime >= firstRealTime) << "\t" << (header->realTime < lastRealTime) << std::endl;	

      Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(solarPhiDeg, recoPhiDeg);
      solarPhiDeg = solarPhiDeg < 0 ? solarPhiDeg + 360 : solarPhiDeg;
      // deltaSolarPhiDeg = deltaSolarPhiDeg < 0 ? deltaSolarPhiDeg + 360 : deltaSolarPhiDeg;	            
      
      Double_t deltaSolarThetaDeg = solarThetaDeg - recoThetaDeg;

      if(TMath::Abs(directionWrtNorth) > maxDirWrtNorth){
	std::cerr << recoPhiDeg << "\t" << pat->heading << std::endl;
      }


      hDeltaSolarPhiDeg->Fill(deltaSolarPhiDeg);
      hDeltaSolarThetaDeg->Fill(deltaSolarThetaDeg);
      hDeltaSolarPhiDegZoom->Fill(solarPhiDeg, deltaSolarPhiDeg);
      hDeltaSolarThetaDegZoom->Fill(solarPhiDeg, deltaSolarThetaDeg);

      hThetaDeg->Fill(recoThetaDeg);

      // std::cerr << solarThetaDeg << "\t" << solarPhiDeg << "\t"
      // 	  << deltaSolarThetaDeg << "\t" << deltaSolarPhiDeg << "\t"
      // 	  << std::endl;

      
      hDeltaSolarThetaDegZoomVsTheta->Fill(solarThetaDeg, deltaSolarThetaDeg);	
      hDeltaSolarThetaDegZoomVsPhi->Fill(solarPhiDeg, deltaSolarThetaDeg);
      hDeltaSolarThetaDegZoomVsTimeOfDay->Fill(pat->timeOfDay/1000, deltaSolarThetaDeg);
      hDeltaSolarPhiDegZoomVsTheta->Fill(solarThetaDeg, deltaSolarPhiDeg);
      hDeltaSolarPhiDegZoomVsPhi->Fill(solarPhiDeg, deltaSolarPhiDeg);
      hDeltaSolarPhiDegZoomVsTimeOfDay->Fill(pat->timeOfDay/1000, deltaSolarPhiDeg);

      hDeltaSolarPhiDegVsDeltaSolarPhiDeg->Fill(deltaSolarPhiDeg, deltaSolarThetaDeg);
      hImagePeakVsDeltaSolarPhiDegVsDeltaSolarPhiDeg->Fill(deltaSolarPhiDeg, deltaSolarThetaDeg, imagePeak);
      

      double deltaPhiCoarseZoom = RootTools::getDeltaAngleDeg(recoPhiDeg0, recoPhiDeg);

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
	

      // reconstruction data quality
      hImagePeakPhi->Fill(recoPhiDeg, imagePeak);
      hImagePeakPhi0->Fill(recoPhiDeg0,
			   eventSummary->peak[pol][0].value);
      hImagePeakDeltaPhi->Fill(deltaPhiCoarseZoom,
			       imagePeak);	
      hPhi0DeltaPhi->Fill(recoPhiDeg0, deltaPhiCoarseZoom);	
      if(deltaPhiCoarseZoom > -4.94 && deltaPhiCoarseZoom < 4.99){
	hPhi0DeltaPhi2->Fill(recoPhiDeg0, deltaPhiCoarseZoom);
      }
      hPhi0Theta0->Fill(recoPhiDeg0,
			recoThetaDeg0+0.0001);

      hHeading->Fill(header->realTime, pat->heading);
    }

    p.inc(entry, maxEntry);
  }
  
  outFile->Write();
  outFile->Close();

  return 0;
}



  
//   RawAnitaHeader* header = NULL;
//   headChain->SetBranchAddress("header", &header);
//   Adu5Pat* pat = NULL;
//   gpsChain->SetBranchAddress("pat", &pat);
//   AnitaEventSummary* eventSummary = NULL;
//   eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);

//   OutputConvention oc(argc, argv);
//   TString outFileName = oc.getOutputFileName();
//   TFile* outFile = new TFile(outFileName, "recreate");
//   if(outFile->IsZombie()){
//     std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
//     return 1;
//   }

//   headChain->BuildIndex("eventNumber");
  
//   headChain->GetEntry(0);
//   UInt_t firstRealTime = header->realTime;
//   headChain->GetEntry(headChain->GetEntries()-2);
//   UInt_t lastRealTime = header->realTime;

//   // std::cerr << firstRealTime << "\t" << lastRealTime << std::endl;

//   Long64_t nEntries = eventSummaryChain->GetEntries();
//   Long64_t maxEntry = 0; //2500;
//   Long64_t startEntry = 0;
//   if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
//   std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
//   ProgressBar p(maxEntry-startEntry);

//   // const int numTrigTypes = 2; // rf=0 and min bias=1

//   TH2D* hPeakDirWrtNorth[NUM_POL];
//   TH2D* hPeakTheta[NUM_POL];
//   TH2D* hImagePeakHilbertPeak[NUM_POL];
//   TH2D* hImagePeakTime[NUM_POL];
//   TH2D* hImagePeakPhi[NUM_POL];
//   TH2D* hImagePeakPhi0[NUM_POL];
//   TH2D* hImagePeakDeltaPhi[NUM_POL];
//   TH2D* hPhi0DeltaPhi[NUM_POL];
//   TH2D* hPhi0DeltaPhi2[NUM_POL];
//   TH2D* hPhi0Theta0[NUM_POL];
//   TH2D* hHilbertPeakTime[NUM_POL];

//   const int numBinsPhi = 420;
//   const int numBinsTheta = numBinsPhi/2;
//   const int numTimeBins = 1024*2;
//   const double maxDirWrtNorth = 180;
//   const double minDirWrtNorth = -180;
//   const double maxTheta = 90;
//   const double minTheta = -90;

//   // const double deltaSolarPhiDegCut = 20; // degrees
//   const double deltaSolarPhiDegCut = 10; // degrees

//   const Int_t numImagePeakBins = 1024;
//   const Int_t numHilbertPeakBins = 1024;
//   // const Double_t maxHilbertPeak = 1000;
//   const Double_t maxHilbertPeak = 2048;  

//   // const Double_t ipCenter = 0.04;
//   // const Double_t hpCenter = 20;
//   // const Double_t ipScale = 0.55 - ipCenter;
//   // const Double_t hpScale = 30 - hpCenter;

//   for(int pol=0; pol<NUM_POL; pol++){

//     TString name = "hPeakDirWrtNorth_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     TString title = " Peak direction w.r.t North vs. time for all min bias events ";
//     title += "; Time; Direction wrt. north (Degrees)";

//     hPeakDirWrtNorth[pol] = new TH2D(name, title,
// 				     numTimeBins,
// 				     firstRealTime, lastRealTime+1,
// 				     numBinsPhi,
// 				     minDirWrtNorth, maxDirWrtNorth);

//     name = "hPeakTheta_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title += "Peak direction #theta vs. time for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; Time; #theta (Degrees)";
	
//     hPeakTheta[pol] = new TH2D(name, title,
// 			       numTimeBins, firstRealTime, lastRealTime+1,
// 			       numBinsTheta, minTheta, maxTheta);

//     name = "hImagePeakHilbertPeak_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "Image peak vs. Hilbert Peak for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; Image peak (no units); Hilbert peak (mV); Number of events";

//     hImagePeakHilbertPeak[pol] = new TH2D(name, title,
// 					  numImagePeakBins, 0, 1,
// 					  numHilbertPeakBins, 0, maxHilbertPeak);


//     name = "hImagePeakPhi_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "Image peak vs. #Phi for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #Phi (Degrees); Image peak (no units); Number of events";

//     hImagePeakPhi[pol] = new TH2D(name, title,
// 				  numBinsPhi, 0, 360,
// 				  numImagePeakBins, 0, 1);

//     name = "hImagePeakPhi0_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "Image peak vs. #Phi for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #Phi (Degrees); Image peak (no units); Number of events";

//     hImagePeakPhi0[pol] = new TH2D(name, title,
// 				   NUM_BINS_PHI*NUM_PHI, 0, 360,
// 				   numImagePeakBins, 0, 1);
	


//     name = "hImagePeakDeltaPhi_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "Image peak vs. #Phi for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #delta#Phi (Degrees); Image peak (no units); Number of events";

//     hImagePeakDeltaPhi[pol] = new TH2D(name, title,
// 				       1024, -6, 6,
// 				       numImagePeakBins, 0, 1);
	

//     name = "hPhi0DeltaPhi_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "#delta#Phi vs. #Phi_{0} for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #Phi_{0} (Degrees); #delta#Phi (Degrees); Number of events";

//     hPhi0DeltaPhi[pol] = new TH2D(name, title,
// 				  NUM_BINS_PHI*NUM_PHI, 0, 360,
// 				  1024, -6, 6);


//     name = "hPhi0DeltaPhi2_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "#delta#Phi vs. #Phi_{0} for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #Phi_{0} (Degrees); #delta#Phi (Degrees); Number of events";

//     hPhi0DeltaPhi2[pol] = new TH2D(name, title,
// 				   NUM_BINS_PHI*NUM_PHI, 0, 360,
// 				   1024, -6, 6);
	
//     name = "hPhi0Theta0_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "#Phi_{0} vs #Theta_{0} for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #Phi_{0} (Degrees); #theta_{0} (Degrees); Number of events";

//     hPhi0Theta0[pol] = new TH2D(name, title,
// 				NUM_BINS_PHI*NUM_PHI, 0, 360,
// 				NUM_BINS_THETA, -75, 75);

//     name = "hImagePeakTime_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "Image peak vs. time for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; Time; Image peak (no units); Number of events";

//     hImagePeakTime[pol] = new TH2D(name, title,
// 				   numTimeBins, firstRealTime, lastRealTime+1,
// 				   numImagePeakBins, 0, 1);
//     name = "hHilbertPeakTime_";
//     // name += trig == 0 ? "" : "goodTime_";
//     // name += sun == 0 ? "awayFromSun_" : "towardsSun_";
//     name += pol == AnitaPol::kHorizontal ? "HPol" : "VPol";

//     title = "Hilbert peak vs. time for all min bias events ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; Time; Hilbert peak (mV); Number of events";

//     hHilbertPeakTime[pol] = new TH2D(name, title,
// 				     numTimeBins, firstRealTime, lastRealTime+1,
// 				     numHilbertPeakBins, 0, maxHilbertPeak);
//   }
  
//   TH2D* hHeading = new TH2D("hHeading", "ANITA heading vs. realTime; realTime; heading (Degrees)",
// 			    1024, firstRealTime, lastRealTime+1,
// 			    numBinsPhi, 0, 360);

//   TH1D* hDeltaSolarPhiDeg[NUM_POL];
//   TH1D* hDeltaSolarThetaDeg[NUM_POL];
//   TH2D* hDeltaSolarPhiDegZoom[NUM_POL];
//   TH2D* hDeltaSolarThetaDegZoom[NUM_POL];
//   TH1D* hThetaDeg[NUM_POL];

//   TH2D* hDeltaSolarThetaDegZoomVsTheta[NUM_POL];
//   TH2D* hDeltaSolarThetaDegZoomVsPhi[NUM_POL];
//   TH2D* hDeltaSolarThetaDegZoomVsTimeOfDay[NUM_POL];
//   TH2D* hDeltaSolarPhiDegZoomVsTheta[NUM_POL];
//   TH2D* hDeltaSolarPhiDegZoomVsPhi[NUM_POL];
//   TH2D* hDeltaSolarPhiDegZoomVsTimeOfDay[NUM_POL];
//   TH2D* hImagePeakSolarTheta[NUM_POL];
  
//   for(int pol=0; pol<NUM_POL; pol++){
//     TString name = pol==AnitaPol::kHorizontal ? "hDeltaSolarPhiDegHPol" : "hDeltaSolarPhiDegVPol";
//     TString title = "#delta#phi_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; #delta#phi_{sun} (Degrees); Number of events / bin";

//     hDeltaSolarPhiDeg[pol] = new TH1D(name, title,
// 				      numBinsPhi, -180, 180);

//     name = pol==AnitaPol::kHorizontal ? "hDeltaSolarThetaDegHPol" : "hDeltaSolarThetaDegVPol";
//     title = "#delta#theta_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; #delta#theta_{sun} (Degrees); Number of events / bin";

//     hDeltaSolarThetaDeg[pol] = new TH1D(name, title,
// 					numBinsPhi, -180, 180);

//     name = pol==AnitaPol::kHorizontal ? "hDeltaSolarPhiDegZoomHPol" : "hDeltaSolarPhiDegZoomVPol";
//     title = "#delta#phi_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; Measured #phi (Degrees); #delta#phi_{sun} (Degrees); Number of events / bin";

//     hDeltaSolarPhiDegZoom[pol] = new TH2D(name, title,
// 					  numBinsPhi, 0, 360,
// 					  numBinsPhi, -10, 10);

//     name = pol==AnitaPol::kHorizontal ? "hDeltaSolarThetaDegZoomHPol" : "hDeltaSolarThetaDegZoomVPol";
//     title = "#delta#theta_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; Measured #phi (Degrees); #delta#theta_{sun} (Degrees); Number of events / bin";

//     hDeltaSolarThetaDegZoom[pol] = new TH2D(name, title,
// 					    numBinsPhi, -0, 360,
// 					    numBinsPhi, -10, 10);

//     name = pol==AnitaPol::kHorizontal ? "hThetaDegHPol" : "hThetaDegVPol";
//     title = "Min bias measured #theta ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "#theta (Degrees); Number of events / bin";

//     hThetaDeg[pol] = new TH1D(name, title,
// 			      numBinsPhi, -90, 90);


//     name = "hDeltaSolarThetaDegZoomVsTheta";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "#delta#theta_{sun} vs. #theta_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; #theta (Degrees); #delta#theta_{sun}; Number of events / bin";

//     hDeltaSolarThetaDegZoomVsTheta[pol] = new TH2D(name, title,
// 						   numBinsPhi, -90, 90,
// 						   numBinsPhi, -10, 10);


//     name = "hDeltaSolarThetaDegZoomVsPhi";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "#delta#theta_{sun} vs. #phi_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; #phi (Degrees); #delta#theta_{sun}; Number of events / bin";

//     hDeltaSolarThetaDegZoomVsPhi[pol] = new TH2D(name, title,
// 						 numBinsPhi, 0, 360,
// 						 numBinsPhi, -10, 10);

//     name = "hDeltaSolarThetaDegZoomVsTimeOfDay";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "#delta#theta_{sun} vs. timeOfDay ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; Time of day; #delta#theta_{sun}; Number of events / bin";

//     hDeltaSolarThetaDegZoomVsTimeOfDay[pol] = new TH2D(name, title,
// 						       1024, 0, 24*60*60,
// 						       numBinsPhi, -10, 10);

//     name = "hDeltaSolarPhiDegZoomVsTheta";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "#delta#phi_{sun} vs. #theta_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; #theta (Degrees); #delta#phi_{sun}; Number of events / bin";

//     hDeltaSolarPhiDegZoomVsTheta[pol] = new TH2D(name, title,
// 						   numBinsPhi, -90, 90,
// 						   numBinsPhi, -10, 10);


//     name = "hDeltaSolarPhiDegZoomVsPhi";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "#delta#phi_{sun} vs. #phi_{sun} ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; #phi (Degrees); #delta#phi_{sun}; Number of events / bin";

//     hDeltaSolarPhiDegZoomVsPhi[pol] = new TH2D(name, title,
// 						 numBinsPhi, 0, 360,
// 						 numBinsPhi, -10, 10);

//     name = "hDeltaSolarPhiDegZoomVsTimeOfDay";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "#delta#phi_{sun} vs. timeOfDay ";
//     title += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title += "; Time of day; #delta#phi_{sun}; Number of events / bin";

//     hDeltaSolarPhiDegZoomVsTimeOfDay[pol] = new TH2D(name, title,
// 						       1024, 0, 24*60*60,
// 						       numBinsPhi, -10, 10);    

//     name = "hImagePeakSolarTheta";
//     name += pol==AnitaPol::kHorizontal ? "HPol" : "VPol";
//     title = "Image peak vs. time of day ";
//     title += pol == AnitaPol::kHorizontal ? "(HPol)" : "(VPol)";
//     title += "; #theta_{sun} (Degrees); Image peak (no units)";
    
//     hImagePeakSolarTheta[pol] = new TH2D(name, title,
// 					 numBinsTheta, minTheta, maxTheta,
// 					 numImagePeakBins, 0, 1);

//   }

//   for(Long64_t entry = startEntry; entry < maxEntry; entry++){
//     eventSummaryChain->GetEntry(entry);
//     Int_t headEntry = headChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
//     // std::cerr << headEntry << "\t" << std::endl;
//     if(headEntry < 0){
//       std::cerr << "Now what!?" << std::endl;      
//     }
//     else{
//       headChain->GetEntry(headEntry);
//       gpsChain->GetEntry(headEntry);

//       Double_t solarPhiDeg, solarThetaDeg;

//       // usefulPat.getSunPosition(solarPhiDeg, solarThetaDeg);
//       solarPhiDeg = eventSummary->sun.phi;
//       solarThetaDeg = eventSummary->sun.theta;
    
//       for(Int_t polInd=0; polInd < NUM_POL; polInd++){

// 	Double_t imagePeak = eventSummary->peak[polInd][1].value;
// 	Double_t hilbertPeak = eventSummary->coherent[polInd][1].peakHilbert;

// 	Double_t recoPhiDeg = eventSummary->peak[polInd][1].phi;
// 	Double_t recoPhiDeg0 = eventSummary->peak[polInd][0].phi;

// 	if(recoPhiDeg < 0) recoPhiDeg += 360;
// 	else if(recoPhiDeg >= 360) recoPhiDeg -= 360;

// 	if(recoPhiDeg0 < 0) recoPhiDeg0 += 360;
// 	else if(recoPhiDeg0 >= 360) recoPhiDeg0 -= 360;


// 	if(recoPhiDeg0 < 0 || recoPhiDeg < 0){
// 	  std::cerr << recoPhiDeg0 << "\t" << recoPhiDeg <<"\t" << std::endl;
// 	}

// 	Double_t recoThetaDeg = eventSummary->peak[polInd][1].theta;
// 	Double_t recoThetaDeg0 = eventSummary->peak[polInd][0].theta;

// 	// std::cerr << recoPhiDeg << "\t" << recoThetaDeg << "\t"
// 	// 	  << recoThetaDeg0 << "\t" << recoPhiDeg << std::endl;

// 	Double_t directionWrtNorth = RootTools::getDeltaAngleDeg(pat->heading, recoPhiDeg);

// 	// std::cerr << "\t" << pat->heading << "\t" << recoPhiDeg << "\t" << directionWrtNorth << std::endl;

// 	directionWrtNorth = directionWrtNorth < -180 ? directionWrtNorth + 360 : directionWrtNorth;
// 	directionWrtNorth = directionWrtNorth > 180 ? directionWrtNorth - 360 : directionWrtNorth;

// 	// std::cerr << "\t" << header->realTime - firstRealTime << "\t" << lastRealTime - header->realTime << "\t" << (header->realTime >= firstRealTime) << "\t" << (header->realTime < lastRealTime) << std::endl;	
	

// 	Double_t deltaSolarPhiDeg = RootTools::getDeltaAngleDeg(solarPhiDeg, recoPhiDeg);
// 	Double_t deltaSolarThetaDeg = solarThetaDeg - recoThetaDeg;


// 	if(TMath::Abs(directionWrtNorth) > maxDirWrtNorth){
// 	  std::cerr << recoPhiDeg << "\t" << pat->heading << std::endl;
// 	}

// 	Int_t sun = TMath::Abs(deltaSolarPhiDeg) > deltaSolarPhiDegCut ? 0 : 1;

// 	hDeltaSolarPhiDeg[polInd]->Fill(deltaSolarPhiDeg);
// 	hDeltaSolarThetaDeg[polInd]->Fill(deltaSolarThetaDeg);
// 	hDeltaSolarPhiDegZoom[polInd]->Fill(solarPhiDeg, deltaSolarPhiDeg);
// 	hDeltaSolarThetaDegZoom[polInd]->Fill(solarPhiDeg, deltaSolarThetaDeg);

// 	hThetaDeg[polInd]->Fill(recoThetaDeg);

// 	if(sun){
// 	  hDeltaSolarThetaDegZoomVsTheta[polInd]->Fill(solarThetaDeg, deltaSolarThetaDeg);	
// 	  hDeltaSolarThetaDegZoomVsPhi[polInd]->Fill(solarPhiDeg, deltaSolarThetaDeg);
// 	  hDeltaSolarThetaDegZoomVsTimeOfDay[polInd]->Fill(pat->timeOfDay/1000, deltaSolarThetaDeg);
// 	  hDeltaSolarPhiDegZoomVsTheta[polInd]->Fill(solarThetaDeg, deltaSolarPhiDeg);
// 	  hDeltaSolarPhiDegZoomVsPhi[polInd]->Fill(solarPhiDeg, deltaSolarPhiDeg);
// 	  hDeltaSolarPhiDegZoomVsTimeOfDay[polInd]->Fill(pat->timeOfDay/1000, deltaSolarPhiDeg);
// 	}

// 	double deltaPhiCoarseZoom = RootTools::getDeltaAngleDeg(recoPhiDeg0, recoPhiDeg);

// 	hPeakDirWrtNorth[polInd]->Fill(header->realTime,
// 				       directionWrtNorth);

// 	hPeakTheta[polInd]->Fill(header->realTime,
// 				 recoThetaDeg);

// 	hImagePeakHilbertPeak[polInd]->Fill(imagePeak,
// 					    hilbertPeak);
// 	hImagePeakPhi[polInd]->Fill(recoPhiDeg, imagePeak);
// 	hImagePeakPhi0[polInd]->Fill(recoPhiDeg0,
// 				     eventSummary->peak[polInd][0].value);
// 	hImagePeakDeltaPhi[polInd]->Fill(deltaPhiCoarseZoom,
// 					 imagePeak);
// 	hPhi0DeltaPhi[polInd]->Fill(recoPhiDeg0, deltaPhiCoarseZoom);

// 	if(deltaPhiCoarseZoom > -4.94 && deltaPhiCoarseZoom < 4.99){
// 	  hPhi0DeltaPhi2[polInd]->Fill(recoPhiDeg0, deltaPhiCoarseZoom);
// 	}
	
// 	hPhi0DeltaPhi[polInd]->Fill(recoPhiDeg0, deltaPhiCoarseZoom);
      
// 	hPhi0Theta0[polInd]->Fill(recoPhiDeg0,
// 				  recoThetaDeg0+0.0001);

// 	hImagePeakTime[polInd]->Fill(header->realTime,
// 				     imagePeak);

// 	hHilbertPeakTime[polInd]->Fill(header->realTime,
// 				       hilbertPeak);

// 	hImagePeakSolarTheta[polInd]->Fill(solarThetaDeg, imagePeak);
//       }

//       hHeading->Fill(header->realTime, pat->heading);
//       if(pat->heading < 0 || pat->heading >= 360){
// 	std::cerr << "bad heading ? " << pat->heading << "\t" << header->run << "\t" << header->eventNumber << std::endl;
//       }
//     }
//     // p++;
//     p.inc(entry, maxEntry);
//   }

//   // grSunPhiDeg->SetName("grSunPhiDeg");
//   // grSunPhiDeg->SetTitle("Predicted sun #phi (Degrees)");
//   // grSunThetaDeg->SetName("grSunThetaDeg");
//   // grSunThetaDeg->SetTitle("Predicted sun #theta (Degrees)");
//   // grSunPhiDeg->Write();
//   // grSunThetaDeg->Write();
  
//   outFile->Write();
//   outFile->Close();

//   return 0;
// }

