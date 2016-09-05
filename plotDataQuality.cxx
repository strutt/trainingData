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
// #include "DataQualityMonitor.h"



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

  TChain* dataQualityChain = new TChain("dataQualityTree");

  TString subFileName = "Decimated";
  // TString subFileName = "Wais";

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    // TString fileName = TString::Format("allMinBiasDistributions260MHzAnd370MHzFiltered/makeDataQualityTreesPlots_%d_*.root", run);
    // TString fileName = TString::Format("finalDataQuality/makeDecimatedDataQualityTreesPlots_%d_*.root", run);
    TString fileName = TString::Format("finalDataQuality/make%sDataQualityTreesPlots_%d_*.root", subFileName.Data(), run);
    // TString fileName = TString::Format("dataQualityTrees/makeDataQualityTreesPlots_%d_*", run);
    // TString fileName = TString::Format("makeDataQualityTreesPlots_%d_*", run);
    dataQualityChain->Add(fileName);
  }

  if(dataQualityChain->GetEntries()==0){
    std::cerr << "Unable to find data quality files!" << std::endl;
    return 1;
  }


  Double_t peakToPeak[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("peakToPeak", peakToPeak);
  Double_t maxVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("maxVolts", maxVolts);
  Double_t minVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("minVolts", minVolts);
  // Double_t power[NUM_POL][NUM_SEAVEYS];
  // dataQualityChain->SetBranchAddress("power", power);

  UInt_t eventNumberDQ;
  dataQualityChain->SetBranchAddress("eventNumber", &eventNumberDQ);


  //  DataQualityMonitor dqm(dataQualityChain);


  OutputConvention oc(argc, argv);
  TString outFileNameTemp = oc.getOutputFileName();
  TString outFileName = TString(outFileNameTemp(0,4)) + subFileName + TString(outFileNameTemp(4, outFileNameTemp.Length()));

  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }


  Long64_t nEntries = dataQualityChain->GetEntries();
  Long64_t maxEntry = 0; //2500;
  Long64_t startEntry = 0;
  if(maxEntry<=0 || maxEntry > nEntries) maxEntry = nEntries;
  std::cout << "Processing " << maxEntry << " of " << nEntries << " entries." << std::endl;
  ProgressBar p(maxEntry-startEntry);

  // const int maxPhiAboveMaxVoltsThresh = 9;

  TH1D* hMaxMaxVolts = new TH1D("hMaxMaxVolts", "Maximum Volts; Maximum Volts (mV); Events per bin", 4096, 0, 4096);
  TH1D* hMinMinVolts = new TH1D("hMinMinVolts", "Minimum Volts; Minimum Volts (mV); Events per bin", 4096, -4096, 0);

  TH1D* hMaxMinVolts = new TH1D("hMaxMinVolts", "Sum of max and min volts; Max Volts + Min Volts", 4096, 0, 4096);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){

    dataQualityChain->GetEntry(entry);

    Double_t absMaxSumVolts = 0;
    Double_t maxMaxVolts = 0;
    Double_t minMinVolts = 0;
    for(int polInd=0; polInd < NUM_POL; polInd++){

      for(int ant=0; ant < NUM_SEAVEYS; ant++){

	if(maxVolts[polInd][ant] > maxMaxVolts){
	  maxMaxVolts = maxVolts[polInd][ant];
	}

	if(minVolts[polInd][ant] < minMinVolts){
	  minMinVolts = minVolts[polInd][ant];
	}

	if(TMath::Abs(minVolts[polInd][ant] + maxVolts[polInd][ant]) > TMath::Abs(absMaxSumVolts)){
	  absMaxSumVolts = minVolts[polInd][ant] + maxVolts[polInd][ant];
	}
      }
    }

    // if(maxMaxVolts > 1.5e3){
    //   TString currentFileName = dataQualityChain->GetCurrentFile()->GetName();
    //   std::cout << std::endl << "max " << eventNumberDQ << "\t" << TString(currentFileName(47, 3)) << std::endl;
    // }
    // else if(minMinVolts < -1e3){
    //   TString currentFileName = dataQualityChain->GetCurrentFile()->GetName();
    //   std::cout << std::endl << "min " << eventNumberDQ << "\t" << TString(currentFileName(47, 3)) << std::endl;
    // }
    // else if(absMaxSumVolts > 500){
    //   TString currentFileName = dataQualityChain->GetCurrentFile()->GetName();
    //   std::cout << std::endl << "diff " << eventNumberDQ << "\t" << TString(currentFileName(47, 3)) << std::endl;
    // }

    hMaxMinVolts->Fill(absMaxSumVolts);
    hMaxMaxVolts->Fill(maxMaxVolts);
    hMinMinVolts->Fill(minMinVolts);

    p.inc(entry, maxEntry);

  }
  outFile->Write();
  outFile->Close();

  return 0;
}
