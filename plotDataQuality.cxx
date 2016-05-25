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

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }

    TString fileName = TString::Format("dataQualityTrees/makeDataQualityTreesPlots_%d_*", run);
    // TString fileName = TString::Format("makeDataQualityTreesPlots_%d_*", run);    
    dataQualityChain->Add(fileName);
  }

  if(dataQualityChain->GetEntries()==0){
    std::cerr << "Unable to find header files!" << std::endl;
    return 1;
  }

  Double_t maxAbsSecondDeriv[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("maxAbsSecondDeriv", &maxAbsSecondDeriv[0][0]);
  Double_t maxVolts[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("maxVolts", &maxVolts[0][0]);
  Int_t numPoints[NUM_POL][NUM_SEAVEYS];
  dataQualityChain->SetBranchAddress("numPoints", &numPoints[0][0]);

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
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

  const double maxVoltsThresh = 400;
  const double saturationVolts = 1500;
  // const int maxPhiAboveMaxVoltsThresh = 9;
  TH2D* hMaxAbsSecondDeriv[NUM_POL];
  TH2D* hMaxVolts[NUM_POL];
  TH1D* hNumChannelsAboveMaxVolts[NUM_POL];
  TH1D* hNumPhiAboveMaxVolts[NUM_POL];

  TString name0 = TString::Format("hAnyChannelsAboveSaturation");
  TString title0 = TString::Format("Any sample > %4.0lf (mV); Sample greater than %4.0lf mV? ; Number of events",
				  saturationVolts, saturationVolts);
  TH1D* hAnyChannelsAboveSaturation = new TH1D(name0, title0, 2, 0, 2);
  hAnyChannelsAboveSaturation->GetXaxis()->SetBinLabel(1, "False");
  hAnyChannelsAboveSaturation->GetXaxis()->SetBinLabel(2, "True");
  hAnyChannelsAboveSaturation->GetXaxis()->SetTitleOffset(0.01);

  name0 = TString::Format("hNumPoints");
  title0 = TString::Format("Number of samples in waveform; Number of samples; Number of channels");
  TH1D* hNumPoints = new TH1D(name0, title0, NUM_SAMP, 0, NUM_SAMP);
  
  const char* polNames[NUM_POL] = {"HPol", "VPol"};
  for(int polInd=0; polInd < NUM_POL; polInd++){
    TString name = TString::Format("hMaxAbsSecondDeriv_%d", polInd);
    TString title = TString::Format("Maximum value of second derivative %s; Antenna Index; abs(2V_{i} - V_{i-1} - V_{i+1}) (mV); Events per bin", polNames[polInd]);
    hMaxAbsSecondDeriv[polInd] = new TH2D(name, title,
					  NUM_SEAVEYS, 0, NUM_SEAVEYS,
					  4096, 0, 4096);

    name = TString::Format("hMaxVolts_%d", polInd);
    title = TString::Format("Maximum volts %s channels; Antenna Index; Volts (mV); Events per bin", polNames[polInd]);
    hMaxVolts[polInd] = new TH2D(name, title,
				 NUM_SEAVEYS, 0, NUM_SEAVEYS,
				 4096, 0, 4096);

    name = TString::Format("hNumChannelsAboveMaxVolts_%d", polInd);
    title = TString::Format("Number of %s channels where maximum volts > above %4.0lf volts; Number of channels; Number of events", polNames[polInd], maxVoltsThresh);
    hNumChannelsAboveMaxVolts[polInd] = new TH1D(name, title,
						 NUM_SEAVEYS+1, 0, NUM_SEAVEYS+1);

    name = TString::Format("hNumPhiAboveMaxVolts_%d", polInd);
    title = TString::Format("Number of %s #phi-sectors where sample > %4.0lf mV; Number of #phi-sectors; Number of events", polNames[polInd], maxVoltsThresh);
    hNumPhiAboveMaxVolts[polInd] = new TH1D(name, title,
						 NUM_PHI+1, 0, NUM_PHI+1);
    
  }

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    dataQualityChain->GetEntry(entry);

    Int_t anyChannelAboveSaturation = 0;
    for(int polInd=0; polInd < NUM_POL; polInd++){
      Int_t numAboveVoltsThresh = 0;
      Int_t phiAboveMaxVoltsThresh[NUM_PHI] = {0};
      for(int phi=0; phi<NUM_PHI; phi++){
	phiAboveMaxVoltsThresh[phi] = 0;
      }
      
      for(int ant=0; ant < NUM_SEAVEYS; ant++){
	hNumPoints->Fill(numPoints[polInd][ant]);
	hMaxVolts[polInd]->Fill(ant, maxVolts[polInd][ant]);
	if(maxVolts[polInd][ant] > maxVoltsThresh){
	  numAboveVoltsThresh++;
	  int phi = ant%NUM_PHI;
	  phiAboveMaxVoltsThresh[phi] = 1;
	}

	if(maxVolts[polInd][ant] > saturationVolts){
	  anyChannelAboveSaturation = 1;
	}
      }

      Int_t numPhiAboveMaxVoltsThresh = 0;
      for(int phi=0; phi<NUM_PHI; phi++){
	if(phiAboveMaxVoltsThresh[phi] > 0){
	  numPhiAboveMaxVoltsThresh++;
	}
      }      
      hNumChannelsAboveMaxVolts[polInd]->Fill(numAboveVoltsThresh);
      hNumPhiAboveMaxVolts[polInd]->Fill(numPhiAboveMaxVoltsThresh);
    }
    hAnyChannelsAboveSaturation->Fill(anyChannelAboveSaturation);
    p++;
  }
  outFile->Write();
  outFile->Close();

  return 0;
}
