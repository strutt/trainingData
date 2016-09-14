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
#include "TGraph2D.h"
#include "TRandom3.h"

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

#include "AnitaClusterer.h"
#include "AnalysisCuts.h"


int main(int argc, char *argv[]){

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

  TChain* eventSummaryChain = new TChain("cutTree");
  // TChain* headChain = new TChain("headTree");
  // TChain* indexedHeadChain = new TChain("headTree");
  // TChain* gpsChain = new TChain("adu5PatTree");

  for(Int_t run=firstRun; run<=lastRun; run++){
    if(run>=257 && run<=263){
      continue;
    }
    if(run == 198 || run == 287){
      continue;
    }

    TString fileName = TString::Format("projectionPlots/projectionDecimatedPlots_%d_*", run);
    eventSummaryChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/headFile%d.root", run, run);
    // TString fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/decimatedHeadFile%d.root", run, run);
    // headChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/indexedBlindHeadFile%d.root", run, run);
    // indexedHeadChain->Add(fileName);

    // fileName = TString::Format("~/UCL/ANITA/flight1415/root/run%d/gpsEvent%d.root", run, run);
    // gpsChain->Add(fileName);
  }

  if(eventSummaryChain->GetEntries()==0){
    std::cerr << "Unable to find eventSummary files!" << std::endl;
    return 1;
  }

  // if(gpsChain->GetEntries()==0){
  //   std::cerr << "Unable to find gps files!" << std::endl;
  //   return 1;
  // }

  // if(indexedHeadChain->GetEntries()==0){
  //   std::cerr << "Unable to find header files!" << std::endl;
  //   return 1;
  // }
  // indexedHeadChain->BuildIndex("eventNumber");

  TFile* fRes = TFile::Open("resolutionPlots.root");
  if(fRes==NULL){
    std::cerr << "Warning! Unable to find resolution plots!" << std::endl;
  }
  TH1D* hResTheta[AnitaPol::kNotAPol];
  TH1D* hResPhi[AnitaPol::kNotAPol];
  hResTheta[AnitaPol::kHorizontal] = (TH1D*) fRes->Get("hDeltaThetaWais2_2__4__4");
  hResPhi[AnitaPol::kHorizontal] = (TH1D*) fRes->Get("hDeltaPhiWais2_2__3__3");
  hResTheta[AnitaPol::kVertical] = (TH1D*) fRes->Get("hDeltaThetaLdb2_2__2__2");
  hResPhi[AnitaPol::kVertical] = (TH1D*) fRes->Get("hDeltaPhiLdb2_2__1__1");
  // std::cout << hResTheta[AnitaPol::kHorizontal] << "\t" << hResPhi[AnitaPol::kHorizontal] << "\t";
  // std::cout << hResTheta[AnitaPol::kVertical] << "\t" << hResPhi[AnitaPol::kVertical] << std::endl;

  // std::cout << hResTheta[AnitaPol::kHorizontal]->GetListOfFunctions()->GetEntries() << "\t"
  // 	    << hResPhi[AnitaPol::kHorizontal]->GetListOfFunctions()->GetEntries() << "\t";
  // std::cout << hResTheta[AnitaPol::kVertical]->GetListOfFunctions()->GetEntries() << "\t"
  // 	    << hResPhi[AnitaPol::kVertical]->GetListOfFunctions()->GetEntries() << std::endl;


  TF1* fThetaFit = (TF1*) (hResTheta[AnitaPol::kHorizontal]->GetListOfFunctions()->At(0));
  TF1* fPhiFit = (TF1*) (hResPhi[AnitaPol::kVertical]->GetListOfFunctions()->At(0));
  // TF1* fThetaFit = NULL;
  // TF1* fPhiFit = NULL;
  // std::cout << hResTheta[AnitaPol::kVertical]->GetListOfFunctions() << "\t"
  // 	    <<  hResPhi[AnitaPol::kHorizontal]->GetListOfFunctions() << std::endl;


  AnitaEventSummary* eventSummary = NULL;
  eventSummaryChain->SetBranchAddress("eventSummary", &eventSummary);
  AnitaPol::AnitaPol_t pol = AnitaPol::kHorizontal;
  eventSummaryChain->SetBranchAddress("pol", &pol);



  // AnalysisCuts::Status_t l3TriggerCut;
  // eventSummaryChain->SetBranchAddress("l3TriggerCut", (int*)&l3TriggerCut);
  // AnalysisCuts::Status_t surfSaturation;
  // eventSummaryChain->SetBranchAddress("surfSaturation", (int*)&surfSaturation);
  // AnalysisCuts::Status_t selfTriggeredBlastCut;
  // eventSummaryChain->SetBranchAddress("selfTriggeredBlastCut", (int*)&selfTriggeredBlastCut);
  // AnalysisCuts::Status_t thermalCut;
  // eventSummaryChain->SetBranchAddress("thermalCut", (int*)&thermalCut);
  // AnalysisCuts::Status_t thetaAngleCut;
  // eventSummaryChain->SetBranchAddress("thetaAngleCut", (int*)&thetaAngleCut);
  // AnalysisCuts::Status_t peakRatioCut;
  // eventSummaryChain->SetBranchAddress("peakRatioCut", (int*)&peakRatioCut);


  Adu5Pat* pat = NULL;
  // gpsChain->SetBranchAddress("pat", &pat);
  eventSummaryChain->SetBranchAddress("pat", &pat);
  // UInt_t eventNumberPat = 0;
  // gpsChain->SetBranchAddress("eventNumber", &eventNumberPat);

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

  const int K = 50;
  const int numIterations = 10;
  AnitaClusterer clusterer(K, numIterations, nEntries);

  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    eventSummaryChain->GetEntry(entry);

    // if(l3TriggerCut==0 && surfSaturation==0 && selfTriggeredBlastCut==0 && thermalCut==0 && thetaAngleCut == 0 && peakRatioCut==0){

      // Long64_t entry2 = indexedHeadChain->GetEntryNumberWithIndex(eventSummary->eventNumber);
      // gpsChain->GetEntry(entry);

      // if(eventSummary->eventNumber!=eventNumberPat){
      //   std::cerr << "fuck" << std::endl;
      // }

      for(Int_t polInd=0; polInd < NUM_POL; polInd++){
	for(int peakInd=0; peakInd < 5; peakInd++){
	  Double_t sourceLat = eventSummary->peak[polInd][peakInd].latitude;
	  Double_t sourceLon = eventSummary->peak[polInd][peakInd].longitude;
	  Double_t sourceAlt = eventSummary->peak[polInd][peakInd].altitude;
	  if(sourceLat > -999 && sourceLon > -999){
	    if(eventSummary->peak[polInd][peakInd].distanceToSource >= 0 && eventSummary->peak[polInd][peakInd].distanceToSource < 1e6){

	      Double_t snr = eventSummary->peak[polInd][peakInd].snr;
	      Double_t sigmaTheta = fThetaFit->Eval(snr);
	      Double_t sigmaPhi = fPhiFit->Eval(snr);

	      clusterer.addPoint(pat, sourceLat,sourceLon,sourceAlt, eventSummary->run, eventSummary->eventNumber, sigmaTheta, sigmaPhi, (AnitaPol::AnitaPol_t)polInd);
	    }
	  }
	}
      }
    // }
    p++;
  }

  // clusterer.kMeansCluster(1);
  // clusterer.llCut = 2500;
  // clusterer.llCut = 2000; //2000;
  // clusterer.llCut = 60; //250;
  clusterer.llCut = 250;
  clusterer.initializeBaseList();
  clusterer.recursivelyAddClusters(0);
  // clusterer.mergeClusters();


  TChain* mcEventSummaryChain = new TChain("eventSummaryTree");
  for(int mcRun=1; mcRun < 30; mcRun++){
    TString fileName = TString::Format("monteCarlo/projectionMonteCarloPlots_%d_*.root", mcRun);
    mcEventSummaryChain->Add(fileName);
    // mcGpsChain->Add(fileName);
  }
  const Long64_t numEntries2 = mcEventSummaryChain->GetEntries();
  std::cerr << "numEntries2 = " << numEntries2 << std::endl;
  ProgressBar prog2(numEntries2);
  AnitaEventSummary* eventSummary2 = NULL;
  Adu5Pat* pat2 = NULL;
  Double_t weight = 0;
  mcEventSummaryChain->SetBranchAddress("weight", &weight);
  mcEventSummaryChain->SetBranchAddress("pat", &pat2);
  mcEventSummaryChain->SetBranchAddress("eventSummary", &eventSummary2);
  for(Long64_t entry2=0; entry2 < mcEventSummaryChain->GetEntries(); entry2++){
    mcEventSummaryChain->GetEntry(entry2);

    for(Int_t polInd=0; polInd < NUM_POL; polInd++){
      for(int peakInd=0; peakInd < 5; peakInd++){
	Double_t sourceLat = eventSummary2->peak[polInd][peakInd].latitude;
	Double_t sourceLon = eventSummary2->peak[polInd][peakInd].longitude;
	Double_t sourceAlt = eventSummary2->peak[polInd][peakInd].altitude;
	if(sourceLat > -999 && sourceLon > -999){
	  // std::cout << (eventSummary->peak[polInd][peakInd].distanceToSource < 1e6 )<< std::endl;
	  // if(eventSummary->peak[polInd][peakInd].theta < 0 && eventSummary->peak[polInd][peakInd].distanceToSource < 1e6){
	  if(eventSummary2->peak[polInd][peakInd].distanceToSource < 1e6){

	    Double_t snr = eventSummary2->peak[polInd][peakInd].snr;
	    // Double_t sigmaTheta = hResTheta[polInd]->GetBinContent(hResPhi[polInd]->FindBin(snr));
	    // Double_t sigmaPhi = hResPhi[polInd]->GetBinContent(hResPhi[polInd]->FindBin(snr));
	    Double_t sigmaTheta = fThetaFit->Eval(snr);
	    Double_t sigmaPhi = fPhiFit->Eval(snr);

	    // std::cout << polInd << "\t" << snr << "\t" << sigmaTheta << "\t" << sigmaPhi << std::endl;

	    clusterer.addMCPoint(pat2, sourceLat, sourceLon, sourceAlt,
				 eventSummary2->run, eventSummary2->eventNumber,
				 sigmaTheta, sigmaPhi, (AnitaPol::AnitaPol_t)polInd, weight);
	    // std::cerr << pat2->altitude << "\t" << pat2->latitude << "\t" << pat->longitude << std::endl;
	    // Int_t n = clusterer.addPoint(sourceLat,sourceLon,sourceAlt);
	    // std::cout << n << std::endl;
	  }
	}
      }
    }
    prog2++;
  }

  clusterer.assignMCPointsToClusters();

  outFile->cd();

  for(int clusterInd=0; clusterInd < clusterer.getNumClusters(); clusterInd++){
    TGraph* gr = clusterer.makeClusterSummaryTGraph(clusterInd);

    if(gr){
      gr->Write();
      delete gr;
    }
    else{
      std::cerr << "??????? " << clusterInd << std::endl;
    }
  }

  // no need to write?
  TTree* clusterTree = clusterer.makeClusterSummaryTree(outFile);
  // clusterTree->BuildIndex("eventNumber");
  std::cout << "Made clusterTree with " << clusterTree->GetEntries() << " entries." << std::endl;

  ClusteredAnitaEvent* clusteredEvent = 0;
  clusterTree->SetBranchAddress("clusteredEvent", &clusteredEvent);
  // eventSummaryChain->BuildIndex("eventNumber");



  TH2D* hPeakSizes = new TH2D("hPeakSizes", "Map peaks; P1; P2", 1024, 0, 1, 1024, 0, 1);
  TH2D* hDeltaPeakVsPeak = new TH2D("hDeltaPeakVsPeak", "Map peaks; P1; P2", 1024, -1, 1, 1024, 0, 1);
  TH1D* hDeltaPeak = new TH1D("hDeltaPeak", "Map peaks; P1-P2", 1024, -1, 1);
  TH1D* hDeltaPeak2 = new TH1D("hDeltaPeak2", "Map peaks; p1/(P1-P2)", 1024, -1, 1);
  for(Long64_t entry = startEntry; entry < maxEntry; entry++){
    clusterTree->GetEntry(entry);

    if(clusteredEvent->inCluster < 0){
      eventSummaryChain->GetEntryWithIndex(clusteredEvent->eventNumber);
      if(clusteredEvent->eventNumber == eventSummary->eventNumber){

	// std::cout << clusteredEvent->eventNumber  << "\t" << eventSummary->eventNumber << "\t";
	Double_t hPeak = eventSummary->peak[0][0].value;
	Double_t vPeak = eventSummary->peak[1][0].value;
	Int_t peakPol = hPeak > vPeak ? 0 : 1;
	// std::cout << eventSummary->peak[0][0].value  << "\t" <<  eventSummary->peak[1][0].value << std::endl;
	// std::cout << eventSummary->peak[peakPol][0].value  << "\t" << eventSummary->peak[peakPol][1].value << std::endl;
	Double_t p1 = eventSummary->peak[peakPol][0].value;
	Double_t p2 = eventSummary->peak[peakPol][1].value;
	hPeakSizes->Fill(p1, p2);
	hDeltaPeakVsPeak->Fill(p1-p2, p2);
	hDeltaPeak->Fill(p1-p2);
	hDeltaPeak2->Fill((p1-p2)/p1);
	// std::cout << eventSummary->peak[1][0].value  << std::endl;
      }
    }
  }

  outFile->Write();
  outFile->Close();

  return 0;
}
