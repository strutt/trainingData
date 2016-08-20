#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TString.h>
#include <TMVA/Factory.h>

#include "OutputConvention.h"

int main(int argc, char *argv[]){

  // Create ouput file, factory object and open the input file

  OutputConvention oc(argc, argv);
  TString outFileName = oc.getOutputFileName();
  TFile* outFile = new TFile(outFileName, "recreate");
  if(outFile->IsZombie()){
    std::cerr << "Error! Unable to open output file " << outFileName.Data() << std::endl;
    return 1;
  }
  
  TMVA::Factory* factory = new TMVA::Factory("discriminationTraining", outFile, "");
  TChain* signalChain = new TChain("signalTree");
  signalChain->Add("makeSignalTreePlots_*.root");
  TChain* backgroundChain = new TChain("thermalTree");
  backgroundChain->Add("makeThermalBackgroundTreePlots_*.root");
  
  // get the TTree objects from the input file

  // int nSig = signalChain->GetEntries();
  // int nBkg = backgroundChain->GetEntries();

  // global event weights (see below for setting event-wise weights)

  double sigWeight = 1.0;
  double bkgWeight = 1.0;
  factory->SetInputTrees(signalChain, backgroundChain, sigWeight, bkgWeight);

  // Define the input variables that shall be used for the MVA training
  // (the variables used in the expression must exist in the original TTree).

  factory->AddVariable("imagePeak", 'D');
  factory->AddVariable("hilbertPeak", 'D');
  // factory->AddVariable("z", 'F');

  // Apply additional cuts on the signal and background sample.
  // for example: TCut mycut = "abs(x)<0.5 && abs(y-0.5)<1";

  TCut mycut = "";

  // Use half of the events for training, half for testing

  // TString splitOpt = "NSigTrain=0:NBkgTrain=0:NSigTest=0:NBkgTest=0";
  // this version no longer works in 5.27 -- replace  7.12.10 GDC
  // factory->PrepareTrainingAndTestTree(mycut, splitOpt);
  factory->PrepareTrainingAndTestTree(mycut, 0, 0, 0, 0);

  // Book MVA methods (see TMVA manual).  

  factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher");   
  // factory->BookMethod(TMVA::Types::kMLP, "MLP", "H:!V:HiddenLayers=3");

  // Train, test and evaluate all methods

  factory->TrainAllMethods();
  factory->TestAllMethods();
  // Following line used to work, causes crash with ROOT 5.20.00, should
  // be fixed in root 5.21, see https://savannah.cern.ch/bugs/?40468
  factory->EvaluateAllMethods();    

  // Save the output and finish up

  outFile->Close();
  std::cout << "==> wrote root file TMVA.root" << std::endl;
  std::cout << "==> TMVAnalysis is done!" << std::endl;

  delete factory;
  return 0;

}
