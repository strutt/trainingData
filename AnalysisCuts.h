/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description: 
             My ANITA-3 analysis cuts
***********************************************************************************************************/

#ifndef ANALYSISCUTS_H
#define ANALYSISCUTS_H

#include "AnitaConventions.h"

#include "CrossCorrelator.h"
#include "RootTools.h"

/**
 * @namespace Analysis Cuts 
 * @brief Put my analysis cuts in one place for multiple scripts.
 * 
 *
 */
namespace AnalysisCuts{

  
  typedef enum EStatus {
    kPass = 0,
    kFail = 1
  } Status_t; // to make the scripts easy to read with comparisons to pass/fail

  // CUT FLOW
  // Aiming for combined reduction of factor O(1e9) for thermal noise events
  // How many signal events will be left?
  // here we go!

  const double ratioCutHigh = 2.8;
  const double ratioCutLow = 1.14;  
  // Step 1: cut self triggered blasts (this is almost a data quality cut)
  Status_t applyBottomToTopRingPeakToPeakRatioCut(AnitaPol::AnitaPol_t pol, Double_t* peakToPeak, Double_t& maxRatio);
  
  

  const int maxAbsDeltaPhiSect = 2;
  Status_t L3TriggerDirectionCut(AnitaPol::AnitaPol_t pol, RawAnitaHeader* header, Double_t recoPhiDeg, Int_t& deltaPhiSect);

  
  const double deltaSolarPhiDegCut = 20;
  Status_t applySunPointingCut(Double_t deltaSolarPhiDeg);



};



#endif //ANALYSISCUTS
