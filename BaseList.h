#ifndef COMNAP_H
#define COMNAP_H

#include "TString.h"
#include <vector>

namespace COMNAP2014{

  class base{
  public:
    base(const TString& theName, double lat, double lon, double alt=0){
      name = theName;
      latitude = lat;
      longitude = lon;
      altitude = alt;
    }
    TString name;
    double latitude;
    double longitude;
    double altitude;
  };

  const base& getBase(UInt_t i);
  size_t getNumBases();

};


#endif // COMNAP_H
