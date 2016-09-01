/* -*- C++ -*-.*********************************************************************************************
 Author: Ben Strutt
 Email: b.strutt.12@ucl.ac.uk

 Description:
             A class to cluster. What were you expecting?
	     It does a kmeans++ cluster algorithm on the positions you feed in.
***********************************************************************************************************/

#ifndef ANITACLUSTERER_H
#define ANITACLUSTERER_H


#include "AnitaGeomTool.h"
#include "UsefulAdu5Pat.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "RootTools.h"
#include "COMNAP.h"

#include "assert.h"


#define nDim 3
//#include "KMeansRex/src/KMeansRexCore.h"





//--------------------------------------------------------------------------------------------------------
/**
 * @class ClusteredAnitaEvent
 * @ Per-event output, perhaps a little gratuitous
 */

class ClusteredAnitaEvent{
public:
  UInt_t eventNumber;
  Int_t run;
  // Double_t eventPosition[nDim]; // Event position cartesian (m)
  Double_t eventLat; // latitude of position of event
  Double_t eventLon; // longitude of position of event
  Double_t eventAlt; // altitude of posiiton of event

  Double_t distanceToClusterCentroid; // distance (km) from point to centroid
  Double_t errorToClusteredCentroid; // distance in some normalized error units TODO

  Int_t inCluster; // ID of cluster
  Int_t numEventsInCluster; // number of events in the cluster containing this event
  // Double_t clusterPosition[nDim]; // centroid of cluster cartesian (m)
  Double_t clusterLat; // latitude of centroid of cluster
  Double_t clusterLon; // longitude of centroid of cluster
  Double_t clusterAlt; // altitude of centroid of cluster

  Int_t numClusters; // Total number of clusters (maybe gratuitous)
  Int_t numIterations; // number of loop iterations (maybe gratuitous)

};



/**
 * @class AnitaCluster
 * @brief A class to cluster locations on the Antarctic ice as reconstructed by ANITA.
 *
*/
class AnitaClusterer{

public:



  //--------------------------------------------------------------------------------------------------------
  // Classes declared inside this class
  //--------------------------------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------------------------------
  /**
   * @class Point
   * @ Position on the contitent in Cartesian Coordinates
   */
  class Point{
  public:
    Double_t centre[nDim];
    Double_t latitude;
    Double_t longitude;
    Double_t altitude;
    Double_t thetaDeg;
    Double_t phiDeg;
    Double_t dTheta; // theta distance to cluster
    Double_t dPhi; // phi distance to cluster
    Double_t sigmaThetaDeg; // resolution associated with this snr?
    Double_t sigmaPhiDeg; // resolution associated with this snr?
    Double_t error; //
    Int_t inCluster; // which cluster am I associated with?

    Point(Adu5Pat* pat, \
	  Double_t lat=0, Double_t lon=0, Double_t alt=0,\
	  Double_t sigmaTheta = 0.5, Double_t sigmaPhi = 1){

      UsefulAdu5Pat usefulPat(pat);
      latitude = lat;
      longitude = lon;
      altitude = alt;
      usefulPat.getThetaAndPhiWave(longitude, latitude, altitude, thetaDeg, phiDeg);


      // convert to degrees
      thetaDeg = -1*thetaDeg*TMath::RadToDeg();
      phiDeg = phiDeg*TMath::RadToDeg();
      // dTheta = 0;
      // dPhi = 0;
      sigmaThetaDeg = sigmaTheta;
      sigmaPhiDeg = sigmaPhi;
      error = 0;
      inCluster = -1;
      AnitaGeomTool* geom = AnitaGeomTool::Instance();
      geom->getCartesianCoords(lat, lon, alt, centre);
      // for(int dim=0; dim < nDim; dim++){
      // 	centre[dim] = 0;
      // }
    }
    Point(){
      latitude = 0;
      longitude = 0;
      altitude = 0;
      thetaDeg = -9999;
      phiDeg = -9999;

      // convert to degrees
      thetaDeg = -9999;
      phiDeg = -9999;
      // dTheta = 0;
      // dPhi = 0;
      sigmaThetaDeg = -9999;
      sigmaPhiDeg = -9999;
      error = 0;
      inCluster = -1;
      for(int dim=0; dim < nDim; dim++){
      	centre[dim] = 0;
      }

      }
    virtual ~Point(){ ;}
    ClassDef(Point, 1)
  };


  //--------------------------------------------------------------------------------------------------------
  /**
   * @class Cluster
   * @ Where the events are clustered
   */
  class Cluster{
  public:
    Cluster(){
      for(int dim=0; dim < nDim; dim++){
      	centre[dim] = 0;
      }
      numEvents = 0;
      totalError = 0;
      latitude = 0;
      longitude = 0;
      altitude = 0;
    }

    explicit Cluster(const Point& seedPoint) {
      latitude = seedPoint.latitude;
      longitude = seedPoint.longitude;
      longitude = seedPoint.longitude;
      for(int dim=0; dim < nDim; dim++){
      	centre[dim] = seedPoint.centre[dim];
      }
      numEvents = 1; // since the seed point will be in this cluster
      totalError = 0;
    }


    explicit Cluster(const COMNAP2014::base& base) {
      latitude = base.latitude;
      longitude = base.longitude;
      longitude = base.longitude;

      AnitaGeomTool* geom = AnitaGeomTool::Instance();
      geom->getCartesianCoords(latitude, longitude, altitude, centre);
      numEvents = 0;
      totalError = 0;
    }

    Double_t centre[nDim];
    Double_t latitude;
    Double_t longitude;
    Double_t altitude;
    Int_t numEvents;
    Double_t totalError;

    virtual ~Cluster(){ ;}
    ClassDef(Cluster, 1)
  };



  AnitaClusterer(Int_t nClusters, Int_t numIterations, Int_t approxNumPoints=0);
  size_t addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg);
  void kMeansCluster(Int_t iterationsPerCout=0);

  TGraph* makeClusterSummaryTGraph(Int_t clusterInd);
  TTree* makeClusterSummaryTree(TFile* fOut);


  void initializeCOMNAP();

private:

  void kMeansPlusPlusInitialize(Int_t seed = 2016);
  void maxDistanceInitialize(Int_t seed = 2016);

  void updateClusterCentres();
  void assignPointsToClosestCluster();
  Double_t assignErrorValues();
  Cluster seedCluster(Point& point);

  Int_t numIter;
  Int_t numClusters;

  const double minimalImprovement = 0.01;
  std::vector<Point> points; // only variables relevant to clustering
  std::vector<Adu5Pat*> pats;
  std::vector<Cluster> clusters; // only variables relevant to clustering
  std::vector<Int_t> seedPoints; // which points seeded which clusters

  std::vector<UInt_t> eventNumbers; // keep track of these separately
  std::vector<Int_t> runs; // keep track of these separately
  std::vector<Double_t> deltaTheta;
  std::vector<Double_t> deltaPhi;
  Bool_t initialized;


  RampdemReader* surfaceModel;
  std::vector<Int_t> pointsInCluster;
  Int_t theMinCluster;
  Double_t sumOfAngularErrorsFromLatLon(const Double_t* latLon);


};


#endif
