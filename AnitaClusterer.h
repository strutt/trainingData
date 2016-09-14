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
#include "BaseList.h"
#include "ProgressBar.h"
#include "AntarcticaMapPlotter.h"

#include "assert.h"
#include "AnitaConventions.h"


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
  AnitaPol::AnitaPol_t pol; // polarization
  // Double_t eventPosition[nDim]; // Event position cartesian (m)
  Double_t eventLat; // latitude of position of event
  Double_t eventLon; // longitude of position of event
  Double_t eventAlt; // altitude of posiiton of event

  Double_t distanceToClusterCentroid; // distance (km) from point to centroid
  Double_t distanceToClusterCentroidSecondBest; // distance (km) from point to centroid
  Double_t minusTwoLogLikelihood; // distance in some normalized error units TODO
  Double_t minusTwoLogLikelihoodSecondBest; // distance in some normalized error units TODO

  Double_t deltaThetaDeg; // angular distance to cluster centre
  Double_t deltaPhiDeg; // angular distance to cluster centre
  Double_t thetaDeg; // reconstruction angle
  Double_t phiDeg; // reconstruction angle

  Int_t inCluster; // ID of cluster
  Int_t secondClosestCluster; // ID of cluster

  Int_t isMC;
  Double_t weight;

  Int_t isBase;
  Int_t numEventsInCluster; // number of events in the cluster containing this event
  // Double_t clusterPosition[nDim]; // centroid of cluster cartesian (m)
  Double_t clusterLat; // latitude of centroid of cluster
  Double_t clusterLon; // longitude of centroid of cluster
  Double_t clusterAlt; // altitude of centroid of cluster

  Double_t anitaLat; // latitude of centroid of cluster
  Double_t anitaLon; // longitude of centroid of cluster
  Double_t anitaAlt; // altitude of centroid of cluster

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
    Double_t errorSecondBest;
    Int_t secondClosestCluster;

    AnitaPol::AnitaPol_t pol; // polarization

    Point(Adu5Pat* pat, \
	  Double_t lat=0, Double_t lon=0, Double_t alt=0,\
	  Double_t sigmaTheta = 0.5, Double_t sigmaPhi = 1,\
	  Int_t polIn=AnitaPol::kVertical){

      pol = (AnitaPol::AnitaPol_t) polIn;
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
      error = DBL_MAX;
      inCluster = -1;
      errorSecondBest = DBL_MAX;
      secondClosestCluster = -1;

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
      error = DBL_MAX;
      inCluster = -1;
      errorSecondBest = DBL_MAX;
      secondClosestCluster = -1;

      for(int dim=0; dim < nDim; dim++){
      	centre[dim] = 0;
      }

      }
    virtual ~Point(){ ;}
    ClassDef(Point, 1)
  };



  class MCPoint : public Point{
  public:
    Double_t weight;
    explicit MCPoint() : Point(){ weight = 0;}
    explicit MCPoint(Adu5Pat* pat,					\
		     Double_t lat=0, Double_t lon=0, Double_t alt=0,	\
		     Double_t sigmaTheta = 0.5, Double_t sigmaPhi = 1,	\
		     Int_t polIn=AnitaPol::kVertical,
		     Double_t theWeight=1) : Point(pat,	lat, lon, alt, sigmaTheta, sigmaPhi, polIn){
      weight = theWeight;
    }
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
      maxDist = 0;
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
      maxDist = 0;
    }


    explicit Cluster(const BaseList::base& base) {
      latitude = base.latitude;
      longitude = base.longitude;
      altitude = base.altitude;

      AnitaGeomTool* geom = AnitaGeomTool::Instance();
      geom->getCartesianCoords(latitude, longitude, altitude, centre);
      numEvents = 0;
      totalError = 0;
      maxDist = 0;
    }

    Double_t centre[nDim];
    Double_t latitude;
    Double_t longitude;
    Double_t altitude;
    Int_t numEvents;
    Double_t totalError;
    Double_t maxDist;

    virtual ~Cluster(){ ;}
    ClassDef(Cluster, 1)
  };



  AnitaClusterer(Int_t nClusters, Int_t numIterations, Int_t approxNumPoints=0);
  size_t addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol);
  size_t addMCPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight);
  void kMeansCluster(Int_t iterationsPerCout=0);

  TGraph* makeClusterSummaryTGraph(Int_t clusterInd);
  TTree* makeClusterSummaryTree(TFile* fOut);


  void initializeBaseList();

  Int_t getNumClusters(){
    return numClusters;
  }

  // Int_t histogramUnclusteredEvents(double& minLat, double& maxLat, double& minLon, double& maxLon);
  Int_t histogramUnclusteredEvents(Int_t& globalMaxBin);
  void recursivelyAddClusters(Int_t minBinContent);
  void assignMCPointsToClusters();

  void mergeClusters();
  double llCut;
  double maxDistCluster;

private:

  void assignSinglePointToCloserCluster(Int_t pointInd, Int_t isMC, Int_t clusterInd);

  Cluster seedCluster(Point& point);

  std::vector<TH2D*> hUnclustereds;

  Int_t numIter;
  Int_t numClusters;

  Int_t numCallsToRecursive;

  std::vector<Point> points; // only variables relevant to clustering
  std::vector<MCPoint> mcpoints; // only variables relevant to clustering
  std::vector<Adu5Pat*> pats;
  std::vector<Adu5Pat*> mcpats;
  std::vector<Cluster> clusters; // only variables relevant to clustering
  std::vector<Int_t> seedPoints; // which points seeded which clusters

  std::vector<UInt_t> eventNumbers; // keep track of these separately
  std::vector<UInt_t> mceventNumbers; // keep track of these separately
  std::vector<Int_t> runs; // keep track of these separately
  std::vector<Int_t> mcruns; // keep track of these separately

  Bool_t initialized;
  Double_t minimalImprovement;


  // ok since things got very slow all of a sudden I need to speed them up
  // since my interupts in lldb seem to fall inside the UsefulAdu5Pat functions
  // I will cache the values.
  // For each point I need the value for each cluster... use pair<pointInd, clusterInd>
  // the figure of merit is the loglikelihood, with a distance d.

  AntarcticaMapPlotter* amp;
  RampdemReader* surfaceModel;
  std::vector<Int_t> pointsInCluster;
  Int_t theMinCluster;
  Double_t sumOfAngularErrorsFromLatLon(const Double_t* latLon);
  std::vector<Int_t> ampBinNumbers;
  std::vector<Int_t> ampBinNumbers2;



};


#endif
