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
  Double_t eventPosition[nDim]; // Event position cartesian (m)
  Double_t eventLat; // latitude of position of event
  Double_t eventLon; // longitude of position of event
  Double_t eventAlt; // altitude of posiiton of event

  Double_t distanceToClusterCentroid; // distance (km) from point to centroid
  Double_t errorToClusteredCentroid; // distance in some normalized error units TODO
 
  Int_t inCluster; // ID of cluster
  Int_t numEventsInCluster; // number of events in the cluster containing this event
  Double_t clusterPosition[nDim]; // centroid of cluster cartesian (m)
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
    Double_t error;
    Int_t inCluster;
    
    Point(Double_t latitude=0, Double_t longitude=0, Double_t altitude=0){
      AnitaGeomTool* geom = AnitaGeomTool::Instance();
      geom->getCartesianCoords(latitude,longitude,altitude, centre);
      inCluster = -1;
      error = 0;
    }
    virtual ~Point(){ ;}
    ClassDef(Point, 1)    
  };


  // //--------------------------------------------------------------------------------------------------------
  // /** 
  //  * @class UncertainPoint
  //  * @ Position on the contitent with an error ellipse in Cartesian Coordinates
  //  */
  // class UncertainPoint : public Point{
  // public:
  //   Double_t posThetaError[nDim];
  //   Double_t negThetaError[nDim];
  //   Double_t posPhiError[nDim];
  //   Double_t negPhiError[nDim];
  //   virtual ~UncertainPoint(){ ;}
  //   ClassDef(Point, 1);
  // }; 

  

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
    }

    explicit Cluster(const Point& seedPoint) {
      numEvents = 1;
      totalError = 0;      
      for(int dim=0; dim < nDim; dim++){
	centre[dim] = seedPoint.centre[dim];
	if(centre[dim] == 0){
	  std::cerr << "???????????" << std::endl;
	}
      }
    }
    
    Double_t centre[nDim];
    Int_t numEvents;
    // Mean "distance" of all points to cluster centre (a cluster figure of merit)
    // this is what the k-means++ algorithm should minimize
    Double_t totalError;
    
    virtual ~Cluster(){ ;}    
    ClassDef(Cluster, 1)
  };


  
  AnitaClusterer(Int_t nClusters, Int_t numIterations, Int_t approxNumPoints=0);
  size_t addPoint(UsefulAdu5Pat usefulPat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber);
  void kMeansCluster(Int_t iterationsPerCout=0);

  TGraph* makeClusterSummaryTGraph(Int_t clusterInd);
  TTree* makeClusterSummaryTree(TFile* fOut);
  
  
private:

  void kMeansPlusPlusInitialize(Int_t seed = 2016); 
  
  void updateClusterCentres();
  void assignPointsToClosestCluster();
  Double_t assignErrorValues();
  Cluster seedCluster(Point& point);  
  
  Int_t numIter;
  Int_t numClusters;

  const double minimalImprovement = 0.01;
  std::vector<Point> points; // only variables relevant to clustering
  std::vector<UsefulAdu5Pat> pats;
  std::vector<Cluster> clusters; // only variables relevant to clustering
  std::vector<Int_t> seedPoints; // which points seeded which clusters

  std::vector<UInt_t> eventNumbers; // keep track of these separately
  std::vector<Int_t> runs; // keep track of these separately
  std::vector<Double_t> deltaTheta;
  std::vector<Double_t> deltaPhi;
  Bool_t initialized;
  
};


#endif
