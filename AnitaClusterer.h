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
#include "TRandom3.h"

#define nDim 3
//#include "KMeansRex/src/KMeansRexCore.h"

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
      for(int dim=0; dim < nDim; dim++){
	// m to km for easy output paring.
	centre[dim]*=1e-3; // just remember to undo this!
      }
      inCluster = -1;      
      error = 0;      
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
    }

    Cluster(Point seedPoint){
      for(int dim=0; dim < nDim; dim++){
	centre[dim] = seedPoint.centre[dim];
      }
      numEvents = 1;
      totalError = 0;
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
  size_t addPoint(Double_t latitude, Double_t longitude, Double_t altitude, UInt_t eventNumber);
  void kMeansCluster(Int_t iterationsPerCout=0);

  Int_t seedVal;
  
private:

  void kMeansPlusPlusInitialize(Int_t seed = 2016); 
  
  void updateClusterCentres();
  void assignPointsToClosestCluster();
  Double_t assignErrorValues();  
  
  Int_t numIter;
  Int_t numClusters;

  const double minimalImprovement = 0.01;
  std::vector<Point> points; // only variables relevant to clustering
  std::vector<Cluster> clusters; // only variables relevant to clustering

  std::vector<UInt_t> eventNumbers; // keep track of these separately 
  Bool_t initialized;
  
};


#endif
