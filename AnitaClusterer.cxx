#include "AnitaClusterer.h"





// utility function hopefully this one gets inlined
inline Double_t getDistSq(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster.centre[dim] - point.centre[dim];
    d2 += d*d;
  }
  return d2;
}




AnitaClusterer::AnitaClusterer(Int_t nClusters, Int_t numIterations, Int_t approxNumPoints){

  numIter = numIterations;
  numClusters = nClusters;
  clusters.reserve(numClusters);
  points.reserve(approxNumPoints);
  initialized = false;
}



size_t AnitaClusterer::addPoint(Double_t latitude, Double_t longitude, Double_t altitude, UInt_t eventNumber){
  points.push_back(Point(latitude, longitude, altitude));
  eventNumbers.push_back(eventNumber);
  return points.size();
}



void AnitaClusterer::kMeansPlusPlusInitialize(Int_t seed){

  // std::cerr << __PRETTY_FUNCTION__ << std::endl;
  
  // https://en.wikipedia.org/wiki/K-means%2B%2B


  // what I'm going to do (from wikipedia):
  // 1. Choose one center uniformly at random from among the data points.
  // 2. For each data point x, compute D(x), the distance between x and the nearest center that has already been chosen.
  // 3. Choose one new data point at random as a new center, using a weighted probability distribution where a point x is chosen with probability proportional to D(x)2.
  // 4. Repeat Steps 2 and 3 until k centers have been chosen.
  // 5. Now that the initial centers have been chosen, proceed using standard k-means clustering.


  // probably want something determinisitic by default so use non-zero seed
  // can specify seed by calling kMeansPlusPlusInitialize explicitly
  TRandom3 rnd(seed);
  

  Int_t numPoints = (Int_t) points.size();
  
  // start cluster with the first point
  Int_t firstPoint = floor(rnd.Uniform(numPoints));
  clusters.push_back(Cluster(points[firstPoint]));

  std::vector<Double_t> dSqIntervals(numPoints, 0);

  // keep doing steps 2 + 3 until we get to the required number of clusters
  while(clusters.size() < (UInt_t) numClusters){

    // calculate distance squared from each point to closest cluster (minDSq)
    for(int pointInd = 0; pointInd < numPoints; pointInd++){
      Double_t minDSq = DBL_MAX;
      
      for(Int_t clusterInd = 0; clusterInd < (Int_t) clusters.size(); clusterInd++){
	Double_t dSq = getDistSq(points[pointInd], clusters[clusterInd]);

	if(dSq < minDSq){
	  minDSq = dSq;
	}
      }


      // want to pick with probability proportional to the minDSq
      // so we're going to roll a uniform dice between 0 and the sum over all minDSq
      // and pick the point that falls in that interval
      // thats probably not very clear and there may be a better way to do this...
      if(pointInd==0){
	dSqIntervals[pointInd] = minDSq;
      }
      else{
	dSqIntervals[pointInd] = dSqIntervals[pointInd-1] + minDSq;
      }
    }

    // Roll uniform dice between 0 and the sum over all dSq
    Double_t val = rnd.Uniform(dSqIntervals[numPoints-1]);


    // pick point where dice roll falls in the dSq intervals
    // this should satisfy point 3 in the explanation above.
    Int_t nextPoint = numPoints - 1;    
    for(int pointInd = 0; pointInd < numPoints-1; pointInd++){
      if(val >= dSqIntervals[pointInd] && val < dSqIntervals[pointInd+1]){
	nextPoint = pointInd;
	break;
      }
    }

    clusters.push_back(Cluster(points[nextPoint]));
    
    // std::cout << "initializing " << clusters.size() << std::endl;
    // std::cout << nextPoint << "\t" << clusters.back().centre[0] << std::endl;
  }

  assignPointsToClosestCluster();

  std::cout << "AnitaClusterer is initialized with " << clusters.size() << " clusters" << std::endl;
  
  initialized = true;
}





void AnitaClusterer::kMeansCluster(Int_t iterationsPerCout){

  // Apparently kMeansPlusPlus is the bees knees
  if(!initialized){
    kMeansPlusPlusInitialize();
  }

  
  
  Double_t lastClusterError = DBL_MAX;
  
  // -----------------------------------------------------------
  // until completion or reach minimal improvement
  for(int i=0; i < numIter; i++){



    updateClusterCentres();

    assignPointsToClosestCluster();
    

    Double_t sumClusterErrors = assignErrorValues();
    
    // -----------------------------------------------------------
    // report status maybe
    if(iterationsPerCout > 0 && (i%iterationsPerCout)==0){
      std::cout << "Iteration " << i << ", total error = " << sumClusterErrors << std::endl;
    }


    if(sumClusterErrors > lastClusterError){
      std::cerr << "Warning in " << __FILE__ << " your sum over errors has increased!" << std::endl;
      std::cerr << "This should be impossible for a well behaved k-means minimization" << std::endl;
      // for more info see this lovey webpage
      // http://nlp.stanford.edu/IR-book/html/htmledition/k-means-1.html
    }
    
    if(TMath::Abs(sumClusterErrors - lastClusterError) < minimalImprovement){
      break;
    }
    
    lastClusterError = sumClusterErrors;
  }
}


void AnitaClusterer::assignPointsToClosestCluster(){
  // -----------------------------------------------------------
  // assign points to closest cluster  
  for(int pointInd=0; pointInd < (Int_t)points.size(); pointInd++){

    // square of distance, to save SQRT
    Double_t minD2 = DBL_MAX;
    Int_t minClusterInd = -1;
    for(int clusterInd=0; clusterInd < numClusters; clusterInd++){

      Double_t d2 = getDistSq(points[pointInd], clusters[clusterInd]);
      
      if(d2 < minD2){
	minD2 = d2;
	minClusterInd = clusterInd;
      }
    }
    points[pointInd].inCluster = minClusterInd;
  }    

}


void AnitaClusterer::updateClusterCentres(){
  // -----------------------------------------------------------
  // update cluster centres...
  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    for(int dim=0; dim < nDim; dim++){
      clusters[clusterInd].centre[dim] = 0;
    }
    clusters[clusterInd].numEvents = 0;
  }
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    Int_t clusterInd = points[pointInd].inCluster;
    for(int dim=0; dim < nDim; dim++){
      clusters[clusterInd].centre[dim] += points[pointInd].centre[dim];
    }
    clusters[clusterInd].numEvents++;
  }
  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    for(int dim=0; dim < nDim; dim++){
      clusters[clusterInd].centre[dim]/=clusters[clusterInd].numEvents;
    }
    clusters[clusterInd].totalError = 0;
  }
}


Double_t AnitaClusterer::assignErrorValues(){
  // -----------------------------------------------------------
  // assign error values
  Double_t sumClusterErrors = 0;    
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    Int_t clusterInd = points[pointInd].inCluster;
    Double_t d2 = getDistSq(points[pointInd], clusters[clusterInd]);
    points[pointInd].error = d2;
    clusters[clusterInd].totalError += d2;
    sumClusterErrors += d2;
  }
  return sumClusterErrors;
}
