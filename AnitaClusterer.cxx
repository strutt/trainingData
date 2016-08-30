#include "AnitaClusterer.h"



// for debugging
inline void prettyPrint(const int n, const double* array){
  for(int i=0; i < n; i++){
    std::cerr << array[i];
    if(i < n - 1){
      std::cerr << ", ";
    }       
  }
  std::cerr << std::endl;  
}



inline void prettyPrintConvert(const int n, const double* array){
  Double_t lat, lon, alt;
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  std::vector<Double_t> temp(array, array+n*sizeof(double));
  geom->getLatLonAltFromCartesian(&temp[0], lat, lon, alt);
  std::cerr << lat << "\t" << lon << "\t" << alt << std::endl;
}



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
  eventNumbers.reserve(approxNumPoints);
  runs.reserve(approxNumPoints);
  initialized = false;
}







size_t AnitaClusterer::addPoint(Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber){
  points.push_back(Point(latitude, longitude, altitude));
  eventNumbers.push_back(eventNumber);
  runs.push_back(run);  
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
  clusters.push_back(Cluster(points.at(firstPoint)));
  seedPoints.push_back(firstPoint);

  std::vector<Double_t> dSqIntervals(numPoints, 0);

  // keep doing steps 2 + 3 until we get to the required number of clusters
  while(clusters.size() < (UInt_t) numClusters){

    // calculate distance squared from each point to closest cluster (minDSq)
    for(int pointInd = 0; pointInd < numPoints; pointInd++){
      Double_t minDSq = DBL_MAX;
      
      for(Int_t clusterInd = 0; clusterInd < (Int_t) clusters.size(); clusterInd++){
	Double_t dSq = getDistSq(points.at(pointInd), clusters.at(clusterInd));

	if(dSq < minDSq){
	  minDSq = dSq;
	}
      }


      // want to pick with probability proportional to the minDSq
      // so we're going to roll a uniform dice between 0 and the sum over all minDSq
      // and pick the point that falls in that interval
      // thats probably not very clear and there may be a better way to do this...
      if(pointInd==0){
	dSqIntervals.at(pointInd) = minDSq;
      }
      else{
	dSqIntervals.at(pointInd) = dSqIntervals[pointInd-1] + minDSq;
      }
    }

    // Roll uniform dice between 0 and the sum over all dSq
    Double_t val = rnd.Uniform(dSqIntervals[numPoints-1]);

    // pick point where dice roll falls in the dSq intervals
    // this should satisfy point 3 in the explanation above.
    Int_t nextPoint = numPoints - 1;
    for(int pointInd = 0; pointInd < numPoints-1; pointInd++){
      if(val >= dSqIntervals.at(pointInd) && val < dSqIntervals[pointInd+1]){
	nextPoint = pointInd;
	break;
      }
    }

    // get index before push back
    const Int_t clusterInd = clusters.size();
    seedPoints.push_back(nextPoint);

    // std::cout << "initializing " << clusterInd << "\t" << nextPoint << "\t"
    // 	      << clusters.size() << "\t" << seedPoints.size() << std::endl;
    

    // initialize the cluster
    Cluster cluster(points.at(nextPoint));
    clusters.push_back(cluster);

    points.at(nextPoint).inCluster = clusterInd;

    prettyPrintConvert(nDim, clusters.at(clusterInd).centre);
    prettyPrintConvert(nDim, points.at(nextPoint).centre);
    
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

    // square of distance, to save on doing Sqrt operation
    Double_t minD2 = DBL_MAX;
    Int_t minClusterInd = -1;
    for(int clusterInd=0; clusterInd < numClusters; clusterInd++){

      // if(pointInd==seedPoints.at(141) && clusterInd==141){
      // 	std::cerr << "in assign..." << std::endl;
      // 	prettyPrint(nDim, clusters.at(clusterInd+1).centre);
      // 	prettyPrintConvert(nDim, clusters.at(clusterInd+1).centre);
      // }
      // if(pointInd==seedPoints.at(142) && clusterInd==142){
      // 	std::cerr << "in assign..." << std::endl;
      // 	prettyPrint(nDim, clusters.at(clusterInd).centre);
      // 	prettyPrintConvert(nDim, clusters.at(clusterInd).centre);
      // }
      // if(pointInd==seedPoints.at(143) && clusterInd==143){
      // 	std::cerr << "in assign..." << std::endl;
      // 	prettyPrint(nDim, clusters.at(clusterInd-1).centre);
      // 	prettyPrintConvert(nDim, clusters.at(clusterInd-1).centre);
      // }

      Double_t d2 = getDistSq(points.at(pointInd), clusters.at(clusterInd));
      
      if(d2 < minD2){
	minD2 = d2;
	minClusterInd = clusterInd;
      }      
    }


    // auto x = std::find(seedPoints.begin(), seedPoints.end(), pointInd);
    // if(x != seedPoints.end()) {    

    //   auto index = std::distance(seedPoints.begin(), x);

    //   std::cerr<< "In assign: at the seed of cluster " <<  index << std::endl;
    //   std::cerr << "It got assigned to cluster " << minClusterInd << std::endl;
    //   std::cerr << "It has minD2 " << minD2 << std::endl;
    //   std::cerr << "Compared to the seed with d2 " << getDistSq(points.at(pointInd), clusters.at(142))
    // 		<< std::endl;
    //   std::cerr << "The point is at " << std::endl;
    //   prettyPrintConvert(nDim, points.at(pointInd).centre);
    //   std::cerr << "142 is at " << std::endl;      
    //   prettyPrintConvert(nDim, clusters.at(142).centre);
    //   std::cerr << minClusterInd << " is at " << std::endl;
    //   prettyPrintConvert(nDim, clusters.at(minClusterInd).centre);      
    // }

    // assert(minClusterInd > -1 && minClusterInd < numClusters);
    points.at(pointInd).inCluster = minClusterInd;
  }    
}













void AnitaClusterer::updateClusterCentres(){
  // -----------------------------------------------------------
  // update cluster centres...
  // std::cerr << "in update (before wipe): " << std::endl;
  // prettyPrint(nDim, clusters.at(142).centre);
  // prettyPrintConvert(nDim, clusters.at(142).centre);      

  
  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    for(int dim=0; dim < nDim; dim++){
      clusters.at(clusterInd).centre[dim] = 0;
    }
    clusters.at(clusterInd).numEvents = 0;
  }
  // Int_t num142 = 0;
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    Int_t clusterInd = points.at(pointInd).inCluster;
    // if(clusterInd==142){
    //   std::cerr << "found " << num142 << ":" << std::endl;
    //   prettyPrintConvert(nDim, points.at(pointInd).centre);
    //   num142++;
    // }
    
    for(int dim=0; dim < nDim; dim++){
      clusters.at(clusterInd).centre[dim] += points.at(pointInd).centre[dim];
    }
    clusters.at(clusterInd).numEvents++;
  }
  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    for(int dim=0; dim < nDim; dim++){
      if(clusters.at(clusterInd).numEvents > 0){
	clusters.at(clusterInd).centre[dim]/=clusters.at(clusterInd).numEvents;
      }
      else{
	// const Point& point = points[seedPoints.at(clusterInd)];
	// const Cluster& cluster = clusters.at(clusterInd);
	std::cerr << "Warning in " << __FILE__ << ", cluster with 0 events" << std::endl;
	// std::cerr << "The cluster has index " << clusterInd << " and centre:" << std::endl;
	// prettyPrint(nDim, cluster.centre);
	// prettyPrintConvert(nDim, cluster.centre);
	// std::cerr << "The seed point has index " << seedPoints.at(clusterInd) << " and centre:" << std::endl;
	// prettyPrint(nDim, point.centre);

	// std::cerr << "Also, the seed point now thinks its closest cluster is " << point.inCluster
	// 	  << ", which is located at " << std::endl;
	// prettyPrint(nDim, clusters[point.inCluster].centre);
      }
    }
    clusters.at(clusterInd).totalError = 0;
  }
}














Double_t AnitaClusterer::assignErrorValues(){
  // -----------------------------------------------------------
  // assign error values
  Double_t sumClusterErrors = 0;    
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    Int_t clusterInd = points.at(pointInd).inCluster;
    Double_t d2 = getDistSq(points.at(pointInd), clusters.at(clusterInd));
    points.at(pointInd).error = d2;
    clusters.at(clusterInd).totalError += d2;
    sumClusterErrors += d2;
  }
  return sumClusterErrors;
}














TGraph* AnitaClusterer::makeClusterSummaryTGraph(Int_t clusterInd){

  TGraph* gr = NULL;
  if(clusterInd >= 0 && clusterInd < numClusters){

    TString name  = TString::Format("grCluster%d", clusterInd);
    TString title  = TString::Format("Cluster %d; Longitude (Degrees); Latitude (Degrees)", clusterInd);
    gr = new TGraph();
    gr->SetName(name);
    gr->SetTitle(title);    
    

    AnitaGeomTool* geom = AnitaGeomTool::Instance();

    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      if(points.at(pointInd).inCluster==clusterInd){

	Double_t lat, lon, alt;
	geom->getLatLonAltFromCartesian(points.at(pointInd).centre, lat, lon, alt);
	gr->SetPoint(gr->GetN(), lon, lat);
      }
    }
  }
  return gr;
}

















TTree* AnitaClusterer::makeClusterSummaryTree(TFile* fOut){

  fOut->cd();

  TTree* clusterTree = new TTree("clusterTree", "Tree of clustered ANITA events");

  ClusteredAnitaEvent* clusteredEvent;
  clusterTree->Branch("clusteredEvent", &clusteredEvent);

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  
  for(Int_t pointInd=0; pointInd < (Int_t)points.size(); pointInd++){

    const Point& point = points.at(pointInd);
    
    clusteredEvent = new ClusteredAnitaEvent();
    clusteredEvent->eventNumber = eventNumbers.at(pointInd);
    clusteredEvent->run = runs.at(pointInd);

    // convert from km back to m for conversion to lat/lon/alt
    for(int dim=0; dim < nDim; dim++){
      clusteredEvent->eventPosition[dim] = point.centre[dim];
    }
    geom->getLatLonAltFromCartesian(clusteredEvent->eventPosition,\
				    clusteredEvent->eventLat,\
				    clusteredEvent->eventLon,\
				    clusteredEvent->eventAlt);

    Int_t clusterInd = points.at(pointInd).inCluster;
    clusteredEvent->inCluster = clusterInd;

    const Cluster& cluster = clusters.at(clusterInd);
    
    for(int dim=0; dim < nDim; dim++){
      clusteredEvent->clusterPosition[dim] = cluster.centre[dim];
    }
    geom->getLatLonAltFromCartesian(clusteredEvent->eventPosition,\
				    clusteredEvent->eventLat,\
				    clusteredEvent->eventLon,\
				    clusteredEvent->eventAlt);

    clusteredEvent->distanceToClusterCentroid = TMath::Sqrt(getDistSq(point, cluster));
    clusteredEvent->errorToClusteredCentroid = points.at(pointInd).error;

    clusteredEvent->numClusters = numClusters;
    clusteredEvent->numIterations = numIter;
    
    clusterTree->Fill();
    delete clusteredEvent;
  }
  
  return clusterTree;
  
}


