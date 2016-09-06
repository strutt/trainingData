#include "AnitaClusterer.h"


AnitaClusterer::AnitaClusterer(Int_t nClusters, Int_t numIterations, Int_t approxNumPoints){

  numIter = numIterations;
  numClusters = nClusters;
  clusters.reserve(numClusters);

  points.reserve(approxNumPoints);
  pats.reserve(approxNumPoints);
  eventNumbers.reserve(approxNumPoints);
  runs.reserve(approxNumPoints);

  initialized = false;
  llCut = 2000;
  logLikelihoodCut = 100;
  numCallsToRecursive = 0;
  minimalImprovement  = 0.01;
  // can use this to move cluster around surface, 3D positions correctly becomes 2D problem
  surfaceModel = RampdemReader::Instance();
  amp = new AntarcticaMapPlotter();
}



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


// utility function hopefully this one gets inlined
inline Double_t getDistSq(AnitaClusterer::Cluster& cluster1, const AnitaClusterer::Cluster& cluster2){
  Double_t d2=0;
  for(int dim=0; dim < nDim; dim++){
    Double_t d = cluster1.centre[dim] - cluster2.centre[dim];
    d2 += d*d;
  }
  return d2;
}


// utility function hopefully this one gets inlined
// inline Double_t getAngDistSq(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, UsefulAdu5Pat& usefulPat){
// inline Double_t getAngDistSq(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, const Adu5Pat* pat){

//   UsefulAdu5Pat usefulPat(pat);
//   Double_t thetaWave, phiWave;
//   usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude,cluster.altitude,
// 			       thetaWave,phiWave);
//   Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
//   Double_t phiDeg = TMath::RadToDeg()*phiWave;

//   Double_t dThetaSq = (thetaDeg - point.thetaDeg)/point.sigmaThetaDeg;
//   Double_t dPhiSq = RootTools::getDeltaAngleDeg(phiDeg, point.phiDeg)/point.sigmaPhiDeg;

//   return dThetaSq*dThetaSq + dPhiSq*dPhiSq;
// }

inline void getDeltaThetaDegDeltaPhiDegCluster(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, const Adu5Pat* pat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  UsefulAdu5Pat usefulPat(pat);

  Double_t thetaWave, phiWave;
  usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude,cluster.altitude,
			       thetaWave,phiWave);
  Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
  Double_t phiDeg = TMath::RadToDeg()*phiWave;

  deltaThetaDeg = (thetaDeg - point.thetaDeg);
  deltaPhiDeg = RootTools::getDeltaAngleDeg(phiDeg, point.phiDeg)/point.sigmaPhiDeg;
}



inline Double_t getAngDistSq(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, const Adu5Pat* pat){

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegCluster(point, cluster, pat, deltaThetaDeg, deltaPhiDeg);

  // normalized
  Double_t dThetaSq = deltaThetaDeg/point.sigmaThetaDeg;
  Double_t dPhiSq = deltaPhiDeg/point.sigmaPhiDeg;

  Double_t angSq =  dThetaSq*dThetaSq + dPhiSq*dPhiSq;

  return angSq;
}




// Int_t AnitaClusterer::addUnclusteredHistogram(){

// }

Int_t AnitaClusterer::histogramUnclusteredEvents(Int_t& globalMaxBin){

  Int_t addClusterIter =  hUnclustereds.size();
  TString name = TString::Format("hUnclustered%d", addClusterIter);
  // const int nBinsLat = 30;
  // const int nBinsLon = 360;
  const int nBins = 128;
  amp->addHistogram(name, name, nBins, nBins);
  hUnclustereds.push_back(amp->getCurrentHistogram());
  // find unclustered points
  ampBinNumbers.clear(); // doesn't clear memory
  ampBinNumbers.resize((Int_t) points.size()); // so this shouldn't reallocate memory
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    const Point& point = points.at(pointInd);
    if(point.inCluster < 0){
      ampBinNumbers.at(pointInd) = amp->Fill(point.latitude, point.longitude);
      // hUnclustereds.back()->Fill(point.longitude, point.latitude);
    }
  }
  globalMaxBin = hUnclustereds.back()->GetMaximumBin();

  Int_t maxEventsInBin = hUnclustereds.back()->GetBinContent(globalMaxBin);

  return maxEventsInBin;


}



void AnitaClusterer::recursivelyAddClusters(Int_t minBinContent){
  Int_t globalMaxBin;
  Int_t maxBinVal = histogramUnclusteredEvents(globalMaxBin);

  // std::cout << minLat << "\t" << maxLat << "\t" << minLon << "\t" << maxLon << std::endl;


  Int_t counter=0;
  Cluster cluster;
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    const Point& point = points.at(pointInd);
    if(point.inCluster < 0){
      if(ampBinNumbers.at(pointInd)==globalMaxBin){
	for(int dim=0; dim < nDim; dim++){
	  cluster.centre[dim] += point.centre[dim];
	}
	counter++;
      }
    }
  }

  if(counter!=maxBinVal){
    std::cerr << "Now what?" << std::endl;
  }
  if(counter > minBinContent){

    for(int dim=0; dim < nDim; dim++){
      cluster.centre[dim]/=counter;
    }
    AnitaGeomTool* geom = AnitaGeomTool::Instance();
    geom->getLatLonAltFromCartesian(cluster.centre, cluster.latitude,\
				    cluster.longitude, cluster.altitude);
    clusters.push_back(cluster);
    numClusters = (Int_t) clusters.size();

    assignPointsToClosestCluster();

    std::cout << counter << "\t" << maxBinVal << "\t" << clusters.back().latitude << "\t"
	      << clusters.back().longitude << "\t" << hUnclustereds.back()->Integral() << std::endl;

    recursivelyAddClusters(minBinContent);
  }





}




void AnitaClusterer::assignPointsToClosestCluster(){
  // -----------------------------------------------------------
  // assign points to closest cluster
  // const double llCut = 30;
  // const double llCut = 1e99;
  const double maxDistCluster = 800e3; // try 800km

  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    clusters.at(clusterInd).numEvents = 0;
  }

  for(int pointInd=0; pointInd < (Int_t)points.size(); pointInd++){

    UsefulAdu5Pat usefulPat(pats.at(pointInd));

    Double_t minD2 = DBL_MAX;
    // Int_t minClusterInd = points.at(pointInd).inCluster;
    Int_t minClusterInd = -1;

    if(minClusterInd >= 0 && minClusterInd < numClusters){
      // minD2 = getDistSq(points.at(pointInd), clusters.at(minClusterInd));
      minD2 = getAngDistSq(points.at(pointInd), clusters.at(minClusterInd), pats.at(pointInd));
    }

    for(int clusterInd=0; clusterInd < numClusters; clusterInd++){

      // if the point and cluster are within 1000km

      Double_t distM = usefulPat.getDistanceFromSource(clusters.at(clusterInd).latitude,
						       clusters.at(clusterInd).longitude,
						       clusters.at(clusterInd).altitude);

      Double_t deltaDistM = 0;
      // Double_t deltaDistM = clusters.at(clusterInd).maxDist;

      if(distM - deltaDistM < maxDistCluster){

	// Double_t d2 = getDistSq(points.at(pointInd), clusters.at(clusterInd));
	Double_t d2 = getAngDistSq(points.at(pointInd), clusters.at(clusterInd), pats.at(pointInd));

	if(hUnclustereds.size() >= 7 && points.at(pointInd).inCluster < 0 && clusterInd==numClusters-1){
	  if(distM < maxDistCluster && points.at(pointInd).sigmaPhiDeg == 0){
	    std::cout << pointInd << "\t" << clusterInd << "\t" << distM/1e3 << "\t" << d2 << "\t"
		      << points.at(pointInd).sigmaPhiDeg << "\t" << points.at(pointInd).sigmaThetaDeg
		      << std::endl;
	  }
	}

	if(llCut <= 0 || d2 < llCut){

	  if(d2 < minD2){
	    minD2 = d2;
	    minClusterInd = clusterInd;
	  }
	}
      }
    }



    // assert(minClusterInd > -1 && minClusterInd < numClusters);
    // const double maxSigmaAway = 5*5;
    // if(minD2 >= maxSigmaAway){
      // minClusterInd = -1;
    // }
    points.at(pointInd).secondClosestCluster = points.at(pointInd).inCluster;
    points.at(pointInd).inCluster = minClusterInd;
    if(minClusterInd < 0 ){
      // std::cerr << "Bollocks \t" << pointInd << std::endl;
      points.at(pointInd).error = -9999;

    }
    else{
      clusters.at(minClusterInd).numEvents++;
      // points.at(pointInd).error = minD2;

      points.at(pointInd).errorSecondBest = points.at(pointInd).error;
      points.at(pointInd).error = minD2;



    }
  }

  // for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
  //   Int_t seedPoint = seedPoints.at(clusterInd);
  //   Double_t dSq = getDistSq(points.at(seedPoint),clusters.at(clusterInd));
  //   Int_t inCluster = points.at(seedPoint).inCluster;
  //   if(inCluster!=clusterInd){
  //     std::cerr << "In assign POST: " << dSq << "\t" << clusterInd << "\t" << seedPoint << "\t" << inCluster << "\t" << std::setprecision(15) << getDistSq(points.at(seedPoint), clusters.at(inCluster)) << std::endl;
  //   }
  // }
}













size_t AnitaClusterer::addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol){

  points.push_back(Point(pat, latitude, longitude, altitude, sigmaThetaDeg, sigmaPhiDeg, pol));
  eventNumbers.push_back(eventNumber);
  runs.push_back(run);
  pats.push_back((Adu5Pat*)pat->Clone());

  return points.size();
}








void AnitaClusterer::initializeBaseList(){

  numClusters = (int) BaseList::getNumBases();
  for(int i=0; i < numClusters; i++){
    const BaseList::base& base = BaseList::getBase(i);
    clusters.push_back(Cluster(base));
  }
  initialized = true;
  assignPointsToClosestCluster();
  assignErrorValues();
}




void AnitaClusterer::maxDistanceInitialize(Int_t seed){

  // std::cerr << __PRETTY_FUNCTION__ << std::endl;
  TRandom3 rnd(seed);

  Int_t numPoints = (Int_t) points.size();
  seedPoints.resize(numPoints, 0);
  // std::vector<Double_t> cumulativeDistanceSquaredToNearestClusters(numPoints, 0);
  std::vector<Double_t> distanceToNearestCluster(numPoints, 0);

  Double_t maxDist = -9999;
  Int_t nextPoint = -1;

  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){


    if(clusterInd==0){

      // pick first point from uniform distribution
      Int_t firstPoint = floor(rnd.Uniform(numPoints));

      // map point to this cluster
      points.at(firstPoint).inCluster = clusterInd;

      // map cluster to piont
      seedPoints.at(clusterInd) = firstPoint;

      // create cluster object
      clusters.push_back(Cluster(points.at(firstPoint)));
      clusters.back().numEvents++;

      if(clusters.size() != (size_t)clusterInd+1){
	std::cerr << "?????????????????" << std::endl;
      }

      // assign all points to this cluster
      for(int pointInd=0; pointInd < numPoints; pointInd++){

	// get linear distance to cluster

	distanceToNearestCluster.at(pointInd) = getDistSq(points.at(pointInd), clusters.at(clusterInd));

	points.at(pointInd).inCluster = clusterInd;

	if(distanceToNearestCluster.at(pointInd) > maxDist){
	  maxDist = distanceToNearestCluster.at(pointInd);
	  nextPoint = pointInd;
	}
      }

    }
    else{

      points.at(nextPoint).inCluster = clusterInd;

      // map cluster to piont
      seedPoints.at(clusterInd) = nextPoint;

      // create cluster object
      clusters.push_back(Cluster(points.at(nextPoint)));
      clusters.back().numEvents++;

      if(clusters.size() != (size_t)clusterInd+1){
	std::cerr << "?????????????????" << std::endl;
      }

      maxDist = -9999;
      for(int pointInd=0; pointInd < numPoints; pointInd++){
	// Double_t dSq = getAngDistSq(points.at(pointInd), clusters.at(clusterInd), pats.at(pointInd));
	Double_t dSq = getDistSq(points.at(pointInd), clusters.back());

	// if new cluster centre is closer then assign point to this cluster
	// no need for else, since this tracks minimum of best produced so far...
	if(dSq < distanceToNearestCluster.at(pointInd)){
	  distanceToNearestCluster.at(pointInd) = dSq;
	  points.at(pointInd).inCluster = clusterInd;
	}

	if(distanceToNearestCluster.at(pointInd) > maxDist){
	  maxDist = distanceToNearestCluster.at(pointInd);
	  nextPoint = pointInd;
	}
      }
    }
    std::cerr << clusterInd << "\t" << TMath::Sqrt(maxDist)/1e3 << std::endl;
    // std::cerr << "Initialized " << clusterInd << std::endl;
  }


  assignPointsToClosestCluster();

  std::cout << "AnitaClusterer is initialized with " << clusters.size() << " clusters" << std::endl;

  initialized = true;
}










void AnitaClusterer::mergeClusters(){

  // const double minClusterSeparation = 100e3;
  const int numBases = BaseList::getNumBases();
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusterInd >= numBases){
      Cluster& cluster1 = clusters.at(clusterInd);
      Double_t thisMinSeparation = 1e99;
      Int_t thisOtherCluster = -1;
      for(int clusterInd2=0; clusterInd2 < numClusters; clusterInd2++){
	if(clusterInd==clusterInd2){
	  continue;
	}
	Cluster& cluster2 = clusters.at(clusterInd2);
	Double_t clusterSeparation = TMath::Sqrt(getDistSq(cluster1, cluster2));

	// std::cout << clusterInd << "\t" << clusterInd2 << "\t" << clusterSeparation/1e3 << std::endl;
	// if(clusterSeparation < minClusterSeparation){
	//   std::cerr << clusterInd << "\t" << clusterInd2 << "\t" << clusterSeparation/1e3 << "\t" << cluster1.numEvents << "\t" << cluster2.numEvents << std::endl;
	//   break;
	// }
	if(clusterSeparation < thisMinSeparation){
	  thisMinSeparation = clusterSeparation;
	  thisOtherCluster = clusterInd2;
	}
      }
      std::cerr << clusterInd << "\t" << thisOtherCluster << "\t" << thisMinSeparation/1e3 << "\t"
		<< cluster1.numEvents << "\t" << clusters.at(thisOtherCluster).numEvents << std::endl;
    }
  }
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
  seedPoints.resize(numClusters, 0);
  std::vector<Double_t> cumulativeDistanceSquaredToNearestClusters(numPoints, 0);


  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){


    if(clusterInd==0){

      // pick first point from uniform distribution
      Int_t firstPoint = floor(rnd.Uniform(numPoints));

      // map point to this cluster
      points.at(firstPoint).inCluster = clusterInd;

      // map cluster to piont
      seedPoints.at(clusterInd) = firstPoint;

      // create cluster object
      clusters.push_back(Cluster(points.at(firstPoint)));
      clusters.back().numEvents++;

      if(clusters.size() != (size_t)clusterInd+1){
	std::cerr << "?????????????????" << std::endl;
      }

      // assign all points to this cluster
      for(int pointInd=0; pointInd < numPoints; pointInd++){
	// points.at(pointInd).error = getDistSq(points.at(pointInd), clusters.at(clusterInd));
	points.at(pointInd).error = getAngDistSq(points.at(pointInd), clusters.at(clusterInd), pats.at(pointInd));
	points.at(pointInd).inCluster = clusterInd;
      }

    }
    else{

      // get cumulative distance to error points squared
      Double_t runningTotalDSq = 0;
      for(int pointInd=0; pointInd < numPoints; pointInd++){
	runningTotalDSq += points.at(pointInd).error; // always filled with distSq
	cumulativeDistanceSquaredToNearestClusters.at(pointInd) = runningTotalDSq;

	// if((clusterInd==65 || clusterInd==142) && pointInd == 218860){
	//   std::cerr << "Updating running total " <<  points.at(pointInd).error << "\t" << runningTotalDSq<< std::endl;
	// }

      }


      Double_t pointInrunningTotal = rnd.Uniform(runningTotalDSq);
      Int_t nextPoint = -1;
      for(int pointInd=0; pointInd < numPoints; pointInd++){
	if(pointInd==0){
	  if(pointInrunningTotal < cumulativeDistanceSquaredToNearestClusters.at(pointInd)){
	    nextPoint = pointInd;
	    break;
	  }
	}
	else{
	  if(pointInrunningTotal >= cumulativeDistanceSquaredToNearestClusters.at(pointInd-1) &&
	     pointInrunningTotal < cumulativeDistanceSquaredToNearestClusters.at(pointInd)){
	    // this one
	    nextPoint = pointInd;

	    break;
	  }
	}
      }

      // if(clusterInd == 142 || clusterInd == 65){
      // 	std::cerr << clusterInd << "\t" << nextPoint << "\t" << pointInrunningTotal << std::endl;
      // 	std::cerr << cumulativeDistanceSquaredToNearestClusters.at(nextPoint) << "\t"
      // 		  << cumulativeDistanceSquaredToNearestClusters.at(nextPoint+1) << std::endl;
      // 	std::cerr << cumulativeDistanceSquaredToNearestClusters.at(nextPoint+1) - cumulativeDistanceSquaredToNearestClusters.at(nextPoint)  << std::endl;
      // 	std::cerr << points.at(nextPoint).error << std::endl;
      // }


      // map point to this cluster
      points.at(nextPoint).inCluster = clusterInd;

      // map cluster to piont
      seedPoints.at(clusterInd) = nextPoint;

      // create cluster object
      clusters.push_back(Cluster(points.at(nextPoint)));
      clusters.back().numEvents++;

      if(clusters.size() != (size_t)clusterInd+1){
	std::cerr << "?????????????????" << std::endl;
      }

      // assign all points to this cluster
      for(int pointInd=0; pointInd < numPoints; pointInd++){
	Double_t dSq = getAngDistSq(points.at(pointInd), clusters.at(clusterInd), pats.at(pointInd));
	// Double_t dSq = getDistSq(points.at(pointInd), clusters.at(clusterInd));

	// if((clusterInd==65 || clusterInd==142) && pointInd == 218860){
	//   std::cerr << "Updating distance squared " << dSq << std::endl;
	// }

	// if new cluster centre is closer then assign point to this cluster
	// no need for else, since this tracks minimum of best produced so far...
	if(dSq < points.at(pointInd).error){
	  points.at(pointInd).error = dSq;
	  points.at(pointInd).inCluster = clusterInd;

	  // if((clusterInd==65 || clusterInd==142) && pointInd == 218860){
	  //   std::cerr << "Updating distance squared " << dSq << std::endl;
	  // }

	}
      }
    }

    // std::cerr << "Initialized " << clusterInd << std::endl;
  }


  assignPointsToClosestCluster();

  std::cout << "AnitaClusterer is initialized with " << clusters.size() << " clusters" << std::endl;

  initialized = true;
}

















void AnitaClusterer::kMeansCluster(Int_t iterationsPerCout){

  // Apparently kMeansPlusPlus is the bees knees
  if(!initialized){
    // kMeansPlusPlusInitialize();
    maxDistanceInitialize();
  }



  Double_t lastClusterError = DBL_MAX;

  // -----------------------------------------------------------
  // until completion or reach minimal improvement
  for(int i=0; i < numIter; i++){



    updateClusterCentres();

    // std::cout << "done cluster centre update" << std::endl;

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



Double_t AnitaClusterer::sumOfAngularErrorsFromLatLon(const Double_t* latLon){

  // probably for a fitter, maybe
  // requires the private variable pointsInCluster to be filled with point indices you wish to loop over
  // I guess latLon should be filled with the current values in the cluster

  Double_t lat = latLon[0];
  Double_t lon = latLon[1];


  // enable wrapping in longitude?
  while(lon >= 180){
    lon -= 360;
  }
  while(lon < -180){
    lon += 360;
  }

  Cluster tempCluster; //
  tempCluster.latitude = lat;
  tempCluster.longitude = lon;
  tempCluster.altitude = surfaceModel->SurfaceAboveGeoid(lon, lat);

  // std::cerr << std::endl;

  Double_t sumOfAngErrors = 0;
  for(int i=0; i < (Int_t)pointsInCluster.size(); i++){
    Int_t pointInd = pointsInCluster.at(i);
    Double_t err = getAngDistSq(points.at(pointInd), tempCluster, pats.at(pointInd));
    sumOfAngErrors += err;
  }


  // std::cerr << theMinCluster << "\t" << lat << "\t" << lon << "\t" << tempCluster.altitude << "\t"
  // 	    << sumOfAngErrors << "\t" << std::endl;

  return sumOfAngErrors;
}


void AnitaClusterer::updateClusterCentres(){
  // -----------------------------------------------------------
  // update cluster centres...

  // pre-check
  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    for(int dim=0; dim < nDim; dim++){
      if(clusters.at(clusterInd).numEvents <= 0){
	std::cerr << std::endl <<  "PRE-CHECK" << std::endl;
	const Point& point = points[seedPoints.at(clusterInd)];
	const Cluster& cluster = clusters.at(clusterInd);
	std::cerr << "Warning in " << __FILE__ << ", cluster with 0 events" << std::endl;
	std::cerr << "The cluster has index " << clusterInd << " and centre:" << std::endl;
	prettyPrint(nDim, cluster.centre);
	prettyPrintConvert(nDim, cluster.centre);
	std::cerr << "The seed point has index " << seedPoints.at(clusterInd) << " and centre:" << std::endl;
	prettyPrint(nDim, point.centre);
	prettyPrint(nDim, point.centre);
	std::cerr << "Seed point has a distance metric " << getDistSq(points.at(seedPoints.at(clusterInd)), cluster) << std::endl;

	std::cerr << "Also, the seed point now thinks its closest cluster is " << point.inCluster
		  << ", which is located at " << std::endl;
	prettyPrint(nDim, clusters[point.inCluster].centre);
	std::cerr << "New point has a distance metric " << getDistSq(point, cluster) << std::endl;
      }
    }
    clusters.at(clusterInd).totalError = 0;
  }

  AnitaGeomTool* geom = AnitaGeomTool::Instance();

  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    Int_t nPoints = clusters.at(clusterInd).numEvents;
    // std::vector<Int_t> pointsInCluster;
    theMinCluster = clusterInd;
    pointsInCluster.clear();
    pointsInCluster.reserve(nPoints);

    // Double_t maxLat = 0;
    // Double_t minLat = 0;

    // Double_t maxLon = 0;
    // Double_t minLon = 0;

    // Double_t maxDeltaLonLo = 0;
    // Double_t maxDeltaLonHi = 0;
    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      if(points.at(pointInd).inCluster==clusterInd){
	pointsInCluster.push_back(pointInd);
      }
    }

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    // set tolerance , etc...
    min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(-1); // SILENCE!!!

    // create funciton wrapper for minmizer
    // a IMultiGenFunction type

    Int_t numVars = 2;
    ROOT::Math::Functor FuncToMin(this, &AnitaClusterer::sumOfAngularErrorsFromLatLon, numVars);

    Double_t stepSize = 1e-3;
    std::vector<Double_t> step = std::vector<Double_t> (numVars, stepSize);

    // starting point
    std::vector<Double_t> variables = std::vector<Double_t> (numVars, 0);
    variables.at(0) = clusters.at(clusterInd).latitude;
    // variables.at(1) = clusters.at(clusterInd).latitude;
    variables.at(1) = clusters.at(clusterInd).longitude;

    min->SetFunction(FuncToMin);

    // min->SetVariable(0, "latitude", variables[0], step[0]);
    // min->SetVariable(1, "longitude", variables[1], step[1]);

    min->SetLimitedVariable(0, "latitude", variables[0], step[0], -90, -60);
    min->SetLimitedVariable(1, "longitude", variables[1], step[1], variables[1]-180, variables[1]+180);

    min->Minimize();

    Int_t s = min->Status();

    if(s!=0){
      std::cerr << "fitter has status " << s << " for " << clusterInd << "\t" << clusters.at(clusterInd).numEvents << "\t" << clusters.at(clusterInd).latitude << "\t" << clusters.at(clusterInd).longitude << std::endl;
    }

    // Time it
    // TStopwatch watch;
    // watch.Start(kTRUE);

    // // do the minimization
    // min->Minimize();

    // // Time!
    // watch.Start(kFALSE);
    // Int_t seconds = Int_t(watch.RealTime());
    // Int_t hours = seconds / 3600;
    // hours = hours < 0 ? 0 : hours;
    // seconds = seconds - hours * 3600;
    // Int_t mins = seconds / 60;
    // mins = mins < 0 ? 0 : mins;
    // seconds = seconds - mins * 60;
    // fprintf(stderr, "Minimization took %02d:%02d:%02d\n", hours, mins, seconds);

    // std::cout << "Minimum = " << min->MinValue() << std::endl;
    const Double_t* xs = min->X();
    clusters.at(clusterInd).latitude = xs[0];
    clusters.at(clusterInd).longitude = xs[1];
    while(clusters.at(clusterInd).longitude >= 180){
      clusters.at(clusterInd).longitude -= 360;
    }
    while(clusters.at(clusterInd).longitude < -180){
      clusters.at(clusterInd).longitude += 360;
    }

    clusters.at(clusterInd).altitude = surfaceModel->SurfaceAboveGeoid(clusters.at(clusterInd).longitude,
								       clusters.at(clusterInd).latitude);

    clusters.at(clusterInd).totalError = 0;


    geom->getCartesianCoords(clusters.at(clusterInd).latitude, clusters.at(clusterInd).longitude,
			     clusters.at(clusterInd).altitude, clusters.at(clusterInd).centre);
  }



  for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
    for(int dim=0; dim < nDim; dim++){
      if(clusters.at(clusterInd).numEvents <= 0){
      // 	clusters.at(clusterInd).centre[dim]/=clusters.at(clusterInd).numEvents;
      // }
      // else{
  	std::cerr << std::endl <<  "POST-CHECK" << std::endl;
  	const Point& point = points[seedPoints.at(clusterInd)];
  	const Cluster& cluster = clusters.at(clusterInd);
  	std::cerr << "Warning in " << __FILE__ << ", cluster with 0 events" << std::endl;
  	std::cerr << "The cluster has index " << clusterInd << " and centre:" << std::endl;
  	prettyPrint(nDim, cluster.centre);
  	prettyPrintConvert(nDim, cluster.centre);
  	std::cerr << "The seed point has index " << seedPoints.at(clusterInd) << " and centre:" << std::endl;
  	prettyPrint(nDim, point.centre);
  	// std::cerr << "Seed point has a distance metric " << getDistSq(points.at(seedPoints.at(clusterInd)), cluster) << std::endl;
  	// std::cerr << "Seed point has a distance metric " << getAngDistSq(points.at(seedPoints.at(clusterInd)), cluster, pats.at(seedPoints.at(clusterInd))) << std::endl;

  	// std::cerr << "Also, the seed point now thinks its closest cluster is " << point.inCluster
  	// 	  << ", which is located at " << std::endl;
  	// prettyPrint(nDim, clusters[point.inCluster].centre);
  	// std::cerr << "New point has a distance metric " << getAngDistSq(point, cluster, pats.at(pointInd)) << std::endl;
      }
    }
  }
  //   AnitaGeomTool* geom = AnitaGeomTool::Instance();
  //   geom->getLatLonAltFromCartesian(clusters.at(clusterInd).centre, clusters.at(clusterInd).latitude, clusters.at(clusterInd).longitude, clusters.at(clusterInd).altitude);

  // clusters.at(clusterInd).totalError = 0;
  // }
}




// void AnitaClusterer::updateClusterCentres(){
//   // -----------------------------------------------------------
//   // update cluster centres...

//   // pre-check
//   for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
//     for(int dim=0; dim < nDim; dim++){
//       if(clusters.at(clusterInd).numEvents <= 0){
// 	std::cerr << std::endl <<  "PRE-CHECK" << std::endl;
// 	const Point& point = points[seedPoints.at(clusterInd)];
// 	const Cluster& cluster = clusters.at(clusterInd);
// 	std::cerr << "Warning in " << __FILE__ << ", cluster with 0 events" << std::endl;
// 	std::cerr << "The cluster has index " << clusterInd << " and centre:" << std::endl;
// 	prettyPrint(nDim, cluster.centre);
// 	prettyPrintConvert(nDim, cluster.centre);
// 	std::cerr << "The seed point has index " << seedPoints.at(clusterInd) << " and centre:" << std::endl;
// 	prettyPrint(nDim, point.centre);
// 	std::cerr << "Seed point has a distance metric " << getDistSq(points.at(seedPoints.at(clusterInd)), cluster) << std::endl;

// 	std::cerr << "Also, the seed point now thinks its closest cluster is " << point.inCluster
// 		  << ", which is located at " << std::endl;
// 	prettyPrint(nDim, clusters[point.inCluster].centre);
// 	std::cerr << "New point has a distance metric " << getDistSq(point, cluster) << std::endl;
//       }
//     }
//     clusters.at(clusterInd).totalError = 0;
//   }

//   for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
//     for(int dim=0; dim < nDim; dim++){
//       clusters.at(clusterInd).centre[dim] = 0;
//     }
//     clusters.at(clusterInd).numEvents = 0;
//   }
//   // Int_t num142 = 0;
//   for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
//     Int_t clusterInd = points.at(pointInd).inCluster;
//     for(int dim=0; dim < nDim; dim++){
//       clusters.at(clusterInd).centre[dim] += points.at(pointInd).centre[dim];
//     }
//     clusters.at(clusterInd).numEvents++;
//   }
//   for(int clusterInd=0; clusterInd < (Int_t) clusters.size(); clusterInd++){
//     for(int dim=0; dim < nDim; dim++){
//       if(clusters.at(clusterInd).numEvents > 0){
// 	clusters.at(clusterInd).centre[dim]/=clusters.at(clusterInd).numEvents;
//       }
//       else{
// 	std::cerr << std::endl <<  "POST-CHECK" << std::endl;
// 	const Point& point = points[seedPoints.at(clusterInd)];
// 	const Cluster& cluster = clusters.at(clusterInd);
// 	std::cerr << "Warning in " << __FILE__ << ", cluster with 0 events" << std::endl;
// 	std::cerr << "The cluster has index " << clusterInd << " and centre:" << std::endl;
// 	prettyPrint(nDim, cluster.centre);
// 	prettyPrintConvert(nDim, cluster.centre);
// 	std::cerr << "The seed point has index " << seedPoints.at(clusterInd) << " and centre:" << std::endl;
// 	prettyPrint(nDim, point.centre);
// 	std::cerr << "Seed point has a distance metric " << getDistSq(points.at(seedPoints.at(clusterInd)), cluster) << std::endl;

// 	std::cerr << "Also, the seed point now thinks its closest cluster is " << point.inCluster
// 		  << ", which is located at " << std::endl;
// 	prettyPrint(nDim, clusters[point.inCluster].centre);
// 	std::cerr << "New point has a distance metric " << getDistSq(point, cluster) << std::endl;
//       }
//     }

//     AnitaGeomTool* geom = AnitaGeomTool::Instance();
//     geom->getLatLonAltFromCartesian(clusters.at(clusterInd).centre, clusters.at(clusterInd).latitude, clusters.at(clusterInd).longitude, clusters.at(clusterInd).altitude);

//     clusters.at(clusterInd).totalError = 0;
//   }
// }






Double_t AnitaClusterer::assignErrorValues(){
  // -----------------------------------------------------------
  // assign error values
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    clusters.at(clusterInd).totalError = 0;
  }

  Double_t sumClusterErrors = 0;
  for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
    Int_t clusterInd = points.at(pointInd).inCluster;
    if(clusterInd >= 0){
      Double_t d2 = getAngDistSq(points.at(pointInd), clusters.at(clusterInd), pats.at(pointInd));
      points.at(pointInd).error = d2;
      clusters.at(clusterInd).totalError += d2;
      sumClusterErrors += d2;

      Double_t deltaDistM = getDistSq(points.at(pointInd), clusters.at(clusterInd));
      if(deltaDistM > clusters.at(clusterInd).maxDist){
	clusters.at(clusterInd).maxDist = deltaDistM;
      }

    }
    else{
      points.at(pointInd).error = -9999;
    }
  }
  return sumClusterErrors;
}









// Double_t AnitaClusterer::assignErrorValues(){
//   // -----------------------------------------------------------
//   // assign error values
//   Double_t sumClusterErrors = 0;
//   for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
//     Int_t clusterInd = points.at(pointInd).inCluster;
//     Double_t d2 = getDistSq(points.at(pointInd), clusters.at(clusterInd));
//     points.at(pointInd).error = d2;
//     clusters.at(clusterInd).totalError += d2;
//     sumClusterErrors += d2;
//   }
//   return sumClusterErrors;
// }














TGraph* AnitaClusterer::makeClusterSummaryTGraph(Int_t clusterInd){

  TGraph* gr = NULL;
  if(clusterInd >= 0 && clusterInd < numClusters){

    TString name  = TString::Format("grCluster%d", clusterInd);
    TString title  = TString::Format("Cluster %d; Longitude (Degrees); Latitude (Degrees)", clusterInd);
    gr = new TGraph();
    gr->SetName(name);
    gr->SetTitle(title);

    // AnitaGeomTool* geom = AnitaGeomTool::Instance();

    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      if(points.at(pointInd).inCluster==clusterInd){

	// Double_t lat, lon, alt;
	// geom->getLatLonAltFromCartesian(points.at(pointInd).centre, lat, lon, alt);
	gr->SetPoint(gr->GetN(), points.at(pointInd).latitude, points.at(pointInd).longitude);
      }
    }
  }
  return gr;
}


















TTree* AnitaClusterer::makeClusterSummaryTree(TFile* fOut){

  std::cout << fOut << std::endl;//->cd();

  TTree* clusterTree = new TTree("clusterTree", "Tree of clustered ANITA events");

  ClusteredAnitaEvent* clusteredEvent;
  clusterTree->Branch("clusteredEvent", &clusteredEvent);

  // AnitaGeomTool* geom = AnitaGeomTool::Instance();
  const int numBases = BaseList::getNumBases();

  for(Int_t pointInd=0; pointInd < (Int_t)points.size(); pointInd++){

    const Point& point = points.at(pointInd);

    clusteredEvent = new ClusteredAnitaEvent();
    clusteredEvent->eventNumber = eventNumbers.at(pointInd);
    clusteredEvent->run = runs.at(pointInd);
    clusteredEvent->pol = points.at(pointInd).pol;

    // convert from km back to m for conversion to lat/lon/alt
    // for(int dim=0; dim < nDim; dim++){
    //   // clusteredEvent->eventPosition[dim] = point.centre[dim];
    //   clusteredEvent->eventPosition[dim] = point.centre[dim];
    // }
    clusteredEvent->eventLat = point.latitude;
    clusteredEvent->eventLon = point.longitude;
    clusteredEvent->eventAlt = point.altitude;


    clusteredEvent->anitaLat = pats.at(pointInd)->latitude;
    clusteredEvent->anitaLon = pats.at(pointInd)->longitude;
    clusteredEvent->anitaAlt = pats.at(pointInd)->altitude;

    clusteredEvent->thetaDeg = point.thetaDeg;
    clusteredEvent->phiDeg = point.phiDeg;

    Int_t clusterInd = points.at(pointInd).inCluster;
    clusteredEvent->inCluster = clusterInd;
    clusteredEvent->isBase = 0;

    if(clusterInd >= 0){
      const Cluster& cluster = clusters.at(clusterInd);

      clusteredEvent->isBase = clusterInd < numBases ? 1 : 0;

      clusteredEvent->clusterLat = cluster.latitude;
      clusteredEvent->clusterLon = cluster.longitude;
      clusteredEvent->clusterAlt = cluster.altitude;

      // Double_t deltaThetaDeg, deltaPhiDeg;
      getDeltaThetaDegDeltaPhiDegCluster(point, cluster, pats.at(pointInd), \
					 clusteredEvent->deltaThetaDeg, \
					 clusteredEvent->deltaPhiDeg);

      // clusteredEvent->distanceToClusterCentroid = TMath::Sqrt(getDistSq(point, cluster));
      // clusteredEvent->distanceToClusterCentroid = get(point, cluster));
      clusteredEvent->minusTwoLogLikelihood = getAngDistSq(point, cluster, pats.at(pointInd));
      UsefulAdu5Pat usefulPat(pats.at(pointInd));
      clusteredEvent->distanceToClusterCentroid = usefulPat.getDistanceFromSource(cluster.latitude, \
										  cluster.longitude, \
										  cluster.altitude);


      if(clusteredEvent->secondClosestCluster >= 0){
	const Cluster& cluster2 = clusters.at(clusteredEvent->secondClosestCluster);


	clusteredEvent->distanceToClusterCentroidSecondBest = usefulPat.getDistanceFromSource(cluster2.latitude, \
											      cluster2.longitude, \
											      cluster2.altitude);
      }

      clusteredEvent->numEventsInCluster = cluster.numEvents;

      clusteredEvent->numClusters = numClusters;
      clusteredEvent->numIterations = numIter;
    }
    clusterTree->Fill();
    delete clusteredEvent;
  }

  return clusterTree;

}
