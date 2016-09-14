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
  llCut = 250;
  maxDistCluster = 800e3; // try 800km
  numCallsToRecursive = 0;
  minimalImprovement  = 0.01;
  // can use this to move cluster around surface, 3D positions correctly becomes 2D problem
  doneDefaultAssignment = false;
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




inline void getDeltaThetaDegDeltaPhiDegCluster(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, UsefulAdu5Pat& usefulPat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  Double_t thetaWave, phiWave;
  usefulPat.getThetaAndPhiWave(cluster.longitude, cluster.latitude,cluster.altitude,
			       thetaWave,phiWave);
  Double_t thetaDeg = -1*TMath::RadToDeg()*thetaWave;
  Double_t phiDeg = TMath::RadToDeg()*phiWave;

  deltaThetaDeg = (thetaDeg - point.thetaDeg);
  deltaPhiDeg = RootTools::getDeltaAngleDeg(phiDeg, point.phiDeg);
}






inline void getDeltaThetaDegDeltaPhiDegCluster(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, const Adu5Pat* pat, Double_t& deltaThetaDeg, Double_t& deltaPhiDeg){

  UsefulAdu5Pat usefulPat(pat);

  getDeltaThetaDegDeltaPhiDegCluster(point, cluster, usefulPat, deltaThetaDeg, deltaPhiDeg);
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

inline Double_t getAngDistSq(const AnitaClusterer::Point& point, const AnitaClusterer::Cluster& cluster, UsefulAdu5Pat& usefulPat){

  Double_t deltaThetaDeg, deltaPhiDeg;
  getDeltaThetaDegDeltaPhiDegCluster(point, cluster, usefulPat, deltaThetaDeg, deltaPhiDeg);

  // normalized
  Double_t dThetaSq = deltaThetaDeg/point.sigmaThetaDeg;
  Double_t dPhiSq = deltaPhiDeg/point.sigmaPhiDeg;

  Double_t angSq =  dThetaSq*dThetaSq + dPhiSq*dPhiSq;

  return angSq;
}











Int_t AnitaClusterer::histogramUnclusteredEvents(Int_t& globalMaxBin){

  Int_t addClusterIter =  hUnclustereds.size();
  TString name = TString::Format("hUnclustered%d", addClusterIter);
  // const int nBinsLat = 30;
  // const int nBinsLon = 360;
  const int nBins = 1024;
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

  if(doneDefaultAssignment==false){
    assignEventsToDefaultClusters();
  }

  numCallsToRecursive++;
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

    const int isMC = 0;
    for(int pointInd=0; pointInd < (Int_t) points.size(); pointInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, numClusters-1);
      // std::cout << pointInd << "\t" << points.at(pointInd).ll << "\t" << std::endl;
    }

    // std::cout << counter << "\t" << maxBinVal << "\t" << clusters.back().latitude << "\t"
    // 	      << clusters.back().longitude << "\t" << hUnclustereds.back()->Integral() << std::endl;


    // if(hUnclustereds.size() < 10){
    recursivelyAddClusters(minBinContent);
    // }
  }
}




void AnitaClusterer::findClosestPointToClustersOfSizeOne(){

  // loop through clusters and find the ones that have a size of one
  Int_t numSinglets = 0;
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    if(clusters.at(clusterInd).numEvents == 1){

      Cluster& cluster = clusters.at(clusterInd);

      // find the point in this cluster
      Int_t isolatedPointInd = -1;
      for(int pointInd=0; pointInd < (int)points.size(); pointInd++){
	if(points.at(pointInd).inCluster==clusterInd){
	  isolatedPointInd = pointInd;
	  break;
	}
      }

      assert(isolatedPointInd!=-1);

      Int_t numNearPoints = 0;

      // now loop through and see how close this isolated cluster is close to any other point...
      for(int pointInd=0; pointInd < (int)points.size(); pointInd++){
	if(pointInd!=isolatedPointInd){
	  UsefulAdu5Pat usefulPat(pats.at(pointInd));
	  Double_t dist = usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

	  if(dist < maxDistCluster){
	    Double_t ll = getAngDistSq(points.at(pointInd), cluster, usefulPat);

	    if(ll < llCut){
	      numNearPoints++;
	      cluster.numPointsWithinMinLL++;
	    }
	  }
	}
      }

      std::cout << runs.at(isolatedPointInd) << "\t" << eventNumbers.at(isolatedPointInd) << "\t"
		<< clusterInd << "\t" << numNearPoints << std::endl;
      numSinglets++;
    }
  }

  // std::cout << "numSinglets = " << numSinglets << std::endl;

}







void AnitaClusterer::assignSinglePointToCloserCluster(Int_t pointInd, Int_t isMC, Int_t clusterInd){
  Point& point = isMC==0 ? points.at(pointInd) : mcpoints.at(pointInd);
  Adu5Pat* pat = isMC==0 ? pats.at(pointInd)   : mcpats.at(pointInd);
  Cluster& cluster = clusters.at(clusterInd);

  UsefulAdu5Pat usefulPat(pat);

  Double_t distM = usefulPat.getDistanceFromSource(cluster.latitude, cluster.longitude, cluster.altitude);

  if(distM < maxDistCluster){ // are we even close?

    Double_t ll = getAngDistSq(point, cluster, usefulPat);

    if(ll < llCut){

      if(ll < point.ll){
	point.secondClosestCluster = point.inCluster;
	point.llSecondBest = point.ll;

	point.ll = ll;
	point.inCluster = clusterInd;

	if(isMC==0){
	  if(point.secondClosestCluster >= 0){
	    clusters.at(point.secondClosestCluster).numEvents--;
	  }
	  cluster.numEvents++;
	}
      }
    }
  }
}




void AnitaClusterer::assignMCPointsToClusters(){

  const int isMC = 1;

  for(int pointInd=0; pointInd < (Int_t)mcpoints.size(); pointInd++){
    for(int clusterInd=0; clusterInd < (Int_t)clusters.size(); clusterInd++){
      assignSinglePointToCloserCluster(pointInd, isMC, clusterInd);
    }
  }
}












size_t AnitaClusterer::addMCPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol, Double_t weight){

  mcpoints.push_back(MCPoint(pat, latitude, longitude, altitude, sigmaThetaDeg, sigmaPhiDeg, pol, weight));
  mceventNumbers.push_back(eventNumber);
  mcruns.push_back(run);
  mcpats.push_back((Adu5Pat*)pat->Clone());

  return points.size();
}






size_t AnitaClusterer::addPoint(Adu5Pat* pat, Double_t latitude, Double_t longitude, Double_t altitude, Int_t run, UInt_t eventNumber, Double_t sigmaThetaDeg, Double_t sigmaPhiDeg, AnitaPol::AnitaPol_t pol){

  points.push_back(Point(pat, latitude, longitude, altitude, sigmaThetaDeg, sigmaPhiDeg, pol));
  eventNumbers.push_back(eventNumber);
  runs.push_back(run);
  pats.push_back((Adu5Pat*)pat->Clone());

  return points.size();
}




void AnitaClusterer::assignEventsToDefaultClusters(){
  ProgressBar p(numClusters);
  const int isMC = 0;
    for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
      for(int pointInd=0; pointInd < (int) points.size(); pointInd++){
	assignSinglePointToCloserCluster(pointInd, isMC, clusterInd);
      }
    p++;
  }
  doneDefaultAssignment = true;
}




void AnitaClusterer::initializeBaseList(){

  numClusters = (int) BaseList::getNumBases();
  std::cout << "Initializing base list..." << std::endl;
  for(int clusterInd=0; clusterInd < numClusters; clusterInd++){
    const BaseList::base& base = BaseList::getBase(clusterInd);

    clusters.push_back(Cluster(base));

  }
  initialized = true;
}



void AnitaClusterer::resetClusters(){

  for(int pointInd=0; pointInd < (int) points.size(); pointInd++){
    points.at(pointInd).ll = DBL_MAX;
    points.at(pointInd).llSecondBest = DBL_MAX;
    points.at(pointInd).inCluster = -1;
    points.at(pointInd).secondClosestCluster = -1;
  }

  const int numBases = BaseList::getNumBases();
  while((int) clusters.size() >= numBases){
    clusters.pop_back();
  }
  numClusters = clusters.size();

  for(int clusterInd=0; clusterInd < (int) numClusters; clusterInd++){
    clusters.at(clusterInd).numEvents = 0;
  }
  doneDefaultAssignment = false;
}


































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

    Int_t clusterInd = points.at(pointInd).inCluster;
    clusteredEvent->inCluster = clusterInd;
    clusteredEvent->isBase = 0;

    if(clusterInd >= 0){
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
      clusteredEvent->isMC = 0;
      clusteredEvent->weight = 1;


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
      clusteredEvent->numEventsInCluster = cluster.numEvents;


      if(clusteredEvent->secondClosestCluster >= 0){
	const Cluster& cluster2 = clusters.at(clusteredEvent->secondClosestCluster);


	clusteredEvent->distanceToClusterCentroidSecondBest = usefulPat.getDistanceFromSource(cluster2.latitude, \
											      cluster2.longitude, \
											      cluster2.altitude);
      }

      clusteredEvent->numIterations = numIter;

      clusteredEvent->numClusters = numClusters;

    }
    else{
      std::cerr << "clusterInd = " << clusterInd << ". ";
      std::cerr << "This shouldn't be possible!" << std::endl;
    }
    clusterTree->Fill();


    delete clusteredEvent;
  }




  for(Int_t pointInd=0; pointInd < (Int_t)mcpoints.size(); pointInd++){

    const MCPoint& point = mcpoints.at(pointInd);

    clusteredEvent = new ClusteredAnitaEvent();
    clusteredEvent->eventNumber = mceventNumbers.at(pointInd);
    clusteredEvent->run = mcruns.at(pointInd);
    clusteredEvent->pol = point.pol;
    clusteredEvent->isMC = 1;
    clusteredEvent->weight = point.weight;

    // convert from km back to m for conversion to lat/lon/alt
    // for(int dim=0; dim < nDim; dim++){
    //   // clusteredEvent->eventPosition[dim] = point.centre[dim];
    //   clusteredEvent->eventPosition[dim] = point.centre[dim];
    // }
    clusteredEvent->eventLat = point.latitude;
    clusteredEvent->eventLon = point.longitude;
    clusteredEvent->eventAlt = point.altitude;


    clusteredEvent->anitaLat = mcpats.at(pointInd)->latitude;
    clusteredEvent->anitaLon = mcpats.at(pointInd)->longitude;
    clusteredEvent->anitaAlt = mcpats.at(pointInd)->altitude;

    clusteredEvent->thetaDeg = point.thetaDeg;
    clusteredEvent->phiDeg = point.phiDeg;

    Int_t clusterInd = mcpoints.at(pointInd).inCluster;
    clusteredEvent->inCluster = clusterInd;
    clusteredEvent->isBase = 0;

    if(clusterInd >= 0){
      const Cluster& cluster = clusters.at(clusterInd);

      clusteredEvent->isBase = clusterInd < numBases ? 1 : 0;

      clusteredEvent->clusterLat = cluster.latitude;
      clusteredEvent->clusterLon = cluster.longitude;
      clusteredEvent->clusterAlt = cluster.altitude;

      // Double_t deltaThetaDeg, deltaPhiDeg;
      getDeltaThetaDegDeltaPhiDegCluster(point, cluster, mcpats.at(pointInd), \
					 clusteredEvent->deltaThetaDeg, \
					 clusteredEvent->deltaPhiDeg);

      // clusteredEvent->distanceToClusterCentroid = TMath::Sqrt(getDistSq(point, cluster));
      // clusteredEvent->distanceToClusterCentroid = get(point, cluster));
      clusteredEvent->minusTwoLogLikelihood = getAngDistSq(point, cluster, mcpats.at(pointInd));
      UsefulAdu5Pat usefulPat(mcpats.at(pointInd));
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
