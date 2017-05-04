/*
 * DDPixelHit.cpp
 *
 *  Created on: Feb 28, 2017
 *      Author: cligtenb
 */

#include "DDKalTest/DDPixelHit.h"
#include <iostream>
#include <iomanip>

#include "DDKalTest/DDPixelMeasVolume.h"

using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::ios;
using std::resetiosflags;

//_________________________________________________________________________
//  --------------------------------
//  Implementation of public methods
//  --------------------------------
//

/** Constructor */

DDPixelHit::DDPixelHit(const TVMeasLayer& measurementVolume, Double_t* x,
		Double_t* dx, Double_t bfield, EVENT::TrackerHit* trkhit):
DDVTrackHit(measurementVolume, x, dx, bfield, 2, trkhit)
			  { /* no op */ }

/** Global to Local coordinates */

TKalMatrix DDPixelHit::XvToMv(const TVector3&, Double_t) const
{
  streamlog_out(ERROR) << "DDPixelHit::XvToMv called without TKalTrackState" << std::endl;
  throw "unexpected call";
}

inline std::ostream& operator<<(std::ostream& o, const TVector3& v) {
	return o<<"("<<v.X()<<","<<v.Y()<<","<<v.Z()<<")";
}

TKalMatrix DDPixelHit::XvToMv(const TVector3 &xv, Double_t, const TKalTrackState& a) const
{
	//get track parameter
	double trackPhi0=fmod(a(1,0) + M_PI/2., 2.*M_PI);
	double trackTanTheta=a(4,0);

	//get hit pos
	TVector3 hitPos(getLCIOTrackerHit()->getPosition());
	TVector3 dist=hitPos-xv;

	double d0=dist.x()*sin(trackPhi0)-dist.y()*cos(trackPhi0);
	double r0=dist.x()*cos(trackPhi0)+dist.y()*sin(trackPhi0); //distance along track in xy to get to point of closest approach in xy
	double z=xv.z()+r0*trackTanTheta; //move an equal step in z

//	std::cout << "called XvToMv( xv = "<<xv<<" , ... )"<<"trackPhi= "<<trackPhi0<<" hitPos= "<<hitPos<<" r0= "<<r0<<" mv becomes ( "<<d0<<", "<<z<<" )"<<endl;

	//calculate measurement vector here
	TKalMatrix mv(2,1);
	mv(0,0)=d0;
	mv(1,0)=z;

//	auto covmat=getLCIOTrackerHit()->getCovMatrix();
//	std::cout<<"pull of d0 = "<<d0/sqrt( (covmat[0]+covmat[2])/2 )<<" dz = "<<sqrt(covmat[5])<<std::endl;

	//check against cylinder hit for debug
   // Calculate hit coordinate information:
   //   mv(0, 0) = r * phi - rhit * phihit
//   double cylinderd0= xv.Perp() * xv.Phi() - hitPos.Perp() * hitPos.Phi() ;
//   std::cout << "cylinderd0 = "<<cylinderd0<<" d0= "<<d0<<" difference= "<<cylinderd0-d0<<endl;
//   mv(0,0)=cylinderd0;

   return mv;
}

/** Print Debug information */

void DDPixelHit::DebugPrint(Option_t *) const
{
  std::cerr << "------------------- Site Info -------------------------" << std::endl;

  for (Int_t i = 0; i < GetDimension(); i++) {
    Double_t x  = (*this)(i, 0);
    Double_t dx = (*this)(i, 1);
    std::cerr  << " x[" << i << "] = " << setw(8) << setprecision(5) << x
    << "    "
    << "dx[" << i << "] = " << setw(6) << setprecision(2) << dx
    << setprecision(7)
    << resetiosflags(ios::showpoint)
    << endl;
  }
  std::cerr  << "-------------------------------------------------------"  << endl;
}
//
//const TVMeasLayer& DDPixelHit::GetMeasLayer() const {
//	DDPixelMeasVolume* pixelVolume=dynamic_cast<DDPixelMeasVolume*>(fMeasLayerPtr);
//	if(!pixelVolume) {
//		streamlog_out(WARNING)<<"DDPixelHit is not on a DDPixelMeasVolume"<<endl;
//	}
//	//todo: set DDPixelVolume to this hit here!
//
//	return *fMeasLayerPtr;
//}
