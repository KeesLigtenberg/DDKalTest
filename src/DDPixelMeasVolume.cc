#include <DDKalTest/DDPixelHit.h>
#include <DDKalTest/DDPixelMeasVolume.h>
#include "TKalTrack.h" 

#include <lcio.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/TrackerHitZCylinder.h>
#include <UTIL/Operators.h>
#include <UTIL/LCTrackerConf.h>

#include "DD4hep/DD4hepUnits.h"
#include "DDSurfaces/Vector3D.h"

#include "aidaTT/trajectory.hh"

#include "streamlog/streamlog.h"

#include "TMath.h"
#include "TCircle.h"

#include <cmath>

#include "DDKalTest/DDGetTPCFieldDescription.h"

namespace{
  /** helper function to restrict the range of the azimuthal angle to ]-pi,pi]*/
  inline double toBaseRange( double phi)  {
    while( phi <= -M_PI ){  phi += 2. * M_PI ; }
    while( phi >   M_PI ){  phi -= 2. * M_PI ; }
    return phi ;
  }

  inline std::ostream& operator<<(std::ostream& os, const TVector3& x) {return os<<"("<<x.x()<<", "<<x.y()<<", "<<x.z()<<")"; }
  inline std::ostream& operator<<(std::ostream& os, const TVector2& x) {return os<<"("<<x.X()<<", "<<x.Y()<<")"; }

}

using namespace UTIL ;

DDPixelMeasVolume::DDPixelMeasVolume(DDSurfaces::ISurface* surf,
					 Double_t   Bz,
					 const Char_t  *name ) :
  DDVMeasLayer(  surf, Bz, name ) ,
  
  TCylinder(  dynamic_cast<DDSurfaces::ICylinder*>(surf)->radius()/dd4hep::mm ,
	      surf->length_along_v()/dd4hep::mm / 2. ,
	      dynamic_cast<DDSurfaces::ICylinder*>(surf)->center().x()/dd4hep::mm,
	      dynamic_cast<DDSurfaces::ICylinder*>(surf)->center().y()/dd4hep::mm ,
	      dynamic_cast<DDSurfaces::ICylinder*>(surf)->center().z()/dd4hep::mm ),
  
  fSortingPolicy(0.),
  
  fMDim( 2 /*d0, z*/ )  {
  
  static double epsilon=1e-4 ;

  UTIL::BitField64 encoder( getDDFieldDescription() ) ;
  encoder.setValue( surf->id() );

  int side = encoder[ UTIL::LCTrackerCellID::side() ] ;

  fSortingPolicy = dynamic_cast<DDSurfaces::ICylinder*>(surf)->radius()/dd4hep::mm + side * epsilon ;

  // assumptions made here: the cylinder runs parallel to z and v ...
  
  streamlog_out(DEBUG1) << "DDPixelMeasVolume created"
			<< " R = " << this->GetR() 
			<< " sorting policy : " << GetSortingPolicy()
			<< " is_active = " << surf->type().isSensitive()  
			<< " CellID = " << UTIL::LCTrackerCellID::valueString( surf->id() )
			<< " layer = " << encoder[ UTIL::LCTrackerCellID::layer() ]
			<< " name = " << this->DDVMeasLayer::GetName()  
			<< std::endl ;

  // for a cylindrical layer we also set the side to 0 in the layerId 
  // ( in the case the cylinder is split between forward and backward )
  encoder[ UTIL::LCTrackerCellID::side() ] = 0;
  encoder[ UTIL::LCTrackerCellID::module() ] = 0;
  encoder[ UTIL::LCTrackerCellID::sensor() ] = 0;
  
  _layerID = encoder.lowWord();

  streamlog_out(DEBUG0)<<"DDPixelMeasVolume with layerID="<<_layerID<<std::endl;

}


/** Global to Local coordinates */

TKalMatrix DDPixelMeasVolume::XvToMv(const TVector3 &xv, const TKalTrackState&) const
{
  
  if(fMDim!=2) {
	  streamlog_out(WARNING)<<"DDPixelMeasVolume::XvToMv fmDim = "<<fMDim<<", but expected 2"<<std::endl;
  }

  TKalMatrix mv( fMDim , 1 );
  
  streamlog_out(ERROR)<<" Did not expect function call DDPixelMeasVolume::XvToMv"<<std::endl;

  //is this function necessary? if I'm correct XvToMv is ONLY called through hits
  
  return mv;
}


/** Local to Global coordinates */

TVector3 DDPixelMeasVolume::HitToXv(const TVTrackHit &vht) const
{
	if(auto hit=dynamic_cast<const DDPixelHit*>(&vht)) {
		TVector3 hitPos( hit->getLCIOTrackerHit()->getPosition() );
		return hitPos;
	} else {
		//not a pixel hit. Try cylinder?
		streamlog_out(WARNING)<<"DDPixelMeasVolume::HitToXv called with hit other than DDPixelHit"<<std::endl;

	   //copied from commented section DDCylindricalMeasLayer
	   Double_t phi = vht(0, 0) / GetR() ;
	   Double_t z   = vht(1, 0);
	   // account for cylinder not centered at x=0.0, y=0.0
	   Double_t x   = GetR() * TMath::Cos(phi) + GetXc().X();
	   Double_t y   = GetR() * TMath::Sin(phi) + GetXc().Y();
	   return TVector3(x, y, z);
	}
}

TMatrixD DDPixelMeasVolume::CalcDSDx(const TVector3 &xx) const
{
//   std::cout<<"DDPixelMeasVolume::CalcDSDx: xx="<<xx<<"is on surface: "<<IsOnSurface(xx)<<" at R="<<GetR()<<std::endl;
   if(IsOnSurface(xx) || xx.Mag2() < 1E-100 /*??*/ ) {
	   //use normal vector if cylinder
	   return TCylinder::CalcDSDx(xx);
   } else {
	   //use same vector if not on surface
//	   std::cout<<"DDPixelMeasVolume::CalcDSDx is not on surface at "<<xx<<", alternative calculation"<<std::endl;
	   TVector3 xxc = (xx - GetXc()).Unit()*GetR();
	   TMatrixD dsdx(1,3);
	   dsdx(0,0) = 2.*xxc.X();
	   dsdx(0,1) = 2.*xxc.Y();
	   dsdx(0,2) = 0.;
	   return dsdx;
   }
}



void DDPixelMeasVolume::CalcDhDa(const TVTrackHit& ht, const TKalTrackState& a,
		const TVector3& xxv, const TKalMatrix& dxphiada, TKalMatrix& H) const {

//	  // Calculate
//	  //    H = (@h/@a) = (@d0/@a, @z/@a)


	  Int_t sdim = H.GetNcols();
	  Int_t hdim = TMath::Max(5, sdim - 1);

	  if(sdim!=5) {
		  streamlog_out(WARNING)<<"DDPixelMeasVolume::CalcDhDa expected H with 5 columns, has "<<sdim<<std::endl;
	  }


//d0 z0

	  TVector3 xxvc = xxv - GetXc();

	  double trackPhi0=a(1,0)+M_PI/2.;
	  double cylinderPhi0=atan2(xxv.Y(), xxv.X());
//	  std::cout<<"Phi track= "<<trackPhi0<<" cylinder= "<<cylinderPhi0<<std::endl;
//	  std::cout<<"dxphiada=\n"<<dxphiada;

//	  Double_t xv   = xxvc.X();
//	  Double_t yv   = xxvc.Y();
//	  Double_t r = xv * xv + yv * yv;

	  // Set H = (@h/@a) = (@d/@a, @z/@a)^t

	  for (Int_t i = 0; i < hdim; i++) {
	    H(0, i)  = -sin(trackPhi0) * dxphiada(0, i) + cos(trackPhi0) * dxphiada(1, i);
	    H(1, i)  = dxphiada(2, i);
	  }

	  if (sdim == 6) {
	    H(0, sdim - 1) = 0.;
	  }

	//cylinder
	  // account for cylinder not centered at x=0.0, y=0.0
//	  TVector3 xxvc = xxv - GetXc();
//
//	  Double_t xv   = xxvc.X();
//	  Double_t yv   = xxvc.Y();
//	  Double_t xxyy = xv * xv + yv * yv;
//
//	  // Set H = (@h/@a) = (@d/@a, @z/@a)^t
//
//	  for (Int_t i = 0; i < hdim; i++) {
//	    H(0, i)  = - (yv / xxyy) * dxphiada(0, i) //y/r2*r=y/r=sin
//	    + (xv / xxyy) * dxphiada(1, i);
//	    H(0, i) *= GetR();
//
//	    H(1, i)  = dxphiada(2, i);
//	  }
//
//	  if (sdim == 6) {
//	    H(0, sdim - 1) = 0.;
//	  }

}

/** Calculate Projector Matrix */
void DDPixelMeasVolume::CalcDhDa(const TVTrackHit &vht, // tracker hit not used here
                                    const TVector3   &xxv,
                                    const TKalMatrix &dxphiada,
                                    TKalMatrix &H) const
{
	streamlog_out(ERROR)<<"DDPixelMeasVolume::CalcDhDa called without TKalTrackState"<<std::endl;
	throw "unexpected call";
}

Int_t DDPixelMeasVolume::calcClosestPointWith(const TVTrack& hel, TVector3& xx,
		Double_t& phi, Double_t eps) const {
	//based on TCylinder::calcXpointWith

//	std::cout<<"DDPixelMeasVolume::calcClosestPointWith(hel, xx="<<xx<<", phi, eps) called."<<std::endl;

   // This assumes nonzero B field.
   //
   // Copy helix parameters to local variables.
   //

   Double_t dr  = hel.GetDrho();
   Double_t fi0 = hel.GetPhi0();
   Double_t cpa = hel.GetKappa();
   Double_t dz  = hel.GetDz();
   Double_t tnl = hel.GetTanLambda();
   TVector3 X0  = hel.GetPivot();


   //
   // Check if charge is nonzero.
   //

   Int_t    chg = (Int_t)TMath::Sign(1.1,cpa/hel.GetPtoR());
   if (!chg) {
	  std::cerr << ">>>> Error >>>> DDPixelMeasVolume::calcClosestPointWith" << std::endl
		   << "      Kappa = 0 is invalid for a helix "          << std::endl;
	  return 0;
   }

   if(TMath::Abs(hel.GetRho())<eps) {
	   std::cerr << ">>>> Error >>>> DDPixelMeasVolume::calcClosestPointWith hel.GetRho()="<<hel.GetRho()<<std::endl;
	   return 0;
   }

   //
   // Project everything to XY plane and calculate crossing points.
   //

   Double_t rho  = hel.GetRho();
   Double_t rdr  = rho + dr;
   Double_t zdz  = X0.Z() + dz;
   Double_t csf0 = TMath::Cos(fi0);
   Double_t snf0 = TMath::Sin(fi0);
   Double_t xc  = X0.X() + rdr*csf0;
   Double_t yc  = X0.Y() + rdr*snf0;
   Double_t zc  = zdz;

   Double_t r   = TMath::Abs(rho);
   TCircle  c(r,xc,yc); //circle representing helix

   Double_t rv  = GetR();
   Double_t xcv = GetXc().X();
   Double_t ycv = GetXc().Y();
   TCircle  cv(rv,xcv,ycv); //circle representing surface

   TVector2 xxp[2];
   Int_t nx = c.CalcXingPointWith(cv, xxp);
   //there is a crossing point somewhere. Hit is probably not on surface, e.g. outside cylinder?
   if(nx) { streamlog_out(DEBUG2)<<"DDPixelMeasVolume::calcClosestPointWith terminated nx="<<nx<<" for volume at radius "<<GetR()<<std::endl; return 0;}
   //calculate closest point here instead of crossing
   nx=1;  //one closest point
   //circle c lies within cv
//   std::cout<<"calculate closest point for c="<<c.GetCenter()<<" & r="<<r<<", cv="<<cv.GetCenter()<<" & rv="<<rv<<std::endl;
   auto centerDistance=c.GetCenter()-cv.GetCenter();
   if(centerDistance.Mod2()<rv*rv) {//helix circle within surface circle
	     xxp[0]=c.GetCenter()+centerDistance.Unit()*r;
   } else if (centerDistance.Mod2()<r*r) { //surface circle within helix circle
	     xxp[0]=c.GetCenter()-centerDistance.Unit()*r;
   } else { //circle far from each other
	     xxp[0]=c.GetCenter()-centerDistance.Unit()*r;
   }


//   std::cout<<"found xxp= "<<xxp[0]<<std::endl;

   static const Double_t kPi     = TMath::Pi();
   static const Double_t kHalfPi = 0.5*TMath::Pi();
   static const Double_t kTwoPi  = 2.0*TMath::Pi();

   phi = 9999.;
   for (Int_t ix=0; ix<nx; ix++) {
	  Double_t x   = xxp[ix].X() - xc;
	  Double_t y   = xxp[ix].Y() - yc;
	  Double_t dfi = TMath::ATan2(y,x) - fi0 - kHalfPi*(1+chg);
	  while (dfi < -kPi) dfi += kTwoPi;
	  while (dfi >= kPi) dfi -= kTwoPi;
	  if (TMath::Abs(dfi) < TMath::Abs(phi)) {
		 phi = dfi;
		 xx.SetXYZ(xxp[ix].X(), xxp[ix].Y(), zc - rho*tnl*phi);
	  }
   }

//   std::cout<<"found xx= "<<xx<<std::endl;
   return -2; //return -2 specifically to indicate that this is a closest point, not a crossing, normally return number of intersections or -1 for error
}

/** Convert LCIO Tracker Hit to an DDPixelHit  */

DDVTrackHit* DDPixelMeasVolume::ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const {
  
  if ( ! trkhit) {
    streamlog_out(ERROR) << "DDPixelMeasVolume::ConvertLCIOTrkHit trkhit pointer is NULL" << std::endl;
    return NULL;
  }

  if(fMDim!=2) {
	  streamlog_out(WARNING)<<"fmDim = "<<fMDim<<", but expected 2"<<std::endl;
  }
  
  Double_t  x[2] { 0 /*d0*/, trkhit->getPosition()[2] /*z*/ };
  auto covMat=trkhit->getCovMatrix();
//  Double_t dx[2] { sqrt( (covMat[0]+covMat[2]) ), sqrt(covMat[5]) } ; //for pads
  Double_t dx[2] { sqrt( (covMat[0]+covMat[2])/2 ), sqrt(covMat[5]) } ; //also good for smearing in both x and y
  
  const TVector3 hit( trkhit->getPosition()[0], trkhit->getPosition()[1], trkhit->getPosition()[2]) ;
  bool hit_on_surface = IsOnSurface(hit);
  
  streamlog_out(DEBUG1) << "DDPixelMeasVolume::ConvertLCIOTrkHit DDPixelHit created"
  << " x = " << trkhit->getPosition()[0]
  << " y = " << trkhit->getPosition()[1]
  << " z = " << trkhit->getPosition()[2]
  << " onSurface = " << hit_on_surface
  << std::endl ;  
  
  // we sometimes get hits which are not on the surface when running with ddsim and the 'old' digitizers
  // for now just return the hit anyways ...
  if( ! hit_on_surface ){
    std::cout<< " DDPixelMeasVolume::ConvertLCIOTrkHit: hit is not on surface: " << *trkhit << std::endl ;
  }
  return new DDPixelHit( *this , x, dx, this->GetBz(), trkhit );
 
  
}

Int_t DDPixelMeasVolume::CalcXingPointWith(
		const TVTrack& hel,
		TVector3& xx,
		Double_t& phi,
		Int_t mode,
		Double_t eps) const {
	auto intersection=TCylinder::CalcXingPointWith(hel, xx, phi, mode, eps);
	streamlog_out(DEBUG2)<<"DDPixelMeasVolume::CalcXingPointWith -> TCylinder::CalcXingPointWith -> "<<intersection<<std::endl;
	return intersection;
}

//DEBUG from DDCylinderMeasLayer
/*
Int_t DDPixelMeasVolume::CalcXingPointWith(const TVTrack  &hel,
					     TVector3 &xx,
					     Double_t &phi,
					     Int_t     mode,
					     Double_t  eps) const {
  streamlog_out(DEBUG4)<<"using DDCylinderMeasLayer::CalcXingPointWith for DDPixelMeasVolume"<<std::endl;

  // check that direction has one of the correct values
  if( !( mode == 0 || mode == 1 || mode == -1) ) return -1 ;

 //fixme: currently the mode is ignored
  //       this might cause issues in track state propagation
  //       for curlers -> to be investigated ...


  // This assumes nonzero B field.
  //
  // Copy helix parameters to local variables.
  //

  Double_t dr     = hel.GetDrho();
  Double_t phi0   = hel.GetPhi0(); //
  Double_t kappa  = hel.GetKappa();
  Double_t rho    = hel.GetRho();
  Double_t omega  = 1.0 / rho;
  Double_t r      = TMath::Abs(rho);
  Double_t z0     = hel.GetDz();
  Double_t tanl   = hel.GetTanLambda();

  TVector3 ref_point = hel.GetPivot();


  // ---  Check if charge is nonzero.
  Int_t    chg = (Int_t)TMath::Sign(1.1,kappa);

  if (!chg) {

    streamlog_out(ERROR) << ">>>> Error >>>> DDCylinderMeasLayer::CalcXingPointWith" << std::endl
			 << "      Kappa = 0 is invalid for a helix "          << std::endl;
    return -1;
  }


  // =============== use code from aidaTT for computing the intersection with the plane =======

  double d0 = - dr ;
  double phi0_lcio =  toBaseRange( phi0 + M_PI/2. );

  // // aidaTT::trackParameters trkParam  ;

  // // trkParam.setTrackParameters( aidaTT::Vector5( omega/dd4hep::mm , tanl, phi0_lcio , d0*dd4hep::mm ,  z0 *dd4hep::mm )  ) ;

  // // //order defined in ./helpers/utilities.cc
  // // //  L3 type: [ Omega, tan(lambda), phi_0, d_0, z_0 ]
  // // trkParam.setReferencePoint( aidaTT::Vector3D( ref_point.X()*dd4hep::mm ,
  // // 						ref_point.Y()*dd4hep::mm,
  // // 						ref_point.Z()*dd4hep::mm ) ) ;

  // // aidaTT::trajectory traj( trkParam , 0 ) ;

  // // double s = 0. ;
  // // aidaTT::Vector3D xxV3 ;
  // // bool foundIntersect = traj._calculateIntersectionWithSurface( _surf , s , ( aidaTT::Vector2D*) 0 , &xxV3 );

  double s = 0. ;
  aidaTT::Vector3D xxV3 ;

  aidaTT::Vector5 hp( omega/dd4hep::mm , tanl, phi0_lcio , d0*dd4hep::mm ,  z0 *dd4hep::mm ) ;

  aidaTT::Vector3D rp( ref_point.X()*dd4hep::mm , ref_point.Y()*dd4hep::mm, ref_point.Z()*dd4hep::mm ) ;

  bool foundIntersect = aidaTT::intersectWithZCylinder( _surf, hp, rp, s, xxV3, mode , true ) ;


  if( foundIntersect ){

    s /= dd4hep::mm ;

    xx.SetXYZ( xxV3[0]/dd4hep::mm , xxV3[1]/dd4hep::mm,  xxV3[2]/dd4hep::mm) ;

    streamlog_out( DEBUG9 ) << " ++++  intersection found for surface : " << DDKalTest::CellIDEncoding::valueString(_surf->id()) << std::endl
     			    << "       at s = " << s
     			    << "       xx   = ( " << xx.X() << ", " << xx.Y() << ", " << xx.Z() << ") " << std::endl
              		    << " track parameters: " <<  aidaTT::trackParameters( hp, rp )
  			    << " mode: " << mode
     			    <<  std::endl ;

    phi = -s * omega ;

  } else {

    streamlog_out( DEBUG4 )<< " ++++ no intersection found for surface : " << DDKalTest::CellIDEncoding::valueString(_surf->id()) << std::endl
			    << " track parameters: " <<  aidaTT::trackParameters( hp, rp )
  			    << " mode : " << mode
  			    << std::endl ;

    return 0 ;
  }


  //=============================================================================================

  streamlog_out(DEBUG4) << "DDCylinderMeasLayer::CalcXingPointWith:on surface:" <<  IsOnSurface(xx)
  			<< "  (chg*phi*mode)<0: " <<  ((chg*phi*mode)<0)
  			<< " x = " << xx.X()
  			<< " y = " << xx.Y()
  			<< " z = " << xx.Z()
  			<< " r = " << xx.Perp()
  			<< " phi = " << xx.Phi()
  			<< " dphi = " <<  phi
  			<< " " << this->TVMeasLayer::GetName()
  			<< std::endl;


  if( mode!=0 && fabs(phi)>1.e-10){ // (+1,-1) = (fwd,bwd)
    if( chg*phi*mode > 0){
      return 0;
    }
  }



  //  return TCylinder::CalcXingPointWith( hel,xx,phi,mode,eps) ;

  if(!IsOnSurface(xx)) streamlog_out(DEBUG4)<<"DDPixelMeasVolume::CalcXingPointWith found point not on surface at "<<xx<<std::endl;
  return (IsOnSurface(xx) ? 1 : 0);

}

//*/


TKalMatrix DDPixelMeasVolume::XvToMv(const TVTrackHit& ht,
		const TVector3& xv) const {
	streamlog_out(ERROR)<<"DDPixelMeasVolume::XvToMv called without TKalTrackState"<<std::endl;
	throw "unexpected call";
}

