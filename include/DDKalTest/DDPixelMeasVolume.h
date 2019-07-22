#ifndef DDPixelMeasVolume_H
#define DDPixelMeasVolume_H

#include "DDVMeasLayer.h"
#include <iostream>
#include <cmath>
#include "streamlog/streamlog.h"

#include "DDSurfaces/ISurface.h"

/** DDPixelMeasVolume provides measurment for a 3d pixel tpc hit
 *  using a DDSurfaces::ISurface.
 *  
 *  @author C Ligtenberg (Nikhef)
 *  @date Feb 2017
 *  @version $Id:$
 */
class DDPixelMeasVolume : public DDVMeasLayer, public TCylinder {
  
public:
  
  /// Constructor: initialize with Surface and B-field
  DDPixelMeasVolume(dd4hep::rec::ISurface* surf,
		      Double_t   Bz,
		      const Char_t    *name = "DDPixelMeasVolume") ;

  // Parent's pure virtual functions that must be implemented
  
  virtual ~DDPixelMeasVolume() {};

  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &ht,
                                  const TVector3   &xv ) const;
  virtual TKalMatrix XvToMv    (const TVector3   &xv,
		    					const TKalTrackState& a)   const;
  
  /** Global to Local coordinates */
  virtual TKalMatrix XvToMv    (const TVTrackHit &,
                                const TVector3   &xv,
							    const TKalTrackState& a)   const
  { return this->XvToMv(xv,a); }

  Int_t CalcXingPointWith(const TVTrack  &hel,
  					     TVector3 &xx,
  					     Double_t &phi,
  					     Int_t     mode,
  					     Double_t  eps) const;

  Int_t calcClosestPointWith(const TVTrack  &hel,
                                           TVector3 &xx,
                                           Double_t &phi,
  				            Double_t  eps = 1.e-8) const override;


  inline virtual       Bool_t     IsOnSurface(const TVector3 &xx) const;

  /** Local to Global coordinates */
  virtual TVector3   HitToXv   (const TVTrackHit &ht)   const;
  
  
  /** override calcDSDx for cases when hit is not on surface */
  TMatrixD CalcDSDx(const TVector3 &xx) const override;

  /** Calculate Projector Matrix */
  virtual void       CalcDhDa  (const TVTrackHit &ht,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)    const;
  
  virtual void       CalcDhDa  (const TVTrackHit &ht,
		  	  	  	  	  	  	const TKalTrackState &a,
                                const TVector3   &xv,
                                const TKalMatrix &dxphiada,
                                TKalMatrix &H)    const;

  /** Convert LCIO Tracker Hit to an DDVolumeHit  */
  virtual DDVTrackHit* ConvertLCIOTrkHit( EVENT::TrackerHit* trkhit) const ;
  
  /** Get the intersection and the CellID, needed for multilayers */
  virtual int getIntersectionAndCellID(const TVTrack  &hel,
                                       TVector3 &xx,
                                       Double_t &phi,
                                       Int_t    &CellID,
                                       Int_t     mode,
                                       Double_t  eps = 1.e-8) const {
    
    CellID = this->getCellIDs()[0]; // not multilayer
    return this->CalcXingPointWith(hel,xx,phi,mode,eps);
  }

  Double_t GetSortingPolicy() const { return fSortingPolicy; }
 
protected:
  Double_t fSortingPolicy;
  unsigned fMDim ;
  
private:
  
};

Bool_t DDPixelMeasVolume::IsOnSurface(const TVector3 &xx) const {

    //fg: leave this code for now - we are restricted to cylinders around the z-axis
    bool z = (xx.Z() >= GetZmin() && xx.Z() <= GetZmax());
    bool r = std::fabs( (xx-this->GetXc()).Perp() - this->GetR() ) < 1.e-3; // for very short, very stiff tracks this can be poorly defined, so we relax this here a bit to 1 micron

	streamlog_out(DEBUG2) << "DDCylinderMeasLayer IsOnSurface for " << this->TVMeasLayer::GetName() << " R =  " << this->GetR()
	<< "  GetZmin() = " << GetZmin() << " GetZmax() = " << GetZmax() << "xx.r = " << xx.Perp()
	<< " dr = " << std::fabs( (xx-this->GetXc()).Perp() - this->GetR() ) << " r = " << r << " z = " << z
	<< std::endl;

    return r && z;
  }

#endif
