/*
 * DDPixelHit.h
 *
 *  Created on: Feb 28, 2017
 *      Author: cligtenb
 */

#ifndef INCLUDE_DDKALTEST_DDPIXELHIT_H_
#define INCLUDE_DDKALTEST_DDPIXELHIT_H_

#include "DDKalTest/DDVTrackHit.h"


class DDPixelHit: public DDVTrackHit {
public:
	  /** Constructor Taking x,y,z coordinates and associated measurement layer, with bfield */
	DDPixelHit(const TVMeasLayer &measurementVolume, Double_t *x, Double_t *dx,
	                 Double_t bfield, EVENT::TrackerHit* trkhit );

//	  virtual const TVMeasLayer & GetMeasLayer() const ;

	  // TVTrackHit's pure virtuals that must be implemented

	  /** Global to Local coordinates */
	  virtual TKalMatrix XvToMv(const TVector3 &xv, Double_t t0) const ;
	  virtual TKalMatrix XvToMv(const TVector3 &xv, Double_t t0, const TKalTrackState& a) const;

	  /** Print Debug information */
	  virtual void DebugPrint(Option_t *opt = "") const;

};

#endif /* INCLUDE_DDKALTEST_DDPIXELHIT_H_ */
