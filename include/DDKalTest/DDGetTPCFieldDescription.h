/*
 * DDGetTPCFieldDescription.h
 *
 *  Created on: May 11, 2017
 *      Author: cligtenb
 */

#ifndef INCLUDE_DDKALTEST_DDGETTPCFIELDDESCRIPTION_H_
#define INCLUDE_DDKALTEST_DDGETTPCFIELDDESCRIPTION_H_

#include "DD4hep/Detector.h"
#include <UTIL/LCTrackerConf.h>

//canoncial is "subdet:5,side:-2,layer:9,module:8,sensor:8" but might be changed!

//inline std::string getDDFieldDescription() { return dd4hep::Detector::getInstance().idSpecification("TPCCollection").fieldDescription(); }
inline std::string getDDFieldDescription(int id) {
	id&=0x7f; //select first seven bits
	if( id==100 || id==36 || id==4 ) //is TPC? (4=no side encoding?)
		return dd4hep::Detector::getInstance().idSpecification("TPCCollection").fieldDescription();
	else
		return UTIL::LCTrackerCellID::encoding_string();
}

#endif /* INCLUDE_DDKALTEST_DDGETTPCFIELDDESCRIPTION_H_ */
