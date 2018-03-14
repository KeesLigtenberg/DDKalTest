/*
 * DDGetTPCFieldDescription.h
 *
 *  Created on: May 11, 2017
 *      Author: cligtenb
 */

#ifndef INCLUDE_DDKALTEST_DDGETTPCFIELDDESCRIPTION_H_
#define INCLUDE_DDKALTEST_DDGETTPCFIELDDESCRIPTION_H_

#include "DD4hep/LCDD.h"

inline std::string getDDFieldDescription() { return DD4hep::Geometry::LCDD::getInstance().idSpecification("TPCCollection").fieldDescription(); }


#endif /* INCLUDE_DDKALTEST_DDGETTPCFIELDDESCRIPTION_H_ */
