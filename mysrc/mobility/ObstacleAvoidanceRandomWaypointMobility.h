/*
 * ObstacleAvoidanceRandomWaypointMobility.h
 *
 *  Created on: Jan 18, 2026
 *      Author: yychen
 */

#ifndef MYSRC_MOBILITY_OBSTACLEAVOIDANCERANDOMWAYPOINTMOBILITY_H_
#define MYSRC_MOBILITY_OBSTACLEAVOIDANCERANDOMWAYPOINTMOBILITY_H_

#include "inet/mobility/single/RandomWaypointMobility.h"
#include "inet/environment/common/PhysicalEnvironment.h"
#include "inet/common/geometry/common/Coord.h"
#include "inet/common/ModuleAccess.h"
#include "inet/common/geometry/base/ShapeBase.h"

using namespace inet;

class ObstacleAvoidanceRandomWaypointMobility : public RandomWaypointMobility {
  private:
    physicalenvironment::PhysicalEnvironment *physicalEnvironment;
    int MAX_RETRIES = 10;
  protected:
    /** @brief override setTargetPosition check obstacle */
    virtual void setTargetPosition() override;

    /** @brief check line */
    bool isPathBlocked(const Coord &start, const Coord &end);

    virtual void initialize(int stage) override;
};


//} // namespace inet



#endif /* MYSRC_MOBILITY_OBSTACLEAVOIDANCERANDOMWAYPOINTMOBILITY_H_ */
