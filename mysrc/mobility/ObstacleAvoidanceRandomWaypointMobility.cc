/*
 * ObstacleAvoidanceRandomWaypointMobility.cc
 *
 *  Created on: Jan 18, 2026
 *      Author: yychen
 */

#include"../../mysrc/mobility/ObstacleAvoidanceRandomWaypointMobility.h"

using namespace inet;

Define_Module(ObstacleAvoidanceRandomWaypointMobility);

void ObstacleAvoidanceRandomWaypointMobility::initialize(int stage)  {
    RandomWaypointMobility::initialize(stage);
}

bool ObstacleAvoidanceRandomWaypointMobility::isPathBlocked(const Coord& start,
                                                           const Coord& end)
{
    if (!physicalEnvironment)
        return false;

    for (int i = 0; i < physicalEnvironment->getNumObjects(); i++) {
        const auto *obj = physicalEnvironment->getObject(i);
        const auto *shape = obj->getShape();

        // --- 1. transfer coord ---
        Quaternion invRot = obj->getOrientation().inverse();

        Coord localStart = invRot.rotate(start - obj->getPosition());
        Coord localEnd   = invRot.rotate(end   - obj->getPosition());

        LineSegment localLine(localStart, localEnd);

        // --- 2. check collision ---
        Coord i1, i2, n1, n2;
        if (shape->computeIntersection(localLine, i1, i2, n1, n2)) {
            return true;   //collision
        }
    }

    return false; // not collision
}


void ObstacleAvoidanceRandomWaypointMobility::setTargetPosition()
{
    // wait
    if (nextMoveIsWait) {
        nextChange = simTime() + waitTimeParameter->doubleValue();
        nextMoveIsWait = false;
        return;
    }

    // get physicalEnvironment
    if (!physicalEnvironment) {
        physicalEnvironment =
            check_and_cast<physicalenvironment::PhysicalEnvironment*>(
                getParentModule()->getParentModule()->getSubmodule("physicalEnvironment"));
    }

    Coord nextPosition;
    bool found = false;

    // try to select proper dist
    for (int retry = 0; retry < MAX_RETRIES; retry++) {
        nextPosition = getRandomPosition();

        //
        if (nextPosition == lastPosition)
            continue;

        if (!isPathBlocked(lastPosition, nextPosition)) {
            found = true;
            break;
        }
    }

    // --- fallback: wait ---
    if (!found) {
        nextPosition = lastPosition;
    }

    // --- same as RandomWaypoint  ---
    targetPosition = nextPosition;

    double speed = speedParameter->doubleValue();
    double distance = lastPosition.distance(targetPosition);
    simtime_t travelTime = found ? distance / speed : 1;

    nextChange = simTime() + travelTime;
    nextMoveIsWait = hasWaitTime;
}


