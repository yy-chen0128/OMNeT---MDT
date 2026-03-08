/*
 * OLSR_RouteData.h
 *
 *  Created on: Mar 7, 2026
 *      Author: yychen
 */

#ifndef MYSRC_OLSR_OLSR_ROUTEDATA_H_
#define MYSRC_OLSR_OLSR_ROUTEDATA_H_

#include <set>
#include <tuple>
#include "inet/networklayer/common/L3Address.h"

namespace mysrc {
namespace olsr {
using namespace inet;
class INET_API OLSR_RouteData : public cObject
{
private:
    unsigned int distance;//hop count
public:
    OLSR_RouteData(){}
    virtual ~OLSR_RouteData(){}

    void setDistance(unsigned int d){distance = d;}
    unsigned int getDistance(){return distance;}
};


}
}
#endif /* MYSRC_OLSR_OLSR_ROUTEDATA_H_ */
