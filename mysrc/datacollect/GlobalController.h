/*
 * GlobalController.h
 *
 *  Created on: Nov 24, 2025
 *      Author: yychen
 */

#ifndef MYSRC_DATACOLLECT_GLOBALCONTROLLER_H_
#define MYSRC_DATACOLLECT_GLOBALCONTROLLER_H_

#include <omnetpp.h>
#include "inet/physicallayer/wireless/common/contract/packetlevel/IRadioMedium.h"
#include "inet/physicallayer/wireless/common/analogmodel/packetlevel/ScalarAnalogModel.h"
#include "inet/physicallayer/wireless/common/pathloss/FreeSpacePathLoss.h"
#include "inet/physicallayer/wireless/common/obstacleloss/DielectricObstacleLoss.h"
#include "inet/physicallayer/wireless/common/backgroundnoise/IsotropicScalarBackgroundNoise.h"

using namespace omnetpp;
using namespace inet;
using namespace inet::physicallayer;

class GlobalController : public cSimpleModule
{
  protected:
    IRadioMedium *radioMedium = nullptr;
    IPathLoss *pathLoss = nullptr;
    IObstacleLoss *obstacleLoss = nullptr;
    IBackgroundNoise *backgroundNoise = nullptr;
    IPropagation *propagation = nullptr;

    long totalTx = 0;
    long totalTheoreticallyRx = 0;

  protected:
    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override { throw cRuntimeError("GlobalController does not accept messages"); }
  public:
    /** Called by apps when they send a packet. */
    void notifyTx(cModule *senderNode, double txPower_dBm, simtime_t startTime, simtime_t endTime, Hz centerFrequency);

    virtual void finish() override;
};



#endif /* MYSRC_DATACOLLECT_GLOBALCONTROLLER_H_ */
