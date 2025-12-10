/*
 * GlobalPrrController.h
 *
 *  Created on: Nov 25, 2025
 *      Author: yychen
 */

#ifndef MYSRC_DATACOLLECT_GLOBALPRRCONTROLLER_H_
#define MYSRC_DATACOLLECT_GLOBALPRRCONTROLLER_H_


#include <omnetpp.h>
#include <vector>
#include <unordered_map>

#include "inet/physicallayer/wireless/common/contract/packetlevel/IRadioMedium.h"
#include "inet/physicallayer/wireless/common/radio/packetlevel/Radio.h"
#include "inet/common/ModuleAccess.h"
#include "inet/common/INETDefs.h"
#include "inet/networklayer/common/L3Address.h"

using namespace omnetpp;
using namespace inet;
using namespace inet::physicallayer;

struct PacketMeta {
    // minimal meta needed; can extend
    double txPower_dBm = NAN;
    Hz centerFrequency = Hz(NaN);
    b bitLength = B(-1);
    double bitrate_bps = NAN;
    simtime_t startTime;
    simtime_t endTime;
};

class GlobalPrrController : public cSimpleModule
{
  public:
    enum Mode { FAST = 0, PHYSICAL = 1 };

  protected:
    std::vector<cModule*> nodes;
    int nodeVectorSize = 0;

    IRadioMedium *radioMedium = nullptr;
    IPathLoss *pathLoss = nullptr;
    IObstacleLoss *obstacleLoss = nullptr;
    IBackgroundNoise *backgroundNoise = nullptr;
    IPropagation *propagation = nullptr;
    IMediumLimitCache *mediumLimitCache = nullptr;

    Mode mode = FAST;

    // stats
    long totalTx = 0;
    long packetsWithAtLeastOneTheoreticalRx = 0;
    long theoreticalTotalReceives = 0;

    // cached list of radios (populated at init)
    std::vector<Radio*> radios;
  protected:
      simsignal_t reachableSignal;
      simsignal_t receiverCountSignal;
      simsignal_t theoreticalPRRSignal;

  protected:
    virtual int numInitStages() const override { return NUM_INIT_STAGES; }
    virtual void initialize(int stage) override;
    virtual void handleMessage(cMessage *msg) override { throw cRuntimeError("GlobalPrrController does not accept messages"); }
    virtual void finish() override;

  public:
    // API for apps to call
    void notifySend(cModule *senderNode, L3Address destAddr, const PacketMeta& meta);

  protected:
    cModule* findNodeByAddress(const L3Address &addr);
    bool canRadioReachRadio_Fast(Radio *txRadio, Radio *rxRadio, const PacketMeta &meta);
    bool canRadioReachRadio_Physical(Radio *txRadio, Radio *rxRadio, const PacketMeta &meta);
    bool canNodeReachNode(cModule *srcNode, cModule *dstNode, const PacketMeta &meta);
    bool isDestinationReachableFromSource(cModule *srcNode, cModule *dstNode, const PacketMeta &meta);

    Radio *findRadioForNode(cModule *node);
    void collectRadios();
    // mode A (fast): compute using pathLoss/obstacle/backgroundNoise (no addTransmission)
    void handleNotifyFast(Radio *txRadio, cModule *senderNode, L3Address destAddr, const PacketMeta &meta);
    // mode B (physical): create transmission, add to medium, ask medium about reception
    void handleNotifyPhysical(Radio *txRadio, cModule *senderNode, L3Address destAddr, const PacketMeta &meta);
};

Define_Module(GlobalPrrController);




#endif /* MYSRC_DATACOLLECT_GLOBALPRRCONTROLLER_H_ */
