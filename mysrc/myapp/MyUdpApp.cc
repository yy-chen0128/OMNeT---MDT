/*
 * MyUdpApp.cc
 *
 *  Created on: Nov 25, 2025
 *      Author: yychen
 */


#include "MyUdpApp.h"
#include "inet/common/ModuleAccess.h"
#include "inet/common/TimeTag_m.h"
#include "inet/networklayer/common/FragmentationTag_m.h"
#include "inet/applications/base/ApplicationPacket_m.h"
#include "inet/common/Simsignals.h"

using namespace inet;

Define_Module(MyUdpApp);

void MyUdpApp::initialize(int stage)
{
    EV_ERROR<<"MyUdpApp::initialize\n";
    UdpBasicApp::initialize(stage); // parent initialize
    if (stage == INITSTAGE_LOCAL) {
        const char *ctrlPath = par("globalControllerModule");
        if (ctrlPath && ctrlPath[0]) {
            EV_ERROR<<"MyUdpApp::initialize: par(globalControllerModule)\n";
            cModule *m = getModuleByPath(ctrlPath);
            controller = check_and_cast<GlobalPrrController *>(m);
        }
        else {
            cModule *m = getModuleByPath("globalController");
            if (m){
                controller = check_and_cast<GlobalPrrController *>(m);
                EV_ERROR<<"MyUdpApp::initialize: getModuleByPath(globalController)\n";
            }
            else{
                controller = nullptr;
                EV_ERROR<<"MyUdpApp::initialize: no controller\n";
            }
        }
    }
}

void MyUdpApp::sendPacket()
{
    // 1) Create packet exactly as UdpBasicApp does
    std::ostringstream str;
    str << packetName << "-" << numSent;
    Packet *packet = new Packet(str.str().c_str());
    if (dontFragment)
        packet->addTag<FragmentationReq>()->setDontFragment(true);

    const auto& payload = makeShared<ApplicationPacket>();
    payload->setChunkLength(B(par("messageLength")));
    payload->setSequenceNumber(numSent);
    payload->addTag<CreationTimeTag>()->setCreationTime(simTime());
    packet->insertAtBack(payload);

    // 2) Choose destination address (reuse base method)
    L3Address destAddr = chooseDestAddr();

    // 3) Build PacketMeta and notify GlobalPrrController (if available)
    if (controller) {
        PacketMeta meta;
        meta.startTime = simTime();
        meta.bitLength = b(packet->getBitLength());

        // from node/radio/params get tx power / freq
        cModule *node = getContainingNode(this);
        if (node && node->hasPar("txPower"))
            meta.txPower_dBm = node->par("txPower").doubleValue();
        else
            meta.txPower_dBm = NAN;

        // TODO centerFrequency from radio/phy
        // meta.centerFrequency = Hz(...);

        // end time
        meta.endTime = meta.startTime + SIMTIME_ZERO; // or compute from bitrate
        EV_ERROR<<"MyUdpApp::sendPacket: controller->notifySend\n";
        controller->notifySend(getContainingNode(this), destAddr, meta);
    }
    else{
        EV_ERROR<<"MyUdpApp::sendPacket: no controller\n";
    }
    // 4) Emit signal and actually send the packet exactly as parent
    emit(packetSentSignal, packet);
    socket.sendTo(packet, destAddr, destPort);

    // 5) Update parent's counters
    numSent++;
}


