/*
 * MyUdpApp.h
 *
 *  Created on: Nov 25, 2025
 *      Author: yychen
 */

#ifndef MYSRC_MYAPP_MYUDPAPP_H_
#define MYSRC_MYAPP_MYUDPAPP_H_

#include "inet/applications/udpapp/UdpBasicApp.h"
#include "../datacollect/GlobalPrrController.h"

using namespace inet;

class MyUdpApp : public UdpBasicApp
{
  protected:
    GlobalPrrController *controller = nullptr;

  protected:
    virtual void initialize(int stage) override;
    virtual void sendPacket() override;
};




#endif /* MYSRC_MYAPP_MYUDPAPP_H_ */
