/*
 * OmnetRpcLogger.h
 *
 *  Created on: Oct 9, 2025
 *      Author: yychen
 */

#ifndef MYSRC_RPC_OMNETRPCLOGGER_H_
#define MYSRC_RPC_OMNETRPCLOGGER_H_

#include <omnetpp.h>
#include <nlohmann/json.hpp>
#include <string>
#include "RpcBridge.h"
using namespace omnetpp;
using json = nlohmann::json;

// NOTE: We intentionally do NOT define `gRpcBridge` in this file. If your project
// needs a global RpcBridge pointer, keep its single definition in the dedicated
// implementation file for RpcBridge (or another globals.cpp). This header may
// include RpcBridge.h elsewhere if needed but does not create/define the global.

// File-local logging helpers (declared here so other translation units can call them)
void initRpcLogFile(const std::string &path = "logs_run1.ndjson", int flushInterval = 200);
void writeJsonToLocalFile(const json &j);
void shutdownRpcLogFile();

class OmnetRpcLogger : public cSimpleModule {
  private:
    std::string logPath;
    int fileFlushInterval = 200;
    RpcBridge *rpc = nullptr;
  protected:
    virtual void initialize() override;
    virtual void handleMessage(cMessage *msg) override;
    virtual void finish() override;
};

// Keep Define_Module here for compatibility with the original layout. If you prefer
// to move it into the .cpp, you can, but leaving it here preserves the original
// compile-time behavior.
Define_Module(OmnetRpcLogger);

#endif /* MYSRC_RPC_OMNETRPCLOGGER_H_ */
