/*
 * RpcBridge.h
 *
 *  Created on: Oct 9, 2025
 *      Author: yychen
 */

#ifndef MYSRC_RPC_RPCBRIDGE_H_
#define MYSRC_RPC_RPCBRIDGE_H_

#include <nlohmann/json.hpp>
#include <memory>
#include <functional>

using json = nlohmann::json;

class RpcBridge {
public:
    RpcBridge(int port = 9002);
    ~RpcBridge();

    void start();
    void stop();
    void enqueueOutgoingLog(const json &logMsg);
    bool tryPopIncomingRequest(json &out);
    void setConnectionHook(std::function<void(size_t)> cb);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

extern RpcBridge* gRpcBridge;

#endif
