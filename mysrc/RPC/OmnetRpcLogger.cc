/*
 * OmnetRpcLogger.cc
 *
 *  Created on: Oct 9, 2025
 *      Author: yychen
 */
#include "OmnetRpcLogger.h"
#include <fstream>
#include <mutex>
#include <iostream>
#include"log2vis.h"

// Include RpcBridge.h only if needed for declarations; we must NOT define gRpcBridge here.
// If RpcBridge.h declares `extern RpcBridge* gRpcBridge;` that's fine; we intentionally
// do not provide a definition in this translation unit.
 // optional: keeps header-level compatibility, but doesn't define gRpcBridge


// File-static globals used by the logging helpers
static std::ofstream g_rpc_log_ofs;
static bool g_rpc_log_file_open = false;
static std::string g_rpc_log_path = "logs_run1.ndjson";
static int g_rpc_unflushed = 0;
static int g_rpc_fileFlushInterval = 200;
static std::mutex g_rpc_file_mutex;

void initRpcLogFile(const std::string &path, int flushInterval) {
    std::lock_guard<std::mutex> lk(g_rpc_file_mutex);
    if (g_rpc_log_file_open) {
        g_rpc_log_ofs.flush();
        g_rpc_log_ofs.close();
        g_rpc_log_file_open = false;
    }


    g_rpc_log_path = path.empty() ? std::string("logs_run1.ndjson") : path;
    g_rpc_fileFlushInterval = flushInterval > 0 ? flushInterval : 200;
    g_rpc_unflushed = 0;


    g_rpc_log_ofs.open(g_rpc_log_path, std::ios::out | std::ios::trunc);
    if (g_rpc_log_ofs.is_open()) {
        g_rpc_log_file_open = true;
    // optionally write a start record
    try {
        json startJ = { {"type", "start"}, {"timestamp", simTime().dbl()} };
        g_rpc_log_ofs << startJ.dump() << '\n';
        ++g_rpc_unflushed;
        if (g_rpc_unflushed >= g_rpc_fileFlushInterval) {
            g_rpc_log_ofs.flush();
            g_rpc_unflushed = 0;
        }
    } catch (...) {
    // if simTime() or json throws (unlikely), ignore and leave file open
    }
    } else {
        g_rpc_log_file_open = false;
        EV_WARN << "OmnetRpcLogger: failed to open log file '" << g_rpc_log_path << "'\n";
    }
}


void writeJsonToLocalFile(const json &j) {
    std::lock_guard<std::mutex> lk(g_rpc_file_mutex);
    if (!g_rpc_log_file_open) return;
    try {
        g_rpc_log_ofs << j.dump() << '\n';
        ++g_rpc_unflushed;
        if (g_rpc_unflushed >= g_rpc_fileFlushInterval) {
            g_rpc_log_ofs.flush();
            g_rpc_unflushed = 0;
        }
    } catch (...) {
    // best-effort: on failure, silently ignore to avoid crashing simulation
    }
}

void shutdownRpcLogFile() {
    std::lock_guard<std::mutex> lk(g_rpc_file_mutex);
    if (!g_rpc_log_file_open) return;
    try {
        json finishJ = { {"type", "logger_shutdown"}, {"timestamp", simTime().dbl()} };
        g_rpc_log_ofs << finishJ.dump() << '\n';
        g_rpc_log_ofs.flush();
    } catch (...) {}
    g_rpc_log_ofs.close();
    g_rpc_log_file_open = false;
}


// ---------------- OmnetRpcLogger implementation ----------------


void OmnetRpcLogger::initialize() {
    // Keep compatibility with existing NED params. If your NED still contains
    // rpcPort / rpcPollInterval they are ignored here (we intentionally do not
    // create a RpcBridge or schedule polling). We support optional parameters
    // for the local logfile path and flush interval.
    if (hasPar("rpcLogPath")) {
        logPath = par("rpcLogPath").stdstringValue();
    } else {
        logPath = "logs_run1.ndjson";
    }
    if (hasPar("rpcFileFlushInterval")) {
        fileFlushInterval = par("rpcFileFlushInterval").intValue();
        if (fileFlushInterval <= 0) fileFlushInterval = 200;
    }
    else {
        fileFlushInterval = 200;
    }
    int port;
    if (hasPar("rpcPort")) {
        port = par("rpcPort");

    } else {
        port = 9002;
        EV_INFO << "rpcPort param not present  RpcBridge not created (local logging only)\n";
    }
    try {
        rpc = new RpcBridge(port);
        gRpcBridge = rpc;

        rpc->start();

        EV_INFO << "RPC Logger started on port " << port << "\n";
    } catch (const std::exception &e) {
        EV_WARN << "Failed to create/start RpcBridge: " << e.what() << "\n";
        if (rpc) {
            delete rpc;
            rpc = nullptr;
            gRpcBridge = nullptr;
        }
    }
    initRpcLogFile(logPath, fileFlushInterval);
    EV_INFO << "OmnetRpcLogger initialized, writing NDJSON to '" << logPath << "'\n";
}

void OmnetRpcLogger::handleMessage(cMessage *msg) {
// This module no longer uses any self-messages (no polling); delete any
// message received to avoid leaks in case something is accidentally sent.
    delete msg;
}

void OmnetRpcLogger::finish() {
    try {
        json finishJ = { {"type", "finish"}, {"timestamp", simTime().dbl()}, {"modulePath", getFullPath()} };
        writeJsonToLocalFile(finishJ);
    } catch (...) {}
    shutdownRpcLogFile();
    if (rpc) {
        rpc->stop();
        delete rpc;
        rpc = nullptr;
        gRpcBridge = nullptr;
    }

    // preserve original behaviour: call logFinish if available
    try { logFinish(this); } catch (...) {}
}
