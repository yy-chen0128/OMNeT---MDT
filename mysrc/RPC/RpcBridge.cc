/*
 * RpcBridge.cc
 *
 *  Created on: Oct 9, 2025
 *      Author: yychen
 */

#include "RpcBridge.h"

#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>
#include <nlohmann/json.hpp>

#include <thread>
#include <mutex>
#include <set>
#include <queue>
#include <atomic>
#include <functional>
#include <iostream>
#include <sstream>

using json = nlohmann::json;
using namespace std;

typedef websocketpp::server<websocketpp::config::asio> WSServer;

// Impl (PIMPL) hides heavy websocketpp/template details from the public header
struct RpcBridge::Impl {
    // websocket server and io
    WSServer server;
    std::thread serverThread;

    // configuration
    int port = 9002;

    // connection management (only touched from IO thread or protected by connMutex)
    std::mutex connMutex;
    std::set<websocketpp::connection_hdl, std::owner_less<websocketpp::connection_hdl>> connections;

    // outgoing queue: pushed from simulation thread, consumed on IO thread
    std::mutex outMutex;
    std::queue<std::string> outgoingQueue;

    // incoming queue: pushed from IO thread, popped by simulation thread
    std::mutex inMutex;
    std::queue<json> incomingQueue;

    // running flag
    std::atomic<bool> running{false};

    // optional connection hook
    std::function<void(size_t)> connHook;

    Impl(int p) : port(p) {
        // initialize server
        server.init_asio();
        server.set_reuse_addr(true);

        // bind handlers
        server.set_open_handler([this](websocketpp::connection_hdl hdl){ on_open(hdl); });
        server.set_close_handler([this](websocketpp::connection_hdl hdl){ on_close(hdl); });
        server.set_message_handler([this](websocketpp::connection_hdl hdl, WSServer::message_ptr msg){ on_message(hdl, msg); });

        // silence websocketpp logging (optional)
        server.clear_access_channels(websocketpp::log::alevel::all);
        server.clear_error_channels(websocketpp::log::elevel::all);

        // start listening (note: actual run happens in separate thread)
        try {
            server.listen(port);
            server.start_accept();
        } catch (const std::exception &e) {
            std::cerr << "RpcBridge::Impl() listen/start_accept exception: " << e.what() << std::endl;
        }
    }

    ~Impl() {
        stop();
    }

    // start server loop
    void start() {
        if (running.exchange(true)) return; // already running
        serverThread = std::thread([this](){
            try {
                server.run();
            } catch (const std::exception &e) {
                std::cerr << "RpcBridge server run exception: " << e.what() << std::endl;
            }
        });
    }

    // stop server and join thread
    void stop() {
        if (!running.exchange(false)) return; // already stopped

        // close clients
        {
            std::lock_guard<std::mutex> lk(connMutex);
            for (auto &hdl : connections) {
                try { server.close(hdl, websocketpp::close::status::going_away, "shutdown"); }
                catch(...) {}
            }
            connections.clear();
        }

        try { server.stop_listening(); } catch(...) {}
        try { server.stop(); } catch(...) {}

        if (serverThread.joinable()) serverThread.join();
    }

    // enqueue outgoing message (thread-safe). Schedules a flush on IO thread.
    void enqueueOutgoingLog(const json &logMsg) {
        std::string s = logMsg.dump();
        {
            std::lock_guard<std::mutex> lk(outMutex);
            outgoingQueue.push(s);
        }
        // schedule flush on io_service so actual sends happen in IO thread
        try {
            server.get_io_service().post([this](){ flushOutgoing(); });
        } catch(...) {
            // fallback: ignore
        }
    }

    // flush outgoing messages (must run on IO thread)
    void flushOutgoing() {
        std::queue<std::string> localQ;
        {
            std::lock_guard<std::mutex> lk(outMutex);
            std::swap(localQ, outgoingQueue);
        }

        std::lock_guard<std::mutex> lk(connMutex);
        while (!localQ.empty()) {
            const std::string &msg = localQ.front();
            for (auto &hdl : connections) {
                try {
                    server.send(hdl, msg, websocketpp::frame::opcode::text);
                } catch (const websocketpp::exception &e) {
                    // ignore send errors; on_close will cleanup
                } catch (...) {}
            }
            localQ.pop();
        }
    }

    // websocket handlers (run on IO thread)
    void on_open(websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(connMutex);
        connections.insert(hdl);
        if (connHook) connHook(connections.size());
    }

    void on_close(websocketpp::connection_hdl hdl) {
        std::lock_guard<std::mutex> lk(connMutex);
        connections.erase(hdl);
        if (connHook) connHook(connections.size());
    }

    void on_message(websocketpp::connection_hdl hdl, WSServer::message_ptr msg) {
        try {
            json j = json::parse(msg->get_payload());
            // push into incoming queue to be processed by simulation thread
            {
                std::lock_guard<std::mutex> lk(inMutex);
                incomingQueue.push(j);
            }
        } catch (const std::exception &e) {
            // reply error immediately (we are on IO thread)
            json resp = { {"type","error"}, {"what", e.what()} };
            try { server.send(hdl, resp.dump(), websocketpp::frame::opcode::text); } catch(...) {}
        }
    }

    // try pop incoming request (called by simulation thread)
    bool tryPopIncomingRequest(json &out) {
        std::lock_guard<std::mutex> lk(inMutex);
        if (incomingQueue.empty()) return false;
        out = incomingQueue.front(); incomingQueue.pop();
        return true;
    }

    void setConnectionHook(std::function<void(size_t)> cb) { connHook = cb; }
};

// --- RpcBridge public forwarding to Impl ---

RpcBridge::RpcBridge(int port) {
    impl_.reset(new Impl(port));
}

RpcBridge::~RpcBridge() {
    if (impl_) {
        impl_->stop();
        impl_.reset();
    }
}

void RpcBridge::start() {
    if (!impl_) return;
    impl_->start();
}

void RpcBridge::stop() {
    if (!impl_) return;
    impl_->stop();
}

void RpcBridge::enqueueOutgoingLog(const json &logMsg) {
    if (!impl_) return;
    impl_->enqueueOutgoingLog(logMsg);
}

bool RpcBridge::tryPopIncomingRequest(json &out) {
    if (!impl_) return false;
    return impl_->tryPopIncomingRequest(out);
}

void RpcBridge::setConnectionHook(std::function<void(size_t)> cb) {
    if (!impl_) return;
    impl_->setConnectionHook(cb);
}

// global pointer is left unchanged; user code may still use gRpcBridge
RpcBridge* gRpcBridge = nullptr;


