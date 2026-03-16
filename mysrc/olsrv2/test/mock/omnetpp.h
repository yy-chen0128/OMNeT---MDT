#ifndef MOCK_OMNETPP_H
#define MOCK_OMNETPP_H

#include <iostream>
#include <string>

namespace omnetpp {

class cObject {
public:
    virtual ~cObject() = default;
};

class cMessage : public cObject {
public:
    cMessage(const char* name = nullptr) {}
    bool isSelfMessage() const { return false; }
};

class cSimpleModule : public cObject {
public:
    virtual void initialize(int stage) {}
    virtual void handleMessage(cMessage *msg) {}
    virtual int numInitStages() const { return 1; }
    
    // Mock parameter access
    virtual double par(const char* name) const { return 0.0; }
    
    // Mock scheduling
    void scheduleAfter(double delay, cMessage *msg) {}
    void cancelAndDelete(cMessage *msg) {}
};

// Define Module macro mock
#define Define_Module(CLASSname)

} // namespace omnetpp

using namespace omnetpp;

#endif // MOCK_OMNETPP_H
