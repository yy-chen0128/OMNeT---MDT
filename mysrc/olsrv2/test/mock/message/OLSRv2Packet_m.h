#ifndef MOCK_OLSRV2PACKET_M_H
#define MOCK_OLSRV2PACKET_M_H

#include "../inet/common/packet/Packet.h"
#include "../inet/networklayer/common/L3Address.h"
#include <vector>

namespace inet {

enum Olsrv2PktType {
    Olsrv2Hello = 1,
    Olsrv2Tc = 2
};

class Olsrv2ControlPacket : public Chunk {
public:
    Olsrv2PktType packetType;
    void setPacketType(Olsrv2PktType type) { packetType = type; }
    void setChunkLength(long len) {}
};

class Olsrv2HelloPacket : public Olsrv2ControlPacket {
public:
    Olsrv2HelloPacket() { packetType = Olsrv2Hello; }
    
    void setOriginator(const L3Address& addr) { originator_ = addr; }
    L3Address getOriginator() const { return originator_; }
    
    void setHTime(double t) { hTime_ = t; }
    // Mock simtime_t as double wrapper
    struct SimTime {
        double t;
        SimTime(double v) : t(v) {}
        double dbl() const { return t; }
    } hTime_{0.0};
    
    SimTime getHTime() const { return hTime_; }
    
    void setWillingness(uint8_t w) { willingness_ = w; }
    uint8_t getWillingness() const { return willingness_; }
    
    void setMsgSeq(uint16_t s) { msgSeq_ = s; }
    uint16_t getMsgSeq() const { return msgSeq_; }
    
    void setNeighAddrsArraySize(size_t s) { neighAddrs_.resize(s); }
    size_t getNeighAddrsArraySize() const { return neighAddrs_.size(); }
    void setNeighAddrs(size_t i, const L3Address& addr) { neighAddrs_[i] = addr; }
    L3Address getNeighAddrs(size_t i) const { return neighAddrs_[i]; }

private:
    L3Address originator_;
    uint8_t willingness_;
    uint16_t msgSeq_;
    std::vector<L3Address> neighAddrs_;
};

class Olsrv2TcGroup : public Olsrv2ControlPacket {
public:
    Olsrv2TcGroup() { packetType = Olsrv2Tc; }

    void setGenTime(double t) { genTime_ = t; }
    
    void setOriginatorsArraySize(size_t s) { originators_.resize(s); }
    size_t getOriginatorsArraySize() const { return originators_.size(); }
    void setOriginators(size_t i, const L3Address& addr) { originators_[i] = addr; }
    L3Address getOriginators(size_t i) const { return originators_[i]; }
    
    void setAnsnsArraySize(size_t s) { ansns_.resize(s); }
    uint16_t getAnsns(size_t i) const { return ansns_[i]; }
    void setAnsns(size_t i, uint16_t v) { ansns_[i] = v; }
    
    void setSeqNumsArraySize(size_t s) { seqNums_.resize(s); }
    void setSeqNums(size_t i, uint16_t v) { seqNums_[i] = v; }
    
    void setHopCountsArraySize(size_t s) { hopCounts_.resize(s); }
    void setHopCounts(size_t i, uint8_t v) { hopCounts_[i] = v; }
    
    void setAdvCountsArraySize(size_t s) { advCounts_.resize(s); }
    int getAdvCounts(size_t i) const { return advCounts_[i]; }
    void setAdvCounts(size_t i, int v) { advCounts_[i] = v; }
    
    void setAdvertisedNeighborsArraySize(size_t s) { advertisedNeighbors_.resize(s); }
    L3Address getAdvertisedNeighbors(size_t i) const { return advertisedNeighbors_[i]; }
    void setAdvertisedNeighbors(size_t i, const L3Address& addr) { advertisedNeighbors_[i] = addr; }

private:
    double genTime_;
    std::vector<L3Address> originators_;
    std::vector<uint16_t> ansns_;
    std::vector<uint16_t> seqNums_;
    std::vector<uint8_t> hopCounts_;
    std::vector<int> advCounts_;
    std::vector<L3Address> advertisedNeighbors_;
};

} // namespace inet

#endif // MOCK_OLSRV2PACKET_M_H
