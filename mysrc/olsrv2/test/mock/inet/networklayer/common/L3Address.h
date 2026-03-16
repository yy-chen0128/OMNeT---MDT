#ifndef MOCK_L3ADDRESS_H
#define MOCK_L3ADDRESS_H

#include <string>
#include <iostream>

namespace inet {

class Ipv4Address {
public:
    static const Ipv4Address UNSPECIFIED_ADDRESS;
    static const Ipv4Address ALLONES_ADDRESS;
    static const Ipv4Address ALLONES_MASK;

    Ipv4Address() : addr_(0) {}
    explicit Ipv4Address(uint32_t addr) : addr_(addr) {}
    explicit Ipv4Address(const std::string& str) {
        // Simple hash to ensure uniqueness for distinct strings
        std::hash<std::string> hasher;
        addr_ = static_cast<uint32_t>(hasher(str));
    }

    uint32_t getInt() const { return addr_; }
    
    bool operator==(const Ipv4Address& other) const { return addr_ == other.addr_; }
    bool operator!=(const Ipv4Address& other) const { return addr_ != other.addr_; }
    bool operator<(const Ipv4Address& other) const { return addr_ < other.addr_; }

private:
    uint32_t addr_;
};

inline const Ipv4Address Ipv4Address::UNSPECIFIED_ADDRESS(0);
inline const Ipv4Address Ipv4Address::ALLONES_ADDRESS(0xFFFFFFFF);
inline const Ipv4Address Ipv4Address::ALLONES_MASK(0xFFFFFFFF);

class L3Address {
public:
    L3Address() : type_(0), v4_() {}
    L3Address(const Ipv4Address& ipv4) : type_(2), v4_(ipv4) {}
    L3Address(const char* str) : type_(2), v4_(str) {} // Simplified
    
    // Add missing constructor for Ipv4Address from string/int if implicit conversion is used
    // But better to be explicit.
    
    Ipv4Address toIpv4() const { return v4_; }
    
    bool operator==(const L3Address& other) const { return v4_ == other.v4_; }
    bool operator!=(const L3Address& other) const { return v4_ != other.v4_; }
    bool operator<(const L3Address& other) const { return v4_ < other.v4_; }

private:
    int type_;
    Ipv4Address v4_;
};

} // namespace inet

#endif // MOCK_L3ADDRESS_H
