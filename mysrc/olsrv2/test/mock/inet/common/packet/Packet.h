#ifndef MOCK_PACKET_H
#define MOCK_PACKET_H

#include <string>
#include <memory>

namespace inet {

class Chunk {
public:
    virtual ~Chunk() = default;
};

template <typename T>
std::shared_ptr<T> makeShared() {
    return std::make_shared<T>();
}

template <typename T>
std::shared_ptr<const T> dynamicPtrCast(std::shared_ptr<const Chunk> chunk) {
    return std::dynamic_pointer_cast<const T>(chunk);
}

class Packet {
public:
    Packet(const char* name = nullptr) {}
    
    template <typename T>
    std::shared_ptr<const T> peekAtFront() const {
        // In real OMNeT++, this returns a chunk.
        // Here we just need to verify type.
        // We'll store a shared_ptr<Chunk> in Packet for mocking.
        return std::dynamic_pointer_cast<const T>(chunk_);
    }

    void insertAtBack(std::shared_ptr<Chunk> chunk) {
        chunk_ = chunk;
    }
    
    template <typename T>
    T* getTag() { return nullptr; } // Mock tag retrieval

private:
    std::shared_ptr<Chunk> chunk_;
};

class L3AddressInd {
public:
    L3Address getSrcAddress() const { return L3Address(); }
};

// Mock B unit
inline long B(long bytes) { return bytes; }

} // namespace inet

#endif // MOCK_PACKET_H
