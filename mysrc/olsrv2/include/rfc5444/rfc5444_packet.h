//
// Created by GuoYudi on 2026/3/11.
//

#ifndef MYSRC_OLSRV2_RFC5444_PACKET_H
#define MYSRC_OLSRV2_RFC5444_PACKET_H

#include <cstdint>
#include <vector>
#include "rfc5444.h"

namespace mysrc::olsrv2::rfc5444 {

    struct PacketHeader {
        uint8_t version = 0;  // RFC5444 version
        uint8_t flags = 0;
        uint16_t packetSeqNo = 0;
    };

    struct DecodedPacket {
        PacketHeader header;
        std::vector<TlvValue> packetTlvs;
        std::vector<DecodedMessage> messages;
    };

} // namespace mysrc::olsrv2::rfc5444

#endif
