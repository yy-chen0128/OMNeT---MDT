#ifndef MYSRC_OLSRV2_RFC5444_ADDRESS_H
#define MYSRC_OLSRV2_RFC5444_ADDRESS_H

#include <cstdint>
#include <vector>

namespace mysrc::olsrv2::rfc5444 {

    struct AddressBlockView {
        uint8_t numAddr = 0;
        uint8_t flags = 0;
        std::vector<uint8_t> head;
        std::vector<uint8_t> tail;
        std::vector<std::vector<uint8_t>> mids;
        std::vector<uint8_t> prefixLengths;
    };

} // namespace inet::olsrv2::rfc5444

#endif
