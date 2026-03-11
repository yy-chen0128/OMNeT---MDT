#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mysrc::olsrv2::rfc5444 {

// RFC5444 parser/serializer return codes used by module internals and callbacks.
enum class Result : int {
    ResultMax = 5,
    DropPacket = 5,
    DropMessage = 4,
    DropMsgButForward = 3,
    DropAddress = 2,
    DropTlv = 1,

    Okay = 0,
    UnsupportedVersion = -1,
    EndOfBuffer = -2,
    BadTlvIdxFlags = -3,
    BadTlvValueFlags = -4,
    BadTlvLength = -5,
    OutOfMemory = -6,
    EmptyAddrBlock = -7,
    BadMsgTailFlags = -8,
    BadMsgPrefixFlags = -9,
    DuplicateTlv = -10,
    OutOfAddrTlvMem = -11,
    MtuTooSmall = -12,
    NoMsgCreator = -13,
    FwMessageTooLong = -14,
    FwBadSize = -15,
    TooLarge = -16,
    FwBadTransform = -17,
    BadAddrHeadLength = -18,
    BadAddrTailLength = -19,
    BadTlvIndex = -20,
    ResultMin = -20,
};

const char *toString(Result result);

inline constexpr int MAX_ADDRLEN = 16;

inline constexpr uint64_t RFC5497_TIMETLV_MIN = 1;
inline constexpr uint64_t RFC5497_TIMETLV_MAX = 3932160000ull;

inline constexpr uint32_t RFC7181_METRIC_MIN = 1;
inline constexpr uint32_t RFC7181_METRIC_MAX = 16776960u;
inline constexpr uint32_t RFC7181_METRIC_INFINITE = 0xffffffffu;

inline constexpr uint8_t PKT_FLAGMASK = 0x0f;
inline constexpr uint8_t PKT_FLAG_SEQNO = 0x08;
inline constexpr uint8_t PKT_FLAG_TLV = 0x04;

inline constexpr uint8_t MSG_FLAG_ORIGINATOR = 0x80;
inline constexpr uint8_t MSG_FLAG_HOPLIMIT = 0x40;
inline constexpr uint8_t MSG_FLAG_HOPCOUNT = 0x20;
inline constexpr uint8_t MSG_FLAG_SEQNO = 0x10;
inline constexpr uint8_t MSG_FLAG_ADDRLENMASK = 0x0f;

inline constexpr uint8_t ADDR_FLAG_HEAD = 0x80;
inline constexpr uint8_t ADDR_FLAG_FULLTAIL = 0x40;
inline constexpr uint8_t ADDR_FLAG_ZEROTAIL = 0x20;
inline constexpr uint8_t ADDR_FLAG_SINGLEPLEN = 0x10;
inline constexpr uint8_t ADDR_FLAG_MULTIPLEN = 0x08;

inline constexpr uint8_t TLV_FLAG_TYPEEXT = 0x80;
inline constexpr uint8_t TLV_FLAG_SINGLE_IDX = 0x40;
inline constexpr uint8_t TLV_FLAG_MULTI_IDX = 0x20;
inline constexpr uint8_t TLV_FLAG_VALUE = 0x10;
inline constexpr uint8_t TLV_FLAG_EXTVALUE = 0x08;
inline constexpr uint8_t TLV_FLAG_MULTIVALUE = 0x04;

// Legacy aliases to reduce migration friction in existing headers/users.
inline constexpr Result RFC5444_RESULT_MAX = Result::ResultMax;
inline constexpr Result RFC5444_DROP_PACKET = Result::DropPacket;
inline constexpr Result RFC5444_DROP_MESSAGE = Result::DropMessage;
inline constexpr Result RFC5444_DROP_MSG_BUT_FORWARD = Result::DropMsgButForward;
inline constexpr Result RFC5444_DROP_ADDRESS = Result::DropAddress;
inline constexpr Result RFC5444_DROP_TLV = Result::DropTlv;
inline constexpr Result RFC5444_OKAY = Result::Okay;
inline constexpr Result RFC5444_UNSUPPORTED_VERSION = Result::UnsupportedVersion;
inline constexpr Result RFC5444_END_OF_BUFFER = Result::EndOfBuffer;
inline constexpr Result RFC5444_BAD_TLV_IDXFLAGS = Result::BadTlvIdxFlags;
inline constexpr Result RFC5444_BAD_TLV_VALUEFLAGS = Result::BadTlvValueFlags;
inline constexpr Result RFC5444_BAD_TLV_LENGTH = Result::BadTlvLength;
inline constexpr Result RFC5444_OUT_OF_MEMORY = Result::OutOfMemory;
inline constexpr Result RFC5444_EMPTY_ADDRBLOCK = Result::EmptyAddrBlock;
inline constexpr Result RFC5444_BAD_MSG_TAILFLAGS = Result::BadMsgTailFlags;
inline constexpr Result RFC5444_BAD_MSG_PREFIXFLAGS = Result::BadMsgPrefixFlags;
inline constexpr Result RFC5444_DUPLICATE_TLV = Result::DuplicateTlv;
inline constexpr Result RFC5444_OUT_OF_ADDRTLV_MEM = Result::OutOfAddrTlvMem;
inline constexpr Result RFC5444_MTU_TOO_SMALL = Result::MtuTooSmall;
inline constexpr Result RFC5444_NO_MSGCREATOR = Result::NoMsgCreator;
inline constexpr Result RFC5444_FW_MESSAGE_TOO_LONG = Result::FwMessageTooLong;
inline constexpr Result RFC5444_FW_BAD_SIZE = Result::FwBadSize;
inline constexpr Result RFC5444_TOO_LARGE = Result::TooLarge;
inline constexpr Result RFC5444_FW_BAD_TRANSFORM = Result::FwBadTransform;
inline constexpr Result RFC5444_BAD_ADDR_HEAD_LENGTH = Result::BadAddrHeadLength;
inline constexpr Result RFC5444_BAD_ADDR_TAIL_LENGTH = Result::BadAddrTailLength;
inline constexpr Result RFC5444_BAD_TLV_INDEX = Result::BadTlvIndex;
inline constexpr Result RFC5444_RESULT_MIN = Result::ResultMin;

enum class CommonIana : int {
    ManetIpProto = 138,
    ManetUdpPort = 269,
};

inline constexpr const char *RFC5444_MANET_IPPROTO_TXT = "138";
inline constexpr const char *RFC5444_MANET_UDP_PORT_TXT = "269";
inline constexpr const char *RFC5444_MANET_MULTICAST_V4_TXT = "224.0.0.109";
inline constexpr const char *RFC5444_MANET_MULTICAST_V6_TXT = "ff02::6d";
inline constexpr const char *RFC7181_WILLINGNESS_DEFAULT_STRING = "7";

enum class MsgTypeIana : uint8_t {
    Hello = 0,
    Tc = 1,
};

enum class PacketTlvIana : uint8_t {
    Icv = 5,
    Timestamp = 6,
};

enum class MsgTlvIana5497 : uint8_t {
    IntervalTime = 0,
    ValidityTime = 1,
};

enum class MsgTlvIana7181 : uint8_t {
    MprWilling = 7,
    ContSeqNum = 8,
};

enum class MsgTlvIana7182 : uint8_t {
    Icv = 5,
    Timestamp = 6,
};

enum class MtIana7722 : uint8_t {
    MprTypes = 7,
    MprTypesExt = 1,
};

enum class SsrIanaDraft : uint8_t {
    Capability = 7,
    CapabilityExt = 2,
};

enum class WillingnessValue : int {
    Undefined = -1,
    Min = 0,
    Never = 0,
    Default = 7,
    Always = 15,
    Max = 15,
    Mask = 0x0f,
    Shift = 4,
};

enum class ContSeqNumExt : uint8_t {
    Complete = 0,
    Incomplete = 1,
    Bitmask = 1,
};

enum class AddrTlvIana5497 : uint8_t {
    IntervalTime = 0,
    ValidityTime = 1,
};

enum class AddrTlvIana6130 : uint8_t {
    LocalIf = 2,
    LinkStatus = 3,
    OtherNeighb = 4,
};

enum class AddrTlvIana7182 : uint8_t {
    Icv = 5,
    Timestamp = 6,
};

enum class AddrTlvIana7181 : uint8_t {
    LinkMetric = 7,
    Mpr = 8,
    NbrAddrType = 9,
    Gateway = 10,
};

enum class GatewayExt7181 : uint8_t {
    DstSpecGateway = 0,
    SrcDstSpecGateway = 1,
    SrcSpecDefGateway = 2,
};

enum class CustomSrcSpecGateway : uint8_t {
    SrcPrefix = 224,
};

enum class LocalIfValue6130 : uint8_t {
    ThisIf = 0,
    OtherIf = 1,
    Bitmask = 1,
};

enum class LinkStatusValue6130 : uint8_t {
    Lost = 0,
    Symmetric = 1,
    Heard = 2,
    Bitmask = 3,
};

enum class OtherNeighValue6130 : uint8_t {
    Symmetric = 1,
    Bitmask = 1,
};

enum class LinkMetricFlags7181 : uint8_t {
    IncomingLink = 1 << 7,
    OutgoingLink = 1 << 6,
    IncomingNeigh = 1 << 5,
    OutgoingNeigh = 1 << 4,
};

enum class MprBitmask7181 : uint8_t {
    Flooding = 1 << 0,
    Routing = 1 << 1,
};

enum class NbrAddrBitmask7181 : uint8_t {
    Originator = 1,
    Routable = 2,
};

enum class IcvExt7182 : uint8_t {
    Generic = 0,
    CryptHash = 1,
    SrcSpecCryptHash = 2,
};

enum class IcvHash7182 : int {
    Unknown = -1,
    Identity = 0,
    Sha1 = 1,
    Sha224 = 2,
    Sha256 = 3,
    Sha384 = 4,
    Sha512 = 5,
    Count,
};

enum class IcvCrypt7182 : int {
    Unknown = -1,
    Identity = 0,
    Rsa = 1,
    Dsa = 2,
    Hmac = 3,
    Des3 = 4,
    Aes = 5,
    Ecdsa = 6,
    Count,
};

enum class TimestampExt7182 : uint8_t {
    Monotonic = 0,
    Unix = 1,
    Ntp = 2,
    Random = 3,
};

const char *rfc7182_get_hash_name(IcvHash7182 hash);
const char *rfc7182_get_crypt_name(IcvCrypt7182 crypt);

const char **rfc7182_get_hashes();
const char **rfc7182_get_crypto();

IcvHash7182 rfc7182_get_hash_id(const char *name);
IcvCrypt7182 rfc7182_get_crypt_id(const char *name);

// Backward-compatible aliases used by existing code.
inline constexpr int RFC5444_MANET_IPPROTO = static_cast<int>(CommonIana::ManetIpProto);
inline constexpr int RFC5444_MANET_UDP_PORT = static_cast<int>(CommonIana::ManetUdpPort);

inline constexpr auto RFC6130_MSGTYPE_HELLO = MsgTypeIana::Hello;
inline constexpr auto RFC7181_MSGTYPE_TC = MsgTypeIana::Tc;

inline constexpr auto RFC7182_PKTTLV_ICV = PacketTlvIana::Icv;
inline constexpr auto RFC7182_PKTTLV_TIMESTAMP = PacketTlvIana::Timestamp;

inline constexpr auto RFC5497_MSGTLV_INTERVAL_TIME = MsgTlvIana5497::IntervalTime;
inline constexpr auto RFC5497_MSGTLV_VALIDITY_TIME = MsgTlvIana5497::ValidityTime;

inline constexpr auto RFC7181_MSGTLV_MPR_WILLING = MsgTlvIana7181::MprWilling;
inline constexpr auto RFC7181_MSGTLV_CONT_SEQ_NUM = MsgTlvIana7181::ContSeqNum;

inline constexpr auto RFC7182_MSGTLV_ICV = MsgTlvIana7182::Icv;
inline constexpr auto RFC7182_MSGTLV_TIMESTAMP = MsgTlvIana7182::Timestamp;

inline constexpr auto RFC7722_MSGTLV_MPR_TYPES = MtIana7722::MprTypes;
inline constexpr auto RFC7722_MSGTLV_MPR_TYPES_EXT = MtIana7722::MprTypesExt;

inline constexpr auto DRAFT_SSR_MSGTLV_CAPABILITY = SsrIanaDraft::Capability;
inline constexpr auto DRAFT_SSR_MSGTLV_CAPABILITY_EXT = SsrIanaDraft::CapabilityExt;

inline constexpr auto RFC7181_WILLINGNESS_UNDEFINED = WillingnessValue::Undefined;
inline constexpr auto RFC7181_WILLINGNESS_MIN = WillingnessValue::Min;
inline constexpr auto RFC7181_WILLINGNESS_NEVER = WillingnessValue::Never;
inline constexpr auto RFC7181_WILLINGNESS_DEFAULT = WillingnessValue::Default;
inline constexpr auto RFC7181_WILLINGNESS_ALWAYS = WillingnessValue::Always;
inline constexpr auto RFC7181_WILLINGNESS_MAX = WillingnessValue::Max;
inline constexpr auto RFC7181_WILLINGNESS_MASK = WillingnessValue::Mask;
inline constexpr auto RFC7181_WILLINGNESS_SHIFT = WillingnessValue::Shift;

inline constexpr auto RFC7181_CONT_SEQ_NUM_COMPLETE = ContSeqNumExt::Complete;
inline constexpr auto RFC7181_CONT_SEQ_NUM_INCOMPLETE = ContSeqNumExt::Incomplete;
inline constexpr auto RFC7181_CONT_SEQ_NUM_BITMASK = ContSeqNumExt::Bitmask;

inline constexpr auto RFC5497_ADDRTLV_INTERVAL_TIME = AddrTlvIana5497::IntervalTime;
inline constexpr auto RFC5497_ADDRTLV_VALIDITY_TIME = AddrTlvIana5497::ValidityTime;

inline constexpr auto RFC6130_ADDRTLV_LOCAL_IF = AddrTlvIana6130::LocalIf;
inline constexpr auto RFC6130_ADDRTLV_LINK_STATUS = AddrTlvIana6130::LinkStatus;
inline constexpr auto RFC6130_ADDRTLV_OTHER_NEIGHB = AddrTlvIana6130::OtherNeighb;

inline constexpr auto RFC7182_ADDRTLV_ICV = AddrTlvIana7182::Icv;
inline constexpr auto RFC7182_ADDRTLV_TIMESTAMP = AddrTlvIana7182::Timestamp;

inline constexpr auto RFC7181_ADDRTLV_LINK_METRIC = AddrTlvIana7181::LinkMetric;
inline constexpr auto RFC7181_ADDRTLV_MPR = AddrTlvIana7181::Mpr;
inline constexpr auto RFC7181_ADDRTLV_NBR_ADDR_TYPE = AddrTlvIana7181::NbrAddrType;
inline constexpr auto RFC7181_ADDRTLV_GATEWAY = AddrTlvIana7181::Gateway;

inline constexpr auto RFC7181_DSTSPEC_GATEWAY = GatewayExt7181::DstSpecGateway;
inline constexpr auto RFC7181_SRCSPEC_GATEWAY = GatewayExt7181::SrcDstSpecGateway;
inline constexpr auto RFC7181_SRCSPEC_DEF_GATEWAY = GatewayExt7181::SrcSpecDefGateway;

inline constexpr auto SRCSPEC_GW_ADDRTLV_SRC_PREFIX = CustomSrcSpecGateway::SrcPrefix;

inline constexpr auto RFC6130_LOCALIF_THIS_IF = LocalIfValue6130::ThisIf;
inline constexpr auto RFC6130_LOCALIF_OTHER_IF = LocalIfValue6130::OtherIf;
inline constexpr auto RFC6130_LOCALIF_BITMASK = LocalIfValue6130::Bitmask;

inline constexpr auto RFC6130_LINKSTATUS_LOST = LinkStatusValue6130::Lost;
inline constexpr auto RFC6130_LINKSTATUS_SYMMETRIC = LinkStatusValue6130::Symmetric;
inline constexpr auto RFC6130_LINKSTATUS_HEARD = LinkStatusValue6130::Heard;
inline constexpr auto RFC6130_LINKSTATUS_BITMASK = LinkStatusValue6130::Bitmask;

inline constexpr auto RFC6130_OTHERNEIGHB_SYMMETRIC = OtherNeighValue6130::Symmetric;
inline constexpr auto RFC6130_OTHERNEIGHB_BITMASK = OtherNeighValue6130::Bitmask;

inline constexpr auto RFC7181_LINKMETRIC_INCOMING_LINK = LinkMetricFlags7181::IncomingLink;
inline constexpr auto RFC7181_LINKMETRIC_OUTGOING_LINK = LinkMetricFlags7181::OutgoingLink;
inline constexpr auto RFC7181_LINKMETRIC_INCOMING_NEIGH = LinkMetricFlags7181::IncomingNeigh;
inline constexpr auto RFC7181_LINKMETRIC_OUTGOING_NEIGH = LinkMetricFlags7181::OutgoingNeigh;

inline constexpr auto RFC7181_MPR_FLOODING = MprBitmask7181::Flooding;
inline constexpr auto RFC7181_MPR_ROUTING = MprBitmask7181::Routing;

inline constexpr auto RFC7181_NBR_ADDR_TYPE_ORIGINATOR = NbrAddrBitmask7181::Originator;
inline constexpr auto RFC7181_NBR_ADDR_TYPE_ROUTABLE = NbrAddrBitmask7181::Routable;

inline constexpr auto RFC7182_ICV_EXT_GENERIC = IcvExt7182::Generic;
inline constexpr auto RFC7182_ICV_EXT_CRYPTHASH = IcvExt7182::CryptHash;
inline constexpr auto RFC7182_ICV_EXT_SRCSPEC_CRYPTHASH = IcvExt7182::SrcSpecCryptHash;

inline constexpr auto RFC7182_ICV_HASH_UNKNOWN = IcvHash7182::Unknown;
inline constexpr auto RFC7182_ICV_HASH_IDENTITY = IcvHash7182::Identity;
inline constexpr auto RFC7182_ICV_HASH_SHA_1 = IcvHash7182::Sha1;
inline constexpr auto RFC7182_ICV_HASH_SHA_224 = IcvHash7182::Sha224;
inline constexpr auto RFC7182_ICV_HASH_SHA_256 = IcvHash7182::Sha256;
inline constexpr auto RFC7182_ICV_HASH_SHA_384 = IcvHash7182::Sha384;
inline constexpr auto RFC7182_ICV_HASH_SHA_512 = IcvHash7182::Sha512;
inline constexpr auto RFC7182_ICV_HASH_COUNT = IcvHash7182::Count;

inline constexpr auto RFC7182_ICV_CRYPT_UNKNOWN = IcvCrypt7182::Unknown;
inline constexpr auto RFC7182_ICV_CRYPT_IDENTITY = IcvCrypt7182::Identity;
inline constexpr auto RFC7182_ICV_CRYPT_RSA = IcvCrypt7182::Rsa;
inline constexpr auto RFC7182_ICV_CRYPT_DSA = IcvCrypt7182::Dsa;
inline constexpr auto RFC7182_ICV_CRYPT_HMAC = IcvCrypt7182::Hmac;
inline constexpr auto RFC7182_ICV_CRYPT_3DES = IcvCrypt7182::Des3;
inline constexpr auto RFC7182_ICV_CRYPT_AES = IcvCrypt7182::Aes;
inline constexpr auto RFC7182_ICV_CRYPT_ECDSA = IcvCrypt7182::Ecdsa;
inline constexpr auto RFC7182_ICV_CRYPT_COUNT = IcvCrypt7182::Count;

inline constexpr auto RFC7182_TIMESTAMP_EXT_MONOTONIC = TimestampExt7182::Monotonic;
inline constexpr auto RFC7182_TIMESTAMP_EXT_UNIX = TimestampExt7182::Unix;
inline constexpr auto RFC7182_TIMESTAMP_EXT_NTP = TimestampExt7182::Ntp;
inline constexpr auto RFC7182_TIMESTAMP_EXT_RANDOM = TimestampExt7182::Random;

struct rfc7181_metric_field {
    uint8_t b[2] = {0, 0};
};

struct TlvValue {
    uint8_t type = 0;
    uint8_t typeExt = 0;
    uint8_t indexStart = 0;
    uint8_t indexEnd = 0;
    std::vector<uint8_t> value;
};

struct DecodedMessage {
    uint8_t type = 0;
    uint8_t flags = 0;
    uint8_t addrLen = 0;
    uint8_t hopCount = 0;
    uint8_t hopLimit = 0;
    uint16_t seqNo = 0;
    std::vector<uint8_t> originator;
    std::vector<TlvValue> messageTlvs;
};

uint8_t rfc5497_timetlv_get_from_vector(const uint8_t *vector, size_t vector_length, uint8_t hopcount);
uint8_t rfc5497_timetlv_encode(uint64_t decoded);
uint64_t rfc5497_timetlv_decode(uint8_t encoded);

int rfc7181_metric_encode(rfc7181_metric_field *encoded, uint32_t decoded);
uint32_t rfc7181_metric_decode(rfc7181_metric_field *encoded);

int rfc5444_seqno_difference(uint16_t seqno1, uint16_t seqno2);

} // namespace mysrc::olsrv2::rfc5444
