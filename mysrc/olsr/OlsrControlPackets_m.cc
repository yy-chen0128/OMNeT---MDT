//
// Generated file, do not edit! Created by opp_msgtool 6.2 from mysrc/OLSR/OlsrControlPackets.msg.
//

// Disable warnings about unused variables, empty switch stmts, etc:
#ifdef _MSC_VER
#  pragma warning(disable:4101)
#  pragma warning(disable:4065)
#endif

#if defined(__clang__)
#  pragma clang diagnostic ignored "-Wshadow"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wunused-parameter"
#  pragma clang diagnostic ignored "-Wc++98-compat"
#  pragma clang diagnostic ignored "-Wunreachable-code-break"
#  pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wshadow"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#  pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#include <iostream>
#include <sstream>
#include <memory>
#include <type_traits>
#include "../olsr/OlsrControlPackets_m.h"

namespace omnetpp {

// Template pack/unpack rules. They are declared *after* a1l type-specific pack functions for multiple reasons.
// They are in the omnetpp namespace, to allow them to be found by argument-dependent lookup via the cCommBuffer argument

// Packing/unpacking an std::vector
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::vector<T,A>& v)
{
    int n = v.size();
    doParsimPacking(buffer, n);
    for (int i = 0; i < n; i++)
        doParsimPacking(buffer, v[i]);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::vector<T,A>& v)
{
    int n;
    doParsimUnpacking(buffer, n);
    v.resize(n);
    for (int i = 0; i < n; i++)
        doParsimUnpacking(buffer, v[i]);
}

// Packing/unpacking an std::list
template<typename T, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::list<T,A>& l)
{
    doParsimPacking(buffer, (int)l.size());
    for (typename std::list<T,A>::const_iterator it = l.begin(); it != l.end(); ++it)
        doParsimPacking(buffer, (T&)*it);
}

template<typename T, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::list<T,A>& l)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        l.push_back(T());
        doParsimUnpacking(buffer, l.back());
    }
}

// Packing/unpacking an std::set
template<typename T, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::set<T,Tr,A>& s)
{
    doParsimPacking(buffer, (int)s.size());
    for (typename std::set<T,Tr,A>::const_iterator it = s.begin(); it != s.end(); ++it)
        doParsimPacking(buffer, *it);
}

template<typename T, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::set<T,Tr,A>& s)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        T x;
        doParsimUnpacking(buffer, x);
        s.insert(x);
    }
}

// Packing/unpacking an std::map
template<typename K, typename V, typename Tr, typename A>
void doParsimPacking(omnetpp::cCommBuffer *buffer, const std::map<K,V,Tr,A>& m)
{
    doParsimPacking(buffer, (int)m.size());
    for (typename std::map<K,V,Tr,A>::const_iterator it = m.begin(); it != m.end(); ++it) {
        doParsimPacking(buffer, it->first);
        doParsimPacking(buffer, it->second);
    }
}

template<typename K, typename V, typename Tr, typename A>
void doParsimUnpacking(omnetpp::cCommBuffer *buffer, std::map<K,V,Tr,A>& m)
{
    int n;
    doParsimUnpacking(buffer, n);
    for (int i = 0; i < n; i++) {
        K k; V v;
        doParsimUnpacking(buffer, k);
        doParsimUnpacking(buffer, v);
        m[k] = v;
    }
}

// Default pack/unpack function for arrays
template<typename T>
void doParsimArrayPacking(omnetpp::cCommBuffer *b, const T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimPacking(b, t[i]);
}

template<typename T>
void doParsimArrayUnpacking(omnetpp::cCommBuffer *b, T *t, int n)
{
    for (int i = 0; i < n; i++)
        doParsimUnpacking(b, t[i]);
}

// Default rule to prevent compiler from choosing base class' doParsimPacking() function
template<typename T>
void doParsimPacking(omnetpp::cCommBuffer *, const T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: No doParsimPacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

template<typename T>
void doParsimUnpacking(omnetpp::cCommBuffer *, T& t)
{
    throw omnetpp::cRuntimeError("Parsim error: No doParsimUnpacking() function for type %s", omnetpp::opp_typename(typeid(t)));
}

}  // namespace omnetpp

namespace inet {

Register_Enum(inet::OlsrPktType, (inet::OlsrPktType::Hello, inet::OlsrPktType::Tc));

Register_Class(OlsrControlPacket)

OlsrControlPacket::OlsrControlPacket() : ::inet::FieldsChunk()
{
}

OlsrControlPacket::OlsrControlPacket(const OlsrControlPacket& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

OlsrControlPacket::~OlsrControlPacket()
{
}

OlsrControlPacket& OlsrControlPacket::operator=(const OlsrControlPacket& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void OlsrControlPacket::copy(const OlsrControlPacket& other)
{
    this->packetType = other.packetType;
}

void OlsrControlPacket::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->packetType);
}

void OlsrControlPacket::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->packetType);
}

OlsrPktType OlsrControlPacket::getPacketType() const
{
    return this->packetType;
}

void OlsrControlPacket::setPacketType(OlsrPktType packetType)
{
    handleChange();
    this->packetType = packetType;
}

class OlsrControlPacketDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_packetType,
    };
  public:
    OlsrControlPacketDescriptor();
    virtual ~OlsrControlPacketDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(OlsrControlPacketDescriptor)

OlsrControlPacketDescriptor::OlsrControlPacketDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::OlsrControlPacket)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

OlsrControlPacketDescriptor::~OlsrControlPacketDescriptor()
{
    delete[] propertyNames;
}

bool OlsrControlPacketDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<OlsrControlPacket *>(obj)!=nullptr;
}

const char **OlsrControlPacketDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *OlsrControlPacketDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int OlsrControlPacketDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int OlsrControlPacketDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_packetType
    };
    return (field >= 0 && field < 1) ? fieldTypeFlags[field] : 0;
}

const char *OlsrControlPacketDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "packetType",
    };
    return (field >= 0 && field < 1) ? fieldNames[field] : nullptr;
}

int OlsrControlPacketDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "packetType") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *OlsrControlPacketDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::OlsrPktType",    // FIELD_packetType
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **OlsrControlPacketDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_packetType: {
            static const char *names[] = { "enum",  nullptr };
            return names;
        }
        default: return nullptr;
    }
}

const char *OlsrControlPacketDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_packetType:
            if (!strcmp(propertyName, "enum")) return "inet::OlsrPktType";
            return nullptr;
        default: return nullptr;
    }
}

int OlsrControlPacketDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void OlsrControlPacketDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'OlsrControlPacket'", field);
    }
}

const char *OlsrControlPacketDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string OlsrControlPacketDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        case FIELD_packetType: return enum2string(pp->getPacketType(), "inet::OlsrPktType");
        default: return "";
    }
}

void OlsrControlPacketDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrControlPacket'", field);
    }
}

omnetpp::cValue OlsrControlPacketDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        case FIELD_packetType: return static_cast<int>(pp->getPacketType());
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'OlsrControlPacket' as cValue -- field index out of range?", field);
    }
}

void OlsrControlPacketDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrControlPacket'", field);
    }
}

const char *OlsrControlPacketDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr OlsrControlPacketDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void OlsrControlPacketDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrControlPacket *pp = omnetpp::fromAnyPtr<OlsrControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrControlPacket'", field);
    }
}

Register_Class(OlsrHello)

OlsrHello::OlsrHello() : ::inet::OlsrControlPacket()
{
}

OlsrHello::OlsrHello(const OlsrHello& other) : ::inet::OlsrControlPacket(other)
{
    copy(other);
}

OlsrHello::~OlsrHello()
{
    delete [] this->linkCodes;
    delete [] this->linkCounts;
    delete [] this->nbIfaceAddrs;
}

OlsrHello& OlsrHello::operator=(const OlsrHello& other)
{
    if (this == &other) return *this;
    ::inet::OlsrControlPacket::operator=(other);
    copy(other);
    return *this;
}

void OlsrHello::copy(const OlsrHello& other)
{
    this->originator = other.originator;
    this->hTime = other.hTime;
    this->willingness = other.willingness;
    this->msgSeq = other.msgSeq;
    delete [] this->linkCodes;
    this->linkCodes = (other.linkCodes_arraysize==0) ? nullptr : new uint8_t[other.linkCodes_arraysize];
    linkCodes_arraysize = other.linkCodes_arraysize;
    for (size_t i = 0; i < linkCodes_arraysize; i++) {
        this->linkCodes[i] = other.linkCodes[i];
    }
    delete [] this->linkCounts;
    this->linkCounts = (other.linkCounts_arraysize==0) ? nullptr : new int[other.linkCounts_arraysize];
    linkCounts_arraysize = other.linkCounts_arraysize;
    for (size_t i = 0; i < linkCounts_arraysize; i++) {
        this->linkCounts[i] = other.linkCounts[i];
    }
    delete [] this->nbIfaceAddrs;
    this->nbIfaceAddrs = (other.nbIfaceAddrs_arraysize==0) ? nullptr : new L3Address[other.nbIfaceAddrs_arraysize];
    nbIfaceAddrs_arraysize = other.nbIfaceAddrs_arraysize;
    for (size_t i = 0; i < nbIfaceAddrs_arraysize; i++) {
        this->nbIfaceAddrs[i] = other.nbIfaceAddrs[i];
    }
}

void OlsrHello::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::OlsrControlPacket::parsimPack(b);
    doParsimPacking(b,this->originator);
    doParsimPacking(b,this->hTime);
    doParsimPacking(b,this->willingness);
    doParsimPacking(b,this->msgSeq);
    b->pack(linkCodes_arraysize);
    doParsimArrayPacking(b,this->linkCodes,linkCodes_arraysize);
    b->pack(linkCounts_arraysize);
    doParsimArrayPacking(b,this->linkCounts,linkCounts_arraysize);
    b->pack(nbIfaceAddrs_arraysize);
    doParsimArrayPacking(b,this->nbIfaceAddrs,nbIfaceAddrs_arraysize);
}

void OlsrHello::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::OlsrControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->originator);
    doParsimUnpacking(b,this->hTime);
    doParsimUnpacking(b,this->willingness);
    doParsimUnpacking(b,this->msgSeq);
    delete [] this->linkCodes;
    b->unpack(linkCodes_arraysize);
    if (linkCodes_arraysize == 0) {
        this->linkCodes = nullptr;
    } else {
        this->linkCodes = new uint8_t[linkCodes_arraysize];
        doParsimArrayUnpacking(b,this->linkCodes,linkCodes_arraysize);
    }
    delete [] this->linkCounts;
    b->unpack(linkCounts_arraysize);
    if (linkCounts_arraysize == 0) {
        this->linkCounts = nullptr;
    } else {
        this->linkCounts = new int[linkCounts_arraysize];
        doParsimArrayUnpacking(b,this->linkCounts,linkCounts_arraysize);
    }
    delete [] this->nbIfaceAddrs;
    b->unpack(nbIfaceAddrs_arraysize);
    if (nbIfaceAddrs_arraysize == 0) {
        this->nbIfaceAddrs = nullptr;
    } else {
        this->nbIfaceAddrs = new L3Address[nbIfaceAddrs_arraysize];
        doParsimArrayUnpacking(b,this->nbIfaceAddrs,nbIfaceAddrs_arraysize);
    }
}

const L3Address& OlsrHello::getOriginator() const
{
    return this->originator;
}

void OlsrHello::setOriginator(const L3Address& originator)
{
    handleChange();
    this->originator = originator;
}

::omnetpp::simtime_t OlsrHello::getHTime() const
{
    return this->hTime;
}

void OlsrHello::setHTime(::omnetpp::simtime_t hTime)
{
    handleChange();
    this->hTime = hTime;
}

uint8_t OlsrHello::getWillingness() const
{
    return this->willingness;
}

void OlsrHello::setWillingness(uint8_t willingness)
{
    handleChange();
    this->willingness = willingness;
}

uint16_t OlsrHello::getMsgSeq() const
{
    return this->msgSeq;
}

void OlsrHello::setMsgSeq(uint16_t msgSeq)
{
    handleChange();
    this->msgSeq = msgSeq;
}

size_t OlsrHello::getLinkCodesArraySize() const
{
    return linkCodes_arraysize;
}

uint8_t OlsrHello::getLinkCodes(size_t k) const
{
    if (k >= linkCodes_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCodes_arraysize, (unsigned long)k);
    return this->linkCodes[k];
}

void OlsrHello::setLinkCodesArraySize(size_t newSize)
{
    handleChange();
    uint8_t *linkCodes2 = (newSize==0) ? nullptr : new uint8_t[newSize];
    size_t minSize = linkCodes_arraysize < newSize ? linkCodes_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        linkCodes2[i] = this->linkCodes[i];
    for (size_t i = minSize; i < newSize; i++)
        linkCodes2[i] = 0;
    delete [] this->linkCodes;
    this->linkCodes = linkCodes2;
    linkCodes_arraysize = newSize;
}

void OlsrHello::setLinkCodes(size_t k, uint8_t linkCodes)
{
    if (k >= linkCodes_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCodes_arraysize, (unsigned long)k);
    handleChange();
    this->linkCodes[k] = linkCodes;
}

void OlsrHello::insertLinkCodes(size_t k, uint8_t linkCodes)
{
    if (k > linkCodes_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCodes_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = linkCodes_arraysize + 1;
    uint8_t *linkCodes2 = new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        linkCodes2[i] = this->linkCodes[i];
    linkCodes2[k] = linkCodes;
    for (i = k + 1; i < newSize; i++)
        linkCodes2[i] = this->linkCodes[i-1];
    delete [] this->linkCodes;
    this->linkCodes = linkCodes2;
    linkCodes_arraysize = newSize;
}

void OlsrHello::appendLinkCodes(uint8_t linkCodes)
{
    insertLinkCodes(linkCodes_arraysize, linkCodes);
}

void OlsrHello::eraseLinkCodes(size_t k)
{
    if (k >= linkCodes_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCodes_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = linkCodes_arraysize - 1;
    uint8_t *linkCodes2 = (newSize == 0) ? nullptr : new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        linkCodes2[i] = this->linkCodes[i];
    for (i = k; i < newSize; i++)
        linkCodes2[i] = this->linkCodes[i+1];
    delete [] this->linkCodes;
    this->linkCodes = linkCodes2;
    linkCodes_arraysize = newSize;
}

size_t OlsrHello::getLinkCountsArraySize() const
{
    return linkCounts_arraysize;
}

int OlsrHello::getLinkCounts(size_t k) const
{
    if (k >= linkCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCounts_arraysize, (unsigned long)k);
    return this->linkCounts[k];
}

void OlsrHello::setLinkCountsArraySize(size_t newSize)
{
    handleChange();
    int *linkCounts2 = (newSize==0) ? nullptr : new int[newSize];
    size_t minSize = linkCounts_arraysize < newSize ? linkCounts_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        linkCounts2[i] = this->linkCounts[i];
    for (size_t i = minSize; i < newSize; i++)
        linkCounts2[i] = 0;
    delete [] this->linkCounts;
    this->linkCounts = linkCounts2;
    linkCounts_arraysize = newSize;
}

void OlsrHello::setLinkCounts(size_t k, int linkCounts)
{
    if (k >= linkCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCounts_arraysize, (unsigned long)k);
    handleChange();
    this->linkCounts[k] = linkCounts;
}

void OlsrHello::insertLinkCounts(size_t k, int linkCounts)
{
    if (k > linkCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = linkCounts_arraysize + 1;
    int *linkCounts2 = new int[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        linkCounts2[i] = this->linkCounts[i];
    linkCounts2[k] = linkCounts;
    for (i = k + 1; i < newSize; i++)
        linkCounts2[i] = this->linkCounts[i-1];
    delete [] this->linkCounts;
    this->linkCounts = linkCounts2;
    linkCounts_arraysize = newSize;
}

void OlsrHello::appendLinkCounts(int linkCounts)
{
    insertLinkCounts(linkCounts_arraysize, linkCounts);
}

void OlsrHello::eraseLinkCounts(size_t k)
{
    if (k >= linkCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = linkCounts_arraysize - 1;
    int *linkCounts2 = (newSize == 0) ? nullptr : new int[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        linkCounts2[i] = this->linkCounts[i];
    for (i = k; i < newSize; i++)
        linkCounts2[i] = this->linkCounts[i+1];
    delete [] this->linkCounts;
    this->linkCounts = linkCounts2;
    linkCounts_arraysize = newSize;
}

size_t OlsrHello::getNbIfaceAddrsArraySize() const
{
    return nbIfaceAddrs_arraysize;
}

const L3Address& OlsrHello::getNbIfaceAddrs(size_t k) const
{
    if (k >= nbIfaceAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)nbIfaceAddrs_arraysize, (unsigned long)k);
    return this->nbIfaceAddrs[k];
}

void OlsrHello::setNbIfaceAddrsArraySize(size_t newSize)
{
    handleChange();
    L3Address *nbIfaceAddrs2 = (newSize==0) ? nullptr : new L3Address[newSize];
    size_t minSize = nbIfaceAddrs_arraysize < newSize ? nbIfaceAddrs_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        nbIfaceAddrs2[i] = this->nbIfaceAddrs[i];
    delete [] this->nbIfaceAddrs;
    this->nbIfaceAddrs = nbIfaceAddrs2;
    nbIfaceAddrs_arraysize = newSize;
}

void OlsrHello::setNbIfaceAddrs(size_t k, const L3Address& nbIfaceAddrs)
{
    if (k >= nbIfaceAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)nbIfaceAddrs_arraysize, (unsigned long)k);
    handleChange();
    this->nbIfaceAddrs[k] = nbIfaceAddrs;
}

void OlsrHello::insertNbIfaceAddrs(size_t k, const L3Address& nbIfaceAddrs)
{
    if (k > nbIfaceAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)nbIfaceAddrs_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = nbIfaceAddrs_arraysize + 1;
    L3Address *nbIfaceAddrs2 = new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        nbIfaceAddrs2[i] = this->nbIfaceAddrs[i];
    nbIfaceAddrs2[k] = nbIfaceAddrs;
    for (i = k + 1; i < newSize; i++)
        nbIfaceAddrs2[i] = this->nbIfaceAddrs[i-1];
    delete [] this->nbIfaceAddrs;
    this->nbIfaceAddrs = nbIfaceAddrs2;
    nbIfaceAddrs_arraysize = newSize;
}

void OlsrHello::appendNbIfaceAddrs(const L3Address& nbIfaceAddrs)
{
    insertNbIfaceAddrs(nbIfaceAddrs_arraysize, nbIfaceAddrs);
}

void OlsrHello::eraseNbIfaceAddrs(size_t k)
{
    if (k >= nbIfaceAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)nbIfaceAddrs_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = nbIfaceAddrs_arraysize - 1;
    L3Address *nbIfaceAddrs2 = (newSize == 0) ? nullptr : new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        nbIfaceAddrs2[i] = this->nbIfaceAddrs[i];
    for (i = k; i < newSize; i++)
        nbIfaceAddrs2[i] = this->nbIfaceAddrs[i+1];
    delete [] this->nbIfaceAddrs;
    this->nbIfaceAddrs = nbIfaceAddrs2;
    nbIfaceAddrs_arraysize = newSize;
}

class OlsrHelloDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_originator,
        FIELD_hTime,
        FIELD_willingness,
        FIELD_msgSeq,
        FIELD_linkCodes,
        FIELD_linkCounts,
        FIELD_nbIfaceAddrs,
    };
  public:
    OlsrHelloDescriptor();
    virtual ~OlsrHelloDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(OlsrHelloDescriptor)

OlsrHelloDescriptor::OlsrHelloDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::OlsrHello)), "inet::OlsrControlPacket")
{
    propertyNames = nullptr;
}

OlsrHelloDescriptor::~OlsrHelloDescriptor()
{
    delete[] propertyNames;
}

bool OlsrHelloDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<OlsrHello *>(obj)!=nullptr;
}

const char **OlsrHelloDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *OlsrHelloDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int OlsrHelloDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 7+base->getFieldCount() : 7;
}

unsigned int OlsrHelloDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_originator
        FD_ISEDITABLE,    // FIELD_hTime
        FD_ISEDITABLE,    // FIELD_willingness
        FD_ISEDITABLE,    // FIELD_msgSeq
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_linkCodes
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_linkCounts
        FD_ISARRAY | FD_ISRESIZABLE,    // FIELD_nbIfaceAddrs
    };
    return (field >= 0 && field < 7) ? fieldTypeFlags[field] : 0;
}

const char *OlsrHelloDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "originator",
        "hTime",
        "willingness",
        "msgSeq",
        "linkCodes",
        "linkCounts",
        "nbIfaceAddrs",
    };
    return (field >= 0 && field < 7) ? fieldNames[field] : nullptr;
}

int OlsrHelloDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "originator") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "hTime") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "willingness") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "msgSeq") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "linkCodes") == 0) return baseIndex + 4;
    if (strcmp(fieldName, "linkCounts") == 0) return baseIndex + 5;
    if (strcmp(fieldName, "nbIfaceAddrs") == 0) return baseIndex + 6;
    return base ? base->findField(fieldName) : -1;
}

const char *OlsrHelloDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_originator
        "omnetpp::simtime_t",    // FIELD_hTime
        "uint8",    // FIELD_willingness
        "uint16",    // FIELD_msgSeq
        "uint8",    // FIELD_linkCodes
        "int",    // FIELD_linkCounts
        "inet::L3Address",    // FIELD_nbIfaceAddrs
    };
    return (field >= 0 && field < 7) ? fieldTypeStrings[field] : nullptr;
}

const char **OlsrHelloDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *OlsrHelloDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int OlsrHelloDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_linkCodes: return pp->getLinkCodesArraySize();
        case FIELD_linkCounts: return pp->getLinkCountsArraySize();
        case FIELD_nbIfaceAddrs: return pp->getNbIfaceAddrsArraySize();
        default: return 0;
    }
}

void OlsrHelloDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_linkCodes: pp->setLinkCodesArraySize(size); break;
        case FIELD_linkCounts: pp->setLinkCountsArraySize(size); break;
        case FIELD_nbIfaceAddrs: pp->setNbIfaceAddrsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'OlsrHello'", field);
    }
}

const char *OlsrHelloDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string OlsrHelloDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_originator: return pp->getOriginator().str();
        case FIELD_hTime: return simtime2string(pp->getHTime());
        case FIELD_willingness: return ulong2string(pp->getWillingness());
        case FIELD_msgSeq: return ulong2string(pp->getMsgSeq());
        case FIELD_linkCodes: return ulong2string(pp->getLinkCodes(i));
        case FIELD_linkCounts: return long2string(pp->getLinkCounts(i));
        case FIELD_nbIfaceAddrs: return pp->getNbIfaceAddrs(i).str();
        default: return "";
    }
}

void OlsrHelloDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_hTime: pp->setHTime(string2simtime(value)); break;
        case FIELD_willingness: pp->setWillingness(string2ulong(value)); break;
        case FIELD_msgSeq: pp->setMsgSeq(string2ulong(value)); break;
        case FIELD_linkCodes: pp->setLinkCodes(i,string2ulong(value)); break;
        case FIELD_linkCounts: pp->setLinkCounts(i,string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrHello'", field);
    }
}

omnetpp::cValue OlsrHelloDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_originator: return omnetpp::toAnyPtr(&pp->getOriginator()); break;
        case FIELD_hTime: return pp->getHTime().dbl();
        case FIELD_willingness: return (omnetpp::intval_t)(pp->getWillingness());
        case FIELD_msgSeq: return (omnetpp::intval_t)(pp->getMsgSeq());
        case FIELD_linkCodes: return (omnetpp::intval_t)(pp->getLinkCodes(i));
        case FIELD_linkCounts: return pp->getLinkCounts(i);
        case FIELD_nbIfaceAddrs: return omnetpp::toAnyPtr(&pp->getNbIfaceAddrs(i)); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'OlsrHello' as cValue -- field index out of range?", field);
    }
}

void OlsrHelloDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_hTime: pp->setHTime(value.doubleValue()); break;
        case FIELD_willingness: pp->setWillingness(omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_msgSeq: pp->setMsgSeq(omnetpp::checked_int_cast<uint16_t>(value.intValue())); break;
        case FIELD_linkCodes: pp->setLinkCodes(i,omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_linkCounts: pp->setLinkCounts(i,omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrHello'", field);
    }
}

const char *OlsrHelloDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr OlsrHelloDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        case FIELD_originator: return omnetpp::toAnyPtr(&pp->getOriginator()); break;
        case FIELD_nbIfaceAddrs: return omnetpp::toAnyPtr(&pp->getNbIfaceAddrs(i)); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void OlsrHelloDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrHello *pp = omnetpp::fromAnyPtr<OlsrHello>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrHello'", field);
    }
}

Register_Class(OlsrTcGroup)

OlsrTcGroup::OlsrTcGroup() : ::inet::OlsrControlPacket()
{
}

OlsrTcGroup::OlsrTcGroup(const OlsrTcGroup& other) : ::inet::OlsrControlPacket(other)
{
    copy(other);
}

OlsrTcGroup::~OlsrTcGroup()
{
    delete [] this->originators;
    delete [] this->ansns;
    delete [] this->seqNums;
    delete [] this->hopCounts;
    delete [] this->advCounts;
    delete [] this->advertisedNeighbors;
}

OlsrTcGroup& OlsrTcGroup::operator=(const OlsrTcGroup& other)
{
    if (this == &other) return *this;
    ::inet::OlsrControlPacket::operator=(other);
    copy(other);
    return *this;
}

void OlsrTcGroup::copy(const OlsrTcGroup& other)
{
    this->genTime = other.genTime;
    delete [] this->originators;
    this->originators = (other.originators_arraysize==0) ? nullptr : new L3Address[other.originators_arraysize];
    originators_arraysize = other.originators_arraysize;
    for (size_t i = 0; i < originators_arraysize; i++) {
        this->originators[i] = other.originators[i];
    }
    delete [] this->ansns;
    this->ansns = (other.ansns_arraysize==0) ? nullptr : new uint16_t[other.ansns_arraysize];
    ansns_arraysize = other.ansns_arraysize;
    for (size_t i = 0; i < ansns_arraysize; i++) {
        this->ansns[i] = other.ansns[i];
    }
    delete [] this->seqNums;
    this->seqNums = (other.seqNums_arraysize==0) ? nullptr : new uint16_t[other.seqNums_arraysize];
    seqNums_arraysize = other.seqNums_arraysize;
    for (size_t i = 0; i < seqNums_arraysize; i++) {
        this->seqNums[i] = other.seqNums[i];
    }
    delete [] this->hopCounts;
    this->hopCounts = (other.hopCounts_arraysize==0) ? nullptr : new uint8_t[other.hopCounts_arraysize];
    hopCounts_arraysize = other.hopCounts_arraysize;
    for (size_t i = 0; i < hopCounts_arraysize; i++) {
        this->hopCounts[i] = other.hopCounts[i];
    }
    delete [] this->advCounts;
    this->advCounts = (other.advCounts_arraysize==0) ? nullptr : new int[other.advCounts_arraysize];
    advCounts_arraysize = other.advCounts_arraysize;
    for (size_t i = 0; i < advCounts_arraysize; i++) {
        this->advCounts[i] = other.advCounts[i];
    }
    delete [] this->advertisedNeighbors;
    this->advertisedNeighbors = (other.advertisedNeighbors_arraysize==0) ? nullptr : new L3Address[other.advertisedNeighbors_arraysize];
    advertisedNeighbors_arraysize = other.advertisedNeighbors_arraysize;
    for (size_t i = 0; i < advertisedNeighbors_arraysize; i++) {
        this->advertisedNeighbors[i] = other.advertisedNeighbors[i];
    }
}

void OlsrTcGroup::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::OlsrControlPacket::parsimPack(b);
    doParsimPacking(b,this->genTime);
    b->pack(originators_arraysize);
    doParsimArrayPacking(b,this->originators,originators_arraysize);
    b->pack(ansns_arraysize);
    doParsimArrayPacking(b,this->ansns,ansns_arraysize);
    b->pack(seqNums_arraysize);
    doParsimArrayPacking(b,this->seqNums,seqNums_arraysize);
    b->pack(hopCounts_arraysize);
    doParsimArrayPacking(b,this->hopCounts,hopCounts_arraysize);
    b->pack(advCounts_arraysize);
    doParsimArrayPacking(b,this->advCounts,advCounts_arraysize);
    b->pack(advertisedNeighbors_arraysize);
    doParsimArrayPacking(b,this->advertisedNeighbors,advertisedNeighbors_arraysize);
}

void OlsrTcGroup::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::OlsrControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->genTime);
    delete [] this->originators;
    b->unpack(originators_arraysize);
    if (originators_arraysize == 0) {
        this->originators = nullptr;
    } else {
        this->originators = new L3Address[originators_arraysize];
        doParsimArrayUnpacking(b,this->originators,originators_arraysize);
    }
    delete [] this->ansns;
    b->unpack(ansns_arraysize);
    if (ansns_arraysize == 0) {
        this->ansns = nullptr;
    } else {
        this->ansns = new uint16_t[ansns_arraysize];
        doParsimArrayUnpacking(b,this->ansns,ansns_arraysize);
    }
    delete [] this->seqNums;
    b->unpack(seqNums_arraysize);
    if (seqNums_arraysize == 0) {
        this->seqNums = nullptr;
    } else {
        this->seqNums = new uint16_t[seqNums_arraysize];
        doParsimArrayUnpacking(b,this->seqNums,seqNums_arraysize);
    }
    delete [] this->hopCounts;
    b->unpack(hopCounts_arraysize);
    if (hopCounts_arraysize == 0) {
        this->hopCounts = nullptr;
    } else {
        this->hopCounts = new uint8_t[hopCounts_arraysize];
        doParsimArrayUnpacking(b,this->hopCounts,hopCounts_arraysize);
    }
    delete [] this->advCounts;
    b->unpack(advCounts_arraysize);
    if (advCounts_arraysize == 0) {
        this->advCounts = nullptr;
    } else {
        this->advCounts = new int[advCounts_arraysize];
        doParsimArrayUnpacking(b,this->advCounts,advCounts_arraysize);
    }
    delete [] this->advertisedNeighbors;
    b->unpack(advertisedNeighbors_arraysize);
    if (advertisedNeighbors_arraysize == 0) {
        this->advertisedNeighbors = nullptr;
    } else {
        this->advertisedNeighbors = new L3Address[advertisedNeighbors_arraysize];
        doParsimArrayUnpacking(b,this->advertisedNeighbors,advertisedNeighbors_arraysize);
    }
}

::omnetpp::simtime_t OlsrTcGroup::getGenTime() const
{
    return this->genTime;
}

void OlsrTcGroup::setGenTime(::omnetpp::simtime_t genTime)
{
    handleChange();
    this->genTime = genTime;
}

size_t OlsrTcGroup::getOriginatorsArraySize() const
{
    return originators_arraysize;
}

const L3Address& OlsrTcGroup::getOriginators(size_t k) const
{
    if (k >= originators_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)originators_arraysize, (unsigned long)k);
    return this->originators[k];
}

void OlsrTcGroup::setOriginatorsArraySize(size_t newSize)
{
    handleChange();
    L3Address *originators2 = (newSize==0) ? nullptr : new L3Address[newSize];
    size_t minSize = originators_arraysize < newSize ? originators_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        originators2[i] = this->originators[i];
    delete [] this->originators;
    this->originators = originators2;
    originators_arraysize = newSize;
}

void OlsrTcGroup::setOriginators(size_t k, const L3Address& originators)
{
    if (k >= originators_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)originators_arraysize, (unsigned long)k);
    handleChange();
    this->originators[k] = originators;
}

void OlsrTcGroup::insertOriginators(size_t k, const L3Address& originators)
{
    if (k > originators_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)originators_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = originators_arraysize + 1;
    L3Address *originators2 = new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        originators2[i] = this->originators[i];
    originators2[k] = originators;
    for (i = k + 1; i < newSize; i++)
        originators2[i] = this->originators[i-1];
    delete [] this->originators;
    this->originators = originators2;
    originators_arraysize = newSize;
}

void OlsrTcGroup::appendOriginators(const L3Address& originators)
{
    insertOriginators(originators_arraysize, originators);
}

void OlsrTcGroup::eraseOriginators(size_t k)
{
    if (k >= originators_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)originators_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = originators_arraysize - 1;
    L3Address *originators2 = (newSize == 0) ? nullptr : new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        originators2[i] = this->originators[i];
    for (i = k; i < newSize; i++)
        originators2[i] = this->originators[i+1];
    delete [] this->originators;
    this->originators = originators2;
    originators_arraysize = newSize;
}

size_t OlsrTcGroup::getAnsnsArraySize() const
{
    return ansns_arraysize;
}

uint16_t OlsrTcGroup::getAnsns(size_t k) const
{
    if (k >= ansns_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)ansns_arraysize, (unsigned long)k);
    return this->ansns[k];
}

void OlsrTcGroup::setAnsnsArraySize(size_t newSize)
{
    handleChange();
    uint16_t *ansns2 = (newSize==0) ? nullptr : new uint16_t[newSize];
    size_t minSize = ansns_arraysize < newSize ? ansns_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        ansns2[i] = this->ansns[i];
    for (size_t i = minSize; i < newSize; i++)
        ansns2[i] = 0;
    delete [] this->ansns;
    this->ansns = ansns2;
    ansns_arraysize = newSize;
}

void OlsrTcGroup::setAnsns(size_t k, uint16_t ansns)
{
    if (k >= ansns_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)ansns_arraysize, (unsigned long)k);
    handleChange();
    this->ansns[k] = ansns;
}

void OlsrTcGroup::insertAnsns(size_t k, uint16_t ansns)
{
    if (k > ansns_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)ansns_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = ansns_arraysize + 1;
    uint16_t *ansns2 = new uint16_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        ansns2[i] = this->ansns[i];
    ansns2[k] = ansns;
    for (i = k + 1; i < newSize; i++)
        ansns2[i] = this->ansns[i-1];
    delete [] this->ansns;
    this->ansns = ansns2;
    ansns_arraysize = newSize;
}

void OlsrTcGroup::appendAnsns(uint16_t ansns)
{
    insertAnsns(ansns_arraysize, ansns);
}

void OlsrTcGroup::eraseAnsns(size_t k)
{
    if (k >= ansns_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)ansns_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = ansns_arraysize - 1;
    uint16_t *ansns2 = (newSize == 0) ? nullptr : new uint16_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        ansns2[i] = this->ansns[i];
    for (i = k; i < newSize; i++)
        ansns2[i] = this->ansns[i+1];
    delete [] this->ansns;
    this->ansns = ansns2;
    ansns_arraysize = newSize;
}

size_t OlsrTcGroup::getSeqNumsArraySize() const
{
    return seqNums_arraysize;
}

uint16_t OlsrTcGroup::getSeqNums(size_t k) const
{
    if (k >= seqNums_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)seqNums_arraysize, (unsigned long)k);
    return this->seqNums[k];
}

void OlsrTcGroup::setSeqNumsArraySize(size_t newSize)
{
    handleChange();
    uint16_t *seqNums2 = (newSize==0) ? nullptr : new uint16_t[newSize];
    size_t minSize = seqNums_arraysize < newSize ? seqNums_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        seqNums2[i] = this->seqNums[i];
    for (size_t i = minSize; i < newSize; i++)
        seqNums2[i] = 0;
    delete [] this->seqNums;
    this->seqNums = seqNums2;
    seqNums_arraysize = newSize;
}

void OlsrTcGroup::setSeqNums(size_t k, uint16_t seqNums)
{
    if (k >= seqNums_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)seqNums_arraysize, (unsigned long)k);
    handleChange();
    this->seqNums[k] = seqNums;
}

void OlsrTcGroup::insertSeqNums(size_t k, uint16_t seqNums)
{
    if (k > seqNums_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)seqNums_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = seqNums_arraysize + 1;
    uint16_t *seqNums2 = new uint16_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        seqNums2[i] = this->seqNums[i];
    seqNums2[k] = seqNums;
    for (i = k + 1; i < newSize; i++)
        seqNums2[i] = this->seqNums[i-1];
    delete [] this->seqNums;
    this->seqNums = seqNums2;
    seqNums_arraysize = newSize;
}

void OlsrTcGroup::appendSeqNums(uint16_t seqNums)
{
    insertSeqNums(seqNums_arraysize, seqNums);
}

void OlsrTcGroup::eraseSeqNums(size_t k)
{
    if (k >= seqNums_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)seqNums_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = seqNums_arraysize - 1;
    uint16_t *seqNums2 = (newSize == 0) ? nullptr : new uint16_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        seqNums2[i] = this->seqNums[i];
    for (i = k; i < newSize; i++)
        seqNums2[i] = this->seqNums[i+1];
    delete [] this->seqNums;
    this->seqNums = seqNums2;
    seqNums_arraysize = newSize;
}

size_t OlsrTcGroup::getHopCountsArraySize() const
{
    return hopCounts_arraysize;
}

uint8_t OlsrTcGroup::getHopCounts(size_t k) const
{
    if (k >= hopCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)hopCounts_arraysize, (unsigned long)k);
    return this->hopCounts[k];
}

void OlsrTcGroup::setHopCountsArraySize(size_t newSize)
{
    handleChange();
    uint8_t *hopCounts2 = (newSize==0) ? nullptr : new uint8_t[newSize];
    size_t minSize = hopCounts_arraysize < newSize ? hopCounts_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        hopCounts2[i] = this->hopCounts[i];
    for (size_t i = minSize; i < newSize; i++)
        hopCounts2[i] = 0;
    delete [] this->hopCounts;
    this->hopCounts = hopCounts2;
    hopCounts_arraysize = newSize;
}

void OlsrTcGroup::setHopCounts(size_t k, uint8_t hopCounts)
{
    if (k >= hopCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)hopCounts_arraysize, (unsigned long)k);
    handleChange();
    this->hopCounts[k] = hopCounts;
}

void OlsrTcGroup::insertHopCounts(size_t k, uint8_t hopCounts)
{
    if (k > hopCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)hopCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = hopCounts_arraysize + 1;
    uint8_t *hopCounts2 = new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        hopCounts2[i] = this->hopCounts[i];
    hopCounts2[k] = hopCounts;
    for (i = k + 1; i < newSize; i++)
        hopCounts2[i] = this->hopCounts[i-1];
    delete [] this->hopCounts;
    this->hopCounts = hopCounts2;
    hopCounts_arraysize = newSize;
}

void OlsrTcGroup::appendHopCounts(uint8_t hopCounts)
{
    insertHopCounts(hopCounts_arraysize, hopCounts);
}

void OlsrTcGroup::eraseHopCounts(size_t k)
{
    if (k >= hopCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)hopCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = hopCounts_arraysize - 1;
    uint8_t *hopCounts2 = (newSize == 0) ? nullptr : new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        hopCounts2[i] = this->hopCounts[i];
    for (i = k; i < newSize; i++)
        hopCounts2[i] = this->hopCounts[i+1];
    delete [] this->hopCounts;
    this->hopCounts = hopCounts2;
    hopCounts_arraysize = newSize;
}

size_t OlsrTcGroup::getAdvCountsArraySize() const
{
    return advCounts_arraysize;
}

int OlsrTcGroup::getAdvCounts(size_t k) const
{
    if (k >= advCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advCounts_arraysize, (unsigned long)k);
    return this->advCounts[k];
}

void OlsrTcGroup::setAdvCountsArraySize(size_t newSize)
{
    handleChange();
    int *advCounts2 = (newSize==0) ? nullptr : new int[newSize];
    size_t minSize = advCounts_arraysize < newSize ? advCounts_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        advCounts2[i] = this->advCounts[i];
    for (size_t i = minSize; i < newSize; i++)
        advCounts2[i] = 0;
    delete [] this->advCounts;
    this->advCounts = advCounts2;
    advCounts_arraysize = newSize;
}

void OlsrTcGroup::setAdvCounts(size_t k, int advCounts)
{
    if (k >= advCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advCounts_arraysize, (unsigned long)k);
    handleChange();
    this->advCounts[k] = advCounts;
}

void OlsrTcGroup::insertAdvCounts(size_t k, int advCounts)
{
    if (k > advCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = advCounts_arraysize + 1;
    int *advCounts2 = new int[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        advCounts2[i] = this->advCounts[i];
    advCounts2[k] = advCounts;
    for (i = k + 1; i < newSize; i++)
        advCounts2[i] = this->advCounts[i-1];
    delete [] this->advCounts;
    this->advCounts = advCounts2;
    advCounts_arraysize = newSize;
}

void OlsrTcGroup::appendAdvCounts(int advCounts)
{
    insertAdvCounts(advCounts_arraysize, advCounts);
}

void OlsrTcGroup::eraseAdvCounts(size_t k)
{
    if (k >= advCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = advCounts_arraysize - 1;
    int *advCounts2 = (newSize == 0) ? nullptr : new int[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        advCounts2[i] = this->advCounts[i];
    for (i = k; i < newSize; i++)
        advCounts2[i] = this->advCounts[i+1];
    delete [] this->advCounts;
    this->advCounts = advCounts2;
    advCounts_arraysize = newSize;
}

size_t OlsrTcGroup::getAdvertisedNeighborsArraySize() const
{
    return advertisedNeighbors_arraysize;
}

const L3Address& OlsrTcGroup::getAdvertisedNeighbors(size_t k) const
{
    if (k >= advertisedNeighbors_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advertisedNeighbors_arraysize, (unsigned long)k);
    return this->advertisedNeighbors[k];
}

void OlsrTcGroup::setAdvertisedNeighborsArraySize(size_t newSize)
{
    handleChange();
    L3Address *advertisedNeighbors2 = (newSize==0) ? nullptr : new L3Address[newSize];
    size_t minSize = advertisedNeighbors_arraysize < newSize ? advertisedNeighbors_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        advertisedNeighbors2[i] = this->advertisedNeighbors[i];
    delete [] this->advertisedNeighbors;
    this->advertisedNeighbors = advertisedNeighbors2;
    advertisedNeighbors_arraysize = newSize;
}

void OlsrTcGroup::setAdvertisedNeighbors(size_t k, const L3Address& advertisedNeighbors)
{
    if (k >= advertisedNeighbors_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advertisedNeighbors_arraysize, (unsigned long)k);
    handleChange();
    this->advertisedNeighbors[k] = advertisedNeighbors;
}

void OlsrTcGroup::insertAdvertisedNeighbors(size_t k, const L3Address& advertisedNeighbors)
{
    if (k > advertisedNeighbors_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advertisedNeighbors_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = advertisedNeighbors_arraysize + 1;
    L3Address *advertisedNeighbors2 = new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        advertisedNeighbors2[i] = this->advertisedNeighbors[i];
    advertisedNeighbors2[k] = advertisedNeighbors;
    for (i = k + 1; i < newSize; i++)
        advertisedNeighbors2[i] = this->advertisedNeighbors[i-1];
    delete [] this->advertisedNeighbors;
    this->advertisedNeighbors = advertisedNeighbors2;
    advertisedNeighbors_arraysize = newSize;
}

void OlsrTcGroup::appendAdvertisedNeighbors(const L3Address& advertisedNeighbors)
{
    insertAdvertisedNeighbors(advertisedNeighbors_arraysize, advertisedNeighbors);
}

void OlsrTcGroup::eraseAdvertisedNeighbors(size_t k)
{
    if (k >= advertisedNeighbors_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advertisedNeighbors_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = advertisedNeighbors_arraysize - 1;
    L3Address *advertisedNeighbors2 = (newSize == 0) ? nullptr : new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        advertisedNeighbors2[i] = this->advertisedNeighbors[i];
    for (i = k; i < newSize; i++)
        advertisedNeighbors2[i] = this->advertisedNeighbors[i+1];
    delete [] this->advertisedNeighbors;
    this->advertisedNeighbors = advertisedNeighbors2;
    advertisedNeighbors_arraysize = newSize;
}

class OlsrTcGroupDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_genTime,
        FIELD_originators,
        FIELD_ansns,
        FIELD_seqNums,
        FIELD_hopCounts,
        FIELD_advCounts,
        FIELD_advertisedNeighbors,
    };
  public:
    OlsrTcGroupDescriptor();
    virtual ~OlsrTcGroupDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(OlsrTcGroupDescriptor)

OlsrTcGroupDescriptor::OlsrTcGroupDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::OlsrTcGroup)), "inet::OlsrControlPacket")
{
    propertyNames = nullptr;
}

OlsrTcGroupDescriptor::~OlsrTcGroupDescriptor()
{
    delete[] propertyNames;
}

bool OlsrTcGroupDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<OlsrTcGroup *>(obj)!=nullptr;
}

const char **OlsrTcGroupDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *OlsrTcGroupDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int OlsrTcGroupDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 7+base->getFieldCount() : 7;
}

unsigned int OlsrTcGroupDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_genTime
        FD_ISARRAY | FD_ISRESIZABLE,    // FIELD_originators
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_ansns
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_seqNums
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_hopCounts
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_advCounts
        FD_ISARRAY | FD_ISRESIZABLE,    // FIELD_advertisedNeighbors
    };
    return (field >= 0 && field < 7) ? fieldTypeFlags[field] : 0;
}

const char *OlsrTcGroupDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "genTime",
        "originators",
        "ansns",
        "seqNums",
        "hopCounts",
        "advCounts",
        "advertisedNeighbors",
    };
    return (field >= 0 && field < 7) ? fieldNames[field] : nullptr;
}

int OlsrTcGroupDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "genTime") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "originators") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "ansns") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "seqNums") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "hopCounts") == 0) return baseIndex + 4;
    if (strcmp(fieldName, "advCounts") == 0) return baseIndex + 5;
    if (strcmp(fieldName, "advertisedNeighbors") == 0) return baseIndex + 6;
    return base ? base->findField(fieldName) : -1;
}

const char *OlsrTcGroupDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "omnetpp::simtime_t",    // FIELD_genTime
        "inet::L3Address",    // FIELD_originators
        "uint16",    // FIELD_ansns
        "uint16",    // FIELD_seqNums
        "uint8",    // FIELD_hopCounts
        "int",    // FIELD_advCounts
        "inet::L3Address",    // FIELD_advertisedNeighbors
    };
    return (field >= 0 && field < 7) ? fieldTypeStrings[field] : nullptr;
}

const char **OlsrTcGroupDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

const char *OlsrTcGroupDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    }
}

int OlsrTcGroupDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_originators: return pp->getOriginatorsArraySize();
        case FIELD_ansns: return pp->getAnsnsArraySize();
        case FIELD_seqNums: return pp->getSeqNumsArraySize();
        case FIELD_hopCounts: return pp->getHopCountsArraySize();
        case FIELD_advCounts: return pp->getAdvCountsArraySize();
        case FIELD_advertisedNeighbors: return pp->getAdvertisedNeighborsArraySize();
        default: return 0;
    }
}

void OlsrTcGroupDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_originators: pp->setOriginatorsArraySize(size); break;
        case FIELD_ansns: pp->setAnsnsArraySize(size); break;
        case FIELD_seqNums: pp->setSeqNumsArraySize(size); break;
        case FIELD_hopCounts: pp->setHopCountsArraySize(size); break;
        case FIELD_advCounts: pp->setAdvCountsArraySize(size); break;
        case FIELD_advertisedNeighbors: pp->setAdvertisedNeighborsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'OlsrTcGroup'", field);
    }
}

const char *OlsrTcGroupDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string OlsrTcGroupDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: return simtime2string(pp->getGenTime());
        case FIELD_originators: return pp->getOriginators(i).str();
        case FIELD_ansns: return ulong2string(pp->getAnsns(i));
        case FIELD_seqNums: return ulong2string(pp->getSeqNums(i));
        case FIELD_hopCounts: return ulong2string(pp->getHopCounts(i));
        case FIELD_advCounts: return long2string(pp->getAdvCounts(i));
        case FIELD_advertisedNeighbors: return pp->getAdvertisedNeighbors(i).str();
        default: return "";
    }
}

void OlsrTcGroupDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: pp->setGenTime(string2simtime(value)); break;
        case FIELD_ansns: pp->setAnsns(i,string2ulong(value)); break;
        case FIELD_seqNums: pp->setSeqNums(i,string2ulong(value)); break;
        case FIELD_hopCounts: pp->setHopCounts(i,string2ulong(value)); break;
        case FIELD_advCounts: pp->setAdvCounts(i,string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrTcGroup'", field);
    }
}

omnetpp::cValue OlsrTcGroupDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: return pp->getGenTime().dbl();
        case FIELD_originators: return omnetpp::toAnyPtr(&pp->getOriginators(i)); break;
        case FIELD_ansns: return (omnetpp::intval_t)(pp->getAnsns(i));
        case FIELD_seqNums: return (omnetpp::intval_t)(pp->getSeqNums(i));
        case FIELD_hopCounts: return (omnetpp::intval_t)(pp->getHopCounts(i));
        case FIELD_advCounts: return pp->getAdvCounts(i);
        case FIELD_advertisedNeighbors: return omnetpp::toAnyPtr(&pp->getAdvertisedNeighbors(i)); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'OlsrTcGroup' as cValue -- field index out of range?", field);
    }
}

void OlsrTcGroupDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: pp->setGenTime(value.doubleValue()); break;
        case FIELD_ansns: pp->setAnsns(i,omnetpp::checked_int_cast<uint16_t>(value.intValue())); break;
        case FIELD_seqNums: pp->setSeqNums(i,omnetpp::checked_int_cast<uint16_t>(value.intValue())); break;
        case FIELD_hopCounts: pp->setHopCounts(i,omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_advCounts: pp->setAdvCounts(i,omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrTcGroup'", field);
    }
}

const char *OlsrTcGroupDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        default: return nullptr;
    };
}

omnetpp::any_ptr OlsrTcGroupDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_originators: return omnetpp::toAnyPtr(&pp->getOriginators(i)); break;
        case FIELD_advertisedNeighbors: return omnetpp::toAnyPtr(&pp->getAdvertisedNeighbors(i)); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void OlsrTcGroupDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrTcGroup *pp = omnetpp::fromAnyPtr<OlsrTcGroup>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrTcGroup'", field);
    }
}

Register_Class(OlsrPacketHolderMessage)

OlsrPacketHolderMessage::OlsrPacketHolderMessage(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
    if (this->ownedPacket != nullptr) take(this->ownedPacket);
}

OlsrPacketHolderMessage::OlsrPacketHolderMessage(const OlsrPacketHolderMessage& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

OlsrPacketHolderMessage::~OlsrPacketHolderMessage()
{
    dropAndDelete(this->ownedPacket);
}

OlsrPacketHolderMessage& OlsrPacketHolderMessage::operator=(const OlsrPacketHolderMessage& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void OlsrPacketHolderMessage::copy(const OlsrPacketHolderMessage& other)
{
    dropAndDelete(this->ownedPacket);
    this->ownedPacket = other.ownedPacket;
    if (this->ownedPacket != nullptr) {
        this->ownedPacket = this->ownedPacket->dup();
        take(this->ownedPacket);
        this->ownedPacket->setName(other.ownedPacket->getName());
    }
}

void OlsrPacketHolderMessage::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->ownedPacket);
}

void OlsrPacketHolderMessage::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->ownedPacket);
}

const Packet * OlsrPacketHolderMessage::getOwnedPacket() const
{
    return this->ownedPacket;
}

void OlsrPacketHolderMessage::setOwnedPacket(Packet * ownedPacket)
{
    if (this->ownedPacket != nullptr) throw omnetpp::cRuntimeError("setOwnedPacket(): a value is already set, remove it first with removeOwnedPacket()");
    this->ownedPacket = ownedPacket;
    if (this->ownedPacket != nullptr) take(this->ownedPacket);
}

Packet * OlsrPacketHolderMessage::removeOwnedPacket()
{
    Packet * retval = this->ownedPacket;
    if (retval != nullptr) drop(retval);
    this->ownedPacket = nullptr;
    return retval;
}

class OlsrPacketHolderMessageDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_ownedPacket,
    };
  public:
    OlsrPacketHolderMessageDescriptor();
    virtual ~OlsrPacketHolderMessageDescriptor();

    virtual bool doesSupport(omnetpp::cObject *obj) const override;
    virtual const char **getPropertyNames() const override;
    virtual const char *getProperty(const char *propertyName) const override;
    virtual int getFieldCount() const override;
    virtual const char *getFieldName(int field) const override;
    virtual int findField(const char *fieldName) const override;
    virtual unsigned int getFieldTypeFlags(int field) const override;
    virtual const char *getFieldTypeString(int field) const override;
    virtual const char **getFieldPropertyNames(int field) const override;
    virtual const char *getFieldProperty(int field, const char *propertyName) const override;
    virtual int getFieldArraySize(omnetpp::any_ptr object, int field) const override;
    virtual void setFieldArraySize(omnetpp::any_ptr object, int field, int size) const override;

    virtual const char *getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const override;
    virtual std::string getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const override;
    virtual omnetpp::cValue getFieldValue(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const override;

    virtual const char *getFieldStructName(int field) const override;
    virtual omnetpp::any_ptr getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const override;
    virtual void setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const override;
};

Register_ClassDescriptor(OlsrPacketHolderMessageDescriptor)

OlsrPacketHolderMessageDescriptor::OlsrPacketHolderMessageDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::OlsrPacketHolderMessage)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

OlsrPacketHolderMessageDescriptor::~OlsrPacketHolderMessageDescriptor()
{
    delete[] propertyNames;
}

bool OlsrPacketHolderMessageDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<OlsrPacketHolderMessage *>(obj)!=nullptr;
}

const char **OlsrPacketHolderMessageDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *OlsrPacketHolderMessageDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int OlsrPacketHolderMessageDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int OlsrPacketHolderMessageDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISCOMPOUND | FD_ISPOINTER | FD_ISCOBJECT | FD_ISCOWNEDOBJECT | FD_ISREPLACEABLE,    // FIELD_ownedPacket
    };
    return (field >= 0 && field < 1) ? fieldTypeFlags[field] : 0;
}

const char *OlsrPacketHolderMessageDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "ownedPacket",
    };
    return (field >= 0 && field < 1) ? fieldNames[field] : nullptr;
}

int OlsrPacketHolderMessageDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "ownedPacket") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *OlsrPacketHolderMessageDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::Packet",    // FIELD_ownedPacket
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **OlsrPacketHolderMessageDescriptor::getFieldPropertyNames(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldPropertyNames(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_ownedPacket: {
            static const char *names[] = { "owned",  nullptr };
            return names;
        }
        default: return nullptr;
    }
}

const char *OlsrPacketHolderMessageDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_ownedPacket:
            if (!strcmp(propertyName, "owned")) return "";
            return nullptr;
        default: return nullptr;
    }
}

int OlsrPacketHolderMessageDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void OlsrPacketHolderMessageDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'OlsrPacketHolderMessage'", field);
    }
}

const char *OlsrPacketHolderMessageDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: { const Packet * value = pp->getOwnedPacket(); return omnetpp::opp_typename(typeid(*value)); }
        default: return nullptr;
    }
}

std::string OlsrPacketHolderMessageDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: { auto obj = pp->getOwnedPacket(); return obj == nullptr ? "" : obj->str(); }
        default: return "";
    }
}

void OlsrPacketHolderMessageDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrPacketHolderMessage'", field);
    }
}

omnetpp::cValue OlsrPacketHolderMessageDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::toAnyPtr(pp->getOwnedPacket()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'OlsrPacketHolderMessage' as cValue -- field index out of range?", field);
    }
}

void OlsrPacketHolderMessageDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: pp->setOwnedPacket(omnetpp::fromAnyPtr<Packet>(value.pointerValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrPacketHolderMessage'", field);
    }
}

const char *OlsrPacketHolderMessageDescriptor::getFieldStructName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructName(field);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::opp_typename(typeid(Packet));
        default: return nullptr;
    };
}

omnetpp::any_ptr OlsrPacketHolderMessageDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::toAnyPtr(pp->getOwnedPacket()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void OlsrPacketHolderMessageDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    OlsrPacketHolderMessage *pp = omnetpp::fromAnyPtr<OlsrPacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: pp->setOwnedPacket(omnetpp::fromAnyPtr<Packet>(ptr)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'OlsrPacketHolderMessage'", field);
    }
}

}  // namespace inet

namespace omnetpp {

}  // namespace omnetpp

