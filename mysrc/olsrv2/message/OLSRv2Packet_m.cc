//
// Generated file, do not edit! Created by opp_msgtool 6.2 from ../src/nhdp/../../message/OLSRv2Packet.msg.
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
#include "OLSRv2Packet_m.h"

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

Register_Enum(inet::Olsrv2PktType, (inet::Olsrv2PktType::Olsrv2Hello, inet::Olsrv2PktType::Olsrv2Tc));

Register_Class(Olsrv2ControlPacket)

Olsrv2ControlPacket::Olsrv2ControlPacket() : ::inet::FieldsChunk()
{
}

Olsrv2ControlPacket::Olsrv2ControlPacket(const Olsrv2ControlPacket& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

Olsrv2ControlPacket::~Olsrv2ControlPacket()
{
}

Olsrv2ControlPacket& Olsrv2ControlPacket::operator=(const Olsrv2ControlPacket& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void Olsrv2ControlPacket::copy(const Olsrv2ControlPacket& other)
{
    this->packetType = other.packetType;
}

void Olsrv2ControlPacket::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->packetType);
}

void Olsrv2ControlPacket::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->packetType);
}

Olsrv2PktType Olsrv2ControlPacket::getPacketType() const
{
    return this->packetType;
}

void Olsrv2ControlPacket::setPacketType(Olsrv2PktType packetType)
{
    handleChange();
    this->packetType = packetType;
}

class Olsrv2ControlPacketDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_packetType,
    };
  public:
    Olsrv2ControlPacketDescriptor();
    virtual ~Olsrv2ControlPacketDescriptor();

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

Register_ClassDescriptor(Olsrv2ControlPacketDescriptor)

Olsrv2ControlPacketDescriptor::Olsrv2ControlPacketDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::Olsrv2ControlPacket)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

Olsrv2ControlPacketDescriptor::~Olsrv2ControlPacketDescriptor()
{
    delete[] propertyNames;
}

bool Olsrv2ControlPacketDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<Olsrv2ControlPacket *>(obj)!=nullptr;
}

const char **Olsrv2ControlPacketDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *Olsrv2ControlPacketDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int Olsrv2ControlPacketDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int Olsrv2ControlPacketDescriptor::getFieldTypeFlags(int field) const
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

const char *Olsrv2ControlPacketDescriptor::getFieldName(int field) const
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

int Olsrv2ControlPacketDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "packetType") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *Olsrv2ControlPacketDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::Olsrv2PktType",    // FIELD_packetType
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **Olsrv2ControlPacketDescriptor::getFieldPropertyNames(int field) const
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

const char *Olsrv2ControlPacketDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_packetType:
            if (!strcmp(propertyName, "enum")) return "inet::Olsrv2PktType";
            return nullptr;
        default: return nullptr;
    }
}

int Olsrv2ControlPacketDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void Olsrv2ControlPacketDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'Olsrv2ControlPacket'", field);
    }
}

const char *Olsrv2ControlPacketDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string Olsrv2ControlPacketDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        case FIELD_packetType: return enum2string(pp->getPacketType(), "inet::Olsrv2PktType");
        default: return "";
    }
}

void Olsrv2ControlPacketDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2ControlPacket'", field);
    }
}

omnetpp::cValue Olsrv2ControlPacketDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        case FIELD_packetType: return static_cast<int>(pp->getPacketType());
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'Olsrv2ControlPacket' as cValue -- field index out of range?", field);
    }
}

void Olsrv2ControlPacketDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2ControlPacket'", field);
    }
}

const char *Olsrv2ControlPacketDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr Olsrv2ControlPacketDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void Olsrv2ControlPacketDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2ControlPacket *pp = omnetpp::fromAnyPtr<Olsrv2ControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2ControlPacket'", field);
    }
}

Register_Class(Olsrv2HelloPacket)

Olsrv2HelloPacket::Olsrv2HelloPacket() : ::inet::Olsrv2ControlPacket()
{
}

Olsrv2HelloPacket::Olsrv2HelloPacket(const Olsrv2HelloPacket& other) : ::inet::Olsrv2ControlPacket(other)
{
    copy(other);
}

Olsrv2HelloPacket::~Olsrv2HelloPacket()
{
    delete [] this->linkStatus;
    delete [] this->neighStatus;
    delete [] this->addrCounts;
    delete [] this->neighAddrs;
}

Olsrv2HelloPacket& Olsrv2HelloPacket::operator=(const Olsrv2HelloPacket& other)
{
    if (this == &other) return *this;
    ::inet::Olsrv2ControlPacket::operator=(other);
    copy(other);
    return *this;
}

void Olsrv2HelloPacket::copy(const Olsrv2HelloPacket& other)
{
    this->originator = other.originator;
    this->hTime = other.hTime;
    this->validityTime = other.validityTime;
    this->willingness = other.willingness;
    this->msgSeq = other.msgSeq;
    delete [] this->linkStatus;
    this->linkStatus = (other.linkStatus_arraysize==0) ? nullptr : new uint8_t[other.linkStatus_arraysize];
    linkStatus_arraysize = other.linkStatus_arraysize;
    for (size_t i = 0; i < linkStatus_arraysize; i++) {
        this->linkStatus[i] = other.linkStatus[i];
    }
    delete [] this->neighStatus;
    this->neighStatus = (other.neighStatus_arraysize==0) ? nullptr : new uint8_t[other.neighStatus_arraysize];
    neighStatus_arraysize = other.neighStatus_arraysize;
    for (size_t i = 0; i < neighStatus_arraysize; i++) {
        this->neighStatus[i] = other.neighStatus[i];
    }
    delete [] this->addrCounts;
    this->addrCounts = (other.addrCounts_arraysize==0) ? nullptr : new int[other.addrCounts_arraysize];
    addrCounts_arraysize = other.addrCounts_arraysize;
    for (size_t i = 0; i < addrCounts_arraysize; i++) {
        this->addrCounts[i] = other.addrCounts[i];
    }
    delete [] this->neighAddrs;
    this->neighAddrs = (other.neighAddrs_arraysize==0) ? nullptr : new L3Address[other.neighAddrs_arraysize];
    neighAddrs_arraysize = other.neighAddrs_arraysize;
    for (size_t i = 0; i < neighAddrs_arraysize; i++) {
        this->neighAddrs[i] = other.neighAddrs[i];
    }
}

void Olsrv2HelloPacket::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::Olsrv2ControlPacket::parsimPack(b);
    doParsimPacking(b,this->originator);
    doParsimPacking(b,this->hTime);
    doParsimPacking(b,this->validityTime);
    doParsimPacking(b,this->willingness);
    doParsimPacking(b,this->msgSeq);
    b->pack(linkStatus_arraysize);
    doParsimArrayPacking(b,this->linkStatus,linkStatus_arraysize);
    b->pack(neighStatus_arraysize);
    doParsimArrayPacking(b,this->neighStatus,neighStatus_arraysize);
    b->pack(addrCounts_arraysize);
    doParsimArrayPacking(b,this->addrCounts,addrCounts_arraysize);
    b->pack(neighAddrs_arraysize);
    doParsimArrayPacking(b,this->neighAddrs,neighAddrs_arraysize);
}

void Olsrv2HelloPacket::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::Olsrv2ControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->originator);
    doParsimUnpacking(b,this->hTime);
    doParsimUnpacking(b,this->validityTime);
    doParsimUnpacking(b,this->willingness);
    doParsimUnpacking(b,this->msgSeq);
    delete [] this->linkStatus;
    b->unpack(linkStatus_arraysize);
    if (linkStatus_arraysize == 0) {
        this->linkStatus = nullptr;
    } else {
        this->linkStatus = new uint8_t[linkStatus_arraysize];
        doParsimArrayUnpacking(b,this->linkStatus,linkStatus_arraysize);
    }
    delete [] this->neighStatus;
    b->unpack(neighStatus_arraysize);
    if (neighStatus_arraysize == 0) {
        this->neighStatus = nullptr;
    } else {
        this->neighStatus = new uint8_t[neighStatus_arraysize];
        doParsimArrayUnpacking(b,this->neighStatus,neighStatus_arraysize);
    }
    delete [] this->addrCounts;
    b->unpack(addrCounts_arraysize);
    if (addrCounts_arraysize == 0) {
        this->addrCounts = nullptr;
    } else {
        this->addrCounts = new int[addrCounts_arraysize];
        doParsimArrayUnpacking(b,this->addrCounts,addrCounts_arraysize);
    }
    delete [] this->neighAddrs;
    b->unpack(neighAddrs_arraysize);
    if (neighAddrs_arraysize == 0) {
        this->neighAddrs = nullptr;
    } else {
        this->neighAddrs = new L3Address[neighAddrs_arraysize];
        doParsimArrayUnpacking(b,this->neighAddrs,neighAddrs_arraysize);
    }
}

const L3Address& Olsrv2HelloPacket::getOriginator() const
{
    return this->originator;
}

void Olsrv2HelloPacket::setOriginator(const L3Address& originator)
{
    handleChange();
    this->originator = originator;
}

::omnetpp::simtime_t Olsrv2HelloPacket::getHTime() const
{
    return this->hTime;
}

void Olsrv2HelloPacket::setHTime(::omnetpp::simtime_t hTime)
{
    handleChange();
    this->hTime = hTime;
}

::omnetpp::simtime_t Olsrv2HelloPacket::getValidityTime() const
{
    return this->validityTime;
}

void Olsrv2HelloPacket::setValidityTime(::omnetpp::simtime_t validityTime)
{
    handleChange();
    this->validityTime = validityTime;
}

uint8_t Olsrv2HelloPacket::getWillingness() const
{
    return this->willingness;
}

void Olsrv2HelloPacket::setWillingness(uint8_t willingness)
{
    handleChange();
    this->willingness = willingness;
}

uint16_t Olsrv2HelloPacket::getMsgSeq() const
{
    return this->msgSeq;
}

void Olsrv2HelloPacket::setMsgSeq(uint16_t msgSeq)
{
    handleChange();
    this->msgSeq = msgSeq;
}

size_t Olsrv2HelloPacket::getLinkStatusArraySize() const
{
    return linkStatus_arraysize;
}

uint8_t Olsrv2HelloPacket::getLinkStatus(size_t k) const
{
    if (k >= linkStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkStatus_arraysize, (unsigned long)k);
    return this->linkStatus[k];
}

void Olsrv2HelloPacket::setLinkStatusArraySize(size_t newSize)
{
    handleChange();
    uint8_t *linkStatus2 = (newSize==0) ? nullptr : new uint8_t[newSize];
    size_t minSize = linkStatus_arraysize < newSize ? linkStatus_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        linkStatus2[i] = this->linkStatus[i];
    for (size_t i = minSize; i < newSize; i++)
        linkStatus2[i] = 0;
    delete [] this->linkStatus;
    this->linkStatus = linkStatus2;
    linkStatus_arraysize = newSize;
}

void Olsrv2HelloPacket::setLinkStatus(size_t k, uint8_t linkStatus)
{
    if (k >= linkStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkStatus_arraysize, (unsigned long)k);
    handleChange();
    this->linkStatus[k] = linkStatus;
}

void Olsrv2HelloPacket::insertLinkStatus(size_t k, uint8_t linkStatus)
{
    if (k > linkStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkStatus_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = linkStatus_arraysize + 1;
    uint8_t *linkStatus2 = new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        linkStatus2[i] = this->linkStatus[i];
    linkStatus2[k] = linkStatus;
    for (i = k + 1; i < newSize; i++)
        linkStatus2[i] = this->linkStatus[i-1];
    delete [] this->linkStatus;
    this->linkStatus = linkStatus2;
    linkStatus_arraysize = newSize;
}

void Olsrv2HelloPacket::appendLinkStatus(uint8_t linkStatus)
{
    insertLinkStatus(linkStatus_arraysize, linkStatus);
}

void Olsrv2HelloPacket::eraseLinkStatus(size_t k)
{
    if (k >= linkStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)linkStatus_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = linkStatus_arraysize - 1;
    uint8_t *linkStatus2 = (newSize == 0) ? nullptr : new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        linkStatus2[i] = this->linkStatus[i];
    for (i = k; i < newSize; i++)
        linkStatus2[i] = this->linkStatus[i+1];
    delete [] this->linkStatus;
    this->linkStatus = linkStatus2;
    linkStatus_arraysize = newSize;
}

size_t Olsrv2HelloPacket::getNeighStatusArraySize() const
{
    return neighStatus_arraysize;
}

uint8_t Olsrv2HelloPacket::getNeighStatus(size_t k) const
{
    if (k >= neighStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighStatus_arraysize, (unsigned long)k);
    return this->neighStatus[k];
}

void Olsrv2HelloPacket::setNeighStatusArraySize(size_t newSize)
{
    handleChange();
    uint8_t *neighStatus2 = (newSize==0) ? nullptr : new uint8_t[newSize];
    size_t minSize = neighStatus_arraysize < newSize ? neighStatus_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        neighStatus2[i] = this->neighStatus[i];
    for (size_t i = minSize; i < newSize; i++)
        neighStatus2[i] = 0;
    delete [] this->neighStatus;
    this->neighStatus = neighStatus2;
    neighStatus_arraysize = newSize;
}

void Olsrv2HelloPacket::setNeighStatus(size_t k, uint8_t neighStatus)
{
    if (k >= neighStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighStatus_arraysize, (unsigned long)k);
    handleChange();
    this->neighStatus[k] = neighStatus;
}

void Olsrv2HelloPacket::insertNeighStatus(size_t k, uint8_t neighStatus)
{
    if (k > neighStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighStatus_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = neighStatus_arraysize + 1;
    uint8_t *neighStatus2 = new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        neighStatus2[i] = this->neighStatus[i];
    neighStatus2[k] = neighStatus;
    for (i = k + 1; i < newSize; i++)
        neighStatus2[i] = this->neighStatus[i-1];
    delete [] this->neighStatus;
    this->neighStatus = neighStatus2;
    neighStatus_arraysize = newSize;
}

void Olsrv2HelloPacket::appendNeighStatus(uint8_t neighStatus)
{
    insertNeighStatus(neighStatus_arraysize, neighStatus);
}

void Olsrv2HelloPacket::eraseNeighStatus(size_t k)
{
    if (k >= neighStatus_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighStatus_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = neighStatus_arraysize - 1;
    uint8_t *neighStatus2 = (newSize == 0) ? nullptr : new uint8_t[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        neighStatus2[i] = this->neighStatus[i];
    for (i = k; i < newSize; i++)
        neighStatus2[i] = this->neighStatus[i+1];
    delete [] this->neighStatus;
    this->neighStatus = neighStatus2;
    neighStatus_arraysize = newSize;
}

size_t Olsrv2HelloPacket::getAddrCountsArraySize() const
{
    return addrCounts_arraysize;
}

int Olsrv2HelloPacket::getAddrCounts(size_t k) const
{
    if (k >= addrCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)addrCounts_arraysize, (unsigned long)k);
    return this->addrCounts[k];
}

void Olsrv2HelloPacket::setAddrCountsArraySize(size_t newSize)
{
    handleChange();
    int *addrCounts2 = (newSize==0) ? nullptr : new int[newSize];
    size_t minSize = addrCounts_arraysize < newSize ? addrCounts_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        addrCounts2[i] = this->addrCounts[i];
    for (size_t i = minSize; i < newSize; i++)
        addrCounts2[i] = 0;
    delete [] this->addrCounts;
    this->addrCounts = addrCounts2;
    addrCounts_arraysize = newSize;
}

void Olsrv2HelloPacket::setAddrCounts(size_t k, int addrCounts)
{
    if (k >= addrCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)addrCounts_arraysize, (unsigned long)k);
    handleChange();
    this->addrCounts[k] = addrCounts;
}

void Olsrv2HelloPacket::insertAddrCounts(size_t k, int addrCounts)
{
    if (k > addrCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)addrCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = addrCounts_arraysize + 1;
    int *addrCounts2 = new int[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        addrCounts2[i] = this->addrCounts[i];
    addrCounts2[k] = addrCounts;
    for (i = k + 1; i < newSize; i++)
        addrCounts2[i] = this->addrCounts[i-1];
    delete [] this->addrCounts;
    this->addrCounts = addrCounts2;
    addrCounts_arraysize = newSize;
}

void Olsrv2HelloPacket::appendAddrCounts(int addrCounts)
{
    insertAddrCounts(addrCounts_arraysize, addrCounts);
}

void Olsrv2HelloPacket::eraseAddrCounts(size_t k)
{
    if (k >= addrCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)addrCounts_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = addrCounts_arraysize - 1;
    int *addrCounts2 = (newSize == 0) ? nullptr : new int[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        addrCounts2[i] = this->addrCounts[i];
    for (i = k; i < newSize; i++)
        addrCounts2[i] = this->addrCounts[i+1];
    delete [] this->addrCounts;
    this->addrCounts = addrCounts2;
    addrCounts_arraysize = newSize;
}

size_t Olsrv2HelloPacket::getNeighAddrsArraySize() const
{
    return neighAddrs_arraysize;
}

const L3Address& Olsrv2HelloPacket::getNeighAddrs(size_t k) const
{
    if (k >= neighAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighAddrs_arraysize, (unsigned long)k);
    return this->neighAddrs[k];
}

void Olsrv2HelloPacket::setNeighAddrsArraySize(size_t newSize)
{
    handleChange();
    L3Address *neighAddrs2 = (newSize==0) ? nullptr : new L3Address[newSize];
    size_t minSize = neighAddrs_arraysize < newSize ? neighAddrs_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        neighAddrs2[i] = this->neighAddrs[i];
    delete [] this->neighAddrs;
    this->neighAddrs = neighAddrs2;
    neighAddrs_arraysize = newSize;
}

void Olsrv2HelloPacket::setNeighAddrs(size_t k, const L3Address& neighAddrs)
{
    if (k >= neighAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighAddrs_arraysize, (unsigned long)k);
    handleChange();
    this->neighAddrs[k] = neighAddrs;
}

void Olsrv2HelloPacket::insertNeighAddrs(size_t k, const L3Address& neighAddrs)
{
    if (k > neighAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighAddrs_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = neighAddrs_arraysize + 1;
    L3Address *neighAddrs2 = new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        neighAddrs2[i] = this->neighAddrs[i];
    neighAddrs2[k] = neighAddrs;
    for (i = k + 1; i < newSize; i++)
        neighAddrs2[i] = this->neighAddrs[i-1];
    delete [] this->neighAddrs;
    this->neighAddrs = neighAddrs2;
    neighAddrs_arraysize = newSize;
}

void Olsrv2HelloPacket::appendNeighAddrs(const L3Address& neighAddrs)
{
    insertNeighAddrs(neighAddrs_arraysize, neighAddrs);
}

void Olsrv2HelloPacket::eraseNeighAddrs(size_t k)
{
    if (k >= neighAddrs_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)neighAddrs_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = neighAddrs_arraysize - 1;
    L3Address *neighAddrs2 = (newSize == 0) ? nullptr : new L3Address[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        neighAddrs2[i] = this->neighAddrs[i];
    for (i = k; i < newSize; i++)
        neighAddrs2[i] = this->neighAddrs[i+1];
    delete [] this->neighAddrs;
    this->neighAddrs = neighAddrs2;
    neighAddrs_arraysize = newSize;
}

class Olsrv2HelloPacketDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_originator,
        FIELD_hTime,
        FIELD_validityTime,
        FIELD_willingness,
        FIELD_msgSeq,
        FIELD_linkStatus,
        FIELD_neighStatus,
        FIELD_addrCounts,
        FIELD_neighAddrs,
    };
  public:
    Olsrv2HelloPacketDescriptor();
    virtual ~Olsrv2HelloPacketDescriptor();

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

Register_ClassDescriptor(Olsrv2HelloPacketDescriptor)

Olsrv2HelloPacketDescriptor::Olsrv2HelloPacketDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::Olsrv2HelloPacket)), "inet::Olsrv2ControlPacket")
{
    propertyNames = nullptr;
}

Olsrv2HelloPacketDescriptor::~Olsrv2HelloPacketDescriptor()
{
    delete[] propertyNames;
}

bool Olsrv2HelloPacketDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<Olsrv2HelloPacket *>(obj)!=nullptr;
}

const char **Olsrv2HelloPacketDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *Olsrv2HelloPacketDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int Olsrv2HelloPacketDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 9+base->getFieldCount() : 9;
}

unsigned int Olsrv2HelloPacketDescriptor::getFieldTypeFlags(int field) const
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
        FD_ISEDITABLE,    // FIELD_validityTime
        FD_ISEDITABLE,    // FIELD_willingness
        FD_ISEDITABLE,    // FIELD_msgSeq
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_linkStatus
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_neighStatus
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_addrCounts
        FD_ISARRAY | FD_ISRESIZABLE,    // FIELD_neighAddrs
    };
    return (field >= 0 && field < 9) ? fieldTypeFlags[field] : 0;
}

const char *Olsrv2HelloPacketDescriptor::getFieldName(int field) const
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
        "validityTime",
        "willingness",
        "msgSeq",
        "linkStatus",
        "neighStatus",
        "addrCounts",
        "neighAddrs",
    };
    return (field >= 0 && field < 9) ? fieldNames[field] : nullptr;
}

int Olsrv2HelloPacketDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "originator") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "hTime") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "validityTime") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "willingness") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "msgSeq") == 0) return baseIndex + 4;
    if (strcmp(fieldName, "linkStatus") == 0) return baseIndex + 5;
    if (strcmp(fieldName, "neighStatus") == 0) return baseIndex + 6;
    if (strcmp(fieldName, "addrCounts") == 0) return baseIndex + 7;
    if (strcmp(fieldName, "neighAddrs") == 0) return baseIndex + 8;
    return base ? base->findField(fieldName) : -1;
}

const char *Olsrv2HelloPacketDescriptor::getFieldTypeString(int field) const
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
        "omnetpp::simtime_t",    // FIELD_validityTime
        "uint8",    // FIELD_willingness
        "uint16",    // FIELD_msgSeq
        "uint8",    // FIELD_linkStatus
        "uint8",    // FIELD_neighStatus
        "int",    // FIELD_addrCounts
        "inet::L3Address",    // FIELD_neighAddrs
    };
    return (field >= 0 && field < 9) ? fieldTypeStrings[field] : nullptr;
}

const char **Olsrv2HelloPacketDescriptor::getFieldPropertyNames(int field) const
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

const char *Olsrv2HelloPacketDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int Olsrv2HelloPacketDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_linkStatus: return pp->getLinkStatusArraySize();
        case FIELD_neighStatus: return pp->getNeighStatusArraySize();
        case FIELD_addrCounts: return pp->getAddrCountsArraySize();
        case FIELD_neighAddrs: return pp->getNeighAddrsArraySize();
        default: return 0;
    }
}

void Olsrv2HelloPacketDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_linkStatus: pp->setLinkStatusArraySize(size); break;
        case FIELD_neighStatus: pp->setNeighStatusArraySize(size); break;
        case FIELD_addrCounts: pp->setAddrCountsArraySize(size); break;
        case FIELD_neighAddrs: pp->setNeighAddrsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'Olsrv2HelloPacket'", field);
    }
}

const char *Olsrv2HelloPacketDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string Olsrv2HelloPacketDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_originator: return pp->getOriginator().str();
        case FIELD_hTime: return simtime2string(pp->getHTime());
        case FIELD_validityTime: return simtime2string(pp->getValidityTime());
        case FIELD_willingness: return ulong2string(pp->getWillingness());
        case FIELD_msgSeq: return ulong2string(pp->getMsgSeq());
        case FIELD_linkStatus: return ulong2string(pp->getLinkStatus(i));
        case FIELD_neighStatus: return ulong2string(pp->getNeighStatus(i));
        case FIELD_addrCounts: return long2string(pp->getAddrCounts(i));
        case FIELD_neighAddrs: return pp->getNeighAddrs(i).str();
        default: return "";
    }
}

void Olsrv2HelloPacketDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_hTime: pp->setHTime(string2simtime(value)); break;
        case FIELD_validityTime: pp->setValidityTime(string2simtime(value)); break;
        case FIELD_willingness: pp->setWillingness(string2ulong(value)); break;
        case FIELD_msgSeq: pp->setMsgSeq(string2ulong(value)); break;
        case FIELD_linkStatus: pp->setLinkStatus(i,string2ulong(value)); break;
        case FIELD_neighStatus: pp->setNeighStatus(i,string2ulong(value)); break;
        case FIELD_addrCounts: pp->setAddrCounts(i,string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2HelloPacket'", field);
    }
}

omnetpp::cValue Olsrv2HelloPacketDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_originator: return omnetpp::toAnyPtr(&pp->getOriginator()); break;
        case FIELD_hTime: return pp->getHTime().dbl();
        case FIELD_validityTime: return pp->getValidityTime().dbl();
        case FIELD_willingness: return (omnetpp::intval_t)(pp->getWillingness());
        case FIELD_msgSeq: return (omnetpp::intval_t)(pp->getMsgSeq());
        case FIELD_linkStatus: return (omnetpp::intval_t)(pp->getLinkStatus(i));
        case FIELD_neighStatus: return (omnetpp::intval_t)(pp->getNeighStatus(i));
        case FIELD_addrCounts: return pp->getAddrCounts(i);
        case FIELD_neighAddrs: return omnetpp::toAnyPtr(&pp->getNeighAddrs(i)); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'Olsrv2HelloPacket' as cValue -- field index out of range?", field);
    }
}

void Olsrv2HelloPacketDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_hTime: pp->setHTime(value.doubleValue()); break;
        case FIELD_validityTime: pp->setValidityTime(value.doubleValue()); break;
        case FIELD_willingness: pp->setWillingness(omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_msgSeq: pp->setMsgSeq(omnetpp::checked_int_cast<uint16_t>(value.intValue())); break;
        case FIELD_linkStatus: pp->setLinkStatus(i,omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_neighStatus: pp->setNeighStatus(i,omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_addrCounts: pp->setAddrCounts(i,omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2HelloPacket'", field);
    }
}

const char *Olsrv2HelloPacketDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr Olsrv2HelloPacketDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        case FIELD_originator: return omnetpp::toAnyPtr(&pp->getOriginator()); break;
        case FIELD_neighAddrs: return omnetpp::toAnyPtr(&pp->getNeighAddrs(i)); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void Olsrv2HelloPacketDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2HelloPacket *pp = omnetpp::fromAnyPtr<Olsrv2HelloPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2HelloPacket'", field);
    }
}

Register_Class(Olsrv2TcGroup)

Olsrv2TcGroup::Olsrv2TcGroup() : ::inet::Olsrv2ControlPacket()
{
}

Olsrv2TcGroup::Olsrv2TcGroup(const Olsrv2TcGroup& other) : ::inet::Olsrv2ControlPacket(other)
{
    copy(other);
}

Olsrv2TcGroup::~Olsrv2TcGroup()
{
    delete [] this->originators;
    delete [] this->ansns;
    delete [] this->seqNums;
    delete [] this->hopCounts;
    delete [] this->advCounts;
    delete [] this->advertisedNeighbors;
}

Olsrv2TcGroup& Olsrv2TcGroup::operator=(const Olsrv2TcGroup& other)
{
    if (this == &other) return *this;
    ::inet::Olsrv2ControlPacket::operator=(other);
    copy(other);
    return *this;
}

void Olsrv2TcGroup::copy(const Olsrv2TcGroup& other)
{
    this->genTime = other.genTime;
    this->validityTime = other.validityTime;
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

void Olsrv2TcGroup::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::Olsrv2ControlPacket::parsimPack(b);
    doParsimPacking(b,this->genTime);
    doParsimPacking(b,this->validityTime);
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

void Olsrv2TcGroup::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::Olsrv2ControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->genTime);
    doParsimUnpacking(b,this->validityTime);
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

::omnetpp::simtime_t Olsrv2TcGroup::getGenTime() const
{
    return this->genTime;
}

void Olsrv2TcGroup::setGenTime(::omnetpp::simtime_t genTime)
{
    handleChange();
    this->genTime = genTime;
}

::omnetpp::simtime_t Olsrv2TcGroup::getValidityTime() const
{
    return this->validityTime;
}

void Olsrv2TcGroup::setValidityTime(::omnetpp::simtime_t validityTime)
{
    handleChange();
    this->validityTime = validityTime;
}

size_t Olsrv2TcGroup::getOriginatorsArraySize() const
{
    return originators_arraysize;
}

const L3Address& Olsrv2TcGroup::getOriginators(size_t k) const
{
    if (k >= originators_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)originators_arraysize, (unsigned long)k);
    return this->originators[k];
}

void Olsrv2TcGroup::setOriginatorsArraySize(size_t newSize)
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

void Olsrv2TcGroup::setOriginators(size_t k, const L3Address& originators)
{
    if (k >= originators_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)originators_arraysize, (unsigned long)k);
    handleChange();
    this->originators[k] = originators;
}

void Olsrv2TcGroup::insertOriginators(size_t k, const L3Address& originators)
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

void Olsrv2TcGroup::appendOriginators(const L3Address& originators)
{
    insertOriginators(originators_arraysize, originators);
}

void Olsrv2TcGroup::eraseOriginators(size_t k)
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

size_t Olsrv2TcGroup::getAnsnsArraySize() const
{
    return ansns_arraysize;
}

uint16_t Olsrv2TcGroup::getAnsns(size_t k) const
{
    if (k >= ansns_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)ansns_arraysize, (unsigned long)k);
    return this->ansns[k];
}

void Olsrv2TcGroup::setAnsnsArraySize(size_t newSize)
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

void Olsrv2TcGroup::setAnsns(size_t k, uint16_t ansns)
{
    if (k >= ansns_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)ansns_arraysize, (unsigned long)k);
    handleChange();
    this->ansns[k] = ansns;
}

void Olsrv2TcGroup::insertAnsns(size_t k, uint16_t ansns)
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

void Olsrv2TcGroup::appendAnsns(uint16_t ansns)
{
    insertAnsns(ansns_arraysize, ansns);
}

void Olsrv2TcGroup::eraseAnsns(size_t k)
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

size_t Olsrv2TcGroup::getSeqNumsArraySize() const
{
    return seqNums_arraysize;
}

uint16_t Olsrv2TcGroup::getSeqNums(size_t k) const
{
    if (k >= seqNums_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)seqNums_arraysize, (unsigned long)k);
    return this->seqNums[k];
}

void Olsrv2TcGroup::setSeqNumsArraySize(size_t newSize)
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

void Olsrv2TcGroup::setSeqNums(size_t k, uint16_t seqNums)
{
    if (k >= seqNums_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)seqNums_arraysize, (unsigned long)k);
    handleChange();
    this->seqNums[k] = seqNums;
}

void Olsrv2TcGroup::insertSeqNums(size_t k, uint16_t seqNums)
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

void Olsrv2TcGroup::appendSeqNums(uint16_t seqNums)
{
    insertSeqNums(seqNums_arraysize, seqNums);
}

void Olsrv2TcGroup::eraseSeqNums(size_t k)
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

size_t Olsrv2TcGroup::getHopCountsArraySize() const
{
    return hopCounts_arraysize;
}

uint8_t Olsrv2TcGroup::getHopCounts(size_t k) const
{
    if (k >= hopCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)hopCounts_arraysize, (unsigned long)k);
    return this->hopCounts[k];
}

void Olsrv2TcGroup::setHopCountsArraySize(size_t newSize)
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

void Olsrv2TcGroup::setHopCounts(size_t k, uint8_t hopCounts)
{
    if (k >= hopCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)hopCounts_arraysize, (unsigned long)k);
    handleChange();
    this->hopCounts[k] = hopCounts;
}

void Olsrv2TcGroup::insertHopCounts(size_t k, uint8_t hopCounts)
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

void Olsrv2TcGroup::appendHopCounts(uint8_t hopCounts)
{
    insertHopCounts(hopCounts_arraysize, hopCounts);
}

void Olsrv2TcGroup::eraseHopCounts(size_t k)
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

size_t Olsrv2TcGroup::getAdvCountsArraySize() const
{
    return advCounts_arraysize;
}

int Olsrv2TcGroup::getAdvCounts(size_t k) const
{
    if (k >= advCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advCounts_arraysize, (unsigned long)k);
    return this->advCounts[k];
}

void Olsrv2TcGroup::setAdvCountsArraySize(size_t newSize)
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

void Olsrv2TcGroup::setAdvCounts(size_t k, int advCounts)
{
    if (k >= advCounts_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advCounts_arraysize, (unsigned long)k);
    handleChange();
    this->advCounts[k] = advCounts;
}

void Olsrv2TcGroup::insertAdvCounts(size_t k, int advCounts)
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

void Olsrv2TcGroup::appendAdvCounts(int advCounts)
{
    insertAdvCounts(advCounts_arraysize, advCounts);
}

void Olsrv2TcGroup::eraseAdvCounts(size_t k)
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

size_t Olsrv2TcGroup::getAdvertisedNeighborsArraySize() const
{
    return advertisedNeighbors_arraysize;
}

const L3Address& Olsrv2TcGroup::getAdvertisedNeighbors(size_t k) const
{
    if (k >= advertisedNeighbors_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advertisedNeighbors_arraysize, (unsigned long)k);
    return this->advertisedNeighbors[k];
}

void Olsrv2TcGroup::setAdvertisedNeighborsArraySize(size_t newSize)
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

void Olsrv2TcGroup::setAdvertisedNeighbors(size_t k, const L3Address& advertisedNeighbors)
{
    if (k >= advertisedNeighbors_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)advertisedNeighbors_arraysize, (unsigned long)k);
    handleChange();
    this->advertisedNeighbors[k] = advertisedNeighbors;
}

void Olsrv2TcGroup::insertAdvertisedNeighbors(size_t k, const L3Address& advertisedNeighbors)
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

void Olsrv2TcGroup::appendAdvertisedNeighbors(const L3Address& advertisedNeighbors)
{
    insertAdvertisedNeighbors(advertisedNeighbors_arraysize, advertisedNeighbors);
}

void Olsrv2TcGroup::eraseAdvertisedNeighbors(size_t k)
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

class Olsrv2TcGroupDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_genTime,
        FIELD_validityTime,
        FIELD_originators,
        FIELD_ansns,
        FIELD_seqNums,
        FIELD_hopCounts,
        FIELD_advCounts,
        FIELD_advertisedNeighbors,
    };
  public:
    Olsrv2TcGroupDescriptor();
    virtual ~Olsrv2TcGroupDescriptor();

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

Register_ClassDescriptor(Olsrv2TcGroupDescriptor)

Olsrv2TcGroupDescriptor::Olsrv2TcGroupDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::Olsrv2TcGroup)), "inet::Olsrv2ControlPacket")
{
    propertyNames = nullptr;
}

Olsrv2TcGroupDescriptor::~Olsrv2TcGroupDescriptor()
{
    delete[] propertyNames;
}

bool Olsrv2TcGroupDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<Olsrv2TcGroup *>(obj)!=nullptr;
}

const char **Olsrv2TcGroupDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *Olsrv2TcGroupDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int Olsrv2TcGroupDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 8+base->getFieldCount() : 8;
}

unsigned int Olsrv2TcGroupDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISEDITABLE,    // FIELD_genTime
        FD_ISEDITABLE,    // FIELD_validityTime
        FD_ISARRAY | FD_ISRESIZABLE,    // FIELD_originators
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_ansns
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_seqNums
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_hopCounts
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_advCounts
        FD_ISARRAY | FD_ISRESIZABLE,    // FIELD_advertisedNeighbors
    };
    return (field >= 0 && field < 8) ? fieldTypeFlags[field] : 0;
}

const char *Olsrv2TcGroupDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "genTime",
        "validityTime",
        "originators",
        "ansns",
        "seqNums",
        "hopCounts",
        "advCounts",
        "advertisedNeighbors",
    };
    return (field >= 0 && field < 8) ? fieldNames[field] : nullptr;
}

int Olsrv2TcGroupDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "genTime") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "validityTime") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "originators") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "ansns") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "seqNums") == 0) return baseIndex + 4;
    if (strcmp(fieldName, "hopCounts") == 0) return baseIndex + 5;
    if (strcmp(fieldName, "advCounts") == 0) return baseIndex + 6;
    if (strcmp(fieldName, "advertisedNeighbors") == 0) return baseIndex + 7;
    return base ? base->findField(fieldName) : -1;
}

const char *Olsrv2TcGroupDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "omnetpp::simtime_t",    // FIELD_genTime
        "omnetpp::simtime_t",    // FIELD_validityTime
        "inet::L3Address",    // FIELD_originators
        "uint16",    // FIELD_ansns
        "uint16",    // FIELD_seqNums
        "uint8",    // FIELD_hopCounts
        "int",    // FIELD_advCounts
        "inet::L3Address",    // FIELD_advertisedNeighbors
    };
    return (field >= 0 && field < 8) ? fieldTypeStrings[field] : nullptr;
}

const char **Olsrv2TcGroupDescriptor::getFieldPropertyNames(int field) const
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

const char *Olsrv2TcGroupDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int Olsrv2TcGroupDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
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

void Olsrv2TcGroupDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_originators: pp->setOriginatorsArraySize(size); break;
        case FIELD_ansns: pp->setAnsnsArraySize(size); break;
        case FIELD_seqNums: pp->setSeqNumsArraySize(size); break;
        case FIELD_hopCounts: pp->setHopCountsArraySize(size); break;
        case FIELD_advCounts: pp->setAdvCountsArraySize(size); break;
        case FIELD_advertisedNeighbors: pp->setAdvertisedNeighborsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'Olsrv2TcGroup'", field);
    }
}

const char *Olsrv2TcGroupDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string Olsrv2TcGroupDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: return simtime2string(pp->getGenTime());
        case FIELD_validityTime: return simtime2string(pp->getValidityTime());
        case FIELD_originators: return pp->getOriginators(i).str();
        case FIELD_ansns: return ulong2string(pp->getAnsns(i));
        case FIELD_seqNums: return ulong2string(pp->getSeqNums(i));
        case FIELD_hopCounts: return ulong2string(pp->getHopCounts(i));
        case FIELD_advCounts: return long2string(pp->getAdvCounts(i));
        case FIELD_advertisedNeighbors: return pp->getAdvertisedNeighbors(i).str();
        default: return "";
    }
}

void Olsrv2TcGroupDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: pp->setGenTime(string2simtime(value)); break;
        case FIELD_validityTime: pp->setValidityTime(string2simtime(value)); break;
        case FIELD_ansns: pp->setAnsns(i,string2ulong(value)); break;
        case FIELD_seqNums: pp->setSeqNums(i,string2ulong(value)); break;
        case FIELD_hopCounts: pp->setHopCounts(i,string2ulong(value)); break;
        case FIELD_advCounts: pp->setAdvCounts(i,string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2TcGroup'", field);
    }
}

omnetpp::cValue Olsrv2TcGroupDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: return pp->getGenTime().dbl();
        case FIELD_validityTime: return pp->getValidityTime().dbl();
        case FIELD_originators: return omnetpp::toAnyPtr(&pp->getOriginators(i)); break;
        case FIELD_ansns: return (omnetpp::intval_t)(pp->getAnsns(i));
        case FIELD_seqNums: return (omnetpp::intval_t)(pp->getSeqNums(i));
        case FIELD_hopCounts: return (omnetpp::intval_t)(pp->getHopCounts(i));
        case FIELD_advCounts: return pp->getAdvCounts(i);
        case FIELD_advertisedNeighbors: return omnetpp::toAnyPtr(&pp->getAdvertisedNeighbors(i)); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'Olsrv2TcGroup' as cValue -- field index out of range?", field);
    }
}

void Olsrv2TcGroupDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_genTime: pp->setGenTime(value.doubleValue()); break;
        case FIELD_validityTime: pp->setValidityTime(value.doubleValue()); break;
        case FIELD_ansns: pp->setAnsns(i,omnetpp::checked_int_cast<uint16_t>(value.intValue())); break;
        case FIELD_seqNums: pp->setSeqNums(i,omnetpp::checked_int_cast<uint16_t>(value.intValue())); break;
        case FIELD_hopCounts: pp->setHopCounts(i,omnetpp::checked_int_cast<uint8_t>(value.intValue())); break;
        case FIELD_advCounts: pp->setAdvCounts(i,omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2TcGroup'", field);
    }
}

const char *Olsrv2TcGroupDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr Olsrv2TcGroupDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        case FIELD_originators: return omnetpp::toAnyPtr(&pp->getOriginators(i)); break;
        case FIELD_advertisedNeighbors: return omnetpp::toAnyPtr(&pp->getAdvertisedNeighbors(i)); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void Olsrv2TcGroupDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    Olsrv2TcGroup *pp = omnetpp::fromAnyPtr<Olsrv2TcGroup>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'Olsrv2TcGroup'", field);
    }
}

}  // namespace inet

namespace omnetpp {

}  // namespace omnetpp

