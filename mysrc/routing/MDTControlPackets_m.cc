//
// Generated file, do not edit! Created by opp_msgtool 6.2 from mysrc/routing/MDTControlPackets.msg.
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
#include "MDTControlPackets_m.h"

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

Register_Enum(inet::MDTControlPacketType, (inet::MDTControlPacketType::JRQ, inet::MDTControlPacketType::JRP, inet::MDTControlPacketType::NSRQ, inet::MDTControlPacketType::NSRP, inet::MDTControlPacketType::NSN, inet::MDTControlPacketType::IT, inet::MDTControlPacketType::LN, inet::MDTControlPacketType::FN, inet::MDTControlPacketType::KA, inet::MDTControlPacketType::CDQ, inet::MDTControlPacketType::CDP, inet::MDTControlPacketType::KAN));

Register_Class(MDTControlPacket)

MDTControlPacket::MDTControlPacket() : ::inet::FieldsChunk()
{
}

MDTControlPacket::MDTControlPacket(const MDTControlPacket& other) : ::inet::FieldsChunk(other)
{
    copy(other);
}

MDTControlPacket::~MDTControlPacket()
{
}

MDTControlPacket& MDTControlPacket::operator=(const MDTControlPacket& other)
{
    if (this == &other) return *this;
    ::inet::FieldsChunk::operator=(other);
    copy(other);
    return *this;
}

void MDTControlPacket::copy(const MDTControlPacket& other)
{
    this->packetType = other.packetType;
}

void MDTControlPacket::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::FieldsChunk::parsimPack(b);
    doParsimPacking(b,this->packetType);
}

void MDTControlPacket::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::FieldsChunk::parsimUnpack(b);
    doParsimUnpacking(b,this->packetType);
}

MDTControlPacketType MDTControlPacket::getPacketType() const
{
    return this->packetType;
}

void MDTControlPacket::setPacketType(MDTControlPacketType packetType)
{
    handleChange();
    this->packetType = packetType;
}

class MDTControlPacketDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_packetType,
    };
  public:
    MDTControlPacketDescriptor();
    virtual ~MDTControlPacketDescriptor();

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

Register_ClassDescriptor(MDTControlPacketDescriptor)

MDTControlPacketDescriptor::MDTControlPacketDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::MDTControlPacket)), "inet::FieldsChunk")
{
    propertyNames = nullptr;
}

MDTControlPacketDescriptor::~MDTControlPacketDescriptor()
{
    delete[] propertyNames;
}

bool MDTControlPacketDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<MDTControlPacket *>(obj)!=nullptr;
}

const char **MDTControlPacketDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *MDTControlPacketDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int MDTControlPacketDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 1+base->getFieldCount() : 1;
}

unsigned int MDTControlPacketDescriptor::getFieldTypeFlags(int field) const
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

const char *MDTControlPacketDescriptor::getFieldName(int field) const
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

int MDTControlPacketDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "packetType") == 0) return baseIndex + 0;
    return base ? base->findField(fieldName) : -1;
}

const char *MDTControlPacketDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::MDTControlPacketType",    // FIELD_packetType
    };
    return (field >= 0 && field < 1) ? fieldTypeStrings[field] : nullptr;
}

const char **MDTControlPacketDescriptor::getFieldPropertyNames(int field) const
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

const char *MDTControlPacketDescriptor::getFieldProperty(int field, const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldProperty(field, propertyName);
        field -= base->getFieldCount();
    }
    switch (field) {
        case FIELD_packetType:
            if (!strcmp(propertyName, "enum")) return "inet::MDTControlPacketType";
            return nullptr;
        default: return nullptr;
    }
}

int MDTControlPacketDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void MDTControlPacketDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'MDTControlPacket'", field);
    }
}

const char *MDTControlPacketDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string MDTControlPacketDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        case FIELD_packetType: return enum2string(pp->getPacketType(), "inet::MDTControlPacketType");
        default: return "";
    }
}

void MDTControlPacketDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MDTControlPacket'", field);
    }
}

omnetpp::cValue MDTControlPacketDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        case FIELD_packetType: return static_cast<int>(pp->getPacketType());
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'MDTControlPacket' as cValue -- field index out of range?", field);
    }
}

void MDTControlPacketDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MDTControlPacket'", field);
    }
}

const char *MDTControlPacketDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr MDTControlPacketDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void MDTControlPacketDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    MDTControlPacket *pp = omnetpp::fromAnyPtr<MDTControlPacket>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'MDTControlPacket'", field);
    }
}

KeepAlive_Base::KeepAlive_Base() : ::inet::MDTControlPacket()
{
}

KeepAlive_Base::KeepAlive_Base(const KeepAlive_Base& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

KeepAlive_Base::~KeepAlive_Base()
{
    delete [] this->coords;
}

KeepAlive_Base& KeepAlive_Base::operator=(const KeepAlive_Base& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void KeepAlive_Base::copy(const KeepAlive_Base& other)
{
    delete [] this->coords;
    this->coords = (other.coords_arraysize==0) ? nullptr : new double[other.coords_arraysize];
    coords_arraysize = other.coords_arraysize;
    for (size_t i = 0; i < coords_arraysize; i++) {
        this->coords[i] = other.coords[i];
    }
    this->attached = other.attached;
}

void KeepAlive_Base::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    b->pack(coords_arraysize);
    doParsimArrayPacking(b,this->coords,coords_arraysize);
    doParsimPacking(b,this->attached);
}

void KeepAlive_Base::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    delete [] this->coords;
    b->unpack(coords_arraysize);
    if (coords_arraysize == 0) {
        this->coords = nullptr;
    } else {
        this->coords = new double[coords_arraysize];
        doParsimArrayUnpacking(b,this->coords,coords_arraysize);
    }
    doParsimUnpacking(b,this->attached);
}

size_t KeepAlive_Base::getCoordsArraySize() const
{
    return coords_arraysize;
}

double KeepAlive_Base::getCoords(size_t k) const
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    return this->coords[k];
}

void KeepAlive_Base::setCoordsArraySize(size_t newSize)
{
    handleChange();
    double *coords2 = (newSize==0) ? nullptr : new double[newSize];
    size_t minSize = coords_arraysize < newSize ? coords_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        coords2[i] = this->coords[i];
    for (size_t i = minSize; i < newSize; i++)
        coords2[i] = 0;
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void KeepAlive_Base::setCoords(size_t k, double coords)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    this->coords[k] = coords;
}

void KeepAlive_Base::insertCoords(size_t k, double coords)
{
    if (k > coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize + 1;
    double *coords2 = new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    coords2[k] = coords;
    for (i = k + 1; i < newSize; i++)
        coords2[i] = this->coords[i-1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void KeepAlive_Base::appendCoords(double coords)
{
    insertCoords(coords_arraysize, coords);
}

void KeepAlive_Base::eraseCoords(size_t k)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize - 1;
    double *coords2 = (newSize == 0) ? nullptr : new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    for (i = k; i < newSize; i++)
        coords2[i] = this->coords[i+1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

bool KeepAlive_Base::getAttached() const
{
    return this->attached;
}

void KeepAlive_Base::setAttached(bool attached)
{
    handleChange();
    this->attached = attached;
}

class KeepAliveDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_coords,
        FIELD_attached,
    };
  public:
    KeepAliveDescriptor();
    virtual ~KeepAliveDescriptor();

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

Register_ClassDescriptor(KeepAliveDescriptor)

KeepAliveDescriptor::KeepAliveDescriptor() : omnetpp::cClassDescriptor("inet::KeepAlive", "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

KeepAliveDescriptor::~KeepAliveDescriptor()
{
    delete[] propertyNames;
}

bool KeepAliveDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<KeepAlive_Base *>(obj)!=nullptr;
}

const char **KeepAliveDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = { "customize",  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *KeepAliveDescriptor::getProperty(const char *propertyName) const
{
    if (!strcmp(propertyName, "customize")) return "true";
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int KeepAliveDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 2+base->getFieldCount() : 2;
}

unsigned int KeepAliveDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_coords
        FD_ISEDITABLE,    // FIELD_attached
    };
    return (field >= 0 && field < 2) ? fieldTypeFlags[field] : 0;
}

const char *KeepAliveDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "coords",
        "attached",
    };
    return (field >= 0 && field < 2) ? fieldNames[field] : nullptr;
}

int KeepAliveDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "coords") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "attached") == 0) return baseIndex + 1;
    return base ? base->findField(fieldName) : -1;
}

const char *KeepAliveDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "double",    // FIELD_coords
        "bool",    // FIELD_attached
    };
    return (field >= 0 && field < 2) ? fieldTypeStrings[field] : nullptr;
}

const char **KeepAliveDescriptor::getFieldPropertyNames(int field) const
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

const char *KeepAliveDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int KeepAliveDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoordsArraySize();
        default: return 0;
    }
}

void KeepAliveDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoordsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'KeepAlive_Base'", field);
    }
}

const char *KeepAliveDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string KeepAliveDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return double2string(pp->getCoords(i));
        case FIELD_attached: return bool2string(pp->getAttached());
        default: return "";
    }
}

void KeepAliveDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,string2double(value)); break;
        case FIELD_attached: pp->setAttached(string2bool(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'KeepAlive_Base'", field);
    }
}

omnetpp::cValue KeepAliveDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoords(i);
        case FIELD_attached: return pp->getAttached();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'KeepAlive_Base' as cValue -- field index out of range?", field);
    }
}

void KeepAliveDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,value.doubleValue()); break;
        case FIELD_attached: pp->setAttached(value.boolValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'KeepAlive_Base'", field);
    }
}

const char *KeepAliveDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr KeepAliveDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void KeepAliveDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAlive_Base *pp = omnetpp::fromAnyPtr<KeepAlive_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'KeepAlive_Base'", field);
    }
}

Register_Class(KeepAliveNotification)

KeepAliveNotification::KeepAliveNotification() : ::inet::MDTControlPacket()
{
}

KeepAliveNotification::KeepAliveNotification(const KeepAliveNotification& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

KeepAliveNotification::~KeepAliveNotification()
{
    delete [] this->coords;
}

KeepAliveNotification& KeepAliveNotification::operator=(const KeepAliveNotification& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void KeepAliveNotification::copy(const KeepAliveNotification& other)
{
    delete [] this->coords;
    this->coords = (other.coords_arraysize==0) ? nullptr : new double[other.coords_arraysize];
    coords_arraysize = other.coords_arraysize;
    for (size_t i = 0; i < coords_arraysize; i++) {
        this->coords[i] = other.coords[i];
    }
    this->attached = other.attached;
}

void KeepAliveNotification::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    b->pack(coords_arraysize);
    doParsimArrayPacking(b,this->coords,coords_arraysize);
    doParsimPacking(b,this->attached);
}

void KeepAliveNotification::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    delete [] this->coords;
    b->unpack(coords_arraysize);
    if (coords_arraysize == 0) {
        this->coords = nullptr;
    } else {
        this->coords = new double[coords_arraysize];
        doParsimArrayUnpacking(b,this->coords,coords_arraysize);
    }
    doParsimUnpacking(b,this->attached);
}

size_t KeepAliveNotification::getCoordsArraySize() const
{
    return coords_arraysize;
}

double KeepAliveNotification::getCoords(size_t k) const
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    return this->coords[k];
}

void KeepAliveNotification::setCoordsArraySize(size_t newSize)
{
    handleChange();
    double *coords2 = (newSize==0) ? nullptr : new double[newSize];
    size_t minSize = coords_arraysize < newSize ? coords_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        coords2[i] = this->coords[i];
    for (size_t i = minSize; i < newSize; i++)
        coords2[i] = 0;
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void KeepAliveNotification::setCoords(size_t k, double coords)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    this->coords[k] = coords;
}

void KeepAliveNotification::insertCoords(size_t k, double coords)
{
    if (k > coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize + 1;
    double *coords2 = new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    coords2[k] = coords;
    for (i = k + 1; i < newSize; i++)
        coords2[i] = this->coords[i-1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void KeepAliveNotification::appendCoords(double coords)
{
    insertCoords(coords_arraysize, coords);
}

void KeepAliveNotification::eraseCoords(size_t k)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize - 1;
    double *coords2 = (newSize == 0) ? nullptr : new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    for (i = k; i < newSize; i++)
        coords2[i] = this->coords[i+1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

bool KeepAliveNotification::getAttached() const
{
    return this->attached;
}

void KeepAliveNotification::setAttached(bool attached)
{
    handleChange();
    this->attached = attached;
}

class KeepAliveNotificationDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_coords,
        FIELD_attached,
    };
  public:
    KeepAliveNotificationDescriptor();
    virtual ~KeepAliveNotificationDescriptor();

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

Register_ClassDescriptor(KeepAliveNotificationDescriptor)

KeepAliveNotificationDescriptor::KeepAliveNotificationDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::KeepAliveNotification)), "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

KeepAliveNotificationDescriptor::~KeepAliveNotificationDescriptor()
{
    delete[] propertyNames;
}

bool KeepAliveNotificationDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<KeepAliveNotification *>(obj)!=nullptr;
}

const char **KeepAliveNotificationDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *KeepAliveNotificationDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int KeepAliveNotificationDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 2+base->getFieldCount() : 2;
}

unsigned int KeepAliveNotificationDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_coords
        FD_ISEDITABLE,    // FIELD_attached
    };
    return (field >= 0 && field < 2) ? fieldTypeFlags[field] : 0;
}

const char *KeepAliveNotificationDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "coords",
        "attached",
    };
    return (field >= 0 && field < 2) ? fieldNames[field] : nullptr;
}

int KeepAliveNotificationDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "coords") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "attached") == 0) return baseIndex + 1;
    return base ? base->findField(fieldName) : -1;
}

const char *KeepAliveNotificationDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "double",    // FIELD_coords
        "bool",    // FIELD_attached
    };
    return (field >= 0 && field < 2) ? fieldTypeStrings[field] : nullptr;
}

const char **KeepAliveNotificationDescriptor::getFieldPropertyNames(int field) const
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

const char *KeepAliveNotificationDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int KeepAliveNotificationDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoordsArraySize();
        default: return 0;
    }
}

void KeepAliveNotificationDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoordsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'KeepAliveNotification'", field);
    }
}

const char *KeepAliveNotificationDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string KeepAliveNotificationDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return double2string(pp->getCoords(i));
        case FIELD_attached: return bool2string(pp->getAttached());
        default: return "";
    }
}

void KeepAliveNotificationDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,string2double(value)); break;
        case FIELD_attached: pp->setAttached(string2bool(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'KeepAliveNotification'", field);
    }
}

omnetpp::cValue KeepAliveNotificationDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoords(i);
        case FIELD_attached: return pp->getAttached();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'KeepAliveNotification' as cValue -- field index out of range?", field);
    }
}

void KeepAliveNotificationDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,value.doubleValue()); break;
        case FIELD_attached: pp->setAttached(value.boolValue()); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'KeepAliveNotification'", field);
    }
}

const char *KeepAliveNotificationDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr KeepAliveNotificationDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        default: return omnetpp::any_ptr(nullptr);
    }
}

void KeepAliveNotificationDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    KeepAliveNotification *pp = omnetpp::fromAnyPtr<KeepAliveNotification>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'KeepAliveNotification'", field);
    }
}

Register_Class(JoinRequest)

JoinRequest::JoinRequest() : ::inet::MDTControlPacket()
{
}

JoinRequest::JoinRequest(const JoinRequest& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

JoinRequest::~JoinRequest()
{
    delete [] this->coords;
}

JoinRequest& JoinRequest::operator=(const JoinRequest& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void JoinRequest::copy(const JoinRequest& other)
{
    delete [] this->coords;
    this->coords = (other.coords_arraysize==0) ? nullptr : new double[other.coords_arraysize];
    coords_arraysize = other.coords_arraysize;
    for (size_t i = 0; i < coords_arraysize; i++) {
        this->coords[i] = other.coords[i];
    }
    this->dim = other.dim;
    this->destAddr = other.destAddr;
    this->sourceAddr = other.sourceAddr;
    this->predAddr = other.predAddr;
}

void JoinRequest::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    b->pack(coords_arraysize);
    doParsimArrayPacking(b,this->coords,coords_arraysize);
    doParsimPacking(b,this->dim);
    doParsimPacking(b,this->destAddr);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->predAddr);
}

void JoinRequest::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    delete [] this->coords;
    b->unpack(coords_arraysize);
    if (coords_arraysize == 0) {
        this->coords = nullptr;
    } else {
        this->coords = new double[coords_arraysize];
        doParsimArrayUnpacking(b,this->coords,coords_arraysize);
    }
    doParsimUnpacking(b,this->dim);
    doParsimUnpacking(b,this->destAddr);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->predAddr);
}

size_t JoinRequest::getCoordsArraySize() const
{
    return coords_arraysize;
}

double JoinRequest::getCoords(size_t k) const
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    return this->coords[k];
}

void JoinRequest::setCoordsArraySize(size_t newSize)
{
    handleChange();
    double *coords2 = (newSize==0) ? nullptr : new double[newSize];
    size_t minSize = coords_arraysize < newSize ? coords_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        coords2[i] = this->coords[i];
    for (size_t i = minSize; i < newSize; i++)
        coords2[i] = 0;
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void JoinRequest::setCoords(size_t k, double coords)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    this->coords[k] = coords;
}

void JoinRequest::insertCoords(size_t k, double coords)
{
    if (k > coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize + 1;
    double *coords2 = new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    coords2[k] = coords;
    for (i = k + 1; i < newSize; i++)
        coords2[i] = this->coords[i-1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void JoinRequest::appendCoords(double coords)
{
    insertCoords(coords_arraysize, coords);
}

void JoinRequest::eraseCoords(size_t k)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize - 1;
    double *coords2 = (newSize == 0) ? nullptr : new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    for (i = k; i < newSize; i++)
        coords2[i] = this->coords[i+1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

int JoinRequest::getDim() const
{
    return this->dim;
}

void JoinRequest::setDim(int dim)
{
    handleChange();
    this->dim = dim;
}

const L3Address& JoinRequest::getDestAddr() const
{
    return this->destAddr;
}

void JoinRequest::setDestAddr(const L3Address& destAddr)
{
    handleChange();
    this->destAddr = destAddr;
}

const L3Address& JoinRequest::getSourceAddr() const
{
    return this->sourceAddr;
}

void JoinRequest::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& JoinRequest::getPredAddr() const
{
    return this->predAddr;
}

void JoinRequest::setPredAddr(const L3Address& predAddr)
{
    handleChange();
    this->predAddr = predAddr;
}

class JoinRequestDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_coords,
        FIELD_dim,
        FIELD_destAddr,
        FIELD_sourceAddr,
        FIELD_predAddr,
    };
  public:
    JoinRequestDescriptor();
    virtual ~JoinRequestDescriptor();

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

Register_ClassDescriptor(JoinRequestDescriptor)

JoinRequestDescriptor::JoinRequestDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::JoinRequest)), "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

JoinRequestDescriptor::~JoinRequestDescriptor()
{
    delete[] propertyNames;
}

bool JoinRequestDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<JoinRequest *>(obj)!=nullptr;
}

const char **JoinRequestDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *JoinRequestDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int JoinRequestDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 5+base->getFieldCount() : 5;
}

unsigned int JoinRequestDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_coords
        FD_ISEDITABLE,    // FIELD_dim
        0,    // FIELD_destAddr
        0,    // FIELD_sourceAddr
        0,    // FIELD_predAddr
    };
    return (field >= 0 && field < 5) ? fieldTypeFlags[field] : 0;
}

const char *JoinRequestDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "coords",
        "dim",
        "destAddr",
        "sourceAddr",
        "predAddr",
    };
    return (field >= 0 && field < 5) ? fieldNames[field] : nullptr;
}

int JoinRequestDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "coords") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "dim") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "predAddr") == 0) return baseIndex + 4;
    return base ? base->findField(fieldName) : -1;
}

const char *JoinRequestDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "double",    // FIELD_coords
        "int",    // FIELD_dim
        "inet::L3Address",    // FIELD_destAddr
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_predAddr
    };
    return (field >= 0 && field < 5) ? fieldTypeStrings[field] : nullptr;
}

const char **JoinRequestDescriptor::getFieldPropertyNames(int field) const
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

const char *JoinRequestDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int JoinRequestDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoordsArraySize();
        default: return 0;
    }
}

void JoinRequestDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoordsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'JoinRequest'", field);
    }
}

const char *JoinRequestDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string JoinRequestDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return double2string(pp->getCoords(i));
        case FIELD_dim: return long2string(pp->getDim());
        case FIELD_destAddr: return pp->getDestAddr().str();
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_predAddr: return pp->getPredAddr().str();
        default: return "";
    }
}

void JoinRequestDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,string2double(value)); break;
        case FIELD_dim: pp->setDim(string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'JoinRequest'", field);
    }
}

omnetpp::cValue JoinRequestDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoords(i);
        case FIELD_dim: return pp->getDim();
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'JoinRequest' as cValue -- field index out of range?", field);
    }
}

void JoinRequestDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,value.doubleValue()); break;
        case FIELD_dim: pp->setDim(omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'JoinRequest'", field);
    }
}

const char *JoinRequestDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr JoinRequestDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void JoinRequestDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinRequest *pp = omnetpp::fromAnyPtr<JoinRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'JoinRequest'", field);
    }
}

JoinReply_Base::JoinReply_Base() : ::inet::MDTControlPacket()
{
}

JoinReply_Base::JoinReply_Base(const JoinReply_Base& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

JoinReply_Base::~JoinReply_Base()
{
}

JoinReply_Base& JoinReply_Base::operator=(const JoinReply_Base& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void JoinReply_Base::copy(const JoinReply_Base& other)
{
    this->sourceAddr = other.sourceAddr;
    this->predAddr = other.predAddr;
    this->destAddr = other.destAddr;
}

void JoinReply_Base::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->predAddr);
    doParsimPacking(b,this->destAddr);
}

void JoinReply_Base::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->predAddr);
    doParsimUnpacking(b,this->destAddr);
}

const L3Address& JoinReply_Base::getSourceAddr() const
{
    return this->sourceAddr;
}

void JoinReply_Base::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& JoinReply_Base::getPredAddr() const
{
    return this->predAddr;
}

void JoinReply_Base::setPredAddr(const L3Address& predAddr)
{
    handleChange();
    this->predAddr = predAddr;
}

const L3Address& JoinReply_Base::getDestAddr() const
{
    return this->destAddr;
}

void JoinReply_Base::setDestAddr(const L3Address& destAddr)
{
    handleChange();
    this->destAddr = destAddr;
}

class JoinReplyDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sourceAddr,
        FIELD_predAddr,
        FIELD_destAddr,
    };
  public:
    JoinReplyDescriptor();
    virtual ~JoinReplyDescriptor();

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

Register_ClassDescriptor(JoinReplyDescriptor)

JoinReplyDescriptor::JoinReplyDescriptor() : omnetpp::cClassDescriptor("inet::JoinReply", "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

JoinReplyDescriptor::~JoinReplyDescriptor()
{
    delete[] propertyNames;
}

bool JoinReplyDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<JoinReply_Base *>(obj)!=nullptr;
}

const char **JoinReplyDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = { "customize",  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *JoinReplyDescriptor::getProperty(const char *propertyName) const
{
    if (!strcmp(propertyName, "customize")) return "true";
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int JoinReplyDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 3+base->getFieldCount() : 3;
}

unsigned int JoinReplyDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_sourceAddr
        0,    // FIELD_predAddr
        0,    // FIELD_destAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeFlags[field] : 0;
}

const char *JoinReplyDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sourceAddr",
        "predAddr",
        "destAddr",
    };
    return (field >= 0 && field < 3) ? fieldNames[field] : nullptr;
}

int JoinReplyDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "predAddr") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 2;
    return base ? base->findField(fieldName) : -1;
}

const char *JoinReplyDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_predAddr
        "inet::L3Address",    // FIELD_destAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeStrings[field] : nullptr;
}

const char **JoinReplyDescriptor::getFieldPropertyNames(int field) const
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

const char *JoinReplyDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int JoinReplyDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void JoinReplyDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'JoinReply_Base'", field);
    }
}

const char *JoinReplyDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string JoinReplyDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_predAddr: return pp->getPredAddr().str();
        case FIELD_destAddr: return pp->getDestAddr().str();
        default: return "";
    }
}

void JoinReplyDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'JoinReply_Base'", field);
    }
}

omnetpp::cValue JoinReplyDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'JoinReply_Base' as cValue -- field index out of range?", field);
    }
}

void JoinReplyDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'JoinReply_Base'", field);
    }
}

const char *JoinReplyDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr JoinReplyDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void JoinReplyDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    JoinReply_Base *pp = omnetpp::fromAnyPtr<JoinReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'JoinReply_Base'", field);
    }
}

Register_Class(NeighborSetRequest)

NeighborSetRequest::NeighborSetRequest() : ::inet::MDTControlPacket()
{
}

NeighborSetRequest::NeighborSetRequest(const NeighborSetRequest& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

NeighborSetRequest::~NeighborSetRequest()
{
    delete [] this->coords;
}

NeighborSetRequest& NeighborSetRequest::operator=(const NeighborSetRequest& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void NeighborSetRequest::copy(const NeighborSetRequest& other)
{
    delete [] this->coords;
    this->coords = (other.coords_arraysize==0) ? nullptr : new double[other.coords_arraysize];
    coords_arraysize = other.coords_arraysize;
    for (size_t i = 0; i < coords_arraysize; i++) {
        this->coords[i] = other.coords[i];
    }
    this->dim = other.dim;
    this->sourceAddr = other.sourceAddr;
    this->destAddr = other.destAddr;
    this->predAddr = other.predAddr;
}

void NeighborSetRequest::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    b->pack(coords_arraysize);
    doParsimArrayPacking(b,this->coords,coords_arraysize);
    doParsimPacking(b,this->dim);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->destAddr);
    doParsimPacking(b,this->predAddr);
}

void NeighborSetRequest::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    delete [] this->coords;
    b->unpack(coords_arraysize);
    if (coords_arraysize == 0) {
        this->coords = nullptr;
    } else {
        this->coords = new double[coords_arraysize];
        doParsimArrayUnpacking(b,this->coords,coords_arraysize);
    }
    doParsimUnpacking(b,this->dim);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->destAddr);
    doParsimUnpacking(b,this->predAddr);
}

size_t NeighborSetRequest::getCoordsArraySize() const
{
    return coords_arraysize;
}

double NeighborSetRequest::getCoords(size_t k) const
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    return this->coords[k];
}

void NeighborSetRequest::setCoordsArraySize(size_t newSize)
{
    handleChange();
    double *coords2 = (newSize==0) ? nullptr : new double[newSize];
    size_t minSize = coords_arraysize < newSize ? coords_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        coords2[i] = this->coords[i];
    for (size_t i = minSize; i < newSize; i++)
        coords2[i] = 0;
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void NeighborSetRequest::setCoords(size_t k, double coords)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    this->coords[k] = coords;
}

void NeighborSetRequest::insertCoords(size_t k, double coords)
{
    if (k > coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize + 1;
    double *coords2 = new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    coords2[k] = coords;
    for (i = k + 1; i < newSize; i++)
        coords2[i] = this->coords[i-1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void NeighborSetRequest::appendCoords(double coords)
{
    insertCoords(coords_arraysize, coords);
}

void NeighborSetRequest::eraseCoords(size_t k)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize - 1;
    double *coords2 = (newSize == 0) ? nullptr : new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    for (i = k; i < newSize; i++)
        coords2[i] = this->coords[i+1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

int NeighborSetRequest::getDim() const
{
    return this->dim;
}

void NeighborSetRequest::setDim(int dim)
{
    handleChange();
    this->dim = dim;
}

const L3Address& NeighborSetRequest::getSourceAddr() const
{
    return this->sourceAddr;
}

void NeighborSetRequest::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& NeighborSetRequest::getDestAddr() const
{
    return this->destAddr;
}

void NeighborSetRequest::setDestAddr(const L3Address& destAddr)
{
    handleChange();
    this->destAddr = destAddr;
}

const L3Address& NeighborSetRequest::getPredAddr() const
{
    return this->predAddr;
}

void NeighborSetRequest::setPredAddr(const L3Address& predAddr)
{
    handleChange();
    this->predAddr = predAddr;
}

class NeighborSetRequestDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_coords,
        FIELD_dim,
        FIELD_sourceAddr,
        FIELD_destAddr,
        FIELD_predAddr,
    };
  public:
    NeighborSetRequestDescriptor();
    virtual ~NeighborSetRequestDescriptor();

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

Register_ClassDescriptor(NeighborSetRequestDescriptor)

NeighborSetRequestDescriptor::NeighborSetRequestDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::NeighborSetRequest)), "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

NeighborSetRequestDescriptor::~NeighborSetRequestDescriptor()
{
    delete[] propertyNames;
}

bool NeighborSetRequestDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<NeighborSetRequest *>(obj)!=nullptr;
}

const char **NeighborSetRequestDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *NeighborSetRequestDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int NeighborSetRequestDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 5+base->getFieldCount() : 5;
}

unsigned int NeighborSetRequestDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_coords
        FD_ISEDITABLE,    // FIELD_dim
        0,    // FIELD_sourceAddr
        0,    // FIELD_destAddr
        0,    // FIELD_predAddr
    };
    return (field >= 0 && field < 5) ? fieldTypeFlags[field] : 0;
}

const char *NeighborSetRequestDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "coords",
        "dim",
        "sourceAddr",
        "destAddr",
        "predAddr",
    };
    return (field >= 0 && field < 5) ? fieldNames[field] : nullptr;
}

int NeighborSetRequestDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "coords") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "dim") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "predAddr") == 0) return baseIndex + 4;
    return base ? base->findField(fieldName) : -1;
}

const char *NeighborSetRequestDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "double",    // FIELD_coords
        "int",    // FIELD_dim
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_destAddr
        "inet::L3Address",    // FIELD_predAddr
    };
    return (field >= 0 && field < 5) ? fieldTypeStrings[field] : nullptr;
}

const char **NeighborSetRequestDescriptor::getFieldPropertyNames(int field) const
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

const char *NeighborSetRequestDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int NeighborSetRequestDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoordsArraySize();
        default: return 0;
    }
}

void NeighborSetRequestDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoordsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'NeighborSetRequest'", field);
    }
}

const char *NeighborSetRequestDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string NeighborSetRequestDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return double2string(pp->getCoords(i));
        case FIELD_dim: return long2string(pp->getDim());
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_destAddr: return pp->getDestAddr().str();
        case FIELD_predAddr: return pp->getPredAddr().str();
        default: return "";
    }
}

void NeighborSetRequestDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,string2double(value)); break;
        case FIELD_dim: pp->setDim(string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetRequest'", field);
    }
}

omnetpp::cValue NeighborSetRequestDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoords(i);
        case FIELD_dim: return pp->getDim();
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'NeighborSetRequest' as cValue -- field index out of range?", field);
    }
}

void NeighborSetRequestDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,value.doubleValue()); break;
        case FIELD_dim: pp->setDim(omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetRequest'", field);
    }
}

const char *NeighborSetRequestDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr NeighborSetRequestDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void NeighborSetRequestDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetRequest *pp = omnetpp::fromAnyPtr<NeighborSetRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetRequest'", field);
    }
}

NeighborSetReply_Base::NeighborSetReply_Base() : ::inet::MDTControlPacket()
{
}

NeighborSetReply_Base::NeighborSetReply_Base(const NeighborSetReply_Base& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

NeighborSetReply_Base::~NeighborSetReply_Base()
{
}

NeighborSetReply_Base& NeighborSetReply_Base::operator=(const NeighborSetReply_Base& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void NeighborSetReply_Base::copy(const NeighborSetReply_Base& other)
{
    this->sourceAddr = other.sourceAddr;
    this->predAddr = other.predAddr;
    this->destAddr = other.destAddr;
}

void NeighborSetReply_Base::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->predAddr);
    doParsimPacking(b,this->destAddr);
}

void NeighborSetReply_Base::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->predAddr);
    doParsimUnpacking(b,this->destAddr);
}

const L3Address& NeighborSetReply_Base::getSourceAddr() const
{
    return this->sourceAddr;
}

void NeighborSetReply_Base::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& NeighborSetReply_Base::getPredAddr() const
{
    return this->predAddr;
}

void NeighborSetReply_Base::setPredAddr(const L3Address& predAddr)
{
    handleChange();
    this->predAddr = predAddr;
}

const L3Address& NeighborSetReply_Base::getDestAddr() const
{
    return this->destAddr;
}

void NeighborSetReply_Base::setDestAddr(const L3Address& destAddr)
{
    handleChange();
    this->destAddr = destAddr;
}

class NeighborSetReplyDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sourceAddr,
        FIELD_predAddr,
        FIELD_destAddr,
    };
  public:
    NeighborSetReplyDescriptor();
    virtual ~NeighborSetReplyDescriptor();

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

Register_ClassDescriptor(NeighborSetReplyDescriptor)

NeighborSetReplyDescriptor::NeighborSetReplyDescriptor() : omnetpp::cClassDescriptor("inet::NeighborSetReply", "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

NeighborSetReplyDescriptor::~NeighborSetReplyDescriptor()
{
    delete[] propertyNames;
}

bool NeighborSetReplyDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<NeighborSetReply_Base *>(obj)!=nullptr;
}

const char **NeighborSetReplyDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = { "customize",  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *NeighborSetReplyDescriptor::getProperty(const char *propertyName) const
{
    if (!strcmp(propertyName, "customize")) return "true";
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int NeighborSetReplyDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 3+base->getFieldCount() : 3;
}

unsigned int NeighborSetReplyDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_sourceAddr
        0,    // FIELD_predAddr
        0,    // FIELD_destAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeFlags[field] : 0;
}

const char *NeighborSetReplyDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sourceAddr",
        "predAddr",
        "destAddr",
    };
    return (field >= 0 && field < 3) ? fieldNames[field] : nullptr;
}

int NeighborSetReplyDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "predAddr") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 2;
    return base ? base->findField(fieldName) : -1;
}

const char *NeighborSetReplyDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_predAddr
        "inet::L3Address",    // FIELD_destAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeStrings[field] : nullptr;
}

const char **NeighborSetReplyDescriptor::getFieldPropertyNames(int field) const
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

const char *NeighborSetReplyDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int NeighborSetReplyDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void NeighborSetReplyDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'NeighborSetReply_Base'", field);
    }
}

const char *NeighborSetReplyDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string NeighborSetReplyDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_predAddr: return pp->getPredAddr().str();
        case FIELD_destAddr: return pp->getDestAddr().str();
        default: return "";
    }
}

void NeighborSetReplyDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetReply_Base'", field);
    }
}

omnetpp::cValue NeighborSetReplyDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'NeighborSetReply_Base' as cValue -- field index out of range?", field);
    }
}

void NeighborSetReplyDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetReply_Base'", field);
    }
}

const char *NeighborSetReplyDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr NeighborSetReplyDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void NeighborSetReplyDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetReply_Base *pp = omnetpp::fromAnyPtr<NeighborSetReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetReply_Base'", field);
    }
}

Register_Class(NeighborSetNotification)

NeighborSetNotification::NeighborSetNotification() : ::inet::MDTControlPacket()
{
}

NeighborSetNotification::NeighborSetNotification(const NeighborSetNotification& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

NeighborSetNotification::~NeighborSetNotification()
{
    delete [] this->coords;
}

NeighborSetNotification& NeighborSetNotification::operator=(const NeighborSetNotification& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void NeighborSetNotification::copy(const NeighborSetNotification& other)
{
    delete [] this->coords;
    this->coords = (other.coords_arraysize==0) ? nullptr : new double[other.coords_arraysize];
    coords_arraysize = other.coords_arraysize;
    for (size_t i = 0; i < coords_arraysize; i++) {
        this->coords[i] = other.coords[i];
    }
    this->dim = other.dim;
    this->sourceAddr = other.sourceAddr;
    this->predAddr = other.predAddr;
    this->destAddr = other.destAddr;
}

void NeighborSetNotification::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    b->pack(coords_arraysize);
    doParsimArrayPacking(b,this->coords,coords_arraysize);
    doParsimPacking(b,this->dim);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->predAddr);
    doParsimPacking(b,this->destAddr);
}

void NeighborSetNotification::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    delete [] this->coords;
    b->unpack(coords_arraysize);
    if (coords_arraysize == 0) {
        this->coords = nullptr;
    } else {
        this->coords = new double[coords_arraysize];
        doParsimArrayUnpacking(b,this->coords,coords_arraysize);
    }
    doParsimUnpacking(b,this->dim);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->predAddr);
    doParsimUnpacking(b,this->destAddr);
}

size_t NeighborSetNotification::getCoordsArraySize() const
{
    return coords_arraysize;
}

double NeighborSetNotification::getCoords(size_t k) const
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    return this->coords[k];
}

void NeighborSetNotification::setCoordsArraySize(size_t newSize)
{
    handleChange();
    double *coords2 = (newSize==0) ? nullptr : new double[newSize];
    size_t minSize = coords_arraysize < newSize ? coords_arraysize : newSize;
    for (size_t i = 0; i < minSize; i++)
        coords2[i] = this->coords[i];
    for (size_t i = minSize; i < newSize; i++)
        coords2[i] = 0;
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void NeighborSetNotification::setCoords(size_t k, double coords)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    this->coords[k] = coords;
}

void NeighborSetNotification::insertCoords(size_t k, double coords)
{
    if (k > coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize + 1;
    double *coords2 = new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    coords2[k] = coords;
    for (i = k + 1; i < newSize; i++)
        coords2[i] = this->coords[i-1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

void NeighborSetNotification::appendCoords(double coords)
{
    insertCoords(coords_arraysize, coords);
}

void NeighborSetNotification::eraseCoords(size_t k)
{
    if (k >= coords_arraysize) throw omnetpp::cRuntimeError("Array of size %lu indexed by %lu", (unsigned long)coords_arraysize, (unsigned long)k);
    handleChange();
    size_t newSize = coords_arraysize - 1;
    double *coords2 = (newSize == 0) ? nullptr : new double[newSize];
    size_t i;
    for (i = 0; i < k; i++)
        coords2[i] = this->coords[i];
    for (i = k; i < newSize; i++)
        coords2[i] = this->coords[i+1];
    delete [] this->coords;
    this->coords = coords2;
    coords_arraysize = newSize;
}

int NeighborSetNotification::getDim() const
{
    return this->dim;
}

void NeighborSetNotification::setDim(int dim)
{
    handleChange();
    this->dim = dim;
}

const L3Address& NeighborSetNotification::getSourceAddr() const
{
    return this->sourceAddr;
}

void NeighborSetNotification::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& NeighborSetNotification::getPredAddr() const
{
    return this->predAddr;
}

void NeighborSetNotification::setPredAddr(const L3Address& predAddr)
{
    handleChange();
    this->predAddr = predAddr;
}

const L3Address& NeighborSetNotification::getDestAddr() const
{
    return this->destAddr;
}

void NeighborSetNotification::setDestAddr(const L3Address& destAddr)
{
    handleChange();
    this->destAddr = destAddr;
}

class NeighborSetNotificationDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_coords,
        FIELD_dim,
        FIELD_sourceAddr,
        FIELD_predAddr,
        FIELD_destAddr,
    };
  public:
    NeighborSetNotificationDescriptor();
    virtual ~NeighborSetNotificationDescriptor();

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

Register_ClassDescriptor(NeighborSetNotificationDescriptor)

NeighborSetNotificationDescriptor::NeighborSetNotificationDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::NeighborSetNotification)), "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

NeighborSetNotificationDescriptor::~NeighborSetNotificationDescriptor()
{
    delete[] propertyNames;
}

bool NeighborSetNotificationDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<NeighborSetNotification *>(obj)!=nullptr;
}

const char **NeighborSetNotificationDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *NeighborSetNotificationDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int NeighborSetNotificationDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 5+base->getFieldCount() : 5;
}

unsigned int NeighborSetNotificationDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISARRAY | FD_ISEDITABLE | FD_ISRESIZABLE,    // FIELD_coords
        FD_ISEDITABLE,    // FIELD_dim
        0,    // FIELD_sourceAddr
        0,    // FIELD_predAddr
        0,    // FIELD_destAddr
    };
    return (field >= 0 && field < 5) ? fieldTypeFlags[field] : 0;
}

const char *NeighborSetNotificationDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "coords",
        "dim",
        "sourceAddr",
        "predAddr",
        "destAddr",
    };
    return (field >= 0 && field < 5) ? fieldNames[field] : nullptr;
}

int NeighborSetNotificationDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "coords") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "dim") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 2;
    if (strcmp(fieldName, "predAddr") == 0) return baseIndex + 3;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 4;
    return base ? base->findField(fieldName) : -1;
}

const char *NeighborSetNotificationDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "double",    // FIELD_coords
        "int",    // FIELD_dim
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_predAddr
        "inet::L3Address",    // FIELD_destAddr
    };
    return (field >= 0 && field < 5) ? fieldTypeStrings[field] : nullptr;
}

const char **NeighborSetNotificationDescriptor::getFieldPropertyNames(int field) const
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

const char *NeighborSetNotificationDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int NeighborSetNotificationDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoordsArraySize();
        default: return 0;
    }
}

void NeighborSetNotificationDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoordsArraySize(size); break;
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'NeighborSetNotification'", field);
    }
}

const char *NeighborSetNotificationDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string NeighborSetNotificationDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return double2string(pp->getCoords(i));
        case FIELD_dim: return long2string(pp->getDim());
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_predAddr: return pp->getPredAddr().str();
        case FIELD_destAddr: return pp->getDestAddr().str();
        default: return "";
    }
}

void NeighborSetNotificationDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,string2double(value)); break;
        case FIELD_dim: pp->setDim(string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetNotification'", field);
    }
}

omnetpp::cValue NeighborSetNotificationDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: return pp->getCoords(i);
        case FIELD_dim: return pp->getDim();
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'NeighborSetNotification' as cValue -- field index out of range?", field);
    }
}

void NeighborSetNotificationDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_coords: pp->setCoords(i,value.doubleValue()); break;
        case FIELD_dim: pp->setDim(omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetNotification'", field);
    }
}

const char *NeighborSetNotificationDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr NeighborSetNotificationDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_predAddr: return omnetpp::toAnyPtr(&pp->getPredAddr()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void NeighborSetNotificationDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    NeighborSetNotification *pp = omnetpp::fromAnyPtr<NeighborSetNotification>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'NeighborSetNotification'", field);
    }
}

Register_Class(CoordsDiscoverRequest)

CoordsDiscoverRequest::CoordsDiscoverRequest() : ::inet::MDTControlPacket()
{
}

CoordsDiscoverRequest::CoordsDiscoverRequest(const CoordsDiscoverRequest& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

CoordsDiscoverRequest::~CoordsDiscoverRequest()
{
}

CoordsDiscoverRequest& CoordsDiscoverRequest::operator=(const CoordsDiscoverRequest& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void CoordsDiscoverRequest::copy(const CoordsDiscoverRequest& other)
{
    this->sourceAddr = other.sourceAddr;
    this->targetAddr = other.targetAddr;
    this->hopcount = other.hopcount;
}

void CoordsDiscoverRequest::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->targetAddr);
    doParsimPacking(b,this->hopcount);
}

void CoordsDiscoverRequest::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->targetAddr);
    doParsimUnpacking(b,this->hopcount);
}

const L3Address& CoordsDiscoverRequest::getSourceAddr() const
{
    return this->sourceAddr;
}

void CoordsDiscoverRequest::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& CoordsDiscoverRequest::getTargetAddr() const
{
    return this->targetAddr;
}

void CoordsDiscoverRequest::setTargetAddr(const L3Address& targetAddr)
{
    handleChange();
    this->targetAddr = targetAddr;
}

int CoordsDiscoverRequest::getHopcount() const
{
    return this->hopcount;
}

void CoordsDiscoverRequest::setHopcount(int hopcount)
{
    handleChange();
    this->hopcount = hopcount;
}

class CoordsDiscoverRequestDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_sourceAddr,
        FIELD_targetAddr,
        FIELD_hopcount,
    };
  public:
    CoordsDiscoverRequestDescriptor();
    virtual ~CoordsDiscoverRequestDescriptor();

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

Register_ClassDescriptor(CoordsDiscoverRequestDescriptor)

CoordsDiscoverRequestDescriptor::CoordsDiscoverRequestDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::CoordsDiscoverRequest)), "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

CoordsDiscoverRequestDescriptor::~CoordsDiscoverRequestDescriptor()
{
    delete[] propertyNames;
}

bool CoordsDiscoverRequestDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<CoordsDiscoverRequest *>(obj)!=nullptr;
}

const char **CoordsDiscoverRequestDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *CoordsDiscoverRequestDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int CoordsDiscoverRequestDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 3+base->getFieldCount() : 3;
}

unsigned int CoordsDiscoverRequestDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_sourceAddr
        0,    // FIELD_targetAddr
        FD_ISEDITABLE,    // FIELD_hopcount
    };
    return (field >= 0 && field < 3) ? fieldTypeFlags[field] : 0;
}

const char *CoordsDiscoverRequestDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "sourceAddr",
        "targetAddr",
        "hopcount",
    };
    return (field >= 0 && field < 3) ? fieldNames[field] : nullptr;
}

int CoordsDiscoverRequestDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "targetAddr") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "hopcount") == 0) return baseIndex + 2;
    return base ? base->findField(fieldName) : -1;
}

const char *CoordsDiscoverRequestDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_targetAddr
        "int",    // FIELD_hopcount
    };
    return (field >= 0 && field < 3) ? fieldTypeStrings[field] : nullptr;
}

const char **CoordsDiscoverRequestDescriptor::getFieldPropertyNames(int field) const
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

const char *CoordsDiscoverRequestDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int CoordsDiscoverRequestDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void CoordsDiscoverRequestDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'CoordsDiscoverRequest'", field);
    }
}

const char *CoordsDiscoverRequestDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string CoordsDiscoverRequestDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_targetAddr: return pp->getTargetAddr().str();
        case FIELD_hopcount: return long2string(pp->getHopcount());
        default: return "";
    }
}

void CoordsDiscoverRequestDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        case FIELD_hopcount: pp->setHopcount(string2long(value)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'CoordsDiscoverRequest'", field);
    }
}

omnetpp::cValue CoordsDiscoverRequestDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_targetAddr: return omnetpp::toAnyPtr(&pp->getTargetAddr()); break;
        case FIELD_hopcount: return pp->getHopcount();
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'CoordsDiscoverRequest' as cValue -- field index out of range?", field);
    }
}

void CoordsDiscoverRequestDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        case FIELD_hopcount: pp->setHopcount(omnetpp::checked_int_cast<int>(value.intValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'CoordsDiscoverRequest'", field);
    }
}

const char *CoordsDiscoverRequestDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr CoordsDiscoverRequestDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_targetAddr: return omnetpp::toAnyPtr(&pp->getTargetAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void CoordsDiscoverRequestDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverRequest *pp = omnetpp::fromAnyPtr<CoordsDiscoverRequest>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'CoordsDiscoverRequest'", field);
    }
}

CoordsDiscoverReply_Base::CoordsDiscoverReply_Base() : ::inet::MDTControlPacket()
{
}

CoordsDiscoverReply_Base::CoordsDiscoverReply_Base(const CoordsDiscoverReply_Base& other) : ::inet::MDTControlPacket(other)
{
    copy(other);
}

CoordsDiscoverReply_Base::~CoordsDiscoverReply_Base()
{
}

CoordsDiscoverReply_Base& CoordsDiscoverReply_Base::operator=(const CoordsDiscoverReply_Base& other)
{
    if (this == &other) return *this;
    ::inet::MDTControlPacket::operator=(other);
    copy(other);
    return *this;
}

void CoordsDiscoverReply_Base::copy(const CoordsDiscoverReply_Base& other)
{
    this->destAddr = other.destAddr;
    this->sourceAddr = other.sourceAddr;
    this->targetAddr = other.targetAddr;
}

void CoordsDiscoverReply_Base::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::inet::MDTControlPacket::parsimPack(b);
    doParsimPacking(b,this->destAddr);
    doParsimPacking(b,this->sourceAddr);
    doParsimPacking(b,this->targetAddr);
}

void CoordsDiscoverReply_Base::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::inet::MDTControlPacket::parsimUnpack(b);
    doParsimUnpacking(b,this->destAddr);
    doParsimUnpacking(b,this->sourceAddr);
    doParsimUnpacking(b,this->targetAddr);
}

const L3Address& CoordsDiscoverReply_Base::getDestAddr() const
{
    return this->destAddr;
}

void CoordsDiscoverReply_Base::setDestAddr(const L3Address& destAddr)
{
    handleChange();
    this->destAddr = destAddr;
}

const L3Address& CoordsDiscoverReply_Base::getSourceAddr() const
{
    return this->sourceAddr;
}

void CoordsDiscoverReply_Base::setSourceAddr(const L3Address& sourceAddr)
{
    handleChange();
    this->sourceAddr = sourceAddr;
}

const L3Address& CoordsDiscoverReply_Base::getTargetAddr() const
{
    return this->targetAddr;
}

void CoordsDiscoverReply_Base::setTargetAddr(const L3Address& targetAddr)
{
    handleChange();
    this->targetAddr = targetAddr;
}

class CoordsDiscoverReplyDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_destAddr,
        FIELD_sourceAddr,
        FIELD_targetAddr,
    };
  public:
    CoordsDiscoverReplyDescriptor();
    virtual ~CoordsDiscoverReplyDescriptor();

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

Register_ClassDescriptor(CoordsDiscoverReplyDescriptor)

CoordsDiscoverReplyDescriptor::CoordsDiscoverReplyDescriptor() : omnetpp::cClassDescriptor("inet::CoordsDiscoverReply", "inet::MDTControlPacket")
{
    propertyNames = nullptr;
}

CoordsDiscoverReplyDescriptor::~CoordsDiscoverReplyDescriptor()
{
    delete[] propertyNames;
}

bool CoordsDiscoverReplyDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<CoordsDiscoverReply_Base *>(obj)!=nullptr;
}

const char **CoordsDiscoverReplyDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = { "customize",  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *CoordsDiscoverReplyDescriptor::getProperty(const char *propertyName) const
{
    if (!strcmp(propertyName, "customize")) return "true";
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int CoordsDiscoverReplyDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 3+base->getFieldCount() : 3;
}

unsigned int CoordsDiscoverReplyDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        0,    // FIELD_destAddr
        0,    // FIELD_sourceAddr
        0,    // FIELD_targetAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeFlags[field] : 0;
}

const char *CoordsDiscoverReplyDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "destAddr",
        "sourceAddr",
        "targetAddr",
    };
    return (field >= 0 && field < 3) ? fieldNames[field] : nullptr;
}

int CoordsDiscoverReplyDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "sourceAddr") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "targetAddr") == 0) return baseIndex + 2;
    return base ? base->findField(fieldName) : -1;
}

const char *CoordsDiscoverReplyDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::L3Address",    // FIELD_destAddr
        "inet::L3Address",    // FIELD_sourceAddr
        "inet::L3Address",    // FIELD_targetAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeStrings[field] : nullptr;
}

const char **CoordsDiscoverReplyDescriptor::getFieldPropertyNames(int field) const
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

const char *CoordsDiscoverReplyDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int CoordsDiscoverReplyDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void CoordsDiscoverReplyDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'CoordsDiscoverReply_Base'", field);
    }
}

const char *CoordsDiscoverReplyDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        default: return nullptr;
    }
}

std::string CoordsDiscoverReplyDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_destAddr: return pp->getDestAddr().str();
        case FIELD_sourceAddr: return pp->getSourceAddr().str();
        case FIELD_targetAddr: return pp->getTargetAddr().str();
        default: return "";
    }
}

void CoordsDiscoverReplyDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'CoordsDiscoverReply_Base'", field);
    }
}

omnetpp::cValue CoordsDiscoverReplyDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_targetAddr: return omnetpp::toAnyPtr(&pp->getTargetAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'CoordsDiscoverReply_Base' as cValue -- field index out of range?", field);
    }
}

void CoordsDiscoverReplyDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'CoordsDiscoverReply_Base'", field);
    }
}

const char *CoordsDiscoverReplyDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr CoordsDiscoverReplyDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        case FIELD_sourceAddr: return omnetpp::toAnyPtr(&pp->getSourceAddr()); break;
        case FIELD_targetAddr: return omnetpp::toAnyPtr(&pp->getTargetAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void CoordsDiscoverReplyDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    CoordsDiscoverReply_Base *pp = omnetpp::fromAnyPtr<CoordsDiscoverReply_Base>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'CoordsDiscoverReply_Base'", field);
    }
}

Register_Class(PacketHolderForNSMessage)

PacketHolderForNSMessage::PacketHolderForNSMessage(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
    if (this->ownedPacket != nullptr) take(this->ownedPacket);
}

PacketHolderForNSMessage::PacketHolderForNSMessage(const PacketHolderForNSMessage& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

PacketHolderForNSMessage::~PacketHolderForNSMessage()
{
    dropAndDelete(this->ownedPacket);
}

PacketHolderForNSMessage& PacketHolderForNSMessage::operator=(const PacketHolderForNSMessage& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void PacketHolderForNSMessage::copy(const PacketHolderForNSMessage& other)
{
    dropAndDelete(this->ownedPacket);
    this->ownedPacket = other.ownedPacket;
    if (this->ownedPacket != nullptr) {
        this->ownedPacket = this->ownedPacket->dup();
        take(this->ownedPacket);
        this->ownedPacket->setName(other.ownedPacket->getName());
    }
    this->nsrqSrcAddr = other.nsrqSrcAddr;
    this->nsrqDestAddr = other.nsrqDestAddr;
}

void PacketHolderForNSMessage::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->ownedPacket);
    doParsimPacking(b,this->nsrqSrcAddr);
    doParsimPacking(b,this->nsrqDestAddr);
}

void PacketHolderForNSMessage::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->ownedPacket);
    doParsimUnpacking(b,this->nsrqSrcAddr);
    doParsimUnpacking(b,this->nsrqDestAddr);
}

const Packet * PacketHolderForNSMessage::getOwnedPacket() const
{
    return this->ownedPacket;
}

void PacketHolderForNSMessage::setOwnedPacket(Packet * ownedPacket)
{
    if (this->ownedPacket != nullptr) throw omnetpp::cRuntimeError("setOwnedPacket(): a value is already set, remove it first with removeOwnedPacket()");
    this->ownedPacket = ownedPacket;
    if (this->ownedPacket != nullptr) take(this->ownedPacket);
}

Packet * PacketHolderForNSMessage::removeOwnedPacket()
{
    Packet * retval = this->ownedPacket;
    if (retval != nullptr) drop(retval);
    this->ownedPacket = nullptr;
    return retval;
}

const L3Address& PacketHolderForNSMessage::getNsrqSrcAddr() const
{
    return this->nsrqSrcAddr;
}

void PacketHolderForNSMessage::setNsrqSrcAddr(const L3Address& nsrqSrcAddr)
{
    this->nsrqSrcAddr = nsrqSrcAddr;
}

const L3Address& PacketHolderForNSMessage::getNsrqDestAddr() const
{
    return this->nsrqDestAddr;
}

void PacketHolderForNSMessage::setNsrqDestAddr(const L3Address& nsrqDestAddr)
{
    this->nsrqDestAddr = nsrqDestAddr;
}

class PacketHolderForNSMessageDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_ownedPacket,
        FIELD_nsrqSrcAddr,
        FIELD_nsrqDestAddr,
    };
  public:
    PacketHolderForNSMessageDescriptor();
    virtual ~PacketHolderForNSMessageDescriptor();

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

Register_ClassDescriptor(PacketHolderForNSMessageDescriptor)

PacketHolderForNSMessageDescriptor::PacketHolderForNSMessageDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::PacketHolderForNSMessage)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

PacketHolderForNSMessageDescriptor::~PacketHolderForNSMessageDescriptor()
{
    delete[] propertyNames;
}

bool PacketHolderForNSMessageDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<PacketHolderForNSMessage *>(obj)!=nullptr;
}

const char **PacketHolderForNSMessageDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *PacketHolderForNSMessageDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int PacketHolderForNSMessageDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 3+base->getFieldCount() : 3;
}

unsigned int PacketHolderForNSMessageDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISCOMPOUND | FD_ISPOINTER | FD_ISCOBJECT | FD_ISCOWNEDOBJECT | FD_ISREPLACEABLE,    // FIELD_ownedPacket
        0,    // FIELD_nsrqSrcAddr
        0,    // FIELD_nsrqDestAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeFlags[field] : 0;
}

const char *PacketHolderForNSMessageDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "ownedPacket",
        "nsrqSrcAddr",
        "nsrqDestAddr",
    };
    return (field >= 0 && field < 3) ? fieldNames[field] : nullptr;
}

int PacketHolderForNSMessageDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "ownedPacket") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "nsrqSrcAddr") == 0) return baseIndex + 1;
    if (strcmp(fieldName, "nsrqDestAddr") == 0) return baseIndex + 2;
    return base ? base->findField(fieldName) : -1;
}

const char *PacketHolderForNSMessageDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::Packet",    // FIELD_ownedPacket
        "inet::L3Address",    // FIELD_nsrqSrcAddr
        "inet::L3Address",    // FIELD_nsrqDestAddr
    };
    return (field >= 0 && field < 3) ? fieldTypeStrings[field] : nullptr;
}

const char **PacketHolderForNSMessageDescriptor::getFieldPropertyNames(int field) const
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

const char *PacketHolderForNSMessageDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int PacketHolderForNSMessageDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void PacketHolderForNSMessageDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'PacketHolderForNSMessage'", field);
    }
}

const char *PacketHolderForNSMessageDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: { const Packet * value = pp->getOwnedPacket(); return omnetpp::opp_typename(typeid(*value)); }
        default: return nullptr;
    }
}

std::string PacketHolderForNSMessageDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: { auto obj = pp->getOwnedPacket(); return obj == nullptr ? "" : obj->str(); }
        case FIELD_nsrqSrcAddr: return pp->getNsrqSrcAddr().str();
        case FIELD_nsrqDestAddr: return pp->getNsrqDestAddr().str();
        default: return "";
    }
}

void PacketHolderForNSMessageDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PacketHolderForNSMessage'", field);
    }
}

omnetpp::cValue PacketHolderForNSMessageDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::toAnyPtr(pp->getOwnedPacket()); break;
        case FIELD_nsrqSrcAddr: return omnetpp::toAnyPtr(&pp->getNsrqSrcAddr()); break;
        case FIELD_nsrqDestAddr: return omnetpp::toAnyPtr(&pp->getNsrqDestAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'PacketHolderForNSMessage' as cValue -- field index out of range?", field);
    }
}

void PacketHolderForNSMessageDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: pp->setOwnedPacket(omnetpp::fromAnyPtr<Packet>(value.pointerValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PacketHolderForNSMessage'", field);
    }
}

const char *PacketHolderForNSMessageDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr PacketHolderForNSMessageDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::toAnyPtr(pp->getOwnedPacket()); break;
        case FIELD_nsrqSrcAddr: return omnetpp::toAnyPtr(&pp->getNsrqSrcAddr()); break;
        case FIELD_nsrqDestAddr: return omnetpp::toAnyPtr(&pp->getNsrqDestAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void PacketHolderForNSMessageDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderForNSMessage *pp = omnetpp::fromAnyPtr<PacketHolderForNSMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: pp->setOwnedPacket(omnetpp::fromAnyPtr<Packet>(ptr)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PacketHolderForNSMessage'", field);
    }
}

Register_Class(PacketHolderMessage)

PacketHolderMessage::PacketHolderMessage(const char *name, short kind) : ::omnetpp::cMessage(name, kind)
{
    if (this->ownedPacket != nullptr) take(this->ownedPacket);
}

PacketHolderMessage::PacketHolderMessage(const PacketHolderMessage& other) : ::omnetpp::cMessage(other)
{
    copy(other);
}

PacketHolderMessage::~PacketHolderMessage()
{
    dropAndDelete(this->ownedPacket);
}

PacketHolderMessage& PacketHolderMessage::operator=(const PacketHolderMessage& other)
{
    if (this == &other) return *this;
    ::omnetpp::cMessage::operator=(other);
    copy(other);
    return *this;
}

void PacketHolderMessage::copy(const PacketHolderMessage& other)
{
    dropAndDelete(this->ownedPacket);
    this->ownedPacket = other.ownedPacket;
    if (this->ownedPacket != nullptr) {
        this->ownedPacket = this->ownedPacket->dup();
        take(this->ownedPacket);
        this->ownedPacket->setName(other.ownedPacket->getName());
    }
    this->destAddr = other.destAddr;
}

void PacketHolderMessage::parsimPack(omnetpp::cCommBuffer *b) const
{
    ::omnetpp::cMessage::parsimPack(b);
    doParsimPacking(b,this->ownedPacket);
    doParsimPacking(b,this->destAddr);
}

void PacketHolderMessage::parsimUnpack(omnetpp::cCommBuffer *b)
{
    ::omnetpp::cMessage::parsimUnpack(b);
    doParsimUnpacking(b,this->ownedPacket);
    doParsimUnpacking(b,this->destAddr);
}

const Packet * PacketHolderMessage::getOwnedPacket() const
{
    return this->ownedPacket;
}

void PacketHolderMessage::setOwnedPacket(Packet * ownedPacket)
{
    if (this->ownedPacket != nullptr) throw omnetpp::cRuntimeError("setOwnedPacket(): a value is already set, remove it first with removeOwnedPacket()");
    this->ownedPacket = ownedPacket;
    if (this->ownedPacket != nullptr) take(this->ownedPacket);
}

Packet * PacketHolderMessage::removeOwnedPacket()
{
    Packet * retval = this->ownedPacket;
    if (retval != nullptr) drop(retval);
    this->ownedPacket = nullptr;
    return retval;
}

const L3Address& PacketHolderMessage::getDestAddr() const
{
    return this->destAddr;
}

void PacketHolderMessage::setDestAddr(const L3Address& destAddr)
{
    this->destAddr = destAddr;
}

class PacketHolderMessageDescriptor : public omnetpp::cClassDescriptor
{
  private:
    mutable const char **propertyNames;
    enum FieldConstants {
        FIELD_ownedPacket,
        FIELD_destAddr,
    };
  public:
    PacketHolderMessageDescriptor();
    virtual ~PacketHolderMessageDescriptor();

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

Register_ClassDescriptor(PacketHolderMessageDescriptor)

PacketHolderMessageDescriptor::PacketHolderMessageDescriptor() : omnetpp::cClassDescriptor(omnetpp::opp_typename(typeid(inet::PacketHolderMessage)), "omnetpp::cMessage")
{
    propertyNames = nullptr;
}

PacketHolderMessageDescriptor::~PacketHolderMessageDescriptor()
{
    delete[] propertyNames;
}

bool PacketHolderMessageDescriptor::doesSupport(omnetpp::cObject *obj) const
{
    return dynamic_cast<PacketHolderMessage *>(obj)!=nullptr;
}

const char **PacketHolderMessageDescriptor::getPropertyNames() const
{
    if (!propertyNames) {
        static const char *names[] = {  nullptr };
        omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
        const char **baseNames = base ? base->getPropertyNames() : nullptr;
        propertyNames = mergeLists(baseNames, names);
    }
    return propertyNames;
}

const char *PacketHolderMessageDescriptor::getProperty(const char *propertyName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? base->getProperty(propertyName) : nullptr;
}

int PacketHolderMessageDescriptor::getFieldCount() const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    return base ? 2+base->getFieldCount() : 2;
}

unsigned int PacketHolderMessageDescriptor::getFieldTypeFlags(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeFlags(field);
        field -= base->getFieldCount();
    }
    static unsigned int fieldTypeFlags[] = {
        FD_ISCOMPOUND | FD_ISPOINTER | FD_ISCOBJECT | FD_ISCOWNEDOBJECT | FD_ISREPLACEABLE,    // FIELD_ownedPacket
        0,    // FIELD_destAddr
    };
    return (field >= 0 && field < 2) ? fieldTypeFlags[field] : 0;
}

const char *PacketHolderMessageDescriptor::getFieldName(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldName(field);
        field -= base->getFieldCount();
    }
    static const char *fieldNames[] = {
        "ownedPacket",
        "destAddr",
    };
    return (field >= 0 && field < 2) ? fieldNames[field] : nullptr;
}

int PacketHolderMessageDescriptor::findField(const char *fieldName) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    int baseIndex = base ? base->getFieldCount() : 0;
    if (strcmp(fieldName, "ownedPacket") == 0) return baseIndex + 0;
    if (strcmp(fieldName, "destAddr") == 0) return baseIndex + 1;
    return base ? base->findField(fieldName) : -1;
}

const char *PacketHolderMessageDescriptor::getFieldTypeString(int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldTypeString(field);
        field -= base->getFieldCount();
    }
    static const char *fieldTypeStrings[] = {
        "inet::Packet",    // FIELD_ownedPacket
        "inet::L3Address",    // FIELD_destAddr
    };
    return (field >= 0 && field < 2) ? fieldTypeStrings[field] : nullptr;
}

const char **PacketHolderMessageDescriptor::getFieldPropertyNames(int field) const
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

const char *PacketHolderMessageDescriptor::getFieldProperty(int field, const char *propertyName) const
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

int PacketHolderMessageDescriptor::getFieldArraySize(omnetpp::any_ptr object, int field) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldArraySize(object, field);
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        default: return 0;
    }
}

void PacketHolderMessageDescriptor::setFieldArraySize(omnetpp::any_ptr object, int field, int size) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldArraySize(object, field, size);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set array size of field %d of class 'PacketHolderMessage'", field);
    }
}

const char *PacketHolderMessageDescriptor::getFieldDynamicTypeString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldDynamicTypeString(object,field,i);
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: { const Packet * value = pp->getOwnedPacket(); return omnetpp::opp_typename(typeid(*value)); }
        default: return nullptr;
    }
}

std::string PacketHolderMessageDescriptor::getFieldValueAsString(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValueAsString(object,field,i);
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: { auto obj = pp->getOwnedPacket(); return obj == nullptr ? "" : obj->str(); }
        case FIELD_destAddr: return pp->getDestAddr().str();
        default: return "";
    }
}

void PacketHolderMessageDescriptor::setFieldValueAsString(omnetpp::any_ptr object, int field, int i, const char *value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValueAsString(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PacketHolderMessage'", field);
    }
}

omnetpp::cValue PacketHolderMessageDescriptor::getFieldValue(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldValue(object,field,i);
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::toAnyPtr(pp->getOwnedPacket()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: throw omnetpp::cRuntimeError("Cannot return field %d of class 'PacketHolderMessage' as cValue -- field index out of range?", field);
    }
}

void PacketHolderMessageDescriptor::setFieldValue(omnetpp::any_ptr object, int field, int i, const omnetpp::cValue& value) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldValue(object, field, i, value);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: pp->setOwnedPacket(omnetpp::fromAnyPtr<Packet>(value.pointerValue())); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PacketHolderMessage'", field);
    }
}

const char *PacketHolderMessageDescriptor::getFieldStructName(int field) const
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

omnetpp::any_ptr PacketHolderMessageDescriptor::getFieldStructValuePointer(omnetpp::any_ptr object, int field, int i) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount())
            return base->getFieldStructValuePointer(object, field, i);
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: return omnetpp::toAnyPtr(pp->getOwnedPacket()); break;
        case FIELD_destAddr: return omnetpp::toAnyPtr(&pp->getDestAddr()); break;
        default: return omnetpp::any_ptr(nullptr);
    }
}

void PacketHolderMessageDescriptor::setFieldStructValuePointer(omnetpp::any_ptr object, int field, int i, omnetpp::any_ptr ptr) const
{
    omnetpp::cClassDescriptor *base = getBaseClassDescriptor();
    if (base) {
        if (field < base->getFieldCount()){
            base->setFieldStructValuePointer(object, field, i, ptr);
            return;
        }
        field -= base->getFieldCount();
    }
    PacketHolderMessage *pp = omnetpp::fromAnyPtr<PacketHolderMessage>(object); (void)pp;
    switch (field) {
        case FIELD_ownedPacket: pp->setOwnedPacket(omnetpp::fromAnyPtr<Packet>(ptr)); break;
        default: throw omnetpp::cRuntimeError("Cannot set field %d of class 'PacketHolderMessage'", field);
    }
}

}  // namespace inet

namespace omnetpp {

}  // namespace omnetpp

