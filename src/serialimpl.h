#ifndef gigide_serialimpl_h
#define gigide_serialimpl_h

#include <src/archive.h>
//#include <src/gigmatrix.h>

#define serial_save(x)        smartptr::serial_call_save(arch, x##_, #x) //arch->serialize(Archive::Write, x##_, #x)
#define serial_save_enum(x) unsigned int tmp = (unsigned int) x##_; arch->serialize_data(Archive::Write, tmp, #x);
#define serial_save_ns(x, ns) smartptr::serial_call_save(arch, x##_, #x, #ns)
#define serial_save_map(x, t1, t2) arch->serialize<t1, t2>(Archive::Write, x##_, #x);

#define macro_heisenbug(x) std::cout << "heisenbug" << #x << " at " << __FILE__ << " and " << __LINE__ << std::endl;

#define serial_load(x)        smartptr::serial_call_load(arch, x##_, #x)
#define serial_load_enum(x) unsigned int tmp; arch->serialize_data(Archive::Read, tmp, #x); arch->enum_cast(tmp, x##_);
#define serial_load_ns(x, ns) smartptr::serial_call_load(arch, x##_, #x, #ns)
#define serial_load_map(x, t1, t2) arch->serialize<t1, t2>(Archive::Read, x##_, #x);

//this has a compile time assertion to make sure the type is correct
#define SetRuntime(x) rtinfo_ = &x##_sc; x* const compile_time_check = this; 

//classname must come before SerialClass<x> to ensure string has been constructed prior
//to static initialization
#define SerialDeclare(x) \
template<> std::string smartptr::SerialClass<x>::classname_ = ""; \
template<> smartptr::SerialRuntime::create_function smartptr::SerialClass<x>::fxn_ = smartptr::build<x>; \
smartptr::SerialClass<x> x##_sc(#x); 

namespace smartptr {

template <class T>
class SerialDecide {

    public:
        static ArchivePtr 
        serialize(
            const ArchivePtr& arch,
            Archive::save_t flag,
            const T& val,
            const std::string& tagname,
            const std::string& nspace
        )
        {
            return arch->serialize(flag, val, tagname, nspace);
        }

        static void
        serialize(
            const ArchivePtr& arch,
            Archive::save_t flag,
            const T& val
        )
        {
            arch->serialize(flag, val);
        }

        static ArchivePtr 
        serialize(
            const ArchivePtr& arch,
            Archive::load_t flag,
            T& val,
            const std::string& tagname,
            const std::string& nspace
        )
        {
            return arch->serialize(flag, val, tagname, nspace);
        }

        static void
        serialize(
            const ArchivePtr& arch,
            Archive::load_t flag,
            T& val
        )
        {
            arch->serialize(flag, val);
        }
};

#define SerialDecideInstance(x) \
template <> class SerialDecide<x> { \
    public: \
        static ArchivePtr  \
        serialize( \
            const ArchivePtr& arch, \
            Archive::save_t flag, \
            const x & val, \
            const std::string& tagname, \
            const std::string& nspace \
        ) \
        { \
            return arch->serialize_data(flag, val, tagname, nspace); \
        } \
        static void  \
        serialize( \
            const ArchivePtr& arch, \
            Archive::save_t flag, \
            const x & val \
        ) \
        { \
            arch->serialize_data(flag, val); \
        } \
        static ArchivePtr  \
        serialize( \
            const ArchivePtr& arch, \
            Archive::load_t flag, \
            x& val, \
            const std::string& tagname, \
            const std::string& nspace \
        ) \
        { \
            return arch->serialize_data(flag, val, tagname, nspace); \
        } \
        static void  \
        serialize( \
            const ArchivePtr& arch, \
            Archive::load_t flag, \
            x & val \
        ) \
        { \
            arch->serialize_data(flag, val); \
        } \
}

SerialDecideInstance(int);
SerialDecideInstance(double);
SerialDecideInstance(unsigned int);
SerialDecideInstance(unsigned long);
SerialDecideInstance(bool);
SerialDecideInstance(std::string);

#define SerialDecideSubptr(x) \
template <> class SerialDecide<x> { \
    public: \
        typedef x::element_type subtype; \
        static ArchivePtr  \
        serialize( \
            const ArchivePtr& arch, \
            Archive::save_t flag, \
            const x & val, \
            const std::string& tagname, \
            const std::string& nspace \
        ) \
        { \
            boost::intrusive_ptr<subtype> tmp = val.get(); \
            return arch->serialize<subtype>(flag, tmp, tagname, nspace); \
        } \
        static ArchivePtr  \
        serialize( \
            const ArchivePtr& arch, \
            Archive::load_t flag, \
            x& val, \
            const std::string& tagname, \
            const std::string& nspace \
        ) \
        { \
            boost::intrusive_ptr<subtype> tmp = val.get(); \
            ArchivePtr node = arch->serialize<subtype>(flag, tmp, tagname, nspace); \
            val = tmp.get(); \
            return node; \ 
        } \
}


template <class T>
ArchivePtr
serial_call_save(
    const ArchivePtr& arch,
    const T& val,
    const std::string& tagname,
    const std::string& nspace = std::string("")
)
{
    return SerialDecide<T>::serialize(arch, Archive::Write, val, tagname, nspace);
}

template <class T>
void
serial_call_save(
    const ArchivePtr& arch,
    const T& val
)
{
    SerialDecide<T>::serialize(arch, Archive::Write, val);
}


template <class T>
ArchivePtr
serial_call_load(
    const ArchivePtr& arch,
    T& val,
    const std::string& tagname,
    const std::string& nspace = std::string("")
)
{
    return SerialDecide<T>::serialize(arch, Archive::Read, val, tagname, nspace);
}

template <class T>
void
serial_call_load(
    const ArchivePtr& arch,
    T& val
)
{
   SerialDecide<T>::serialize(arch, Archive::Read, val);
}

class SerialRuntime {

    public:
        typedef Serializable*(*create_function)(const ArchivePtr&);

        SerialRuntime(refsmartstr str, create_function fxnptr);

        static SerializablePtr getObject(const ArchivePtr& arch, refsmartstr classname);

        static void free();

        virtual std::string classname() const = 0;

    private:
        static std::map<std::string, create_function> *classlist_;



};


template <class T>
Serializable*
build(const ArchivePtr& arch)
{
    return new T(arch);
}

template <class T> 
class SerialClass : public SerialRuntime {

    private:
        static std::string classname_;
        static create_function fxn_;

    public:

        SerialClass(refsmartstr name) :
            SerialRuntime(name, fxn_)
        {
            classname_ = name;
        }

        std::string
        classname() const
        {
            return classname_;
        }

};

} //end namespace

#endif
