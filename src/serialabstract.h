#ifndef _gigide_serialabstract_h_
#define _gigide_serialabstract_h_

#include <string>
#include <src/ref.h>
#include <src/set.h>
#include <iostream>
#include <map>

//#define SerialClassnameGet(node) node->getAttribute("classtype");

#define heisenbug std::cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << std::endl

#define SerializablePtr boost::intrusive_ptr<smartptr::Serializable>  
#define ConstSerializablePtr boost::intrusive_ptr<const smartptr::Serializable>  

typedef unsigned long serid_t;
#define serid_format_str "%ld"

namespace smartptr {

template <bool t>
struct AssertNotHere {
  enum { N = 1 - 2 * int(!t) };
      // 1 if t is true, -1 if t is false.
  static char A[N];
};

#define smartstr std::string
#define refsmartstr const smartstr&


#define ArchivePtr boost::intrusive_ptr<smartptr::Archive>  
class Archive;
class SerialRuntime;
class Serializable : public Countable {
    
    protected:
        SerialRuntime* rtinfo_;

    private:
        serid_t id_;

    public:
        serid_t id() const;

        virtual void serialize(const ArchivePtr& archive) const;

        static void xmlCommit(refsmartstr filename, const ArchivePtr& archive);

        const SerialRuntime* runtime_info() const;

    protected:
        Serializable();

        Serializable(const ArchivePtr& archive);

        ~Serializable();


};

template <class T>
ArchivePtr
serial_call_load(
    const ArchivePtr& arch,
    T& val,
    const std::string& tagname,
    const std::string& nspace = std::string("")
);

template <class T>
void
serial_call_load(
    const ArchivePtr& arch,
    T& val
);

template <class T>
ArchivePtr
serial_call_save(
    const ArchivePtr& arch,
    const T& val,
    const std::string& tagname,
    const std::string& nspace = std::string("")
);

template <class T>
void
serial_call_save(
    const ArchivePtr& arch,
    const T& val
);

} //end namespace

#endif
