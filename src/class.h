#ifndef yeti_class_h
#define yeti_class_h

/** Macros for making class operations easier when using smart pointers */

#include <sstream>
#include <vector>
#include <cmath>
#include <src/smartptr/src/serialize.h>
#include <src/smartptr/src/set.h>

#define make(name, cls, ...) cls##Ptr name (new cls(__VA_ARGS__))
#define polymake(name, super, sub, ...) super##Ptr name (new sub(__VA_ARGS__))
#define asgn(name, cls, ...) name = new cls(__VA_ARGS__)
#define make_ptr(name, cls, ...) cls * name (new cls(__VA_ARGS__))

#define ptrcast(cls, x) reinterpret_cast<cls*>(x)

#define LOG_ZERO -100

#define heisenbug std::cout << __FILE__ << " " << __LINE__ << std::endl

#define findbug(x) std::cout << __FILE__ << " " << x << std::endl

#define DESTRUCTOR_PRINT 0

#define YETI_SANITY_CHECK 1

#define NMALLOC_BLOCKS_DEBUG 1000000

namespace yeti {

typedef unsigned short int usi;
typedef unsigned long int uli;

class DefaultCTFailure {};

template <bool t, typename msg = DefaultCTFailure>
struct FailCompileTime {
    enum { N = 1 - 2 * int(t) };
    static char A[N];
};

/** if t true, send true to fail compile check */
template <class T, bool t> class AssertPointer       { public: AssertPointer(){ FailCompileTime<t> check; } };
/** if t true, send false to fail compile check */
template <class T, bool t> class AssertPointer<T*,t> { public: AssertPointer(){ FailCompileTime<!t> check; } };

/** if t true, send true to fail compile check */
template <class T, bool t> class AssertSmartPtr       { public: AssertSmartPtr(){ FailCompileTime<t> check; } };
/** if t true, send false to fail compile check */
template <class T, bool t> class AssertSmartPtr< boost::intrusive_ptr<T>, t> { public: AssertSmartPtr(){ FailCompileTime<!t> check; } };

/** if t true, send true to fail compile check */
template <class T, bool t> class AssertVector       { public: AssertVector(){ FailCompileTime<t> check; } };
/** if t true, send false to fail compile check */
template <class T, bool t> class AssertVector< boost::intrusive_ptr<T>, t> { public: AssertVector(){ FailCompileTime<!t> check; } };


template <class T>
const char*
yeti_classname(T* x)
{
    return typeid(x).name();
}

#define class_status(x) std::cout << #x << " " << yeti_classname(this) << " " << (void*) this << " at " << __FILE__ << " " << __LINE__ << endl;

template <class T>
class ConstClass
{
    public:
        typedef const T const_type;
        typedef T non_const_type;
};

template <class T>
class ConstClass<const T>
{
    public:
        typedef const T const_type;
        typedef T non_const_type;
};

template <typename t1, typename t2>
struct TypeCheck {

    const static bool mismatch = true;
    const static bool match = false;

};

template <typename t>
struct TypeCheck<t,t> {

    const static bool mismatch = false;
    const static bool match = true;

};

/** if t true, send true to fail compile check */
template <bool T=true> class TypeNotImplemented  { public: TypeNotImplemented(){ FailCompileTime<T> check; } };

struct TemplateInfo {

    /**
        Integers should be cause up to doubles or floats.
    */
    typedef enum { integer_type = 0, double_type = 1, float_type = 2 } type_t;

    static const usi ntypes = 3;
};


template <class T>
class TypeInfo {

    public:
        TypeNotImplemented<> typecheck;
};

template <>
class TypeInfo<double> {

    public:
        static const char* name(){return "double";}
        static TemplateInfo::type_t type(){return TemplateInfo::double_type;}
        static const char* printf_str;

};

template <>
class TypeInfo<float> {


    public:
        static const char* name(){return "float";}
        static TemplateInfo::type_t type(){return TemplateInfo::float_type;}
        static const char* printf_str;
};

template <>
class TypeInfo<int> {

    public:
        static const char* name(){return "int";}
        static TemplateInfo::type_t type(){return TemplateInfo::integer_type;}
        static const char* printf_str;
};

//Default string output for unit tests
template <class T> class ClassOutput {

    public:
        static std::string
        str(T t)
        {
            //AssertPointer<T,false> should_not_be_pointer;
            std::stringstream sstr;
            sstr << t;
            return sstr.str();
        }

        static std::string
        str(size_t n, T t)
        {
            AssertPointer<T,true> should_be_pointer;
            std::stringstream sstr;
            sstr << "{";
            for (size_t i=0; i < n; ++i)
                sstr << " " << t[i];
            sstr << " }";
            return sstr.str();
        }

};

template <> class ClassOutput<char *> {

    public:
        static std::string
        str(char* ptr)
        {
            void* vptr = (void *) ptr;
            std::stringstream sstr;
            sstr << vptr;
            return sstr.str();
        }

};

template <class T> class ClassOutput< boost::intrusive_ptr<T> > {

    public:
        static std::string str(const boost::intrusive_ptr<T>& test)
        {
            std::stringstream sstr;
            test->print(sstr);
            return sstr.str();
        }

};

template <class T> class ClassOutput< std::vector<T> > {

    public:
        static std::string str(const std::vector<T>& test)
        {
            std::stringstream sstr;
            sstr << "{";
            typename std::vector<T>::const_iterator itest(test.begin());
            for ( ; itest != test.end(); ++itest)
            {
                sstr << " " << ClassOutput<T>::str(*itest);
            }
            sstr << " }";
            return sstr.str();
        }
};

template <class T, class U> class ClassOutput< std::map<T, U> > {

    public:
        static std::string str(const std::map<T, U>& test)
        {
            std::stringstream sstr;
            sstr << "map {" << std::endl;
            typename std::map<T, U>::const_iterator itest(test.begin());
            for ( ; itest != test.end(); ++itest)
            {
                sstr << ClassOutput<T>::str(itest->first)
                     << "->"
                     << ClassOutput<U>::str(itest->second)
                     << std::endl;
            }
            sstr << "   }";
            return sstr.str();
        }

};

#define FLOAT_CUTOFF 1e-10
#define DOUBLE_CUTOFF 1e-12

template <class T>
inline bool arrays_equal(size_t n, const T* test, const T* right);

template <class T> class TestEquals {

    public:
        static bool equals(size_t n, const T* test, const T* right)
        {
            return arrays_equal<T>(n, test, right);
        }

        static bool equals(T test, T right)
        {
            //AssertPointer<T,false> should_not_be_pointer;
            AssertSmartPtr<T,false> should_not_be_smartptr;
            AssertVector<T,false> should_not_be_vector;
            return test == right;
        }

};

template <> class TestEquals<float> {

    public:
        static float cutoff;

        static bool equals(float test, float right)
        {
            return fabs(test - right) < cutoff;
        }

        static bool equals(size_t n, const float* test, const float* right)
        {
            return arrays_equal<float>(n, test, right);
        }


};



template <> class TestEquals<double> {

    public:
        static double cutoff;

        static bool equals(double test, double right)
        {
            return fabs(test - right) < cutoff;
        }


};

template <class T> class TestEquals< std::vector<T> > {

    public:
        static bool equals(const std::vector<T>& test, const std::vector<T>& right)
        {
            AssertVector< std::vector<T>, false> fail;

            if (test.size() != right.size())
                return false;

            typename std::vector<T>::const_iterator itest(test.begin());
            typename std::vector<T>::const_iterator iright(right.begin());
            for ( ; itest != test.end(); ++itest, ++iright)
            {
                if ( !(TestEquals<T>::equals(*itest, *right)) )
                    return false;
            }
            return true;
        }
};

template <class T> class TestEquals< boost::intrusive_ptr<T> > {

    public:
        static bool equals(const boost::intrusive_ptr<T>& test, const boost::intrusive_ptr<T>& right)
        {
            return right->equals(test);
        }

};

template <class T>
inline bool arrays_equal(size_t n, const T* test, const T* right)
{
    for (size_t i=0; i < n; ++i)
    {
        if ( !(TestEquals<T>::equals(test[i], right[i])) )
            return false;
    }
    return true;
}

template <class T>
std::string
output_str(const T& obj)
{
    return ClassOutput<T>::str(obj);
}


}

#endif


