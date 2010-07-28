#include <src/serialabstract.h>
#include <src/serialimpl.h>
#include <src/archive.h>
#include <src/printstream.h>
#include <src/exception.h>
#include <src/timer.h>
#include <src/archive.h>
#include <sstream>

using namespace smartptr;
using namespace std;

#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

map<string, SerialRuntime::create_function>* SerialRuntime::classlist_ = NULL;

Serializable::Serializable()
    : Countable()
{
    id_ = (serid_t) this;
}

Serializable::Serializable(const ArchivePtr& arch)
    : Countable()
{ 
    serial_load(id);
    if (arch.get() == NULL)
    {
        cerr << "Serializable class received null node in constructor" << endl;
        abort();
    }

    arch->registry()->add(id_, this);
}

void
Serializable::serialize(const ArchivePtr& arch) const
{
    serial_save(id);
}

serid_t
Serializable::id() const
{
    return id_;
}

Serializable::~Serializable() //delete from the registry
{
}

const SerialRuntime* 
Serializable::runtime_info() const
{
    return rtinfo_;
}

void 
Serializable::xmlCommit(refsmartstr filename, const ArchivePtr& arch)
{
    arch->xml()->toFile("gigide.arch");
}

#if 0
SerializablePtr
Serializable::getObject(const ArchivePtr& arch)
{
    serid_t id; arch->xml()->getValue(id);
    return getObject(arch, id);
}

SerializablePtr
Serializable::getObject(const ArchivePtr& arch, serid_t id)
{
    gigtimer::Timer::start("build");
    ConstSerializablePtr obj = arch->registry()->get(id);
    if (obj != NULL)
    {
        SerializablePtr ptr = boost::const_pointer_cast<Serializable, const Serializable>(obj);
        gigtimer::Timer::stop("build");
        return ptr;
    }

    gigtimer::Timer::start("get branch");
    std::string nodeid = stream_printf("node" serid_format_str, id);
    ArchivePtr objnode = arch->getObjectBranch()->getBranch(nodeid);
    if (objnode.get() == NULL)
    {
        cerr << "Cannot find object " << id << " in object registry" << endl;
        cerr << "Node text: " << arch->xml()->getText() << endl;
        abort();
    }
    gigtimer::Timer::stop("get branch");
    gigtimer::Timer::stop("build");

    return SerialRuntime::getObject(objnode);
}
#endif

SerialRuntime::SerialRuntime(refsmartstr name, create_function fxn)
{
    if (classlist_ == NULL)
        classlist_ = new map<string, create_function>();

    cout << "Registering " << name << endl;
    classlist_->insert(pair<string,create_function>(name,fxn));
}

boost::intrusive_ptr<Serializable>
SerialRuntime::getObject(const ArchivePtr& arch, refsmartstr classname)
{
    map<string, create_function>::const_iterator it(classlist_->find(classname));
    if (it == classlist_->end())
    {
        cerr << stream_printf("%s is not a valid serializable class name", classname.c_str()) << endl;
        abort();
    }
    return (it->second)(arch);
}

void
SerialRuntime::free()
{
    delete classlist_;
}

