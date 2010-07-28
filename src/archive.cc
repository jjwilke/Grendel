#include <src/archive.h>
#include <src/serialimpl.h>
#include <src/serialabstract.h>
#include <src/printstream.h>
#include <iostream>

#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

using namespace smartptr;
using namespace pyxml;
using namespace std;


SerialMap::SerialMap()
{
}

bool
SerialMap::has(serid_t id) const
{
    registry_map::const_iterator it(registry_.find(id));
    return it != registry_.end();
}

SerializablePtr
SerialMap::get(serid_t id) const
{
    registry_map::const_iterator it(registry_.find(id));
    if (it == registry_.end())
        return NULL;

    return boost::const_pointer_cast<Serializable>(it->second);
}

void
SerialMap::add(serid_t id, const ConstSerializablePtr& obj)
{
    registry_[id] = obj;
}

Archive::Archive()
    : topnode_(new pyxml::PyXMLDomParser()), 
      node_(topnode_),
      registry_(new SerialMap)
{
      objnode_ = new Archive(topnode_, topnode_->getNode("objects"), registry_);
}

Archive::Archive(refxmlstr filename)
    : topnode_(new pyxml::PyXMLDomParser(filename)),
      node_(topnode_),
      registry_(new SerialMap)
{
      objnode_ = new Archive(topnode_, topnode_->getNode("objects"), registry_);
}

Archive::Archive(
    const PyXMLDomParserPtr& topnode,
    const PyXMLDomParserPtr& node,
    const SerialMapPtr& registry,
    const ArchivePtr& objnode
)
    : topnode_(topnode), node_(node), registry_(registry), objnode_(objnode)
{
}

Archive::Archive(
    const PyXMLDomParserPtr& topnode,
    const PyXMLDomParserPtr& node,
    const SerialMapPtr& registry
)
    : topnode_(topnode), node_(node), registry_(registry)
{
    //I am the objectnode
    Archive* ptr = const_cast<Archive*>(this);
    objnode_ = ptr;
}

Archive::~Archive()
{
}

ArchivePtr
Archive::getTop() const
{
    ArchivePtr node(new Archive(topnode_, topnode_, registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::fetchBranch(refxmlstr tagname, refxmlstr nspace) const
{
    XMLParser node = node_->fetchNode(tagname, nspace);
    if (node.get() == NULL)
        return NULL;

    ArchivePtr arch(new Archive(topnode_, node, registry_, objnode_));
    return arch;
}

ArchivePtr 
Archive::getBranch(pyxml::refxmlstr tagname, pyxml::refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->getNode(tagname, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::makeBranch(pyxml::refxmlstr tagname, pyxml::refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->createNode(tagname, nspace), registry_, objnode_));
    return node;
}

void
Archive::fetchBranches(std::vector<ArchivePtr>& nodes, pyxml::refxmlstr tagname, pyxml::refxmlstr nspace) const
{
    _verify();
    vector<XMLParser> xmls;
    node_->getElementsByTagName(xmls, tagname, nspace);

    for (int i=0; i < xmls.size(); ++i)
    {
        ArchivePtr node(new Archive(topnode_, xmls[i], registry_, objnode_));
        nodes.push_back(node);
    }
}

void 
Archive::getAttribute(int& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::getAttribute(long& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::getAttribute(size_t& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::getAttribute(double& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::getAttribute(unsigned int& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::getAttribute(bool& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::getAttribute(xmlstr& val, refxmlstr attrname) const
{
    node_->getAttribute(val, attrname);
}

void 
Archive::setAttribute(int val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void 
Archive::setAttribute(long val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void 
Archive::setAttribute(size_t val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void 
Archive::setAttribute(double val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void 
Archive::setAttribute(unsigned int val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void 
Archive::setAttribute(bool val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void 
Archive::setAttribute(refxmlstr val, refxmlstr attrname)
{
    node_->setAttribute(val, attrname);
}

void
Archive::getValue(int& val)
{
    node_->getValue(val);
}

void
Archive::getValue(long& val)
{
    node_->getValue(val);
}

void
Archive::getValue(unsigned int& val)
{
    node_->getValue(val);
}

void
Archive::getValue(unsigned long& val)
{
    node_->getValue(val);
}

void
Archive::getValue(double& val)
{
    node_->getValue(val);
}

void
Archive::getValue(bool& val)
{
    node_->getValue(val);
}

void
Archive::getValue(xmlstr& val)
{
    node_->getValue(val);
}

ArchivePtr 
Archive::getValue(int& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::getValue(long& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::getValue(unsigned int& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::getValue(unsigned long& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::getValue(double& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::getValue(bool& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::getValue(xmlstr& val, refxmlstr name, refxmlstr nspace) const
{
    ArchivePtr node(new Archive(topnode_, node_->getValue(val, name, nspace), registry_, objnode_));
    return node;
}

void
Archive::setValue(int val)
{
    node_->writeValue(val);
}

void
Archive::setValue(long val)
{
    node_->writeValue(val);
}

void
Archive::setValue(unsigned int val)
{
    node_->writeValue(val);
}

void
Archive::setValue(unsigned long val)
{
    node_->writeValue(val);
}

void
Archive::setValue(double val)
{
    node_->writeValue(val);
}

void
Archive::setValue(bool val) 
{
    node_->writeValue(val);
}

void
Archive::setValue(const pyxml::xmlstr& val)
{
    node_->writeValue(val);
}

ArchivePtr 
Archive::setValue(int val, refxmlstr name, refxmlstr nspace) 
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::setValue(long val, refxmlstr name, refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::setValue(unsigned int val, refxmlstr name, refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::setValue(unsigned long val, refxmlstr name, refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::setValue(double val, refxmlstr name, refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::setValue(bool val, refxmlstr name, refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::setValue(refxmlstr val, refxmlstr name, refxmlstr nspace)
{
    ArchivePtr node(new Archive(topnode_, node_->writeValue(val, name, nspace), registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::writeBinary(const void* vals, size_t size, pyxml::refxmlstr name, pyxml::refxmlstr nspace)
{
    XMLWriter xml = node_->writeBinary(vals, size, name, nspace);
    if (xml.get() == NULL)
        return NULL; 

    ArchivePtr node(new Archive(topnode_, xml, registry_, objnode_));
    return node;
}

ArchivePtr 
Archive::loadBinary(void** vals, size_t& size, pyxml::refxmlstr name, pyxml::refxmlstr nspace) const 
{
    XMLParser xml = node_->loadBinary(vals, size, name, nspace);
    if (xml.get() == NULL)
        return NULL; 

    ArchivePtr node(new Archive(topnode_, xml, registry_, objnode_));
    return node;
}

ArchivePtr
Archive::getObjectBranch()
{
    if (objnode_.get() != NULL)
    {
        return objnode_;
    }

    objnode_ = new Archive(topnode_, topnode_->getNode("objects"), registry_);
    return objnode_;
}

XMLParser
Archive::xml() const
{
    return node_;
}

void
Archive::_verify() const
{
    if (node_.get() == NULL)
    {
        cerr << "Node is null" << endl;
        abort();
    }
}

ConstSerialMapPtr 
Archive::registry() const
{
    return registry_;
}

SerialMapPtr 
Archive::registry()
{
    return registry_;
}

SerializablePtr
Archive::getObject(serid_t id)
{
    SerializablePtr obj = registry_->get(id);
    if (obj.get() != NULL)
    {
        return obj;
    }
    else 
    {
        ArchivePtr objnode = fetchObjectBranch(id);
        std::string classname; objnode->getAttribute(classname, "classname");
        obj = SerialRuntime::getObject(objnode, classname);
        registry_->add(id, obj);
        return obj;
    }
}

void
Archive::setObject(const ConstSerializablePtr& obj)
{
    serid_t id = obj->id();
    if (registry_->has(id))
        return;

    registry_->add(id, obj);
    ArchivePtr node = makeObjectBranch(id);
    const SerialRuntime* rtinfo = obj->runtime_info();
    std::string classname = rtinfo->classname(); 
    node->setAttribute(classname, "classname");
    obj->serialize(node);
}


ArchivePtr 
Archive::serialize_data(load_t flag, int& val, refsmartstr name, refsmartstr nspace)
{
    return getValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(load_t flag, double& val, refsmartstr name, refsmartstr nspace)
{
    return getValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(load_t flag, unsigned int& val, refsmartstr name, refsmartstr nspace)
{
    return getValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(load_t flag, unsigned long& val, refsmartstr name, refsmartstr nspace)
{
    return getValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(load_t flag, std::string& val, refsmartstr name, refsmartstr nspace)
{
    return getValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(load_t flag, bool& val, refsmartstr name, refsmartstr nspace)
{
    return getValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(save_t flag, int val, refsmartstr name, refsmartstr nspace)
{
    return setValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(save_t flag, double val, refsmartstr name, refsmartstr nspace)
{
    return setValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(save_t flag, unsigned int val, refsmartstr name, refsmartstr nspace)
{
    return setValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(save_t flag, unsigned long val, refsmartstr name, refsmartstr nspace)
{
    return setValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(save_t flag, std::string val, refsmartstr name, refsmartstr nspace)
{
    return setValue(val, name, nspace);
}

ArchivePtr 
Archive::serialize_data(save_t flag, bool val, refsmartstr name, refsmartstr nspace)
{
    ArchivePtr node = setValue(val, name, nspace);
    return node;
}

void
Archive::serialize_data(load_t flag, int& val)
{
    getValue(val);
}

void
Archive::serialize_data(load_t flag, double& val)
{
    getValue(val);
}

void
Archive::serialize_data(load_t flag, unsigned int& val)
{
    getValue(val);
}

void
Archive::serialize_data(load_t flag, unsigned long& val)
{
    getValue(val);
}

void
Archive::serialize_data(load_t flag, std::string& val)
{
    getValue(val);
}

void
Archive::serialize_data(load_t flag, bool& val)
{
    getValue(val);
}

void
Archive::serialize_data(save_t flag, int val)
{
    setValue(val);
}

void
Archive::serialize_data(save_t flag, double val)
{
    setValue(val);
}

void
Archive::serialize_data(save_t flag, unsigned int val)
{
    setValue(val);
}

void
Archive::serialize_data(save_t flag, unsigned long val)
{
    setValue(val);
}

void
Archive::serialize_data(save_t flag, std::string val)
{
    setValue(val);
}

void
Archive::serialize_data(save_t flag, bool val)
{
    setValue(val);
}

std::string
Archive::getIDTag(serid_t id)
{
    return stream_printf("node%ld", id);
}

ArchivePtr 
Archive::fetchObjectBranch(serid_t id)
{
    std::string tagname = getIDTag(id);
    return objnode_->fetchBranch(tagname);
}

ArchivePtr 
Archive::makeObjectBranch(serid_t id)
{
    std::string tagname = getIDTag(id);
    return objnode_->makeBranch(tagname);
}

