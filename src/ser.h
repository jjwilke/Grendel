#ifndef _gigide_ser_h_
#define _gigide_ser_h_

#include <src/serialize.h>

namespace smartptr {

class Archive;
class Serializable;

template <class T>
Serializable*
serialize(const ArchivePtr& arch)
{
    if (arch.get() == NULL)
    {
        std::cerr << "Cannot build object. Null xml node received" << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    Serializable* obj = new T(arch);
    return obj;
}

template <class T>
boost::intrusive_ptr<T>
getObject(const ArchivePtr& arch)
{
    if (arch.get() == NULL)
    {
        std::cerr << "Cannot build object. Null xml node received." << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    SerializablePtr parent(Serializable::getObject(arch));
    boost::intrusive_ptr<T> obj(boost::dynamic_pointer_cast<T, Serializable>(parent));
    return obj;
}

template <class T>
ArchivePtr
getObject(boost::intrusive_ptr<T>& item, const ArchivePtr& arch, refsmartstr tagname)
{
    ArchivePtr node = arch->fetchBranch(tagname);
    if (node.get() == NULL)
        return NULL;

    item = getObject<T>(node);
    return node;
}

template <class T>
ArchivePtr
getObject(boost::intrusive_ptr<T>& item, const ArchivePtr& arch, refsmartstr tagname, refsmartstr nspace)
{
    ArchivePtr node = arch->fetchBranch(tagname, nspace);
    if (node.get() == NULL)
    {
        return NULL;
    }

    item = getObject<T>(node);
    return node;
}

template <class T>
void
getObjects(std::vector<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr tagname)
{  
    if (arch.get() == NULL)
    {
        std::cerr << "Cannot build objects for tagname " << tagname << ". Null xml node received" << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    std::vector<ArchivePtr> nodes; arch->fetchBranches(nodes, tagname);

    std::vector<ArchivePtr>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); ++it)
    {
        ArchivePtr node = *it;
        items.push_back(getObject<T>(node)); 
    }
}

template <class T>
ArchivePtr
getObjects(std::vector<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr topname, refsmartstr tagname)
{
    if (arch.get() == NULL) //wtf man
    {
        std::cerr << "Cannot build objects for topname: " << topname << " tagname: " << tagname << ". Null xml node received" << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    ArchivePtr topnode = arch->fetchBranch(topname);
    if (topnode.get() == NULL) //wtf man
    {
        std::cerr << "Tagname " << topname << " yielded null node in getObjects." << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    getObjects(items, topnode, tagname);
    return topnode;
}

template <class T>
ArchivePtr
getObjects(std::vector<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr topname, refsmartstr tagname, refsmartstr name_space)
{
    if (arch.get() == NULL) //wtf man
    {
        std::cerr << "Cannot build objects for topname: " << topname << " tagname: " << tagname << ". Null xml node received" << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    ArchivePtr topnode = arch->fetchBranch(topname, name_space);
    if (topnode.get() == NULL) //wtf man
    {
        std::cerr << "Tagname " << topname << " yielded null node in getObjects." << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    getObjects(items, topnode, tagname);
    return topnode;
}

template <class T>
void
getObjects(gigide::Set<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr tagname)
{
    if (arch.get() == NULL)
    {
        std::cerr << "Cannot build objects for tagname " << tagname << ". Null xml node received." << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }

    std::vector<ArchivePtr> nodes; arch->fetchBranches(nodes, tagname);

    std::vector<ArchivePtr>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); ++it)
    {
        ArchivePtr node = *it;

        items.append(getObject<T>(node));
    }
}

template <class T>
ArchivePtr
getObjects(gigide::Set<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr topname, refsmartstr tagname)
{
    if (arch.get() == NULL)
    {
        std::cerr << "Cannot build objects for topname " << topname << " on tagname " << tagname << ". Null xml node received." << std::endl;
        std::cerr << "Line: " << __LINE__ << std::endl;
        std::cerr << "File: " << __FILE__ << std::endl;
        abort();
    }
    ArchivePtr topnode = arch->fetchBranch(topname);
    getObjects(items, topnode, tagname);

    return topnode;
}

template <class Key_t, class Val_t>
ArchivePtr
getObjects(std::map<Key_t, boost::intrusive_ptr<Val_t> >& keymap, const ArchivePtr& arch, refsmartstr tagname)
{
    //build a map of the ids
    std::map<Key_t, serid_t> ids;
    ArchivePtr topnode = arch->getMap(ids, tagname);
    if (topnode.get() == NULL)
        return NULL;

    typename std::map<Key_t, serid_t>::const_iterator it(ids.begin());
    for ( ; it != ids.end(); ++it)
    {
        SerializablePtr parent = Serializable::getObject(arch, it->second);
        //add the entry to the map
        boost::intrusive_ptr<Val_t> obj(boost::dynamic_pointer_cast<Val_t, Serializable>(parent));
        keymap[it->first] = obj;
    }
    return topnode;
}

template <class T>
ArchivePtr
setObject(const boost::intrusive_ptr<T>& obj, const ArchivePtr& arch, refsmartstr tagname)
{
    if (obj.get() == NULL)
        return NULL;

    return obj->write(arch, tagname);
}

template <class T>
ArchivePtr
setObject(const boost::intrusive_ptr<T>& obj, const ArchivePtr& arch, refsmartstr tagname, refsmartstr name_space)
{
    if (obj.get() == NULL)
        return NULL;

    return obj->write(arch, tagname, name_space);
}

template <class T>
ArchivePtr
setObjects(const std::vector<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr group_name, refsmartstr node_name, refsmartstr name_space)
{
    ArchivePtr topnode = arch->getBranch(group_name, name_space);
    for (int i=0; i < items.size(); ++i)
        items[i]->write(topnode, node_name);
    return topnode;
}

template <class T>
ArchivePtr
setObjects(const std::vector<boost::intrusive_ptr<T> >& items, const ArchivePtr& arch, refsmartstr group_name, refsmartstr node_name)
{
    ArchivePtr topnode = arch->getBranch(group_name);
    for (int i=0; i < items.size(); ++i)
        items[i]->write(topnode, node_name);
    return topnode;
}

template <class T>
ArchivePtr
setObjects(const gigide::Set<boost::intrusive_ptr<T> >& items,  const ArchivePtr& arch, refsmartstr group_name, refsmartstr node_name)
{
    ArchivePtr topnode = arch->getBranch(group_name);
    for (int i=0; i < items.size(); ++i)
        items[i]->write(topnode, node_name);
    return topnode;
}

template <class Key_t, class Val_t>
ArchivePtr
setObjects(const std::map<Key_t, boost::intrusive_ptr<Val_t> >& keymap, const ArchivePtr& arch, refsmartstr tagname)
{
    //build a map of the ids
    std::map<Key_t, serid_t> ids;
    typename std::map<Key_t, boost::intrusive_ptr<Val_t> >::const_iterator it(keymap.begin());
    for ( ; it != keymap.end(); ++it)
    {
        //add the id to the map
        ids[it->first] = it->second->id();
        //write the objects to the archive
        it->second->writeObjectNode(arch);
    }

    ArchivePtr topnode = arch->setMap(ids, tagname);
    return topnode;
}

template <class T>
ArchivePtr
value(const ArchivePtr& arch, const boost::intrusive_ptr<const T>& obj, T val, refsmartstr name, refsmartstr nspace = blankstr)
{
    return arch->setValue(val, name, nspace);
}

template <class T>
ArchivePtr
value(const ArchivePtr& arch, const boost::intrusive_ptr<T>& obj, T& val, refsmartstr name, refsmartstr nspace = blankstr)
{
    return arch->getValue(val, name, nspace);
}

template <class T>
ArchivePtr
object(const ArchivePtr& arch, const boost::intrusive_ptr<const T>& obj, T val, refsmartstr name, refsmartstr nspace = blankstr)
{
    return arch->setValue(val, name, nspace);
}

template <class T>
ArchivePtr
object(const ArchivePtr& arch, const boost::intrusive_ptr<T>& obj, T& val, refsmartstr name, refsmartstr nspace = blankstr)
{
    return arch->getValue(val, name, nspace);
}

} //end namespace smarptr

#endif
