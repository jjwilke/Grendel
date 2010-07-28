#ifndef gigide_archive_h_
#define gigide_archive_h_

#include <src/pyxml/pyxml.h>
#include <src/ref.h>
#include <src/set.h>
#include <vector>
#include <src/serialabstract.h>

#define ArchivePtr boost::intrusive_ptr<smartptr::Archive>  
#define ConstArchivePtr boost::intrusive_ptr<const smartptr::Archive>  

#define SerialMapPtr boost::intrusive_ptr<smartptr::SerialMap>  
#define ConstSerialMapPtr boost::intrusive_ptr<const smartptr::SerialMap>  

#define vector_item_tagname "vector_item"
#define set_item_tagname "set_item"
#define map_item_tagname "map_item"
#define map_key_tagname "key"
#define map_value_tagname "value"

#define blankstr std::string("")

namespace smartptr {

class SerialMap : public Countable {

    private:
        typedef std::map<serid_t, boost::intrusive_ptr<const Serializable> > registry_map;

        registry_map registry_;

    public:
        SerialMap();

        boost::intrusive_ptr<Serializable> get(serid_t id) const;

        bool has(serid_t id) const;
        
        void add(serid_t id, const ConstSerializablePtr& ptr);

};

class Archive : public Countable {
    
    public:
        typedef enum { Read } load_t;
        typedef enum { Write } save_t;

    private:

        PyXMLDomParserPtr topnode_;

        PyXMLDomParserPtr node_;

        ArchivePtr objnode_;
        
        void _verify() const;

        SerialMapPtr registry_;

        Archive(
            const PyXMLDomParserPtr& topnode,
            const PyXMLDomParserPtr& node,
            const SerialMapPtr& registry
        );

        ArchivePtr fetchObjectBranch(serid_t id);

        ArchivePtr makeObjectBranch(serid_t id);

        SerializablePtr getObject(serid_t id);

        void setObject(const ConstSerializablePtr& obj);

        static std::string getIDTag(serid_t id);

    public:

        Archive();

        Archive(refsmartstr filename);

        Archive(
            const PyXMLDomParserPtr& topnode,
            const PyXMLDomParserPtr& node,
            const SerialMapPtr& registry,
            const ArchivePtr& objnode
        );

        ~Archive();

        PyXMLDomParserPtr xml() const;

        ArchivePtr getTop() const;

        ArchivePtr getObjectBranch();

        SerialMapPtr registry();

        ConstSerialMapPtr registry() const;

        ArchivePtr fetchBranch(refsmartstr tagname, refsmartstr nspace = blankstr) const;

        ArchivePtr getBranch(refsmartstr tagname, refsmartstr nspace = blankstr);

        void fetchBranches(std::vector<ArchivePtr>& nodes, refsmartstr tagname, refsmartstr nspace = blankstr) const;

        ArchivePtr makeBranch(refsmartstr tagname, refsmartstr nspace = blankstr);

        bool hasAttribute(refsmartstr attrname) const;
            
        void getAttribute(int& val, refsmartstr attrname) const;

        void getAttribute(long& val, refsmartstr attrname) const;

        void getAttribute(size_t& val, refsmartstr attrname) const;

        void getAttribute(double& val, refsmartstr attrname) const;

        void getAttribute(unsigned int& val, refsmartstr attrname) const;

        void getAttribute(bool& val, refsmartstr attrname) const;

        void getAttribute(pyxml::xmlstr& val, refsmartstr attrname) const;

        void setAttribute(int val, refsmartstr attrname);

        void setAttribute(long val, refsmartstr attrname);

        void setAttribute(size_t val, refsmartstr attrname);

        void setAttribute(double val, refsmartstr attrname);

        void setAttribute(unsigned int val, refsmartstr attrname);

        void setAttribute(bool val, refsmartstr attrname);

        void setAttribute(refsmartstr val, refsmartstr attrname);

        //These are not const as the the get function may need to scan through a binary file, changing the file
        void getValue(void*& val);

        void getValue(int& val);

        void getValue(long& val);

        void getValue(unsigned int& val);

        void getValue(unsigned long& val);

        void getValue(double& val);

        void getValue(bool& val);

        void getValue(pyxml::xmlstr& value);

        void setValue(void* val);

        void setValue(int val);

        void setValue(long val);

        void setValue(unsigned int val);

        void setValue(unsigned long val);

        void setValue(double val);

        void setValue(bool val);

        void setValue(const pyxml::xmlstr& value);

        ArchivePtr getValue(void*& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(int& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(long& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(unsigned int& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(unsigned long& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(double& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(bool& val, refsmartstr val, refsmartstr nspace = blankstr) const;

        ArchivePtr getValue(pyxml::xmlstr& value, refsmartstr value, refsmartstr nspace = blankstr) const;

        ArchivePtr setValue(void* val, refsmartstr val, refsmartstr nspace = blankstr);

        ArchivePtr setValue(int val, refsmartstr name, refsmartstr nspace = blankstr);

        ArchivePtr setValue(long val, refsmartstr name, refsmartstr nspace = blankstr);

        ArchivePtr setValue(unsigned int val, refsmartstr val, refsmartstr nspace = blankstr);

        ArchivePtr setValue(unsigned long val, refsmartstr val, refsmartstr nspace = blankstr);

        ArchivePtr setValue(double val, refsmartstr name, refsmartstr nspace = blankstr);

        ArchivePtr setValue(bool val, refsmartstr val, refsmartstr nspace = blankstr);

        ArchivePtr setValue(refsmartstr val, refsmartstr name, refsmartstr nspace = blankstr);

        ArchivePtr writeBinary(const void* vals, size_t size, refsmartstr name, refsmartstr nspace = blankstr);

        ArchivePtr loadBinary(void** vals, size_t& size, refsmartstr name, refsmartstr nspace = blankstr) const; 
        
        template <class T>
        ArchivePtr 
        serialize(load_t flag, boost::intrusive_ptr<T>& obj, refsmartstr name, refsmartstr nspace = blankstr)
        {
            ArchivePtr node = fetchBranch(name, nspace);
            if (node.get() == NULL)
            {
                obj = NULL;
                return NULL;
            }
            node->serialize(flag, obj);
            return node;
        }

        template <class T>
        void
        serialize(load_t flag, boost::intrusive_ptr<T>& obj)
        {
            serid_t id; 
            getValue(id);
            obj = boost::dynamic_pointer_cast<T, Serializable>(getObject(id));
        }

        template <class T>
        ArchivePtr 
        serialize(save_t flag, const boost::intrusive_ptr<T>& obj, refsmartstr name, refsmartstr nspace = blankstr)
        {
            if (obj.get() == NULL)
                return NULL;

            //AssertNotHere<false> test;
            ArchivePtr node = makeBranch(name, nspace);
            node->serialize(flag, obj);
            return node;
        }

        template <class T>
        void
        serialize(save_t flag, const boost::intrusive_ptr<T>& obj)
        {
            serid_t id = obj->id(); 
            setValue(id);
            setObject(obj);
        }

        void serialize_data(load_t flag, int& val);
        void serialize_data(load_t flag, double& val);
        void serialize_data(load_t flag, unsigned int& val);
        void serialize_data(load_t flag, unsigned long& val);
        void serialize_data(load_t flag, std::string& val);
        void serialize_data(load_t flag, bool& val);

        void serialize_data(save_t flag, int val);
        void serialize_data(save_t flag, double val);
        void serialize_data(save_t flag, unsigned int val);
        void serialize_data(save_t flag, unsigned long val);
        void serialize_data(save_t flag, std::string val);
        void serialize_data(save_t flag, bool val);

        ArchivePtr serialize_data(load_t flag, int& val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(load_t flag, double& val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(load_t flag, unsigned int& val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(load_t flag, unsigned long& val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(load_t flag, std::string& val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(load_t flag, bool& val, refsmartstr name, refsmartstr nspace = blankstr);

        ArchivePtr serialize_data(save_t flag, int val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(save_t flag, double val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(save_t flag, unsigned int val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(save_t flag, unsigned long val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(save_t flag, std::string val, refsmartstr name, refsmartstr nspace = blankstr);
        ArchivePtr serialize_data(save_t flag, bool val, refsmartstr name, refsmartstr nspace = blankstr);

#if 1
        template <class T>
        ArchivePtr 
        serialize(load_t flag, std::vector<T>& vals, refsmartstr tagname, refsmartstr nspace = blankstr)
        {
            ArchivePtr topnode = getBranch(tagname, nspace);
            if (topnode.get() == NULL)
                return NULL;

            std::vector<ArchivePtr> nodes; topnode->fetchBranches(nodes, vector_item_tagname, nspace);
            std::vector<ArchivePtr>::iterator it;
            for (it = nodes.begin(); it != nodes.end(); ++it)
            {
                ArchivePtr node = *it;
                T val; 
                serial_call_load(node, val);
                vals.push_back(val);
            }
            return topnode;
        }

        template <class T>
        ArchivePtr 
        serialize(save_t flag, const std::vector<T>& vals, refsmartstr tagname, refsmartstr nspace = blankstr)
        {
            ArchivePtr topnode = getBranch(tagname, nspace);
            typename std::vector<T>::const_iterator it(vals.begin());
            for ( ; it != vals.end(); ++it)
            {
                T val = *it;
                serial_call_save(topnode, val, vector_item_tagname);
            }
            return topnode;
        }
#endif

        template <class T>
        ArchivePtr
        serialize(save_t flag, const gigide::Set<T>& items, refsmartstr tagname, refsmartstr nspace = blankstr)
        {
            ArchivePtr topnode = getBranch(tagname, nspace);
            for (int i=0; i < items.size(); ++i)
                serial_call_save(topnode, items[i], set_item_tagname);
            return topnode;
        }

        template <class T>
        ArchivePtr
        serialize(load_t flag, gigide::Set<T>& items, refsmartstr tagname, refsmartstr nspace = blankstr)
        {
            ArchivePtr topnode = getBranch(tagname, nspace);
            std::vector<ArchivePtr> nodes; topnode->fetchBranches(nodes, set_item_tagname);
            for (std::vector<ArchivePtr>::iterator it(nodes.begin()); it != nodes.end(); ++it)
            {
                ArchivePtr node = *it;
                T val;
                serial_call_load(node, val);
                items.append(val);
            }
            return topnode;
        }

        template <class Key_t, class Val_t>
        ArchivePtr
        serialize(save_t flag, const std::map<Key_t, Val_t>& keymap, refsmartstr tagname, refsmartstr nspace = blankstr)
        {
            ArchivePtr topnode = getBranch(tagname, nspace);
            typename std::map<Key_t, Val_t>::const_iterator it(keymap.begin());
            for ( ; it != keymap.end(); ++it)
            {
                ArchivePtr node = topnode->makeBranch(map_item_tagname);
                serial_call_save(node, it->first, map_key_tagname);
                serial_call_save(node, it->second, map_value_tagname);
            }
            return topnode;
        }

        template <class Key_t, class Val_t>
        ArchivePtr
        serialize(load_t flag, std::map<Key_t, Val_t>& keymap, refsmartstr tagname, refsmartstr nspace = blankstr)
        {
            ArchivePtr topnode = getBranch(tagname, nspace);
            std::vector<ArchivePtr> nodes; topnode->fetchBranches(nodes, map_item_tagname);
            for (std::vector<ArchivePtr>::iterator it(nodes.begin()); it != nodes.end(); ++it)
            {
                ArchivePtr node = *it;
                Key_t key; Val_t val;
                serial_call_load(node, key, map_key_tagname);
                serial_call_load(node, val, map_value_tagname);
                keymap[key] = val;
            }
            return topnode;
        }

        template <class T>
        void
        enum_cast(unsigned int tmp, T& val)
        {
            val = (T) tmp;
        }

};

#undef blankstr

}

#endif
