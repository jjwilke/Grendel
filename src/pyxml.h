#ifndef _gigide_pyxml_h_
#define _gigide_pyxml_h_

#include <Python.h>
#include <vector>
#include <map>
#include <src/deftypes.h>
#include <string>

namespace pyxml {

typedef std::string xmlstr;
typedef const xmlstr& refxmlstr;
#define blankstr std::string("")

#define ConstPyXMLDomParserPtr boost::intrusive_ptr<const PyXMLDomParser>  
#define PyXMLDomParserPtr boost::intrusive_ptr<pyxml::PyXMLDomParser>  
#define XMLParser boost::intrusive_ptr<pyxml::PyXMLDomParser>  
#define XMLWriter boost::intrusive_ptr<pyxml::PyXMLDomParser>  

class PyXMLDomParser : public smartptr::Countable {
    
    private:
        PyObject* parser_;

        void _getElementsByTagName(std::vector<PyXMLDomParserPtr>& nodes, refxmlstr tagname) const;

        PyXMLDomParserPtr _createNode(refxmlstr name);

    public:
        PyXMLDomParser();

        PyXMLDomParser(refxmlstr filename);

        PyXMLDomParser(PyObject* parser);

        void getElementsByTagName(std::vector<PyXMLDomParserPtr>& nodes, refxmlstr tagname, refxmlstr nspace = blankstr) const;

        std::string getAttribute(refxmlstr attrname) const;

        template <class T>
        void getAttribute(T& val, refxmlstr attrname) const;

        void getAttribute(int& val, refxmlstr attrname) const;

        void getAttribute(long& val, refxmlstr attrname) const;

        void getAttribute(size_t& val, refxmlstr attrname) const;

        void getAttribute(double& val, refxmlstr attrname) const;

        void getAttribute(unsigned int& val, refxmlstr attrname) const;

        void getAttribute(bool& val, refxmlstr attrname) const;

        bool hasAttribute(refxmlstr attrname) const;

        std::string getText() const;

        template <class T>
        void getValue(T& val) const;

        void getValue(void* val) const;

        void getValue(double& val) const;

        void getValue(std::string& val) const;

        void getValue(int& val) const;

        void getValue(long& val) const;

        void getValue(unsigned long& val) const;

        void getValue(unsigned int& val) const;

        void getValue(bool& val) const;

        PyXMLDomParserPtr getTopnode() const;

        PyXMLDomParserPtr getValue(void*& value, refxmlstr name, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr getValue(int& value, refxmlstr name, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr getValue(double& value, refxmlstr name, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr getValue(long& value, refxmlstr name, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr getValue(unsigned long& val, refxmlstr name, refxmlstr nspace) const;

        PyXMLDomParserPtr getValue(unsigned int& val, refxmlstr name, refxmlstr nspace) const;

        PyXMLDomParserPtr getValue(bool& val, refxmlstr name, refxmlstr nspace) const;

        PyXMLDomParserPtr getValue(xmlstr& value, refxmlstr name, refxmlstr nspace = blankstr) const;

        template <class T>
        PyXMLDomParserPtr readList(std::vector<T>& list, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr readList(std::vector<double>&, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr readList(std::vector<int>&, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr readList(std::vector<std::string>& list, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = blankstr) const;

        template <typename T, typename U>
        PyXMLDomParserPtr readMap(std::map<T, U>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr readMap(std::map<int, int>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr readMap(std::map<int, unsigned long>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr readMap(std::map<std::string, int>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr getNode(refxmlstr tagname, refxmlstr name_space = blankstr);

        PyXMLDomParserPtr fetchNode(refxmlstr tagname, refxmlstr name_space = blankstr) const;

        PyXMLDomParserPtr getChildElement(refxmlstr tagname, refxmlstr nspace = blankstr) const;

        PyXMLDomParserPtr getChildByAttribute(refxmlstr tagname, refxmlstr attrname, refxmlstr attrval) const;

        PyXMLDomParserPtr createNode(refxmlstr name, refxmlstr nspace = blankstr);

        void toFile(refxmlstr filename) const;

        void setAttribute(refxmlstr val, refxmlstr name);

        void setAttribute(const char* val, refxmlstr name);

        void setAttribute(long val, refxmlstr name);

        void setAttribute(int val, refxmlstr name);

        void setAttribute(unsigned int val, refxmlstr name);

        void setAttribute(size_t val, refxmlstr name);

        void setAttribute(double val, refxmlstr name);

        void setAttribute(bool val, refxmlstr name);

        void writeValue(void* val);

        void writeValue(double val);

        void writeValue(int val);

        void writeValue(long val);

        void writeValue(unsigned long val);

        void writeValue(unsigned int val);

        void writeValue(bool val);

        void writeValue(refxmlstr val);

        PyXMLDomParserPtr writeValue(void* val, refxmlstr name);

        PyXMLDomParserPtr writeValue(int val, refxmlstr name);

        PyXMLDomParserPtr writeValue(long val, refxmlstr name);

        PyXMLDomParserPtr writeValue(double val, refxmlstr name);

        PyXMLDomParserPtr writeValue(refxmlstr val, refxmlstr name);

        PyXMLDomParserPtr writeValue(int val, refxmlstr name, refxmlstr nspace);

        PyXMLDomParserPtr writeValue(long val, refxmlstr name, refxmlstr nspace);

        PyXMLDomParserPtr writeValue(double val, refxmlstr name, refxmlstr nspace);

        PyXMLDomParserPtr writeValue(unsigned long val, refxmlstr name, refxmlstr nspace);

        PyXMLDomParserPtr writeValue(unsigned int val, refxmlstr name, refxmlstr nspace);

        PyXMLDomParserPtr writeValue(bool val, refxmlstr name, refxmlstr nspace);

        PyXMLDomParserPtr writeValue(refxmlstr val, refxmlstr name, refxmlstr nspace);

        template <class T>
        PyXMLDomParserPtr writeList(const std::vector<T>& list, refxmlstr tagname, refxmlstr itemname, refxmlstr nspace);

        PyXMLDomParserPtr writeList(const std::vector<double>&, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr writeList(const std::vector<int>&, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr writeList(const std::vector<std::string>&, refxmlstr listname, refxmlstr itemname, refxmlstr nspace = std::string(""));

        template <class T, class U>
        PyXMLDomParserPtr writeMap(const std::map<T, U>& vals, refxmlstr tagname, refxmlstr nspace);

        PyXMLDomParserPtr writeMap(const std::map<int, int>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr writeMap(const std::map<int, unsigned long>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        PyXMLDomParserPtr writeMap(const std::map<std::string, int>& vals, refxmlstr tagname, refxmlstr nspace = blankstr);

        XMLWriter writeBinary(const void* vals, size_t size, refxmlstr name, refxmlstr nspace = blankstr);

        XMLParser loadBinary(void** vals, size_t& size, refxmlstr name, refxmlstr nspace = blankstr) const; 

        void setText(refxmlstr text);

        void addBinary(PyObject* obj);

        void getBinary(void** vals, size_t& size) const; 

        ~PyXMLDomParser();
};

}

#endif
