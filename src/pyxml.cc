#include <src/pyxml.h>
#include <src/exception.h>
#include <src/pyregexp.h>
#include <src/printstream.h>

#define heisenbug cout << "Heisenbug: " << __FILE__ << " " << __LINE__ << endl

#define refcount(x) cout << stream_printf("%s = %ld", #x, x->ob_refcnt) << endl;
#define pointer(x) cout << stream_printf("%s = %p", #x, x) << endl;

using namespace pyxml;
using namespace gigide;
using namespace std;

//using pyxml::refxmlstr;

void
PY_CHECK(PyObject* obj, const char* msg)
{
    if (obj == NULL)
    {
        except(msg);
    }
}

void
ASSERT_NONNULL(XMLParser node, refxmlstr name)
{
    if (node.get() == NULL)
    {
        cerr << "Null node obtained for tagname " << name << endl;
        abort();
    }
}

PyXMLDomParser::PyXMLDomParser(refxmlstr filename)
{
    PyObject *pmod = NULL, *pargs = NULL, *pfunc = NULL, *preturn = NULL, *test = NULL;
    pmod = PyImport_ImportModule("PyXML"); //I could validate, but this should never fail
	if (pmod == NULL)
	{
		except("Could not find module PyXML\n");
	}


    pfunc = PyObject_GetAttrString(pmod, "construct_parser");
	if (pfunc == NULL)
	{
		except("Could not find python function construct_parser\n");
	}

    pargs = Py_BuildValue("(s)", filename.c_str());
    preturn = PyEval_CallObject(pfunc, pargs);
    if (preturn == NULL)
    {
		except(stream_printf("XML error on file was not returned %s", filename.c_str()));
    }


    PyArg_Parse(preturn, "O", &parser_);
    Py_INCREF(parser_);

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);
    Py_DECREF(pmod);
}

PyXMLDomParser::PyXMLDomParser(PyObject* parser)
    : parser_(parser)
{
    Py_INCREF(parser_);
}

PyXMLDomParser::PyXMLDomParser()
{
    PyObject *pmod = NULL, *pargs = NULL, *pfunc = NULL, *preturn = NULL;
    pmod = PyImport_ImportModule("PyXML"); //I could validate, but this should never fail
	if (pmod == NULL)
	{
		except("Could not find module PyXML\n");
	}

    pfunc = PyObject_GetAttrString(pmod, "construct_writer");
	if (pfunc == NULL)
	{
		except("Could not find python function construct_writer\n");
	}

    pargs = Py_BuildValue("()");
    preturn = PyEval_CallObject(pfunc, pargs);
    if (preturn == NULL)
    {
		except(stream_printf("Error creating XML writer"));
    }

    PyArg_Parse(preturn, "O", &parser_);
    Py_INCREF(parser_);

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);
    Py_DECREF(pmod);
}

PyXMLDomParserPtr
PyXMLDomParser::getTopnode() const
{
    PyObject *pargs = NULL, *pfunc = NULL, *preturn = NULL, *test = NULL, *obj = NULL;
    pfunc = PyObject_GetAttrString(parser_, "getTopnode");
	if (pfunc == NULL)
	{
		except("Could not find python function getTopnode\n");
	}
    pargs = Py_BuildValue("()");

    preturn = PyEval_CallObject(pfunc, pargs);
    PY_CHECK(preturn, "top node is null for xml parser");
    PyArg_Parse(preturn, "O", &obj);

    Py_DECREF(pargs);
    Py_DECREF(pfunc);
    //must do this first to increment refcount
    PyXMLDomParserPtr parser(new PyXMLDomParser(obj));
    Py_DECREF(preturn);
    return parser;
}

template <typename T, typename U>
PyXMLDomParserPtr
PyXMLDomParser::readMap(
    std::map<T, U>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    XMLParser node = fetchNode(tagname);
    if (node.get() == NULL)
        return NULL;

    vector<XMLParser> nodes; node->getElementsByTagName(nodes, "map_entry");
    vector<XMLParser>::iterator it;
    for (it = nodes.begin(); it != nodes.end(); ++it)
    {
        XMLParser pair_node = *it;
        U value; 
        pair_node->getValue(value);
        T key; 
        pair_node->getAttribute(key, "key");
        vals[key] = value;
    }
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::readMap(
    std::map<int, int>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    return readMap<int,int>(vals, tagname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::readMap(
    std::map<int, unsigned long>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    return readMap<int,unsigned long>(vals, tagname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::readMap(
    std::map<std::string, int>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    return readMap<std::string,int>(vals, tagname, nspace);
}

std::string
PyXMLDomParser::getText() const
{
    PyObject* pfunc, *pargs, *preturn;

    pfunc = PyObject_GetAttrString(parser_, "getText");
    if (pfunc == NULL)
    {
        cout << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("()");

    preturn = PyEval_CallObject(pfunc, pargs);
    char* val;
    PyArg_Parse(preturn, "s", &val);

    std::string str = val;

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);

    return str;
}

template <class T>
void
PyXMLDomParser::getValue(T& val) const
{
    stringstream sstr(getText());
    sstr >> val;
}

void
PyXMLDomParser::getValue(std::string& val) const
{
    val = getText();
}

void
PyXMLDomParser::getValue(long& val) const
{
    getValue<long>(val);
}

void
PyXMLDomParser::getValue(unsigned int& val) const
{
    getValue<unsigned int>(val);
}

void
PyXMLDomParser::getValue(unsigned long& val) const
{
    getValue<unsigned long>(val);
}

void
PyXMLDomParser::getValue(double& val) const
{
    getValue<double>(val);
}

void
PyXMLDomParser::getValue(bool& val) const
{
    std::string text = getText();
    if (text == "true")
        val = true;
    else
        val = false;
}

void
PyXMLDomParser::getValue(int& val) const
{
    getValue<int>(val);
}


PyXMLDomParserPtr
PyXMLDomParser::getValue(int& val, refxmlstr name, refxmlstr nspace) const
{

    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue<int>(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::getValue(double& val, refxmlstr name, refxmlstr nspace) const
{
    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue<double>(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::getValue(unsigned int& val, refxmlstr name, refxmlstr nspace) const
{

    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue<unsigned int>(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::getValue(unsigned long& val, refxmlstr name, refxmlstr nspace) const
{

    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue<unsigned long>(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::getValue(long& val, refxmlstr name, refxmlstr nspace) const
{

    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue<long>(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::getValue(bool& val, refxmlstr name, refxmlstr nspace) const
{

    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::getValue(xmlstr& val, refxmlstr name, refxmlstr nspace) const
{

    XMLParser node = fetchNode(name);
    ASSERT_NONNULL(node, name);
    node->getValue<xmlstr>(val);
    return node;
}

void
PyXMLDomParser::getAttribute(bool& val, refxmlstr name) const
{
    std::string attr = getAttribute(name);
    if (attr == "true")
        val = true;
    else
        val = false;
}

void
PyXMLDomParser::getAttribute(unsigned int& val, refxmlstr attrname) const
{
    getAttribute<unsigned int>(val, attrname);
}

void
PyXMLDomParser::getAttribute(long& val, refxmlstr attrname) const
{
    getAttribute<long>(val, attrname);
}

void
PyXMLDomParser::getAttribute(int& val, refxmlstr attrname) const
{
    getAttribute<int>(val, attrname);
}

void
PyXMLDomParser::getAttribute(size_t& val, refxmlstr attrname) const
{
    getAttribute<size_t>(val, attrname);
}

void
PyXMLDomParser::getAttribute(double& val, refxmlstr attrname) const
{
    getAttribute<double>(val, attrname);
}

template <class T>
void
PyXMLDomParser::getAttribute(T& val, refxmlstr attrname) const
{
    stringstream sstr(getAttribute(attrname));
    sstr >> val;
}

std::string
PyXMLDomParser::getAttribute(refxmlstr attrname) const
{
    if (!hasAttribute(attrname))
    {
        cerr << "node does not have attribute " << attrname << endl;
        abort();
    }

    PyObject* pfunc, *pargs, *preturn;

    pfunc = PyObject_GetAttrString(parser_, "getAttribute");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(s)", attrname.c_str());

    preturn = PyEval_CallObject(pfunc, pargs);
    char* val;
    PyArg_Parse(preturn, "s", &val);
    std::string str = val;

    if (str == "")
    {
        //cout << "Warning! Attribute " << attrname << " has empty value.  Is this right?" << endl;
    }

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);

    return str;
}

bool
PyXMLDomParser::hasAttribute(refxmlstr attrname) const
{
    PyObject* pfunc = NULL, *pargs = NULL, *preturn = NULL;

    pfunc = PyObject_GetAttrString(parser_, "hasAttribute");
    if (pfunc == NULL)
    {
        cout << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(s)", attrname.c_str());

    preturn = PyEval_CallObject(pfunc, pargs);

    int hasattr;
    PyArg_Parse(preturn, "i", &hasattr);

    Py_DECREF(preturn);
    Py_DECREF(pfunc);
    Py_DECREF(pargs);

    return hasattr;
}

PyXMLDomParserPtr
PyXMLDomParser::getNode(refxmlstr tagname, refxmlstr nspace)
{
    PyXMLDomParserPtr node = fetchNode(tagname, nspace);
    if (node.get() == NULL)
        return createNode(tagname, nspace);
    else
        return node;
}

PyXMLDomParserPtr
PyXMLDomParser::fetchNode(refxmlstr tagname, refxmlstr nspace) const
{
    return getChildElement(tagname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::getChildElement(refxmlstr tagname, refxmlstr nspace) const
{
    vector<PyXMLDomParserPtr> children;
    getElementsByTagName(children, tagname, nspace);
    int nchild = children.size();
    if (nchild > 1)
    {
        cerr << stream_printf("getChildElement expected to find only 1 node with tagname %s"
                             ", but found %d nodes", tagname.c_str(), nchild) << endl;
    }
    else if (nchild == 0)
    {
        //cerr << "0 nodes found for " << tagname << " " << name_space << endl;
        return NULL; //send back null
    }

    return children[0];
}

PyXMLDomParserPtr
PyXMLDomParser::getChildByAttribute(
    refxmlstr tagname,
    refxmlstr attrname,
    refxmlstr attrval
) const
{
    vector<PyXMLDomParserPtr> children;
    getElementsByTagName(children, tagname);
    int nchild = children.size();

    vector<PyXMLDomParserPtr>::iterator it;
    for (it = children.begin(); it != children.end(); ++it)
    {
        PyXMLDomParserPtr node = *it;
        if (!node->hasAttribute(attrname))
            continue;

        if (node->getAttribute(attrname) == attrval)
            return node;
    }

    //not a valid node
    return 0;
}

template <class T>
PyXMLDomParserPtr
PyXMLDomParser::readList(
    vector<T>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
) const
{
    PyXMLDomParserPtr topnode = getChildElement(listname, nspace);
    if (topnode.get() == NULL)
        return NULL;

    vector<XMLParser> nodes; topnode->getElementsByTagName(nodes, itemname);
    for (int i=0; i < nodes.size(); ++i)
    {
        T next;
        nodes[i]->getValue(next);
        values.push_back(next);
    }
    return topnode;
}

PyXMLDomParserPtr
PyXMLDomParser::readList(
    vector<int>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
) const
{
    return readList<int>(values, listname, itemname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::readList(
    vector<double>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
) const
{
    return readList<double>(values, listname, itemname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::readList(
    vector<std::string>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
) const
{
    return readList<std::string>(values, listname, itemname, nspace);
}

void
PyXMLDomParser::getElementsByTagName(
    vector<PyXMLDomParserPtr>& children,
    refxmlstr tagname,
    refxmlstr nspace
) const
{
    if (nspace == "")
    {
        _getElementsByTagName(children, tagname);
        return;
    }

    PyObject* pfunc, *pargs;

    pfunc = PyObject_GetAttrString(parser_, "getElementsByTagNameInNamespace");
    if (pfunc == NULL)
    {
        cout << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(ss)", nspace.c_str(), tagname.c_str());

    PyObject* nodes;
    nodes = PyEval_CallObject(pfunc, pargs);
    if (nodes == NULL)
    {
        cout << "null tuple returned" << endl;
        abort();
    }

    int size = PyList_Size(nodes);
    for (int i=0; i < size; ++i)
    {
        children.push_back(new PyXMLDomParser(PyList_GetItem(nodes, i)));
    }

    Py_DECREF(pargs);
    Py_DECREF(nodes);
    Py_DECREF(pfunc);
}

void
PyXMLDomParser::_getElementsByTagName(
    vector<PyXMLDomParserPtr>& children,
    refxmlstr tagname
) const
{
    PyObject* pfunc, *pargs;

    pfunc = PyObject_GetAttrString(parser_, "getElementsByTagName");
    if (pfunc == NULL)
    {
        cout << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(s)", tagname.c_str());

    PyObject* nodes;
    nodes = PyEval_CallObject(pfunc, pargs);
    if (nodes == NULL)
    {
        cout << "null tuple returned" << endl;
        abort();
    }

    int size = PyList_Size(nodes);
    for (int i=0; i < size; ++i)
    {
        children.push_back(new PyXMLDomParser(PyList_GetItem(nodes, i)));
    }

    Py_DECREF(pargs);
    Py_DECREF(nodes);
    Py_DECREF(pfunc);
}

PyXMLDomParser::~PyXMLDomParser()
{
    Py_DECREF(parser_);
}


PyXMLDomParserPtr
PyXMLDomParser::createNode(refxmlstr name, refxmlstr nspace)
{
    if (nspace == "") //no namespae specified
        return _createNode(name);

    PyObject *pargs = NULL, *pfunc = NULL, *preturn = NULL;
    pfunc = PyObject_GetAttrString(parser_, "createNodeInNamespace");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(ss)", nspace.c_str(), name.c_str());

    preturn = PyEval_CallObject(pfunc, pargs);

    PyObject* val;
    PyArg_Parse(preturn, "O", &val);

    PyXMLDomParserPtr writer(new PyXMLDomParser(val));

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);
    return writer;
}

PyXMLDomParserPtr
PyXMLDomParser::_createNode(refxmlstr name)
{
    PyObject *pargs = NULL, *pfunc = NULL, *preturn = NULL;

    pfunc = PyObject_GetAttrString(parser_, "createNode");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(s)", name.c_str());

    preturn = PyEval_CallObject(pfunc, pargs);

    PyObject* val;

    PyArg_Parse(preturn, "O", &val);


    //this must go first to increment refcount
    PyXMLDomParserPtr writer(new PyXMLDomParser(val));
    //Py_DECREF(preturn);
    Py_INCREF(val);
    //abort();

    Py_DECREF(pargs);
    Py_DECREF(pfunc);

    return writer;
}

void
PyXMLDomParser::setAttribute(long val, refxmlstr name)
{
    setAttribute(stream_printf("%ld", val), name);
}

void
PyXMLDomParser::setAttribute(int val, refxmlstr name)
{
    setAttribute(stream_printf("%d", val), name);
}

void
PyXMLDomParser::setAttribute(unsigned int val, refxmlstr name)
{
    setAttribute(stream_printf("%d", val), name);
}

void
PyXMLDomParser::setAttribute(bool val, refxmlstr name)
{
    if (val)
        setAttribute("true", name);
    else
        setAttribute("false", name);
}

void
PyXMLDomParser::setAttribute(size_t val, refxmlstr name)
{
    setAttribute(stream_printf("%ld", val), name);
}

void
PyXMLDomParser::setAttribute(double val, refxmlstr name)
{
    setAttribute(stream_printf("%20.14f", val), name);
}

void
PyXMLDomParser::setAttribute(const char* val, refxmlstr name)
{
    std::string str(val);
    setAttribute(str, name);
}

void
PyXMLDomParser::setAttribute(refxmlstr value, refxmlstr name)
{
    PyObject* pfunc, *pargs;

    pfunc = PyObject_GetAttrString(parser_, "setAttribute");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(ss)", name.c_str(), value.c_str());

    PyEval_CallObject(pfunc, pargs);

    Py_DECREF(pargs);
    Py_DECREF(pfunc);
}

void
PyXMLDomParser::setText(refxmlstr text)
{
    PyObject* pfunc, *pargs;

    pfunc = PyObject_GetAttrString(parser_, "setText");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(s)", text.c_str());

    PyEval_CallObject(pfunc, pargs);

    Py_DECREF(pargs);
    Py_DECREF(pfunc);
}

void
PyXMLDomParser::addBinary(PyObject* obj)
{
    PyObject* pfunc = NULL, *pargs = NULL, *preturn = NULL;

    pfunc = PyObject_GetAttrString(parser_, "addBinary");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }

    pargs = Py_BuildValue("(O)", obj);

    preturn = PyEval_CallObject(pfunc, pargs);

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(refxmlstr value, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->setText(value);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(int val, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(long val, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(unsigned int val, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(unsigned long val, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(double val, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(bool val, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(refxmlstr value, refxmlstr name)
{
    XMLWriter node = createNode(name);
    node->setText(value);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(int val, refxmlstr name)
{
    XMLWriter node = createNode(name);
    node->writeValue(val);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeValue(double val, refxmlstr name)
{
    XMLWriter node = createNode(name);
    node->writeValue(val);
    return node;
}

void
PyXMLDomParser::writeValue(refxmlstr value)
{
    setText(value);
}

void
PyXMLDomParser::writeValue(double val)
{
    setText(stream_printf("%20.14f", val));
}

void
PyXMLDomParser::writeValue(long val)
{
    setText(stream_printf("%ld", val));
}

void
PyXMLDomParser::writeValue(unsigned int val)
{
    setText(stream_printf("%ld", val));
}

void
PyXMLDomParser::writeValue(unsigned long val)
{
    setText(stream_printf("%ld", val));
}

void
PyXMLDomParser::writeValue(bool val)
{
    if (val)
        setText("true");
    else
        setText("false");
}

void
PyXMLDomParser::writeValue(int val)
{
    setText(stream_printf("%d", val));
}

template <class T>
PyXMLDomParserPtr
PyXMLDomParser::writeList(
    const vector<T>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
)
{
    PyXMLDomParserPtr node = createNode(listname, nspace);

    if (values.size() == 0)
        return node;

    node->setAttribute(stream_printf("%d", values.size()), "n");

    for (int i=0; i < values.size(); ++i)
        node->writeValue(values[i], itemname);
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeList(
    const vector<std::string>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
)
{
    return writeList<std::string>(values, listname, itemname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::writeList(
    const vector<int>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
)
{
    return writeList<int>(values, listname, itemname, nspace);
}


PyXMLDomParserPtr
PyXMLDomParser::writeList(
    const vector<double>& values,
    refxmlstr listname,
    refxmlstr itemname,
    refxmlstr nspace
)
{
    return writeList<double>(values, listname, itemname, nspace);
}

template <class T, class U>
PyXMLDomParserPtr
PyXMLDomParser::writeMap(
    const std::map<T, U>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    XMLWriter node = createNode(tagname, nspace);
    typename map<T,U>::const_iterator it;
    for (it = vals.begin(); it != vals.end(); ++it)
    {
        XMLWriter pair_node = node->createNode("map_entry");
        pair_node->setAttribute(it->first, "key");
        pair_node->writeValue(it->second);
    }
    return node;
}

PyXMLDomParserPtr
PyXMLDomParser::writeMap(
    const std::map<int, int>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    return writeMap<int,int>(vals, tagname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::writeMap(
    const std::map<int, unsigned long>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    return writeMap<int,unsigned long>(vals, tagname, nspace);
}

PyXMLDomParserPtr
PyXMLDomParser::writeMap(
    const std::map<std::string, int>& vals,
    refxmlstr tagname,
    refxmlstr nspace
)
{
    return writeMap<std::string,int>(vals, tagname, nspace);
}

void
PyXMLDomParser::toFile(refxmlstr filename) const
{
    PyObject* pfunc, *pargs;

    pfunc = PyObject_GetAttrString(parser_, "toFile");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("(s)", filename.c_str());

    PyEval_CallObject(pfunc, pargs);

    Py_DECREF(pargs);
    Py_DECREF(pfunc);
}

XMLWriter
PyXMLDomParser::writeBinary(const void* vals, size_t size, refxmlstr name, refxmlstr nspace)
{
    XMLWriter node = createNode(name, nspace);
    node->setAttribute(size, "size");

    PyObject* obj = PyString_FromStringAndSize((const char*) vals, size);
    node->addBinary(obj);
    return node;
}


void
PyXMLDomParser::getBinary(void** vals, size_t& size) const
{
    getAttribute(size, "size");

    char* data = new char[size];
    (*vals) = (void *) data;

    PyObject* pfunc = NULL, *pargs = NULL, *preturn = NULL;
    pfunc = PyObject_GetAttrString(parser_, "getBinary");
    if (pfunc == NULL)
    {
        cerr << "null function returned" << endl;
        abort();
    }
    pargs = Py_BuildValue("()");
    preturn = PyEval_CallObject(pfunc, pargs); 

    int pysize;
    char* tmp;
    PyArg_Parse(preturn, "s#", &tmp, &pysize);

    memcpy(data, tmp, size);

#if 0
    Py_buffer buf;
    int check = PyObject_GetBuffer(obj, &buf, PyBUF_C_CONTIGUOUS);
    if (check != 0)
    {
        cerr << "buffer build failed" << endl;
        abort();
    }

    memcpy(data, buf.buf, size);
    PyBuffer_Release(&buf);
#endif

    Py_DECREF(pargs);
    Py_DECREF(preturn);
    Py_DECREF(pfunc);

}

XMLParser
PyXMLDomParser::loadBinary(void** vals, size_t& size, refxmlstr name, refxmlstr nspace) const
{
    XMLWriter node = fetchNode(name, nspace);
    if (node.get() == NULL)
        return NULL;

    node->getBinary(vals, size);
    return node;
}

