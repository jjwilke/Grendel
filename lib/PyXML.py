
teststr = "hello"

def traceback(error=None):
    import sys
    import traceback
    if error:
        tb_list = traceback.format_tb(sys.exc_info()[2])
        return "%s\n%s" % (error, "\n".join(tb_list))
    else:
        import StringIO
        output = StringIO.StringIO()
        traceback.print_stack(file=output)
        contents = output.getvalue()
        return contents


class XMLParser(object):

    def __init__(self, filename = None, node = None):
        try:
            import xml.dom.minidom
            if filename:
                self.parser = xml.dom.minidom.parse(filename)
                self.firstChild = XMLParser(node=self.parser.firstChild)
            elif isinstance(node, xml.dom.minidom.Node):
                self.parser = node
            else:
                impl = xml.dom.getDOMImplementation()
                doc = impl.createDocument("gigide", "document", None)
                self.parser = doc
                self.firstChild = XMLParser(node=self.parser.firstChild)

        except Exception, error:
            import sys
            import time
            time.sleep(10)
            print error
            print traceback(error)
            sys.stderr.write("%s\n" % traceback(error))

    def getTopnode(self):
        try:
            return XMLParser(node=self.parser.ownerDocument.firstChild)
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def toFile(self, filename):
        try:
            import codecs
            text = self.parser.toxml()
            fileobj = codecs.open(filename, "w", "utf-8")
            fileobj.write(text)
            fileobj.close()
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def setAttribute(self, name, value):
        try:
            self.parser.setAttribute(name, value)
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % error)

    def getNextSibling(self, name):
        parser = None
        check = None
        node = self.parser
        try:
            while node:
                node = node.nextSibling
                if not node:
                    return None
                check = node.nodeName
                if node.__class__ == self.parser.__class__ and name == check: #this is an element
                    parser = XMLParser(node=node)
                    return parser
            return None
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))
            return None

    def createNode(self, name):
        try:
            import xml.dom.minidom
            node = self.parser.ownerDocument.createElement(name)
            self.parser.appendChild(node)
            return XMLParser(node=node)
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def createNodeInNamespace(self, namespace, name):
        try:
            import xml.dom.minidom
            node = self.parser.ownerDocument.createElementNS(namespace, name)
            #python is retarded and won't generate namespaces
            node.setAttribute("xmlns", node.namespaceURI)
            self.parser.appendChild(node)
            return XMLParser(node=node)
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def addBinary(self, data):
        try:
            import base64
            ascii = base64.encodestring(data)
            self.setText(ascii)
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def setText(self, text):
        try:
            textNode = self.parser.ownerDocument.createTextNode("text")
            #clear everything
            self.parser.childNodes = []
            self.parser.appendChild(textNode)
            textNode.nodeValue = text
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def getElementsByTagName(self, tag):
        try:
            nodes = self.parser.getElementsByTagName(tag)

            parsers = []
            for node in nodes:
                parsers.append(XMLParser(node =  node))

            return parsers
        except Exception, error:
            return []

    def getElementsByTagNameInNamespace(self, namespace, tag):
        try:
            nodes = self.parser.getElementsByTagNameNS(namespace, tag)

            parsers = []
            for node in nodes:
                parsers.append(XMLParser(node =  node))

            return parsers
        except Exception, error:
            return []

    def getBinary(self):
        try:
            import base64
            text = self.getText()
            data = base64.decodestring(text)
            return data
        except Exception, error:
            return ""

    def getText(self):
        try:
            node = self.parser.firstChild
            if node.nodeType == node.TEXT_NODE:
                return node.data
            else:
                return ""
        except Exception, error:
            return ""

    def hasAttribute(self, attrname):
        try:
            return self.parser.hasAttribute(attrname)
        except Exception, error:
            return False

    def getAttribute(self, attrname):
        try:
            attr = self.parser.getAttribute(attrname)
            return str(attr) #no unicode for now, just in case
        except Exception, error:
            return ""

    def getUnicodeString(self):
        x = u"\xf1\xf2\xf3\xf4"
        return x

    def tostr(self):
        try:
            return self.parser.toprettyxml()
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

    def testUnicode(self, uni):
        try:
            print "testing unicode"
            print uni.__class__
            print len(uni)
            for i in range(len(uni)):
                print i,
                print uni[i].encode("utf-8")
            print uni.encode("utf-8")
        except Exception, error:
            import sys
            sys.stderr.write("%s\n" % traceback(error))

def construct_parser_from_node(node):
    try:
        parser = XMLParser(node = node)
        return parser
    except Exception, error:
        import sys
        sys.stderr.write("%s\n" % traceback(error))
        return None #send back null

def construct_writer():
    try:
        writer = XMLParser()
        return writer
    except Exception, error:
        import sys
        sys.stderr.write("%s\n" % traceback(error))
        return None

def construct_parser(filename):
    try:
        parser = XMLParser(filename = filename)
        return parser
    except Exception, error:
        import sys
        sys.stderr.write("%s\n" % traceback(error))
        return None #send back null

if __name__ == "__main__":
    writer = construct_writer()


