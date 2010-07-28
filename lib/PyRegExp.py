IgnoreCaseShift = 0
CapitalizeShift = 1
FindAllShift = 2
LowerCaseShift = 3
IncludeLineBreakShift = 4
StripWhitespace = 5

def setpath(path):
    import sys
    if not path in sys.path:
        sys.path.append(path)

def setup_regexp(text, flags):
    try:
        import re
        lowercase = (flags >> LowerCaseShift) % 2
        capitalize = (flags >> CapitalizeShift) % 2
        ignorecase = (flags >> IgnoreCaseShift) % 2
        includelinebreak = (flags >> IncludeLineBreakShift) % 2
        stripwhitespace = (flags >> StripWhitespace) % 2

        reflags = 0
        if ignorecase: reflags = reflags | re.IGNORECASE
        if includelinebreak: reflags = reflags | re.DOTALL

        searchtext = text
        if      lowercase: searchtext = text.lower()
        elif    capitalize: searchtext = text.upper()

        if    stripwhitespace: searchtext = searchtext.strip()

        return reflags, searchtext

    except Exception, error:
        print "Error in regexp setup", error

def get_num_matches_from_regexp(regexp, text, flags):
    try:
        matches = get_all_from_regexp(regexp, text, flags)
        return len(matches)
    except Exception, error:
        print "Error in num matches:", error

def has_match(regexp, text, flags):
    try:
        import re
        reflags, searchtext = setup_regexp(text, flags)
        expObj = re.compile(regexp, reflags)
        match = expObj.search(searchtext)
        if match: return 1
        else: return 0
    except Exception, error:
        print "Error in has match: ", error

def get_value_from_regexp(regexp, text, flags = 0):
    try:
        reflags, searchtext = setup_regexp(text, flags)
        import re
        expObj = re.compile(regexp, reflags)
        match = expObj.search(searchtext)
        val = match.groups()
        return val
    except Exception, error:
        print "Error in get value:", error
        print "Regular expression:", regexp
        print "Text:\n", text

def get_all_from_regexp(regexp, text, flags = 0):
    try:
        reflags, searchtext = setup_regexp(text, flags)
        import re
        expObj = re.compile(regexp, reflags)
        match = expObj.findall(searchtext)
        return match
    except Exception, error:
        print "Error in get all:", error
        print "Regular expression:", regexp
        print "Text:\n", text

def get_double_from_regexp(regexp, text, flags = 0):
    try:
        val = get_value_from_regexp(regexp, text, flags)[0].strip().strip("+")
        return float(val)
    except Exception, error:
        print "Error in get double:", error
    
def get_int_from_regexp(regexp, text, flags = 0):
    try:
        val = int(get_value_from_regexp(regexp, text, flags)[0])
        return val
    except Exception, error:
        print "Error in get int:", error

def get_string_from_regexp(regexp, text, flags = 0):
    try:
        val  = get_value_from_regexp(regexp, text, flags)[0]
        return val
    except Exception, error:
        print "Error in get string:", error

def get_double_array_from_regexp(regexp, text, flags = 0):
    try:
        dbl_arr = []
        findall = (flags >> FindAllShift) % 2
        if findall:
            #this will be a 2D array
            allvals = get_all_from_regexp(regexp, text, flags)
            for i in range(0, len(allvals)):
                match = allvals[i]
                if isinstance(match, list) or isinstance(match, tuple):
                    for value in match:
                        dbl_arr.append( eval(value) )
                else:
                    value = eval(match) 
                    dbl_arr.append(value)

        else: #simple match
            vals = get_value_from_regexp(regexp, text, flags)
            for value in vals:
                dbl_arr.append( eval(value) )

        return dbl_arr;
    except Exception, error:
        print "Error in get double array:", error

def get_int_array_from_regexp(regexp, text, flags = 0):

    try:
        int_arr = []
        findall = (flags >> FindAllShift) % 2
        if findall:
            #this will be a 2D array
            allvals = get_all_from_regexp(regexp, text, flags)
            for i in range(0, len(allvals)):
                match = allvals[i]
                if isinstance(match, list) or isinstance(match, tuple):
                    for value in match:
                        int_arr.append( eval(value) )
                else:
                    value = int(match) 
                    int_arr.append(value)

        else: #simple match
            vals = get_value_from_regexp(regexp, text, flags)
            for value in vals:
                int_arr.append( eval(value) )

        return int_arr;

    except Exception, error:
        print "Error in get int array:", error

def build_match_from_regexp(regexp, text, pydict, flags = 0):
    try:
        match = get_value_from_regexp(regexp, text, flags)
        for i in range(0, len(match)):
            pydict["%d" % i] = match[i]
    except Exception, error:
        print "Error in build match:", error

def build_findall_from_regexp(regexp, text, pydict, flags = 0):
    try:
        allvals = get_all_from_regexp(regexp, text, flags)
        matchnum=0 ; groupnum = 0
        for match in allvals:
            if isinstance(match, tuple):
                groupnum = 0
                for group in match:
                    key = "%d,%d" % (matchnum, groupnum)
                    pydict[key] = allvals[matchnum][groupnum]
                    groupnum += 1
                matchnum += 1
            else:
                key = "%d,0" % matchnum
                pydict[key] = match
                matchnum += 1

        #return the number of matches
        return len(allvals)
    except Exception, error:
        print "Error in find all:", error
        


