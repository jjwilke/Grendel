rep_list = {
    "RefSCMatrix" : "RectMatrixPtr",
    "RefSymmSCMatrix" : "SymmMatrixPtr",
    "RefSCVector" : "VectorPtr",
}

import glob, sys
headers = glob.glob("*.h")
source = glob.glob("*.cc")
filelist = headers + source

for file in filelist:
    if file in ['deftypes.h']:
        continue

    text = open(file).read()
    for entry in rep_list:
        text = text.replace(entry, rep_list[entry])

    fileobj = open(file, 'w')
    fileobj.write(text)
    fileobj.close()



