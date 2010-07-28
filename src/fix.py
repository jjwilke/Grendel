
import glob, sys
headers = glob.glob("*.h")
source = glob.glob("*.cc")
filelist = headers + source

for file in filelist:
    text = open(file).read().splitlines()
    newtext = [""]
    for line in text:
        if line == newtext[-1]:
            continue #repeat line, don't add
        newtext.append(line)
    fileobj = open(file, 'w')
    fileobj.write("\n".join(newtext))
    fileobj.close()
