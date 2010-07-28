
text = open("deftypes.h").read()

str_arr = []
str_arr.extend(text.splitlines())

for line in text.splitlines():
    if not "intrusive_ptr" in line:
        continue

    newline = line.replace("define ", "define Const").replace("ptr<", "ptr<const ")
    str_arr.append(newline)

fileobj = open("deftypes.h", "w")
fileobj.write("\n".join(str_arr))
fileobj.close()

