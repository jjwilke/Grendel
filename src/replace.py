rep_list = [
    #"InternalCoordinate",
    #"SimpleInternalCoordinate",
    #"SymmetryInternalCoordinate",
    #"BondAngle",
    #"BongLength",
    #"Torsion",
    #"LinX",
    #"LinY",
    #"Lin1",
    #"Displacement",
    #"DisplacementMapping",
    #"DisplacementIterator",
    #"Derivative",
    #"DerivativeIterator",
    #"Fit",
    #"FitPoint",
    #"FitFactory",
    #"XYZPoint",
    #"Midpoint",
    #"Axis",
    #"EmpiricalHessianTerm",
    #"ConstantTerm",
    #"BondBond",
    #"BendBend",
    #"BondBend",
    #"TorsTors",
    #"InputFile",
    #"GigideInputFile",
    #"KeywordValue",
    #"KeywordSet",
    #"KeywordIterator",
    #"Atom",
    #"Molecule",
    #"PermutationGenerator",
    #"PyXMLDomParser",
    #"ForceField",
    #"Molecule",
    #"XYZPoint",
    #"SymmetryOperation",
    #"ImproperRotation",
    #"Rotation",
    #"Reflection",
    #"Inversion",
    #"IdentityElement",
    #"PointGroupClass",
    #"PointGroup",
    "TaylorTerm",
]

import glob, sys
headers = glob.glob("*.h")
source = glob.glob("*.cc")
filelist = headers + source

old_list = ["RefPtr<%s", "Ref<%s", "intrusive_ptr<%s", "%sPtr>"]

for file in filelist:
    if file in ['deftypes.h']:
        continue

    text = open(file).read().splitlines()
    newtext = []
    for line in text:
        if "typedef" in line:
            newtext.append(line) #skip it
            continue
        newline = line
        for entry in rep_list:
            for old in old_list:
                old = old % entry
                new = "%sPtr" % entry
                newline = newline.replace(old, new)
        newtext.append(newline)
    fileobj = open(file, 'w')
    fileobj.write("\n".join(newtext))
    fileobj.close()



