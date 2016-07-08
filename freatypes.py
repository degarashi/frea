from dumper import *

Fr_ElementName = ["x", "y", "z", "w"]
Fr_ElementIndex = []
for idx in range(64):
    Fr_ElementIndex.append("[" + str(idx) + "]")
def Fr_GetElementIndex(size):
    return Fr_ElementName if size<=len(Fr_ElementName) else Fr_ElementIndex

def fr_DumpType(b_integral, bitwidth, b_align):
    typv = "A" if b_align else ""
    if b_integral:
        typv += "I"
    else:
        typv += "" if bitwidth <= 32 else "D"
    return typv
def fr_DumpValues(elem, size, typ):
    s = "["
    comma = False
    for idx in range(size):
        if comma:
            s += ", "
        e = elem[idx]
        if typ == float:
            s += "{0:4.6f}".format(float(e))
        elif typ != None:
            s += str(typ(e))
        else:
            s += str(e)
        comma = True
    s += "]"
    return s

def fr_DumpVector(d, elem, size, typ):
    d.putValue(fr_DumpValues(elem, size, typ))
    d.putNumChild(size)

    if d.isExpanded():
        index = Fr_GetElementIndex(size)
        with Children(d):
            for idx in range(size):
                d.putSubItem(index[idx], elem[idx])

# ------------------------------------------------------------
def fr_DumpData(d, value, name):
    bInt = value["is_integral"]
    d.putType(fr_DumpType(
        bInt,
        value["bit_width"],
        value["align"]) + name + str(value["size"]))
    fr_DumpVector(d, value["m"], int(value["size"]), float if bInt==False else int)

def qdump__frea__Data(d, value):
    fr_DumpData(d, value, "Data")
def qdump__frea__VecT(d, value):
    fr_DumpData(d, value, "VecT")
def qdump__frea__VecT_spec(d, value):
    fr_DumpData(d, value, "VecT(s)")

# ------------------------------------------------------------
def fr_DumpWrap(d, value, name):
    elem = value["m"]
    size = value["size"]
    d.putType(name + str(size))
    if value["is_integral"] == True:
        d.putArrayData(elem.address, 4, d.lookupType("unsigned int"))
    else:
        d.putValue(fr_DumpValues(elem, size, float))
        d.putNumChild(size)
        if d.isExpanded():
            index = Fr_ElementName if size<=len(Fr_ElementName) else Fr_ElementIndex
            with Children(d):
                for idx in range(size):
                    d.putSubItem(index[idx], elem[idx])

def qdump__frea__wrap(d, value):
    fr_DumpWrap(d, value, "Wrap")

def qdump__frea__wrap_spec(d, value):
    fr_DumpWrap(d, value, "Wrap(s)")

# ------------------------------------------------------------
def fr_DumpWrapM(d, value, name):
    m = value["dim_m"]
    n = value["dim_n"]
    elem = value["v"]
    s = "["
    comma = False
    for idx in range(m):
        if comma:
            s += ", "
        comma = True
        s += fr_DumpValues(elem[idx]["m"], n, float if value["is_integral"]==False else None)
    s += "]"
    d.putValue(s)
    d.putType(name + "[" + str(m) + "][" + str(n) + "]")
    d.putNumChild(m)

    if d.isExpanded():
        index = Fr_ElementIndex
        with Children(d):
            for idx in range(m):
                d.putSubItem(index[idx], elem[idx])

def qdump__frea__wrapM(d, value):
    fr_DumpWrapM(d, value, "WrapM")
def qdump__frea__wrapM_spec(d, value):
    fr_DumpWrapM(d, value, "WrapM(s)")

# ------------------------------------------------------------
def fr_DumpMData(d, value, name):
    m = value["dim_m"]
    n = value["dim_n"]
    elem = value["m"]
    d.putType(name + "[" + str(m) + "][" + str(n) + "]")
    s = "["
    comma = False
    for idx in range(m):
        if comma:
            s += ", "
        comma = True
        s += fr_DumpValues(elem[idx]["m"], n, None)
    s += "]"
    d.putValue(s)
    d.putNumChild(m)
    if d. isExpanded():
        index = Fr_ElementIndex
        with Children(d):
            for idx in range(m):
                d.putSubItem(index[idx], elem[idx])

def qdump__frea__DataM(d, value):
    fr_DumpMData(d, value, "DataM")
def qdump__frea__MatT_spec(d, value):
    fr_DumpMData(d, value, "MatT(s)")
def qdump__frea__MatT(d, value):
    fr_DumpMData(d, value, "MatT")

# ------------------------------------------------------------
def qdump__frea__Range(d, value):
    d.putValue("[" + str(value["from"]) + " | " + str(value["to"]) + "]")
    d.putType("Range")
    d.putNumChild(0)
