# coding: utf-8
# QtCreator Debugger用のpretty-printer
from dumper import *

Fr_ElementName = ["x", "y", "z", "w"]
Fr_ElementIndex = []
for idx in range(64):
    Fr_ElementIndex.append("[" + str(idx) + "]")

def fr_GetElementIndex(size):
    u''' 要素数が4つまでのベクトルはXYZW表記、それ以上の場合は数値インデックスのリストを返す '''
    return Fr_ElementName if size<=len(Fr_ElementName) else Fr_ElementIndex

def fr_TypePrefix(b_integral, bitwidth, b_align):
    # Alignedなら先頭にA(aligned)を付与
    typv = "A" if b_align else ""
    if b_integral:
        # 整数型にはI(integral)を付与
        typv += "I"
    else:
        # 内部表現に使うサイズが32bitを超える場合はD(double)を付与
        typv += "" if bitwidth <= 32 else "D"
    return typv

def fr_ArrayString(elem, size, typ):
    s = ""
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
    return s
def fr_ArrayStringB(elem, size, typ):
    return Embrace(fr_ArrayString(elem, size, typ))
def Embrace(s):
    return "[" + s + "]"

def fr_DumpVector(d, elem, size, typ):
    d.putValue(fr_ArrayStringB(elem, size, typ))
    d.putNumChild(size)

    if d.isExpanded():
        index = fr_GetElementIndex(size)
        with Children(d):
            for idx in range(size):
                d.putSubItem(index[idx], elem[idx])

# ------------------------------------------------------------
# Vector
def fr_DumpData(d, value, name):
    bInt = bool(value["is_integral"])
    d.putType(
        fr_TypePrefix(
            bInt,
            value["bit_width"],
            value["align"]
        )
        + name
        + str(value["size"])
    )
    fr_DumpVector(d, value["m"], int(value["size"]), int if bInt else float)

def qdump__frea__Data(d, value):
    fr_DumpData(d, value, "Data")
def qdump__frea__VecT(d, value):
    fr_DumpData(d, value, "Vec")
def qdump__frea__VecT_spec(d, value):
    fr_DumpData(d, value, "Vec(s)")

# ------------------------------------------------------------
# Vector wrapper
def fr_DumpWrap(d, value, name):
    elem = value["m"]
    size = int(value["size"])
    bw = int(value["bit_width"])
    bInt = bool(value["is_integral"])
    d.putType(fr_TypePrefix(bInt, bw, False) + name + str(size))
    if bInt:
        d.putArrayData(elem.address, size, d.lookupType("unsigned int"))
    else:
        d.putValue(fr_ArrayStringB(elem, size, float))
        d.putNumChild(size)
        if d.isExpanded():
            index = fr_GetElementIndex(size)
            with Children(d):
                for idx in range(size):
                    d.putSubItem(index[idx], elem[idx])

def qdump__frea__wrap(d, value):
    fr_DumpWrap(d, value, "Wrap")
def qdump__frea__wrap_spec(d, value):
    fr_DumpWrap(d, value, "Wrap(s)")

def fr_TupString(value):
    data = value["data"]
    wcap = int(value["w_capacity"])
    bInt = bool(value["is_integral"])
    size = int(value["size"])
    numChild = int(size/wcap)
    modChild = size%wcap
    if modChild != 0:
        numChild += 1
    s = "["
    if not bInt:
        comma = False
        for i in range(numChild-1):
            if comma:
                s += ", "
            s += fr_ArrayString(data[i]["m"], wcap, float)
            comma = True
        s += ", "
        s += fr_ArrayString(
                data[numChild-1]["m"],
                modChild if modChild>0 else wcap,
                float
             )
    s += "]"
    return s

def fr_DumpTup(d, value, name):
    bw = int(value["bit_width"])
    data = value["data"]
    wcap = int(value["w_capacity"])
    bInt = bool(value["is_integral"])
    size = int(value["size"])
    numChild = int(size/wcap)
    modChild = size%wcap
    if modChild != 0:
        numChild += 1
    # 整数型は今の所未対応
    if not bInt:
        d.putValue(fr_TupString(value))
        d.putType(fr_TypePrefix(bInt, bw, False) + name + str(size))
    d.putNumChild(numChild)
    if d.isExpanded():
        with Children(d):
            for idx in range(numChild):
                d.putSubItem(idx, data[idx])

def qdump__frea__tup(d, value):
    fr_DumpTup(d, value, "Wrap")
def qdump__frea__tup_spec(d, value):
    fr_DumpTup(d, value, "Wrap(s)")

# ------------------------------------------------------------
# Matrix wrapper
def fr_DumpWrapM(d, value, name):
    m = int(value["dim_m"])
    n = int(value["dim_n"])
    bInt = bool(value["is_integral"])
    elem = value["v"]

    # 整数型は今の所未対応
    if not bInt:
        s = "["
        comma = False
        for idx in range(m):
            if comma:
                s += ", "
            comma = True
            el = elem[idx]
            try:
                s += fr_ArrayStringB(el["m"], n, None if bInt else float)
            except RuntimeError:
                s += fr_TupString(el)
        s += "]"
        d.putValue(s)
        d.putType(name + "[" + str(m) + "][" + str(n) + "]")
    d.putNumChild(m)
    if d.isExpanded():
        with Children(d):
            for idx in range(m):
                d.putSubItem(idx, elem[idx])

def qdump__frea__wrapM(d, value):
    fr_DumpWrapM(d, value, "WrapM")
def qdump__frea__wrapM_spec(d, value):
    fr_DumpWrapM(d, value, "WrapM(s)")

# ------------------------------------------------------------
# Matrix
def fr_DumpMData(d, value, name):
    m = int(value["dim_m"])
    n = int(value["dim_n"])
    bInt = bool(value["is_integral"])
    bw = int(value["bit_width"])
    align = bool(value["align"])
    elem = value["m"]
    d.putType(
        fr_TypePrefix(bInt, bw, align)
        + name + "[" + str(m) + "][" + str(n) + "]"
    )
    s = "["
    comma = False
    for idx in range(m):
        if comma:
            s += ", "
        comma = True
        s += fr_ArrayStringB(elem[idx]["m"], n, None if bInt else float)
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
# Range
def qdump__frea__Range(d, value):
    d.putValue("[" + str(value["from"]) + " | " + str(value["to"]) + "]")
    d.putType("Range")
    d.putNumChild(0)
