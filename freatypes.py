# coding: utf-8
# QtCreator Debugger用のpretty-printer
from dumper import *
import math
from lubee.lubeetypes import *

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
            s += "{0:4.6f}".format(e.value())
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

def fr_AxisString(axis):
    fmt = "{0:3.3f}"
    str_axis = []
    for a in axis:
        str_axis.append(fmt.format(a))
    s = "Axis=[" + ",".join(str_axis) + "]"
    return s
def fr_RotationString(axis, angle):
    fmt = "{0:3.3f}"
    s = fr_AxisString(axis) + " Deg=[" + fmt.format(angle) + "]"
    return s
def fr_DumpChildren(d, node, size, indexIsNumber=False):
    d.putNumChild(size)
    if d.isExpanded():
        if indexIsNumber:
            index = Fr_ElementIndex
        else:
            index = fr_GetElementIndex(size)
        with Children(d):
            for idx in range(size):
                d.putSubItem(index[idx], node[idx])

def fr_DumpVector(d, elem, size, typ):
    d.putValue(fr_ArrayStringB(elem, size, typ))
    fr_DumpChildren(d, elem, size)

# ------------------------------------------------------------
# Vector
def fr_DumpData(d, value, name):
    nv = value.native
    bInt = bool(nv["is_integral"])
    d.putType(
        fr_TypePrefix(
            bInt,
            nv["bit_width"],
            nv["align"]
        )
        + name
        + str(nv["size"])
    )
    fr_DumpVector(d, value["m"], int(nv["size"]), int if bInt else float)

def qdump__frea__Data(d, value):
    fr_DumpData(d, value, "Data")
def qdump__frea__VecT(d, value):
    fr_DumpData(d, value, "Vec")
def qdump__frea__VecT_spec(d, value):
    fr_DumpData(d, value, "Vec(s)")

# ------------------------------------------------------------
# Vector wrapper
def fr_DumpWrap(d, value, name):
    nv = value.native
    elem = value["m"]
    size = int(nv["size"])
    bw = int(nv["bit_width"])
    bInt = bool(nv["is_integral"])
    d.putType(fr_TypePrefix(bInt, bw, False) + name + str(size))
    if bInt:
        d.putArrayData(elem.address, size, d.lookupType("unsigned int"))
    else:
        d.putValue(fr_ArrayStringB(elem, size, float))
        fr_DumpChildren(d, elem, size)

def qdump__frea__wrap(d, value):
    fr_DumpWrap(d, value, "Wrap")
def qdump__frea__wrap_spec(d, value):
    fr_DumpWrap(d, value, "Wrap(s)")

def fr_TupString(value):
    nv = value.native
    data = value["data"]
    wcap = int(nv["w_capacity"])
    bInt = bool(nv["is_integral"])
    size = int(nv["size"])
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
    nv = value.native
    bw = int(nv["bit_width"])
    data = value["data"]
    wcap = int(nv["w_capacity"])
    bInt = bool(nv["is_integral"])
    size = int(nv["size"])
    numChild = int(size/wcap)
    modChild = size%wcap
    if modChild != 0:
        numChild += 1
    # 整数型は今の所未対応
    if not bInt:
        d.putValue(fr_TupString(value))
        d.putType(fr_TypePrefix(bInt, bw, False) + name + str(size))
    fr_DumpChildren(d, data, numChild, True)

def qdump__frea__tup(d, value):
    fr_DumpTup(d, value, "Wrap")
def qdump__frea__tup_spec(d, value):
    fr_DumpTup(d, value, "Wrap(s)")

# ------------------------------------------------------------
# Matrix wrapper
def fr_DumpWrapM(d, value, name):
    nv = value.native
    m = int(nv["dim_m"])
    n = int(nv["dim_n"])
    bInt = bool(nv["is_integral"])
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
    fr_DumpChildren(d, elem, m, True)

def qdump__frea__wrapM(d, value):
    fr_DumpWrapM(d, value, "WrapM")
def qdump__frea__wrapM_spec(d, value):
    fr_DumpWrapM(d, value, "WrapM(s)")

# ------------------------------------------------------------
# Matrix
def fr_DumpMData(d, value, name):
    nv = value.native
    m = int(nv["dim_m"])
    n = int(nv["dim_n"])
    bInt = bool(nv["is_integral"])
    bw = int(nv["bit_width"])
    align = bool(nv["align"])
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
    fr_DumpChildren(d, elem, m, True)

def qdump__frea__DataM(d, value):
    fr_DumpMData(d, value, "DataM")
def qdump__frea__MatT_spec(d, value):
    fr_DumpMData(d, value, "MatT(s)")
def qdump__frea__MatT(d, value):
    fr_DumpMData(d, value, "MatT")

# ------------------------------------------------------------
# Quaternion
def qdump__frea__QuatT(d, value):
    val = value["m"]
    s_theta = math.sqrt(1.0 - val[3].value()*val[3].value())
    if s_theta < 1e-6:
        # 無回転クォータニオンとして扱う
        axis = [0,0,0]
        angle = 0
    else:
        s_theta = 1.0 / s_theta
        axis = []
        for i in range(3):
            axis.append(float(val[i].value() * s_theta))
        angle = math.acos(val[3].value())*2 / math.pi * 180

    d.putValue("Quat:" + fr_RotationString(axis, angle))
    nv = value.native
    d.putType(fr_TypePrefix(nv["is_integral"], nv["bit_width"], nv["align"]) + "Quat")
    fr_DumpChildren(d, val, 4)

# ------------------------------------------------------------
# ExpQuaternion
def qdump__frea__ExpQuatT(d, value):
    val = value["m"]
    theta = math.sqrt(val[0].value()*val[0].value() + val[1].value()*val[1].value() + val[2].value()*val[2].value())
    if theta < 1e-6:
        # 無回転クォータニオンとして扱う
        axis = [0,0,0]
        angle = 0
    else:
        axis = []
        for i in range(3):
            axis.append(float(val[i].value() / theta))
        angle = theta*2 / math.pi * 180
    d.putValue("ExpQuat: " + fr_RotationString(axis, angle))
    fr_DumpChildren(d, val, 3)

# ------------------------------------------------------------
# Plane
def qdump__frea__PlaneT(d, value):
    val = value["m"]
    axis = [val[0].value(), val[1].value(), val[2].value()]
    fmt = "{0:3.3f}"
    d.putValue("Plane: " + fr_AxisString(axis) + " D=" + fmt.format(val[3].value()))
    fr_DumpChildren(d, val, 4)
