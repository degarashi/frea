# coding: utf-8
# GDB用のpretty-printer
import math
import lubee.printer
import re
import gdb

def buildPrinter(obj):
    obj.pretty_printers.append(Lookup)
    lubee.printer.buildPrinter(obj)

Re_Frea = re.compile("^frea::(.+)$")
def Lookup(val):
    name = val.type.strip_typedefs().name
    if name == None:
        return None
    obj = Re_Frea.match(name)
    if obj:
        return LookupFrea(obj.group(1), val)
    return None

Re_Vector = re.compile("^Vec.+$")
Re_Matrix = re.compile("^Mat.+$")
Re_VectorWrap = re.compile("^wrap[^M].+$")
Re_MatrixWrap = re.compile("^wrapM.+$")
Re_Quat = re.compile("^Quat.+$")
Re_EQuat = re.compile("^ExpQuat.+$")
Re_Plane = re.compile("^Plane.+$")
Re_TupleVec = re.compile("^tup.+$")
def LookupFrea(name, val):
    if Re_Vector.match(name) or Re_VectorWrap.match(name):
        return Vector(val)
    if Re_TupleVec.match(name):
        return TupleVec(val)
    if Re_Matrix.match(name):
        return Matrix(val, "m")
    if Re_MatrixWrap.match(name):
        return Matrix(val, "v")
    if Re_Quat.match(name):
        return Quat(val)
    if Re_EQuat.match(name):
        return ExpQuat(val)
    if Re_Plane.match(name):
        return Plane(val)
    return None

class Iterator:
    def next(self):
        return self.__next__()

ElementName = ["x", "y", "z", "w"]
ElementIndex = []
for idx in range(64):
    ElementIndex.append("[" + str(idx) + "]")
def GetElementIndex(size):
    u''' 要素数が4つまでのベクトルはXYZW表記、それ以上の場合は数値インデックスのリストを返す '''
    return ElementName if size<=len(ElementName) else ElementIndex

def Embrace(s):
    return "[" + s + "]"
def MakeTag(align, integral, name, size0=0, size1=0):
    return "%s%s%s%s%s" % (
        "A" if align else "",
        "I" if integral else "",
        name,
        "" if size0==0 else str(size0),
        "" if size1==0 else str(size1)
    );

class Vector:
    def __init__(self, val, size=-1):
        self._ar = val["m"]
        try:
            self._ar = self._ar["m"]
        except gdb.error:
            pass
        self._size = size if size>=0 else int(val["size"])
        try:
            self._align = bool(val["align"])
        except:
            self._align = True

        code = self._ar[0].type.strip_typedefs().code
        self._integral = code == gdb.TYPE_CODE_INT
        if self._integral:
            self._makeEnt = self._makeIntEntry
        else:
            self._makeEnt = self._makeFloatEntry

    def _makeIntEntry(self, val):
        return "%d" % (int(val))
    def _makeFloatEntry(self, val):
        return "%.5g" % (float(val))
    def _raw_string(self):
        a = []
        for i in range(self._size):
            a.append(self._makeEnt(self._ar[i]))
        return Embrace(", ".join(a))
    def _makeTag(self):
        return MakeTag(self._align, self._integral, "Vec", self._size)
    def to_string(self):
        return self._makeTag() + self._raw_string()
    def isIntegral(self):
        return self._integral

class TupleVec:
    class _iterator(Iterator):
        def __init__(self, src):
            self._src = src
            self._cur = 0
            self._bEnd = False
        def __iter__(self):
            return self
        def __next__(self):
            src = self._src
            cur = self._cur
            self._cur += 1
            if cur >= self._src._numChild:
                if self._modChild==0 or self._bEnd:
                    raise StopIteration
                self._bEnd = True
                return ("Tail(%d)" % self._src._modChild, self._src._ar[cur])
            return ("[%d]" % cur, self._src._ar[cur])

    def __init__(self, val):
        self._ar = val["data"]
        self._wcap = int(val["w_capacity"])
        self._size = int(val["size"])
        self._numChild = int(self._size / self._wcap)
        self._modChild = self._size % self._wcap
        if self._modChild != 0:
            self._numChild += 1

    def to_string(self):
        return "(tuple)"
    def children(self):
        return self._iterator(self)
    def display_hint(self):
        return "array"

class Matrix:
    class _iterator(Iterator):
        def __init__(self, ar, size):
            self._ar = ar
            self._size = size
            self._cur  = 0
        def __iter__(self):
            return self
        def __next__(self):
            if self._cur == self._size:
                raise StopIteration
            cur = self._cur
            self._cur += 1
            return ("[%d]" % cur, self._ar[cur])

    def __init__(self, val, ar_name):
        self._ar = val[ar_name]
        self._width = val["dim_n"]
        self._height = val["dim_m"]
        self._integral = Vector(self._ar[0]).isIntegral()
        try:
            self._align = bool(val["align"])
        except gdb.error:
            self._align = True

    def children(self):
        return self._iterator(self._ar, self._height)
    def to_string(self):
        return MakeTag(self._align, self._integral, "Mat", self._height, self._width)
    def display_hint(self):
        return "array"

class Quat(Vector, object):
    def __init__(self, val):
        super(Quat, self).__init__(val)
    def to_string(self):
        ar = self._ar
        ar3 = float(ar[3])
        s_theta = math.sqrt(1.0 - min(1.0, ar3*ar3))
        if s_theta < 1e-6:
            # 無回転クォータニオンとして扱う
            axis = [0,0,0]
            angle = 0
        else:
            s_theta = 1.0 / s_theta
            axis = []
            for i in range(3):
                axis.append(float(ar[i]) * s_theta)
            angle = math.acos(float(ar[3]))*2 / math.pi * 180
        return MakeTag(self._align, self._integral, "Quat", 0, 0)\
                + "{axis=[%.3f, %.3f, %.3f], deg=%.3f}" % (axis[0], axis[1], axis[2], angle)

class ExpQuat(Vector, object):
    def __init__(self, val):
        super(ExpQuat, self).__init__(val)
    def to_string(self):
        ar = self._ar
        theta = math.sqrt(
            float(ar[0])*float(ar[0])
            + float(ar[1])*float(ar[1])
            + float(ar[2])*float(ar[2])
        )
        if theta < 1e-6:
            # 無回転クォータニオンとして扱う
            axis = [0,0,0]
            angle = 0
        else:
            axis = []
            for i in range(3):
                axis.append(float(ar[i]) / theta)
            angle = theta*2 / math.pi * 180
        return MakeTag(self._align, self._integral, "EQuat", 0, 0)\
                + "{axis=[%.3f, %.3f, %.3f], deg=%.3f}" % (axis[0], axis[1], axis[2], angle)

class Plane(Vector, object):
    def __init__(self, val):
        super(Plane, self).__init__(val)
    def to_string(self):
        ar = self._ar
        axis = [float(ar[0]), float(ar[1]), float(ar[2])]
        angle = float(ar[3])
        return MakeTag(self._align, self._integral, "Plane", 0, 0)\
                + "{axis=[%.3f, %.3f, %.3f], dist=%.3f}" % (axis[0], axis[1], axis[2], angle)

