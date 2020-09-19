import numpy as np
from interpolant1d import Interpolant, Spline
from math import sqrt

# %%


class Point3d:
    def __init__(self, x=0., y=0., z=0.):
        self.X = x
        self.Y = y
        self.Z = z

    def toList(self):
        return [self.X, self.Y, self.Z]

    def copy(self):
        return Point3d(self.X, self.Y, self.Z)

# %%


class Vector3d:
    def __init__(self, pt1: Point3d, pt2: Point3d):
        self.vX = pt2.X - pt1.X
        self.vY = pt2.Y - pt1.Y
        self.vZ = pt2.Z - pt1.Z
        self.norma = sqrt(self.vX ** 2 + self.vY ** 2 + self.vZ ** 2)

    def toList(self):
        return [self.vX, self.vY, self.vZ]

    def toNumpy(self):
        return np.array([self.vX, self.vY, self.vZ])

# %%


class Line2d:
    def __init__(self, pt1: Point3d, pt2: Point3d):
        self.start = pt1
        self.end = pt2
        self.V = Vector3d(pt1, pt2)
        self.A = pt1.Y - pt2.Y
        self.Az = pt1.Z - pt2.Z
        self.B = pt2.X - pt1.X
        self.C = pt1.X * pt2.Y - pt1.Y * pt2.X
        self.Cz = pt1.X * pt2.Z - pt1.Z * pt2.X
        self.L = sqrt(self.A * self.A + self.B * self.B)

    def interpolate(self, x):
        if self.B == 0.:
            return np.inf
        return (-self.A * x - self.C) / self.B

    def __interpolateY(self, x: float):
        if self.B == 0.:
            return np.inf
        return (-self.A * x - self.C) / self.B

    def __interpolateZ(self, x: float):
        if self.B == 0.:
            return np.inf
        return (-self.Az * x - self.Cz) / self.B

    def interpolateEps(self, x):
        return Point3d(self.__interpolateY(x), self.__interpolateZ(x), x)
# %%


class Lspline:

    @classmethod
    def build(cls, x: list, y: list):
        res = Interpolant(x, y)
        res.typ = 'L'

        n = len(x) - 1
        if n < 1:
            print('ERROR! Длина входных массивов должна быть не менее 2.')
            return

        for i in range(0, len(x) - 1):
            pt1: Point3d = Point3d(x[i], y[i])
            pt2: Point3d = Point3d(x[i + 1], y[i + 1])
            ln = Line2d(pt1, pt2)
            res.a.append(ln.A)
            res.b.append(ln.B)
            res.c.append(ln.C)
        res.isBuild = True
        return res
