
# %%


class Spline(object):
    typ: str

    def interpolate(self, val):
        pass
# %%


class Interpolant(Spline):

    def __init__(self, x: list, y: list):
        if len(x) != len(y):
            print('ERROR! Входные массивы различной длины. Проверьте исходные данные.')
            return
        self.x = x
        self.y = y
        self.a = []
        self.b = []
        self.c = []
        self.d = []
        self.isBuild = False

    def interpolate(self, val):
        v = list(val)
        res = []
        #Линейная интерполяция
        if self.isBuild and self.typ == 'L':
            for item in v:
                if item < self.x[0]:
                    res.append((-self.a[0] * item - self.c[0]) / self.b[0])
                if item >= self.x[len(self.x) - 1]:
                    n = len(self.x) - 2
                    res.append((-self.a[n] * item - self.c[n]) / self.b[n])
                for i in range(0, len(self.x)-1):
                    if item >= self.x[i] and item < self.x[i + 1]:
                        res.append((-self.a[i] * item - self.c[i]) / self.b[i])
        #Интерполяция сплайном Акимы
        elif self.isBuild and self.typ == 'A':
            for item in v:
                if item == self.x[len(self.x) - 1]:
                    n = len(self.x) - 1
                    res.append(self.y[n])
                for i in range(0, len(self.x) - 1):
                    if item >= self.x[i] and item < self.x[i + 1]:
                        res.append((self.a[i] + self.b[i] * (item - self.x[i]) + self.c[i] *
                                    (item - self.x[i]) ** 2 + self.d[i] * (item - self.x[i]) ** 3))
        #Параметрическая интероляция сплайном Эрмита
        elif self.isBuild and self.typ == 'H' and type(self.a[0]) == tuple:
            resx = []
            resy = []
            for item in v:
                if item == self.x[len(self.x) - 1]:
                    n = len(self.x) - 1
                    resx.append(self.x[n])
                    resy.append(self.y[n])
                for i in range(0, len(self.x)-1):
                    if item >= i and item < i + 1:
                        w = item - i
                        resx.append(self.a[i][0] + self.b[i][0] * w +
                                    self.c[i][0] * w ** 2 + self.d[i][0] * w ** 3)
                        resy.append(self.a[i][1] + self.b[i][1] * w +
                                    self.c[i][1] * w ** 2 + self.d[i][1] * w ** 3)
            res = (resx, resy)
        #Интероляция сплайном Эрмита
        elif self.isBuild and self.typ == 'H' and (type(self.a[0]) == float or type(self.a[0]) == int):
            for item in v:
                p = self.a
                q = self.b
                if item == self.x[len(self.x) - 1]:
                    n = len(self.x) - 1
                    res.append(self.y[n])
                    # res.append(self.y[n])
                for i in range(0, len(self.x) - 1):
                    if item >= self.x[i] and item < self.x[i + 1]:
                        d = self.x[i + 1] - self.x[i]
                        w = (item - self.x[i]) / d
                        res.append((2 * w**3 - 3 * w**2 + 1) * p[i] + (w**3 - 2 * w**2 + w) * d * q[i] + (
                            -2 * w ** 3 + 3 * w ** 2) * p[i + 1] + (w**3 - w**2) * d * q[i+1])
        #Интерполяция кубическим сплайном
        elif self.isBuild and self.typ == 'С':
            for item in v:
                if item == self.x[len(self.x) - 1]:
                    n = len(self.x) - 1
                    res.append(self.y[n])
                for i in range(0, len(self.x) - 1):
                    if item >= self.x[i] and item < self.x[i + 1]:
                        res.append((self.a[i] + self.b[i] * (item - self.x[i]) + self.c[i] *
                                    (item - self.x[i]) ** 2 + self.d[i] * (item - self.x[i]) ** 3))
        if len(res) > 1:
            return res
        elif len(res) == 1:
            return res[0]
        else:
            return False
