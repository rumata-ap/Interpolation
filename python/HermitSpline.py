from interpolant1d import Interpolant, Spline

# %%


class Hspline:

    @staticmethod
    def __diffp1(x: list, y: list):
        qx = []
        qy = []

        for i in range(1, len(x) - 1):
            qx.append(0.5 * (x[i] - x[i - 1]) + 0.5 * (x[i + 1] - x[i]))
            qy.append(0.5 * (y[i] - y[i - 1]) + 0.5 * (y[i + 1] - y[i]))

        qx.insert(0, 0)
        qx.insert(len(qx), 0)
        qx[0] = 2 * (x[1] - x[0]) - qx[1]
        m = len(x)-1
        qx[m] = 2 * (x[m] - x[m - 1]) - qx[m - 1]

        qy.insert(0, 0)
        qy.insert(len(qy), 0)
        qy[0] = 2 * (y[1] - y[0]) - qy[1]
        m = len(y)-1
        qy[m] = 2 * (y[m] - y[m - 1]) - qy[m - 1]

        return qx, qy

    @staticmethod
    def __diffp2(x: list, y: list):
        qx = []
        qy = []

        for i in range(1, len(x) - 1):
            sxl = abs((x[i] - x[i - 1]))
            syl = abs((y[i] - y[i - 1]))
            sxr = abs((x[i + 1] - x[i]))
            syr = abs((y[i + 1] - y[i]))
            qx.append(sxr * ((x[i] - x[i - 1]) / (sxr + sxl)) +
                      sxl * ((x[i + 1] - x[i]) / (sxr + sxl)))
            if syl + syr > 0:
                qy.append(syr * ((y[i] - y[i - 1]) / (syr + syl)) +
                          syl * ((y[i + 1] - y[i]) / (syr + syl)))
            else:
                qy.append(0)

        qx.insert(0, 0)
        qx.insert(len(qx), 0)
        qx[0] = 2 * (x[1] - x[0]) - qx[1]
        m = len(x)-1
        qx[m] = 2 * (x[m] - x[m - 1]) - qx[m - 1]

        qy.insert(0, 0)
        qy.insert(len(qy), 0)
        qy[0] = 2 * (y[1] - y[0]) - qy[1]
        m = len(y)-1
        qy[m] = 2 * (y[m] - y[m - 1]) - qy[m - 1]

        return qx, qy

    @staticmethod
    def __diff1d1(x: list, y: list):
        q = []
        for i in range(1, len(x) - 1):
            q.append(0.5*((y[i] - y[i - 1]) / (x[i] - x[i - 1]) +
                          (y[i + 1] - y[i]) / (x[i + 1] - x[i])))

        m = len(x) - 1
        q.insert(0, 0)
        q.insert(len(q), 0)
        q[0] = 2 * ((y[1] - y[0]) / (x[1] - x[0])) - q[1]
        q[m] = 2 * ((y[m] - y[m - 1]) / (x[m] - x[m - 1])) - q[m - 1]

        return q

    @staticmethod
    def __diff1d2(x: list, y: list):
        q = []
        for i in range(1, len(x) - 1):
            sl = abs((y[i] - y[i - 1])/(x[i] - x[i - 1]))
            sr = abs((y[i + 1] - y[i]) / (x[i + 1] - x[i]))
            if sl + sr > 0:
                q.append(sr * (((y[i] - y[i - 1])/(x[i] - x[i - 1])) / (sr + sl)) +
                        sl * (((y[i + 1] - y[i]) / (x[i + 1] - x[i])) / (sr + sl)))
            else:
                q.append(0)
        m = len(x) - 1

        q.insert(0, 0)
        q.insert(len(q), 0)
        q[0] = 2 * ((y[1] - y[0]) / (x[1] - x[0])) - q[1]
        q[m] = 2 * ((y[m] - y[m - 1]) / (x[m] - x[m - 1])) - q[m - 1]

        return q

    @classmethod
    def build(cls, x: list, y: list, diff=1, isParametric=False):
        """
        Построение сплайна Эрмита
        """
        res: Interpolant = Interpolant(x, y)
        res.typ = 'H'

        n = len(x) - 1
        if n < 1:
            print('ERROR! Длина входных массивов должна быть не менее 2.')
            return
        if isParametric:
            qx: list
            qy: list
            if diff == 1:
                qx, qy = cls.__diffp1(x, y)
            else:
                qx, qy = cls.__diffp2(x, y)

            for i in range(0, len(x) - 1):
                res.a.append((x[i], y[i]))
                res.b.append((qx[i], qy[i]))
                res.c.append((-3 * x[i] + 3 * x[i + 1] - 2 * qx[i] - qx[i + 1],
                              -3 * y[i] + 3 * y[i + 1] - 2 * qy[i] - qy[i + 1]))
                res.d.append((2 * x[i] - 2 * x[i + 1] + qx[i] + qx[i + 1],
                              2 * y[i] - 2 * y[i + 1] + qy[i] + qy[i + 1]))
        else:
            dy: list
            if diff == 1:
                dy = cls.__diff1d1(x, y)
            else:
                dy = cls.__diff1d2(x, y)

            for i in range(0, len(x)):
                res.a.append(y[i])
                res.b.append(dy[i])

        res.isBuild = True
        return res
