import numpy as np
from interpolant1d import Interpolant, Spline

# %%


class Aspline:

    @staticmethod
    def __akimaBC2(x1, x2, x3, y1, y2, y3):
        """Вычисление граничных условий 1-го типа для сплайна Акимы.

        Args:
            x1 (float): аргумент точечно заданной функции в 1-й точке
            x2 (float): аргумент точечно заданной функции во 2-й точке
            x3 (float): аргумент точечно заданной функции в 3-й точке
            y1 (float): значение точечно заданной функции в 1-й точке
            y2 (float): значение точечно заданной функции во 2-й точке
            y3 (float): значение точечно заданной функции в 3-й точке

        Returns:
            tuple: (производная в 4-й точке, производная в 5-й точке)
        """          
        sx = x3 - x1
        x4 = sx + x2
        x5 = sx + x3
        sy = (y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)
        y4 = (sy + (y3 - y2) / (x3 - x2)) * (x4 - x3) + y3
        y5 = (sy + (y4 - y3) / (x4 - x3)) * (x5 - x4) + y4
        return (y4 - y3) / (x4 - x3), (y5 - y4) / (x5 - x4)

    @staticmethod
    def __akimaBC1(x2, x3, y2, y3):
        """
        Вычисление граничных условий сплайна Акимы 1-го типа
        """
        x4 = x3 + (x3 - x2)
        x5 = x4 + (x3 - x2)
        y4 = y3 + (y3 - y2)
        y5 = y4 + (y3 - y2)
        return (y4 - y3) / (x4 - x3), (y5 - y4) / (x5 - x4)

    @classmethod
    def build(cls, x: list, y: list, bc=2):
        """Построение сплайна Акимы с заданными граничными условиями.
      
        Args:
            x (list): список координат аргументов точечно заданной функции(len(x) >= 5)
            y (list): список координат значений точечно заданной функции(len(y) == len(x) >= 5)
            bc (int, optional): тип граничных условий. Defaults to 2.

        Returns:
            Interpolant: интеполирующая функция
        """        

        res: Interpolant = Interpolant(x, y)
        res.typ = 'A'

        if len(x) - 1 < 4:
            print('ERROR! Длина входных массивов должна быть не менее 5.')
            return
        n = len(x) - 1
        m = []

        for i in range(0, n):
            m.append((y[i + 1] - y[i]) / (x[i + 1] - x[i]))
        if bc == 1:
            m1n, m2n = cls.__akimaBC2(
                x[n - 2], x[n - 1], x[n], y[n - 2], y[n - 1], y[n])
            m1m, m2m = cls.__akimaBC2(x[2], x[1], x[0], y[2], y[1], y[0])
            m.insert(0, m1m)
            m.insert(0, m2m)
            m.insert(n+2, m1n)
            m.insert(n + 3, m2n)
        elif bc == 2:
            m1m = 2*m[0] - 2*m[1]
            m.insert(0, m1m)
            m2m = 3*m[0] - 2*m[1]
            m.insert(0, m2m)

            m1n = 2*m[n-1+2] - 2*m[n-2+2]
            m.append(m1n)
            m2n = 3*m[n-1+3] - 2*m[n-2+3]
            m.append(m2n)
        elif bc == 3:
            m1m = 2*m[n-1]
            m.insert(0, m1m)
            m2m = 3*m[n-2+1]
            m.insert(0, m2m)

            m1n = 2*m[1]
            m.append(m1n)
            m2n = 3*m[0]
            m.append(m2n)
        else:
            m1n, m2n = cls.__akimaBC1(x[n - 1], x[n], y[n - 1], y[n])
            m1m, m2m = cls.__akimaBC1(x[1], x[0], y[1], y[0])

            m.insert(0, m1m)
            m.insert(0, m2m)
            m.insert(n+2, m1n)
            m.insert(n + 3, m2n)

        tL = []
        tR = []
        for i in range(2, n + 3):
            NE = abs(m[i + 1] - m[i]) + abs(m[i - 1] - m[i - 2])
            if NE > 0:
                t = (abs(m[i + 1] - m[i]) * m[i - 1] +
                     abs(m[i - 1] - m[i - 2]) * m[i]) / NE
                tL.append(t)
            else:
                t = (m[i - 1] + m[i]) / 2
                tL.append(t)

        for i in range(0, n):
            tR.append(tL[i+1])

        for i in range(0, n):
            res.a.append(y[i])
            res.b.append(tL[i])
            h = x[i + 1] - x[i]
            res.c.append((3 * (y[i + 1] - y[i]) /
                          (x[i + 1] - x[i]) - 2 * tL[i] - tR[i]) / h)
            res.d.append((tL[i] + tR[i] - 2 * (y[i + 1] - y[i]) /
                          (x[i + 1] - x[i])) / (h ** 2))

        res.isBuild = True
        return res
