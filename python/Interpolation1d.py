import numpy as np
import matplotlib.pyplot as plt
from interpolant1d import Interpolant
from LinearSpline import Lspline
from AkimaSpline import Aspline
from HermitSpline import Hspline
from CubicSpline import Cspline

# %%


class Spline1d:
    @staticmethod
    def create(typeSpline: str):
        if typeSpline == 'L':
            return Lspline()
        if typeSpline == 'A':
            return Aspline()
        if typeSpline == 'H':
            return Hspline()
        if typeSpline == 'C':
            return Cspline()


# %%

if __name__ == "__main__":
    xm = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    ym = [10, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85]
    # xm = [0, 1, 2, 3, 4, 5]
    # ym = [0.1, 0.25, 0.12, 0.62, 0.38, 0.31]
    #ym = [0.1, 0.25, 0.12, 0.62, 0.38, 0.25]
    s = Spline1d.create('L')
    lspline = s.build(xm, ym)
    #aspline1 = Aspline.build(xm, ym, 1)
    #s = Spline1d.create('A')
    s = Aspline()
    aspline = s.build(xm, ym , 2)
    s = Spline1d.create('H')
    hspline = s.build(xm, ym, diff=2)
    hsplinep = s.build(xm, ym, diff=1, isParametric=True)
    s = Spline1d.create('C')
    сspline = s.build(xm, ym)
    vals = np.linspace(xm[0], xm[len(xm) - 1], num=101, endpoint=True)
    valt = np.linspace(0, len(xm) - 1, num=101, endpoint=True)
    resL = lspline.interpolate(vals)
    resH = hspline.interpolate(vals)
    resHp = hsplinep.interpolate(valt)
    #resA1 = aspline1.interpolate(vals)
    resA = aspline.interpolate(vals)
    resС = сspline.interpolate(vals)
    plt.figure('Одномерная интерполяция')
    plt.plot(xm, ym, 'bo', label='data')
    plt.plot(vals, resL, 'r--', label='linear')
    #plt.plot(vals, resA1, label='akima1')
    plt.plot(vals, resA, label='akima')
    plt.plot(vals, resH, label='hermit')
    # plt.plot(resHp[0], resHp[1], label='hermitp')
    plt.plot(vals, resС, label='cubic')
    # plt.legend(['data', 'linear', 'akima', 'hermit',
    #             'hermit_p', 'cubic'], loc='best')
    plt.show()
