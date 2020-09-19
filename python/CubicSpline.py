from interpolant1d import Interpolant
from tridiagonalSolver import SolverTriD

# %%


class Cspline:

    @classmethod
    def build(cls, x: list, y: list):
        """Построение кубического сплайна с простыми граничными условиями.

        Args:
            x (list): список координат аргументов точечно заданной функции(len(x) >= 3)
            y (list): список координат значений точечно заданной функции(len(y) == len(x) >= 3)

        Returns:
            Interpolant: интеполирующая функция
        """
        res: Interpolant = Interpolant(x, y)
        res.typ = 'С'

        n = len(x) - 1

        if n < 2:
            print('ERROR! Длина входных массивов должна быть не менее 3.')
            return

        h = []
        matrix = []
        for i in range(n):
            res.a.append(y[i])
            h.append(x[i + 1] - x[i])
            matrix.append([])
            for j in range(n):
                matrix[i].append(0)
                
        matrix = []
        for i in range(n-1):
            matrix.append([])
            for j in range(n-1):
                matrix[i].append(0)

        matrix[0][0] = 2 * (h[0] + h[1])
        matrix[0][1] = h[1]

        for i in range(1, n - 2):
            matrix[i][i - 1] = h[i]
            matrix[i][i] = 2 * (h[i] + h[i+1])
            matrix[i][i + 1] = h[i + 1]

        matrix[n - 2][n - 3] = h[n - 2]
        matrix[n - 2][n - 2] = 2 * (h[n - 2] + h[n - 1])
        
        f = [0]*(n-1)
        for i in range(1, n):
            f[i - 1] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
            
        tmp = SolverTriD.solve(matrix, f)
        
        c = [0]*(n+1)
        for i in range(1, n):
            c[i] = tmp[i - 1]

        d = [0]*(n)
        for i in range(n):
            d[i] = (c[i + 1] - c[i]) / (3 * h[i])

        b = [0]*(n)
        for i in range(n):
            b[i] = ((y[i + 1] - y[i]) / h[i]) - \
                (h[i] / 3) * (c[i + 1] + 2 * c[i])

        for i in range(n):
            res.b.append(b[i])
            res.c.append(c[i])
            res.d.append(d[i])

        res.isBuild = True
        return res
