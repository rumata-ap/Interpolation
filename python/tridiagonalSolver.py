class SolverTriD:
    @staticmethod
    def solve(matrix, b):
        eps = [0]*len(b)
        et = [0]*len(b)
        n = len(b) - 1
        eps[0] = -matrix[0][1] / matrix[0][0]
        et[0] = b[0] / matrix[0][0]

        for i in range(1, n):
            z = matrix[i][i] + matrix[i][i - 1] * eps[i - 1]
            eps[i] = -matrix[i][i + 1] / z
            et[i] = (b[i] - matrix[i][i - 1] * et[i - 1]) / z

        res = list(range(n + 1))
        res[n] = (b[n] - matrix[n][n - 1] * et[n - 1]) / \
            (matrix[n][n] + matrix[n][n - 1] * eps[n - 1])

        rew = list(range(n))
        rew.reverse()
        for i in rew:
            res[i] = eps[i] * res[i + 1] + et[i]

        return res


if __name__ == "__main__":
    matrix = [
        [2.0, 1.0, 0],
        [-3.0, -1.0, 2.0],
        [0, 1.0, 2.0]]

    b = [8, -11, -3]

    print(SolverTriD.solve(matrix, b))
