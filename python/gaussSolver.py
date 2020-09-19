import math

def get_roots(matrix): 
    for i in range(0, len(matrix) - 1):
        max_row = get_index_of_max_in_column_by_abs(matrix, i)
        swap_row(matrix, max_row, i)
        make_all_below_rows_to_zero(matrix, i)
    x = solve_system(matrix)
    return x


def get_index_of_max_in_column_by_abs(matrix, column):
    max_element = matrix[column][column]
    number_of_max_row = 0
    for i in range(column + 1, len(matrix)):
        if math.fabs(matrix[i][column]) > math.fabs(max_element):
            max_element = matrix[i][column]
            number_of_max_row = i
    return number_of_max_row


def swap_row(matrix, n, m):
    temp = matrix[n]
    matrix[n] = matrix[m]
    matrix[m] = temp


def make_all_below_rows_to_zero(matrix, i):
    for k in range(i + 1, len(matrix)):
        c = -matrix[k][i] / matrix[i][i]
        for j in range(i, len(matrix) + 1):
            if i == j:
                matrix[k][j] = 0
            else:
                matrix[k][j] += c * matrix[i][j]


def solve_system(matrix):
    x = [0] * len(matrix)
    n = len(matrix)
    for i in range(n - 1, -1, -1):
        x[i] = (matrix[i][n] / matrix[i][i])
        for k in range(i - 1, -1, -1):
            matrix[k][n] -= matrix[k][i] * x[i]
    return x


def print_matrix(matrix):
    for string in matrix:
        for column in string:
            print("\t" + str(column), end=" ")
        print()
    print("-------------------------------------------")
