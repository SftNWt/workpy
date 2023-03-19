import math
from numpy import *
import scipy
from functions import *

# Выставляем количество знаков после запятой в матрице
set_printoptions(precision=3)

# Получаем вариант
N = int(input("Ваш вариант: "))

MATRIX_A = array([
    [12 * N, N, N + 2],
    [N, 9 * N, 3 * N],
    [2 * N, 3 * N, 11 * N],
])
MATRIX_B = array([
    [2 * N], [N], [3 * N]
])

ROOTS         = gauss (MATRIX_A, MATRIX_B)
[roots, L, U] = crout (MATRIX_A, MATRIX_B)
JACOBI        = jacobi(MATRIX_A, MATRIX_B)
GAUSS_SAIDEL = gauss_saidel(MATRIX_A, MATRIX_B)

print(ROOTS)