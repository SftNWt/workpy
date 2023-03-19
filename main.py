from numpy import *
import scipy 

set_printoptions(precision=2)

# N = int(input("Ваш вариант: "))
N = 21

# Находим матрицу
MATRIX = matrix([
    [N, 4*N, N+2, N, N, N],
    [0, N+1, N, N, N+1, 2],
    [2*N, N, N+2, 2*N, 3*N, N+1],
    [N, 7*N, N+1, N, 4*N, N],
    [N, N+1, N, N, N+1, 2],
    [2*N, N, N+2, 2*N, N, N+1]
])

# Находим союзную матрицу и множитель
M_UNION = f"{N} * {MATRIX / N}"

# Находим транспонированную матрицу
M_TRANS = transpose(MATRIX)

# Выполняем сложение обычной и транспонированной матриц
M_SUMOT = MATRIX + M_TRANS

# Выполняем умножение обычной и транспонированной матриц
M_MULOT = MATRIX * M_TRANS

# Находим определитель
DET = int(linalg.det(MATRIX))

# Находим ранг матрицы
M_RANK = linalg.matrix_rank(MATRIX)

# Находим обратную матрицу
M_INV = linalg.inv(MATRIX)

# Выполняем умножение обычной и обратной матриц
M_MULOI = MATRIX + M_INV

# Находим собственные значения и собственные вектора матрицы
M_W, M_V = linalg.eig(MATRIX)

# Выполняем LU-разложение матрицы
P, M_L, M_U = scipy.linalg.lu(MATRIX)

# Выполняем перемножение LU-матриц
M_LU = M_L @ M_U

print(M_LU)