import math
from numpy import *
import scipy

def forward_elimination(A, b, n):
    for row in range(0, n-1):
        for i in range(row+1, n):
            factor = A[i,row] / A[row,row]
            for j in range(row, n):
                A[i,j] = A[i,j] - factor * A[row,j]

            b[i] = b[i] - factor * b[row]
        str1 = 'A = '
        str2 = ' и B = '
        string = ''
        for i in range(0, len(A)):
            if i == 0 or i == 2:
                string += '    '
                string += str(A[i])
                string += '       '
                string += str(b[i])
                string += '\n'
            elif i == 1:
                string += str1
                string += str(A[i])
                string += str2
                string += str(b[i])
                string += '\n'
        # print(string)
    return A, b

def back_substitution(a, b, n):
    x = zeros((n,1))
    x[n-1] = b[n-1] / a[n-1, n-1]
    for row in range(n-2, -1, -1):
        sums = b[row]
        for j in range(row+1, n):
            sums = sums - a[row,j] * x[j]
        x[row] = sums / a[row,row]
    return x

def gauss(A, b):
    n = A.shape[0]
    if any(diag(A)==0):
        raise ZeroDivisionError(('Division by zero will occur; '
                                  'pivoting currently not supported'))

    A, b = forward_elimination(A, b, n)
    return back_substitution(A, b, n)

def crout(a, b):
    cout = 0
    m, n = a.shape
    if (m !=n ):
        print("Crout cannot be used.")
    else:
        l = zeros((n,n))
        u = zeros((n,n))
        s1 = 0
        s2 = 0

        for i in range(n):
           l[i][0] = a[i][0]
           u[i][i] = 1
        for j in range(1, n):
            u[0][j] = a[0][j] / l[0][0]
        for k in range(1, n):
            for i in range(k, n):
                for r in range(k): s1 += l[i][r] * u[r][k]
                l[i][k] = a[i][k] - s1
                s1 = 0
            for j in range(k+1, n):
                for r in range(k): s2 += l[k][r] * u[r][j]
                u[k][j] = (a[k][j] - s2) / l[k][k]
                s2 = 0

        y = zeros(n)
        s3 = 0
        y[0] = b[0] / l[0][0]
        for k in range(1, n):
            for r in range(k):
                s3 += l[k][r] * y[r]
            y[k] = (b[k]-s3) / l[k][k]
            s3 = 0

        x = zeros(n)
        s4 = 0
        x[n-1] = y[n-1]
        for k in range(n-2, -1, -1):
            for r in range(k+1, n):
                s4 += u[k][r] * x[r]
            x[k] = y[k] - s4
            s4 = 0

        roots = {'x': list(x)}
        for i in range(n):
            roots["x" + str(i + 1)] = x[i]
        return [roots, l, u]

def jacobi(A,b,N=2,x=None):
    if x is None:
        x = zeros(len(A[0]))

    D = diag(A)
    R = A - diagflat(D)

    for i in range(N):
        x = (b - dot(R,x)) / D
    return x

def gauss_saidel(A,b, N=2, esp=1e-6):
        n = len(A)

        btol = linalg.norm(b) * esp

        x0   = zeros(n)
        k    = 0 ;
        isActive = False
        x1   = empty(n)

        while not(isActive) and k < N:
            print (f'Начало цикла при k = {k}')
            x1 = zeros(n)
            for i in range(n):          # rows of A
               x1[i] = ( b[i] - dot(A[i,0:i], x1[0:i]) - dot(A[i,i+1:n], x0[i+1:n]) ) / A[i,i]
               print("x1  = ", x1)

            r = b - dot(A,x1)
            isActive = (linalg.norm(r) < btol) and (linalg.norm(x1-x0) < esp)
            print(f'\n')
            print("x0  = ", x0)
            print("btol = %e;\nla.norm(r) = %e;\ntol = %e;\nla.norm(x1-x0) = %e;\nisActive = %s " % (btol, linalg.norm(r), esp, linalg.norm(x1-x0), isActive))
            x0   = x1
            print("x0  = ", x0, end='\n')
            print("Окончание цикла \n\n")
            k    = k + 1

        if not(isActive): # or if k >= N
            print(f'Не сходится за {N} итерации.')

        return x1
