import sys
import scipy
import scipy.linalg as la
import numpy as np

n = int(sys.argv[1])

A = np.random.random((n,n))
b = np.random.random(n)
x = la.solve(A, b)
