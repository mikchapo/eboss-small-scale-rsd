import numpy as np
import sys

par = np.loadtxt(sys.argv[1])
fib = np.loadtxt(sys.argv[2])
norms = np.loadtxt(sys.argv[3], skiprows=1)

par[:, -1] = par[:, -1] / norms[0]
fib[:, -1] = fib[:, -1] / norms[1]

weight = np.copy(par)
weight[:, -1] = weight[:, -1] / fib[:, -1]

np.savetxt(sys.argv[4], weight)
