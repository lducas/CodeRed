from numpy import zeros, matrix, array, random
from middleware import CodeRedLib
from time import time
from math import floor, ceil, sqrt, log

from multiprocessing import Pool
from random import randint
from weights import *

def dGV(n, k):
    d = 0
    aux = 2**(n-k)
    b = 1
    while aux >= 0:
        aux -= b
        d += 1
        b *= (n-d+1)
        b /= d
    return d 

def one_experiment(par):
    n, k, seed = par
    
    dgv = dGV(n, k)
    goal = ceil(1.05 * dgv)
    # print()
    #print("length = ", n,  "dgv = ", dgv, "goal hw = ", goal)

    data = zeros(2)
    G = random.randint(0,2, size=(k, n), dtype="bool")
    
    # make sure its invertible
    for i in range(k):
        G[i, i] = 1
        for j in range(i):
            G[i, j] = 0
            G[j, i] = 0

    red = CodeRedLib(G, seed=seed)

    # Defining the skip parameter to only expore 
    # a portion of the space and get relevant data in reasonable time
    # This speeds up things by (1+skip)^{w2-1} where w2=3 is the LB/LBB parameter

    T0 = time()
    red.Sieve(goal)
    # print(time() - T0)
    data[0] = time() - T0
    
    T0 = time()
    red.Prange(goal)
    # print(time() - T0)
    data[1] = time() - T0
    return data
    
    # L = []

    # for loop in range(2**n):
    #     red.Randomize()
    #     red.Systematize()
    #     L += [(loop, sum(x), x) for x in red.B]
    #     # print(L)
    #     N = len(L)
    #     for i in range(N):
    #         birth, hx, x = L[i]
    #         if birth < loop - 1:
    #             continue

    #         for j in range(i):
    #             _, hy, y = L[j]
    #             z = x ^ y
    #             hz = sum(z)
    #             if hz == 0:
    #                 continue
    #             if hz <= h:
    #                 data[0] = time() - T0
    #                 data[1] = N
    #                 data[2] = loop
    #                 return data
    #             if hz < hy:
    #                 L[j] = loop, hz, z
    #             elif hz < hx:
    #                 L[i] = loop, hz, z
    #                 break





def experiment(n, k, samples, cores=1):
    if cores == 1:
        res = [one_experiment((n, k, randint(0,2**63))) for i in range(samples)]
    else:
        p = Pool(cores)
        res = p.map(one_experiment, [(n, k, randint(0,2**63)) for i in range(samples)])
            
    return sum(res)/samples

for n in range(20, 200, 10):
    print(n, experiment(n, n//2, 12, cores=4))