import warnings
import operator as op
import numpy as np 
from functools import reduce

import sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


# n choose r
def comb(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

# convolution of two distributions/measures with support {0, 1, ..., n}
def convol(A, B):
    C = [0 for i in range(len(A)+len(B))]  
    for x in range(len(A)):
        for y in range(len(B)):
            C[x+y] += A[x]*B[y]
    return C

# returns the weight enumerator of a ball of radius r in F_2^n
def weights_ball(n, r):
    return [comb(n, int(i)) for i in range(r+1)]

def volume_ball(n,r):
    return sum( weights_ball(n,r) )

# returns the weight enumerator of a fundamental ball of length n
def weights_fundamental_ball(n):
    if n%2 > 0:
        L = [comb(n, i) for i in range((n+1)//2)]
    else:
        L = [comb(n, i) for i in range((n+1)//2)] + [comb(n, (n+1)//2)//2]
    return L

# returns the weight enumerator of the SizeRed fundamental domain
def weights_fundamental_domain(profile):
    fundamental_balls = [weights_fundamental_ball(int(l)) for l in profile]
    n = sum(profile)
    k = len(profile)
    C = fundamental_balls[0]
    for i in range(1, k):
        C = convol(C, fundamental_balls[i])
    return C

# return the expected weights of visited codewords by Lee-Brickell
def weights_LB(n, k, w2):
    a,b = weights_LB_absolute(n, k, w2)
    return [x/float(b) for x in a]

# returns integer list a and integer b such that a/b = weights_LB(n,k,w2)
# prevents overflow
def weights_LB_absolute(n, k, w2):
    A = weights_ball(n-k, n-k)
    B = weights_ball(k, w2)
    return convol(A, B), 2**int(n-k)

# returns the expected weights of visited codewords by Lee-Brickell-Babai
def weights_LBB(profile, k1, w2):
    a,b = weights_LBB_absolute(profile, k1, w2)
    return [x/float(b) for x in a]

# returns integer list a and integer b such that a/b = weights_LBB(profile,k1,w2)
# prevents overflow
def weights_LBB_absolute(profile, k1, w2):
    n = sum(profile)
    k = len(profile)
    A = weights_fundamental_domain(profile[:k1])
    B = weights_ball(k - k1, w2)
    return convol(A, B), 2**int(sum(profile[:k1])-k1)

