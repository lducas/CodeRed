from numpy import zeros, matrix, array, random
from middleware import CodeRedLib
from time import time
from math import floor, ceil, sqrt, log

from multiprocessing import Pool
from random import randint
from weights import *

def one_experiment(par):
    n, k, seed = par
    
    h = int(ceil(0.115 * n))

    data = zeros(9)
    G = random.randint(0,2, size=(k, n), dtype="bool")
    red = CodeRedLib(G, seed=seed)

    # Defining the skip parameter to only expore 
    # a portion of the space and get relevant data in reasonable time
    # This speeds up things by (1+skip)^{w2-1} where w2=3 is the LB/LBB parameter

    skip = floor(sqrt(k/256.))
    red.set_skip(skip)

    T0 = time()
    red.Randomize()
    red.Systematize()
    red.LB(3)
    T1 = time()
    predicted_distr, denom = weights_LB_absolute(n, k, 3)
    data[0] += log(sum(predicted_distr[:h+1]), 2) - log(denom,2)
    data[1] += T1 - T0

    T0 = time()
    red.Randomize()
    red.Systematize()
    red.EpiSort()
    red.LLL()
    red.KillTwos()
    k1 = red.SemiSystematize()
    red.LBB(k1, 3)

    T1 = time()
    predicted_distr, denom = weights_LBB_absolute(red.l, k1, 3)
    data[2] += log(sum(predicted_distr[:h+1]), 2) - log(denom,2)
    data[3] += T1 - T0
    data[7] += k1
    return data


def experiment(n, k, samples, cores=1):
    if cores == 1:
        res = [one_experiment((n, k, randint(0,2**63))) for i in range(samples)]
    else:
        p = Pool(cores)
        res = p.map(one_experiment, [(n, k, randint(0,2**63)) for i in range(samples)])
            
    return sum(res)/samples

print("n, \t LB_logP, \tLB_time,  \tLBB_logP, \tLBB_time,  \t ProbaGain,  \tTimeLoss,  \tGain,  \t\tk1,   \tCGain")
for n in [128, 192, 256, 384, 512, 768, 1024, 1280, 1536, 2048, 3072, 4096, 6144, 8192, 12288, 16384]:
    k = n//2
    C = experiment(n, k, 12, cores=4)
    C[4] =  2**(C[2]-C[0])  # Proba gain LBB/LB
    C[5] =  C[3]/C[1]       # time loss LBB/LB
    C[6] =  C[4]/C[5]       # overall gain LBB/LB
    C[8] =  C[4]/C[7]       # overall Corrected gain LBB/LB
    print(n, ",\t %.1f,  \t%.3f,  \t%.1f,   \t%.3f,  \t%.3f,  \t%.3f,  \t%.3f,  \t%.1f,  \t%.3f"%tuple(C))
