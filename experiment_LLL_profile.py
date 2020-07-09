from numpy import zeros, matrix, array, random
from middleware import CodeRedLib
from time import time
from math import floor

def experiment(n, k, samples):
    profiles = zeros((6, k))
    G = random.randint(0,2, size=(k, n), dtype="bool")
    red = CodeRedLib(G)

    for s in range(samples):

        red.Randomize()
        red.LLL()
        profiles[0] += array(sorted(list(red.l), reverse=True)) / (1.*samples)


        red.Randomize()
        red.Systematize()
        profiles[1] += array(sorted(list(red.l), reverse=True)) / (1.*samples)
        red.LLL()
        profiles[2] += array(sorted(list(red.l), reverse=True)) / (1.*samples)

        red.Randomize()
        red.Systematize()
        red.EpiSort()
        profiles[3] += array(sorted(list(red.l), reverse=True)) / (1.*samples)
        red.LLL()
        profiles[4] += array(sorted(list(red.l), reverse=True)) / (1.*samples)

        red.SizeRedBasis()
        red.KillTwos()
        profiles[5] += array(sorted(list(red.l), reverse=True)) / (1.*samples)

    return profiles

n = 1280
k = int(n/2)
samples = 100
M = experiment(n, k, samples)
C = M.transpose()

print("index,       pLLL_raw,      pSys,    pLLL_Sys,      pSort,     pLLL_Sort,           pLLL_Sort_K2")
for i in range(25):
    print(i+1, ", \t   %.2f,\t   %.2f,\t   %.2f,\t   %.2f,\t   %.2f, \t   %.2f \t"%tuple(C[i]))

