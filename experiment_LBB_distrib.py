from numpy import zeros, matrix, array, random
from middleware import CodeRedLib
from time import time
from math import floor
from weights import *

def experiment(n, k, samples):
    curves = zeros((4, n+1))
    G = random.randint(0,2, size=(k, n), dtype="bool")
    red = CodeRedLib(G)

    for s in range(samples):

        T0 = time()
        red.Randomize()
        red.Systematize()
        T1 = time()
        print("Produced Systematic form in %.4fs"%(T1 - T0))

        T0 = time()
        stats = red.LB(3, stats = True)
        curves[0] += stats
        T1 = time()
        print("Ran LB(3) in %.4fs"%(T1 -T0))
        print("Visited %.4f codewords"%sum(stats))

        stats = weights_LB(n, k, 3) 
        print("Predicted %.4f codewords"%sum(stats))
        curves[1] += array(stats+(n+1)*[0])[:n+1]

        print()
        T0 = time()
        red.Randomize()
        red.Systematize()
        red.EpiSort()
        red.LLL()
        red.KillTwos()
        k1 = red.SemiSystematize()

        T1 = time()
        print("Produced Semisystematic form in %.4fs"%(T1 - T0))
        print("k1 = %d"%k1)
        print("profile = ", red.l[:k1+1])
        assert(sum(red.l[:k1]) == k+k1)
        assert(min(red.l[k1:]) == 1)
        assert(max(red.l[k1:]) == 1)

        T0 = time()
        stats = red.LBB(k1, 3, stats=True)
        curves[2] += stats
        T1 = time()
        print("Ran LBB(k1, 3) in %.4fs"%(T1 -T0))

        print("Visited %.4f codewords"%sum(stats))
        stats = weights_LBB(red.l, k1, 3)
        stats = array(stats +(n+1)*[0])
        print("Predicted %.4f codewords"%sum(stats))
        curves[3] += stats[:n+1]

    return curves/samples


n = 1280
k = 640

curves = experiment(n, k, 1)

print("w,       LB_exp,     LB_pred,    LBB_exp,    LBB_pred")
C = curves.transpose()
for i in range(230, 340):
    print(i, "\t %.3e \t%.3e \t%.3e \t%.3e \t"%tuple(C[i]))
