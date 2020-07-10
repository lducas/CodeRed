from numpy import zeros, array, random
from middleware import CodeRedLib
from time import time


def experiment(n, k, samples):
    times = zeros(5)
    G = random.randint(0,2, size=(k, n), dtype="bool")
    red = CodeRedLib(G)

    for s in range(samples):


        T0 = time()
        red.Randomize()
        red.LLL()
        T1 = time()


        red.Randomize()
        red.Systematize()
        T2 = time()
        red.LLL()
        T3 = time()

        red.Randomize()
        red.Systematize()
        red.EpiSort()
        T4 = time()
        red.LLL()

        T5 = time()
        times += array([T1 - T0, T3 - T1, T3 - T2, T5 - T3, T5 - T4])

    return times/samples


print("n,       tLLL_raw,      tLLL_Sys,     tLLL_afterSys,     tLLL_Sort,    tLLL_afterSort")
for n in [128, 192, 256, 384, 512, 768, 1024, 1280, 1536, 2048, 3072, 4096, 6144, 8192, 12288, 16384]:
    k = int(n/2)
    v = experiment(n, k, 10)
    print (n, ",\t %.4e,\t %.4e,\t %.4e,\t %.4e,\t %.4e"%tuple(v))