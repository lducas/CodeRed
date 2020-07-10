from numpy import array, random
from middleware import CodeRedLib

B = random.randint(0,2, size=(5, 16), dtype="bool")  # Create a random Basis for a [5,16]-code
red = CodeRedLib(B)                                  # Load it into a fresh CodeRedLib object  

def niceprint(B):
    for v in B:
        print("".join(["1" if x else "." for x in v]))
    print()

  
niceprint(red.B)    # Print current basis
niceprint(red.E)    # Print current Epipodal matrix
print(red.l)        # Print current Profile

red.LLL()           # Apply LLL

niceprint(red.B)    # Print current basis
niceprint(red.E)    # Print current Epipodal matrix
print(red.l)        # Print current Profile
