# ![CodeRed](CodeRedSmall.png)  CodeRed
Basis Reduction Algorithms for Codes (LLL and more)

---

This small library is research Artefact, paired with the article  
  
**An Algorithmic Reduction Theory for Binary Codes: LLL and more**  
*Thomas Debris--Alazard, Léo Ducas and Wessel P.J. van Woerden*

The article is available on this [repository](paper.pdf) and soon on the [IACR eprint](https://eprint.iacr.org/).

It provides a core in C++ for reduction algorithms for codes (LLL, Size-Reduction, Lee-Brickell, Lee-Brickell-Babai), and a Python interface.

#### License:   
MIT License see [LICENSE.txt](LICENSE.txt).

#### How to cite:   

```
@misc{DDW20,
    author = {Thomas {Debris-Alazard} and Léo Ducas and Wessel P.J. van Woerden},
    title = {An Algorithmic Reduction Theory for Binary Codes: LLL and more},
    howpublished = {Cryptology ePrint Archive, Report 2020/XXX},
    year = {2020},
    note = {\url{https://eprint.iacr.org/2020/XXX}},
}
```

---
## Installation

#### Requirement  
gcc  
python  
numpy  

#### Compilation

To compile the C++ core before reproducing the experiments simply run
`bash compile_cpp_core.sh`.  

A maximum dimension (n, multiple of 64) is hardcoded at compile time, and the above creates binaries for `n=256,384,512,768...` If you use this library for other purpose than reproducing experiment, please adjust `compile_cpp_core.sh` to your needs.

---
## Reproducing Experiments from the paired article

Figure 4: `python experiment_LLL_time.py`  
Figure 5: `python experiment_LLL_profile.py`  
Figure 6: `python experiment_LBB_distrib.py`  
Figure 7: `python experiment_LBB_perf.py`  

To reach very large dimensions, the last experiment may be run using several cores by editing the python script.

---
## Usage from Python

All reduction algorithms are available through the class `CodeRedLib` and are applied on the internally stored basis. All input/outputs are `numpy` arrays.  
  
#### Example:  

``` python
>>> from numpy import array, random
>>> from middleware import CodeRedLib
>>>  
>>> B = random.randint(0,2, size=(5, 16), dtype="bool")  # Create a random Basis for a [5,16]-code
>>> red = CodeRedLib(B)                                  # Load it into a fresh CodeRedLib object  
>>> # The above fails in the unlucky case of the code C(B) not of full rank or not of full length.  
>>>  
>>> def niceprint(B):
...     for v in B:
...         print("".join(["1" if x else "." for x in v]))
...     print
... 
>>>   
>>> niceprint(red.B)    # Print current basis
1..1.1.1...1.1..
11..1..1111..1.1
11.11.1.1...111.
.111.11.1....1.1
11....1.1.1.1...
>>> niceprint(red.E)    # Print current Epipodal matrix
1..1.1.1...1.1..
.1..1...111....1
......1.....1.1.
..1.............
................
>>> print(red.l)        # Print current Profile
[6 6 3 1 0]
>>> 
>>> red.LLL()           # Apply LLL
>>> 
>>> niceprint(red.B)    # Print current basis
1..1.1.1...1.1..
..1....1..111..1
.111.11.1....1.1
1.111111.11.....
11.11.1.1...111.
>>> niceprint(red.E)    # Print current Epipodal matrix
1..1.1.1...1.1..
..1.......1.1..1
.1....1.1.......
....1....1......
..............1.
>>> print(red.l)        # Print current Profile
[6 4 3 2 1]
```

#### API:  
  
Functions names match with the one of the paper.  

Basis Processing functions:
``` python
Randomize()
Systematize()
LLL()
EpiSort()
SizeRedBasis()
SemiSystematize()
KillTwos()
```

Functions for finding short words in the code C (or in the coset t+C):
``` python
SizeRed(t)
LB(w2, goal_w=None, t=None, stats=False)
LBB(k1, w2, goal_w=None, t=None, stats=False)
```
Parameters:  
- The default value t=None is interpreted as 0 (*i.e.* LB/LBB searches for a short codeword rather than a close codeword).   
- Parameters `k1` and `w2` are integers, whose role is described in the paper.   
- Leave default value `goal_w=None` to get the shortest visited codeword as output. Set `goal_w` to an integer to return as soon as a coderword of at most that length is found (and returns `None` if goal not met).  
- Set `Stats=True` to instead get as return value the counts on visited codeword of each length.

#### Fundamental Domain and Probabilities  
  
Function for computing distributions of words in the fundamental domain W(l) for a given profile l are provided in [weights.py](weights.py).

---
## Extending the C++ Core

The "plumbing" from C++ to python is done via the `ctypes` library for `python`. This require boilerplates both in the `middleware.py` file defining the `CodeRedLib` class, and within the 
``` c++
extern "C" {
...
}
``` 
statement at the end of the `coderedlib.cpp` file.
