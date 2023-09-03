This project is a C++ implementation of my paper "Improved Algorithms for Integer Complexity"[1] for computing the integer complexity $f(n)$ of a single integer $n$.  
For more context, see [OEIS A005245](https://oeis.org/A005245).

[[1]](https://arxiv.org/pdf/2308.10301.pdf) Qizheng He, Improved Algorithms for Integer Complexity, arXiv:2308.10301 [cs.DS], 2023.


# Complexity
The theoretical running time of the single-target algorithm is $O(n^{0.6154})$. For our applications (searching for counterexamples for various number-theoretical conjectures) it is usually much faster, because of practical pruning strategies.


# Conjectures
1. $f(2^i)=2i$?
2. $f(2^i3^j5^k)=2i+3j+5k$ ($k\leq 5$)?
3. $f(2^i+1)=2i+1$? (except 3 and 9)
4. $f(p^i)=i\cdot f(p)$, for $p=577,811,109,433,163,487,2$?


# Implementations
For integers within long long ($\sim 2^{63}$): see [integer_complexity_single_2i3j5k.cpp](https://github.com/hqztrue/integer_complexity/blob/main/integer_complexity_single_2i3j5k.cpp).  

For larger integers, see [single_new.cpp](https://github.com/hqztrue/integer_complexity/blob/main/C-Quadratic-Sieve_C%2B%2B/single_new.cpp). It used [C-Quadratic-Sieve](https://github.com/michel-leonard/C-Quadratic-Sieve) for factoring.  
Run [compile.bat]() to compile.

# License
This project is shared under the terms of the **GNU General Public License**.
