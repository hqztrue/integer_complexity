This project is a C++ implementation of my paper "Improved Algorithms for Integer Complexity"[1] for computing the integer complexity $f(n)$ of a single integer $n$.  
For more context, see [OEIS A005245](https://oeis.org/A005245).

[[1]](https://arxiv.org/pdf/2308.10301.pdf) Qizheng He, Improved Algorithms for Integer Complexity, arXiv:2308.10301 [cs.DS], 2023.


# Conjectures
1. $f(2^i)=2i$?
2. $f(2^i3^j5^k)=2i+3j+5k$ ($k\leq 5$)?
3. $f(2^i+1)=2i+1$? (except 3 and 9)
4. $f(p^i)=i\cdot f(p)$, for $p=109,433,163,487,2$?


# Single-Target Algorithm

## Complexity
The theoretical running time of the single-target algorithm is $O(n^{0.6154})$. For our applications (searching for counterexamples for various number-theoretical conjectures) it is usually much faster, because of practical pruning strategies.

## Implementations
For integers within long long ($\sim 2^{63}$): see [single64.cpp](https://github.com/hqztrue/integer_complexity/blob/main/old/single64.cpp).

For larger integers, see [single.cpp](https://github.com/hqztrue/integer_complexity/blob/main/integer-complexity-128/single.cpp). It used [C-RHO](https://github.com/michel-leonard/C-RHO) and [C-Quadratic-Sieve](https://github.com/michel-leonard/C-Quadratic-Sieve) for factoring. It also supports printing the formulas.

Run [compile.bat]() to compile.


# Upper Bound
We designed a more efficient algorithm that can compute an upper bound for the integer complexity, which holds for a set of numbers with natural density 1. The result is correct with high probability.
[upper_bound_slow.cpp](https://github.com/hqztrue/integer_complexity/blob/main/upper_bound/upper_bound_slow.cpp): $O(n^3/w)$ per sample.

[upper_bound.cpp](https://github.com/hqztrue/integer_complexity/blob/main/upper_bound/upper_bound.cpp): $O(n^2)$ per sample.

# TODO
1. Use a better factorization library.
2. Use parallelization.
3. Improve the lower bound computation.
4. Need u256 and more memory for a significantly larger range.


# License
This project is shared under the terms of the **GNU General Public License**.

