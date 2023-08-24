This project is a C++ implementation of my paper "[Improved Algorithms for Integer Complexity](https://arxiv.org/pdf/2308.10301.pdf)" for computing the integer complexity $f(n)$ of a single integer $n$.

# Complexity
The theoretical running time of the single-target algorithm is $O(n^{0.6154})$. For our applications (searching for counterexamples for various number-theoretical conjectures) it is usually much faster, because of pruning.


# Conjectures
1. $f(2^i)=2i$?
2. $f(2^i3^j5^k)=2i+3j+5k$ ($k\leq 5$)?
3. $f(p^i)=i\cdot f(p)$, for $p=577,811,109,433,163,487,2$?


# Implementations
For integers within long long ($\sim 2^{63}$): see [integer_complexity_single_2i3j5k.cpp](https://github.com/hqztrue/integer_complexity/blob/main/integer_complexity_single_2i3j5k.cpp).  
For larger integers, see [single_new.cpp](https://github.com/hqztrue/integer_complexity/blob/main/C-Quadratic-Sieve_C%2B%2B/single_new.cpp). It used [C-Quadratic-Sieve](https://github.com/michel-leonard/C-Quadratic-Sieve) for factoring.  


# License
This project is shared under the terms of the **GNU General Public License**.

