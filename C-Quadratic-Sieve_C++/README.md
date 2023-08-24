# C Factorization using Quadratic Sieve

Pure C factorizer (*Kraitchik family algorithm*) using self-initialising  **Quadratic Sieve**.

Your ~2500 lines GitHub project :

- is immediately compatible with Microsoft Windows, Linux (no dependancy)
- is a C99 **command line** factorizer from 0 to 300 bits (330 bits were factored in the lab)
- is built so that you can easily use and test the software
- use its own "big num" library named **cint**
- use **[AVL trees](https://en.wikipedia.org/wiki/AVL_tree)** to organize information
- use **[Lanczos Block](https://en.wikipedia.org/wiki/Lanczos_algorithm)**, a pure C iterative matrix eigenvalues finder algorithm
- use **[Pollard's Rho](https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm)** algorithm to answer under 64 bits
- try to benefit from the single [large prime variant](https://www.google.com/search?q=%22Factoring+with+two+large+primes%22+AK+Lenstra) of the quadratic sieve

# GNU General Public License

- As an undergraduate student, this project is part of my computer science + maths training.
- This software proposition is from Michel Leonard (student at Université de Franche-Comté, Mon, 11 Jul 2022)
- There is no guarantee on the software
- C code is shared under the terms of the **GNU General Public License**
- The **main mathematical and logical inspiration source** is located at :
  - [http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp](http://web.mit.edu/sage/export/flintqs-0.0.20070817/QS.cpp) - **GNU General Public License**

This software implementation would have been impossible without "FLINT: Fast Library for Number Theory" maintained by [William Hart](https://github.com/wbhart).

# Usage
If you don't know how to get an executable, try to follow this procedure :

- Ubuntu provide you a C compiler by the command `sudo apt install build-essential`
- On Windows you can install [MinGW](https://winlibs.com/), it will, like as Ubuntu, provide you a **GCC** compiler

With **Terminal** on Ubuntu or **Powershell** on Windows you can go to the downloaded directory ([cd](https://en.wikipedia.org/wiki/Cd_(command)) command) :

- Execute `gcc -Wall -pedantic -O3 main.c -o qs`, GCC will create the right executable for your device

The compilation took a few seconds, you can use the software :

- Execute `./qs 1389767379868103466550952369852270458062404848980444227682523391`

The software will quickly show you the answer :

```
C:\qs.exe 1389767379868103466550952369852270458062404848980444227682523391
(1062949 ^ 2) * 19221289015581393571881329 * 63993335829453341906998279
```
					
# Options

Some options designed for software development haven't been removed, they are not intrusive.

|Option example| Action
|--|--|
| `-h`  | prints the help
| `-l=270`  | causes the quadratic sieve to reject inputs larger than 270-bit, default to 220-bit
| `-t=150`  | execute a minute of factorization test with random 150-bit odd numbers 
| `-m=3`  | causes the quadratic sieve to choose the multiplier 3 for its input N
| `-r=1`  | causes the quadratic sieve to choose the seed 1 for its random number generator
| `-s`  | enables the silent mode, it disables the printing of progress

# RSA factorization

Small and larger RSA numbers have been factored by the software, such as the 100 decimal digit number [**RSA-100**](https://en.wikipedia.org/wiki/RSA_numbers#RSA-100).

|Bits| RSA Number | Took
|--|--|--|
| 130 | `982374584994591973035454918323883152991`  | 0.1 s
| 150 | `1179676342138800493769972880071793674490763037`  | 0.2 s
| 170 | `1107814594796361407351721529681773419585636720020897`  | 1 s
| 190 | `995209482127497644492758995962031762505687007535997155961`  | 3 s
| 210 | `1143938601578848425045957187554857460827103223569512762813742971`  | 15 s
| 230 | `1268631359685752304166485771456206118633849920528997030903788793364933`  | 45 s
| 250 | `1401811817899460116600945074728583412740519573015376930481203561750251051823`  | 4 min
| 270 | `1415606447884291776606783262139201189953436249643759632827004228713595295320953939`  | 13 min

The initial software goal was to **factor 200-bit RSA in 30 seconds**, after which there were fewer situations tested.

- RSA numbers greater than 250 bits have been tested, they have been factorized
- During a test that lasted 2 hours, the software factored a **300-bit RSA** number
- The software factored the **321-bit** RSA number relating to the "[bank card case](https://www.enseignement.polytechnique.fr/profs/informatique/Eric.Goubault/Cours09/qs.pdf)"

# Fermat numbers factorization
|F| Fermat Number | Took  |
|--|--|--|
| 7 | `340282366920938463463374607431768211457`  | 150 ms
| 8 | `115792089237316195423570985008687907853269984665640564039457584007913129639937`  | 3 min

- These tests like software development were made by laptop Honor MagicBook on Windows (64-bit)
- One of the largest number factored during development was the 79 digits 8th Fermat Number

# Mersenne numbers factorization

|Number| Decimal digits | Fully factored in  | Using option |
|--|--|--|--|
| 2^259 - 1 | 78  | 1 min 40 s| **-limit**=250
| 2^360 - 1 | 109  | 3 s|
| 2^468 - 1 | 141  | 40 s| **-limit**=230
| 2^504 - 1 | 152  | 6 min 50 s | **-limit**=270
| 2^630 - 1 | 190  | 55 min | **-limit**=290

Mersenne numbers were factored by the trial division algorithm and then completed by the quadratic sieve.

# Other factorizations

The software, although modest, is designed to be a general purpose factoring solution.

# Random

Thousands of random factorization tests have been performed using :

- a single-core Debian system (Linux) on [OVH-cloud-amd64](https://www.ovhcloud.com/fr/vps/)
- PHP7.4 which took for time measurements and compared data
- [GMP](https://www.php.net/manual/en/book.gmp.php) who provided inputs and checked the factorizations validity
- the numbers provided in text files, among others

# Testing

A basic ~100 line test feature is available for developer convenience :

- the goal is to provide input when testing new configurations
- test durations are between 30 seconds (130-bit) and 3 minutes (200-bit)
- the `./qs test=1` test offers a 2 minute crescendo up to 200+ bits
- the `./qs test=160` test offers a minute of 160-bit random odd numbers

The time measurements provided in this web page are indicative in 2022, not all were taken by the same device.

# Primes

- the usual answers contain prime numbers
- quotes around a number means the software know it's not a prime

# Memory

Memory allocations are reasonably sized, so this project passes pointers to **assert**.

- the program will stop if the memory is refused, showing you an error message
- it would be recommended to restart your device if you see this kind of message

Technical : [valgrind](https://valgrind.org/) usually shows around **10MB** and **80MB** allocated depending on the SIQS input size.

# cint

You have access to the source code of **cint**, the lite "big num" library **cint** is designed to :

- take input from a regular integer
- take input from an "arbitrary" long string in base 2 to 62
- perform basic math operations, including nth_root, modular_inverse, is_prime
- output its content as a string in base from 2 to 62

To perform intermediates computations **cint** does not use global variable, excepting rand it's thread-safe.

# AVL trees

You have access to the source code of **AVL**, a fast self-balancing binary search tree utility  :

- used to store relations and answers
- uses a **height field** and a **parent pointer** for each tree node
- tree can store and retrieve keys in worst case `1.44 * log2 (number of entries)`
- this software implementation use **cint** as tree entries (or keys)
- tested and fast with many more keys than quadratic sieve uses

**cint** and **AVL** are initially separate projects from the quadratic sieve.

# Improvements

Developers can for example :

- replace **cint** with GMP
- replace **AVL** with a hashmap
- better reconfigure the software
- criticize the `preparation_part_3_michel` function
- report any inaccuracies they notice, please

During development you may like to compile your code faster by disabling GCC compiler optimizations :

```
gcc -Wall -pedantic -O0 main.c -o qs
```

- GCC `-O3` optimization improves the speed of the executable compared to `-O0`, it's just slower to compile
- the current source code does not have assembly or pre-processing tricks because it's designed for beginners

# The file main.c

This file only contains a few lines so you can easily perform your refactoring :

```c
int main(void) {
	fac_params config = {0};
	config.silent = 1;
	cint N; // initialize a number N (max 1024-bit) from a string in base 10 to be factored.
	cint_init_by_string(&N, 1024, "1179676342138800493769972880071793674490763037", 10);
	fac_cint **answer = c_factor(&N, &config); // get a factorization of N.
	char *str = fac_answer_to_string(answer); // format the answer as a string.
	puts(str);
	free(str); // release the string.
	free(answer); // release the answer memory.
	free(N.mem); // release N memory.
}
```

# The file fac_utils.c

The file contains utilities that are not specifically intended for a quadratic sieve :

- `c_factor`, takes on a factorization request and coordinates the work
- `pollard_rho`
- `is_prime_1062961`
- `log_computation`
- `multiplication_modulo`
- `power_modulo`
- `kronecker_symbol`
- `tonelli_shanks`
- `modular_inverse`
- `answer_to_string`, translates a factorization result into a string
- `mem_aligned`

# The file fac_quadratic.c

The quadratic sieve file structure is as follows:

- the ~40 lines function that allows to see the algorithm **structure**
- the 2 important loop **conditions**
- the algorithm **parameters**
- the **functions** approximately in the order they are called

### preparation_part_2 .. 3

The input **N is** duplicated and called "**kN**" after this function complete :

- multiply **N** by a prime number to reach 120 bits
- apply a multiplier to **N** intended to optimize runtime

The `preparation_part_3_michel` function provide a multiplier that make the algorithm completes **faster** on average. 

|Variable| Information|
|--|--|
| N | Prime factors are removed from **N** and **N** is updated until `N = 1`  |
| kN | Algorithm computes with **kN** which always remains a constant |

See how multipliers affect the execution of `./qs 144311040110679000656983950712977705460301891839699388207` :

| Option | Makes the algorithm really factoring | Comment on the multiplier | Took |
|:-:|---|:-:|:-:|
|`-m=47`|`kN = 47 * N`|One of the best|5.8 s|
|`-m=41`|`kN = 41 * N`|Current proposition|7.8 s|
|`-m=1`|`kN = 1 * N`|No multiplier|12.5 s|
|`-m=11`|`kN = 11 * N`|One of the worst|27.3 s|

My proposal don't always pick the very best multiplier, so "on average" is fairer. A significant **speed difference** can be noticed between good and no multiplier, so the question ([Knuth-Schroeppel analysis](https://www.google.com/search?q=Knuth+Schroeppel+analysis)) deserves some attention. Your `preparation_part_3` function or concept would be gladly compared and why not integrated to this project.

### qs_parametrize, preparation_part_4

*Good parameters can improve the speed*.

- define the algorithm parameters
- allocates a block of memory for the quadratic sieve computations
- prepare constants, variables, buffers, data arrays ...

There is a struct inside the **qs_sheet** (or manager) called **mem** :

- it holds the  **base** entry point of the allocated memory
- it holds a **now** void* pointer which represent the current available memory

**qs_sheet** holds 3 AVL tree manager :

- one to store the regular relations
- one to store the relations that wait to be paired (partials)
- one to store the known divisors of N

With small precautions you are supposed to be able to store anything in **now**, then to update **now** accordingly to what you stored. **now** is always supposed to contain only zero until its end. This is the main memory management technique used by this software. `mem_align` aims to provide aligned pointers, wasting a few bits if necessary.

### preparation_part_5 .. 6
Fill the manager's **base** array with prime numbers provided [by](https://stackoverflow.com/a/61895974/18765627) a constant expression :

- verify that kN is a square mod prime, ignore it otherwise
- associates the prime with square root of kN mod prime
- associates the prime with its size (log2)
- computes invariants like **D** used to generate the **A** polynomial coefficient

### get_started_iteration

- can restore the relations previously saved by Lanczos algorithm
- can be used to perform analysis, save/restore factorization to file or other periodic actions

### Polynomial `AX^2 + 2BX + C` coefficients :

| coefficient | is a constant after  | is a constant until|
|--|--|--|
| A | `iteration_part_1` |`inner_continuation_condition` completion
| B | `iteration_part_4` |for loop final expression
| C | `iteration_part_6` |for loop final expression

- Polynomial and related data are computed until `iteration_part_6` complete
- The algorithm prepares data that will help generate polynomial values that will be multiples of the factor base
- `iteration_part_7` and `iteration_part_8` are used for sieving

### Search sieve for relations
`register_relations` searches sieve for relations, it reads the sieve. When the function finds interesting to define the variable **X** , calculates the value of the polynomial in **X** then divides the value with the factor base, it thus tries to establish relations.

Buffered knowledge is structured by `register_relation_kind_1` and `register_relation_kind_2` :

- informations potentially useful are saved into a struct **qs_relation***
- `register_relation_kind_1` immediately build a matrix of regular full **qs_relation*** using AVL tree
- `register_relation_kind_2` combine (single large primes variation) **qs_relation*** together using AVL tree

The AVL trees are used to identify duplicates and retrieve their data.

### `register_relation_kind_2` : Large prime variation

*In short : given two partial relation having same large prime, the software merge them to obtain a full relation.*

This step is not necessary to complete a factorization, so it could easily be skipped completely, but it leads to a runtime speedup. To disable this feature it is sufficient to never call the `register_relation_kind_2` function. So, during the execution of `register_relations`, after the trial division by the factor base primes (`cint_remove`), if the quotient is still greater than 1 but less than the large prime bound, then the software register this partial relation with `register_relation_kind_2` which stores a new `qs_relation` in a different location than full relations. After that, if the large prime number collected is equal to a prime number previously collected in another **partial relation** (the AVL tree recognizes it), the software **combines these two relations**, then calls the `register_relation_kind_1` function (which register all regular full relations), and at this point there is no difference between this newly created relation and any other one. A third partial relation sharing the same prime number as the previous two would simply be ignored, as it is more complicated (we must not provide linearly dependent rows to `lanczos_block`). There are solutions to add even more relations (double large primes variation), but I didn't go that far mathematically, this software only try to benefit from the **single large primes variation** of the quadratic sieve. I think it's classic, the operation consumes memory for an expected gain in execution time ; as shown in the table below, this allow us to gather our amount of needed relations faster. Because of this operation the percentage of progression of the quadratic sieve does not evolve linearly, but accelerates progressively.

|N bits| Relations added immediately | Relations added by single large primes variation | Relations needs
|:-:|:-:|:-:|:-:|
| 150 | 2551 | 3| 2554
| 170 | 3168 | 659 | 3827
| 190 | 3212 | 1786 | 4998
| 210 | 3409 | 3012 | 6421
| 230 | 5659 | 4686 | 10345
| 250 | 10051 | 9793 | 19844
| 270 | 14026 | 15484 | 29510
| 290 | 12820 | 20201 | 33021

### inner_continuation_condition
Decides if sieving shoud continue or break, the relation counter is usually the condition.

### lanczos_block

This algorithm find a subset of all exponent vectors such that the sum of their exponent vectors is the zero vector :

- the process need memory, **fac_lanczos.c** have its own array builder
- all memory taken in **mem.now** is zeroed and reusable after calculations
- finding the matrix eigenvalues is usually fast, it can still take 400+ iterations with large N
- the algorithm is probabilistic with high success rate, generally it is lucky enough
- several tries with and without **reduce_matrix** are performed before giving up

The `reduce_matrix` function helps `lanczos_block` and could be called unconditionally, but for now it still first tries not to call it. This was an experience that helped me understand why we should sometimes discard a few relations.

### finalization_part_1 .. 2 .. 3

The `finalization_part_1` function performs the **square root + GCD** step based on the `lanczos_block` answer :

It often leads to a factorization of N but in certain cases such as :
```
./qs 51460938795049063955433175628971167803839994111348342302522016010379
```
when N is not fully factored, `finalization_part_2` completes the factorization with GCDs and perfect power checks.

### outer_continuation_condition
This function is used to decide if the algorithm should return the control, or return to sieving.

- normal case is to return the control, answer was found
- unusual case is when N isn't fully factored (maybe parameters was wrong)

So before giving up, the algorithm searches 10%, 25% then 50% more relations.

# See also

Here are some useful online tools I've used :

- Dario Alpern [Integer factorization calculator](https://www.alpertron.com.ar/ECM.HTM)
- numberempire.com [Factoring Calculator](https://www.numberempire.com/factoringcalculator.php)
- bigprimes.org [RSA numbers generator](https://bigprimes.org/RSA-challenge)

# Thank you
There are many people to thank :

- Carl Pomerance
- Jason Papadopoulos
- William Hart
- my professors at University of Franche-Comté ❤️
- GitHub and SourceForge users reporting issues
