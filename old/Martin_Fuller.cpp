/*
martin_n_fuller@btinternet.com, 01 Feb 2008

Program to calculate A005245, A005520, A005421
  A005245(1) = 1
  A005245(n) = MIN {A005245(x)+A005245(y) : x+y = n or x*y = n}

  A005520(n) = least m with A005245(m) = n
  A005421(n) = number of m with A005245(m) = n



Conjectures to check for counterexamples:

David Wilson [2]: MIN {A005245(x)+A005245(y) : x+y = n} = A005245(n-1)+1

(Counterexample n = 21080618 was suggested by CD_Eater at Russian forum lib.mexmat.ru and verified by Max Alekseyev.  This program confirms it is the smallest.)

David Wilson [3]: A005245(n) = MIN {A005245(n-1)+1} u {A005245(x)+A005245(y) : x*y = n}

(This program finds smallest counterexamples n = 353942783 and n = 516743639)

UPINT section F26: A005245(p) = A005245(p-1)+1 for p prime.

(Smallest counterexample to [3] is also prime, so p = 353942783 is counterexample)



Program to calculate A005245:

1. Create an array large enough to store the sequence values
   and fill it with a fixed upper bound.

2. Set A(1)=1

3. Loop from n=2
  a. The loop starts with all seqeunce terms < n confirmed, 
     and all products where both multiplicands < n checked as new minimum in the array.

  b. Check all additions A(m)+A(n-m). This confirms term n of the sequence.

  c. Check all multiplications from n*2 to n^2, and store any improved minimum in the array.

  d. Now ready to continue to the next loop.


Optimisation of step 3b:
  Want to find A(m)+A(n-m) < best so far.
  Set target = best so far - 1.

  We know that m > A000792(k) implies A(m) > k.
  So if A000792(k)+A000792(target-k) < n then there is no solution with A(m)=k.
  The function A000792(k)+A000792(target-k) is minimised for k = target/2 and increases as k->1.
  After finding the largest k <= target/2 with a potential solution, check all m <= A000792(k).

  In this program we want to find counterexamples to David Wilson's conjecture [2],
  so we relax the target from (best so far)-1 to known value of A(n-1).

  There is no need to check any m with A(m) = A(m-1)+1 as they cannot be better than A(m-1).



Future extension:
  This program is memory intensive, which restricts its range.
  But there is no need to store all the sequence simultaneously because of the optimisation in step 3b.
  All that is required is enough recent terms for the additions, and then n/2, n/3, etc for multiplications.
  
  Solution 1:
  The recent terms could be held in a small array which is sieved as necessary for multiplications.
  This would give an immediate doubling of the range.

  Solution 2:
  This idea could be extended, using more small arrays for n/2, n/3, etc until the quotient is small enough.
  The program will slow down as terms are recalculated, but I expect 10^10 to be an overnight run.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

typedef long long ll;
typedef unsigned char A005245_value_t;
typedef struct { unsigned size; A005245_value_t *array; } A005245_array_t;

#define MAX_A005245_VALUE 127 /* Keep a safety factor of 2 to avoid overflow */


unsigned
A000792(A005245_value_t n)
{
  unsigned result = 1;
  while (n >= 5 || n == 3)
  {
    result *= 3;
    n -= 3;
  }
  return result << (n/2);
}

void
A005245_free(A005245_array_t *a)
{
  a->size = 0;
  free(a->array);
  a->array = 0;
}

int
A005245_init(A005245_array_t *a, unsigned size)
{
  unsigned i;

  //A005245_free(a);
  a->array = (A005245_value_t*)malloc(size * sizeof(A005245_value_t));
  if (a->array)
  {
    a->size = size;
    a->array[1] = 1;
    for (i = 2; i < size; i++)
      a->array[i] = MAX_A005245_VALUE;
  }
  return (a->array != 0);
}

ll S=0;
void
A005245_additions_to_n(A005245_array_t *a, unsigned n)
{
  unsigned limit_m, m;
  A005245_value_t target, k;

  if (a->array[n] > a->array[n-1] + 1)
    a->array[n] = a->array[n-1] + 1;

  target = a->array[n-1];
  /*k = target / 2;
  while (A000792(k) + A000792(target - k) < n)
    k--;*/
  k = 1;
  while (A000792(k) + A000792(target - k) >= n && k<target / 2)
    k++;
  //S+=(target/2-k);
  //S+=k;
  limit_m = A000792(k);
  S+=limit_m;

  /* Already used m=1 earlier, and don't need m=2..5 as they cannot be better than m=1 */
  for (m = 6; m <= limit_m; m++)
  {
    if (a->array[n] > a->array[m] + a->array[n-m])
    {
      printf("Counterexample to [3]: A(%u) + A(%u) = %u < conjecture(%u) = %u\n",
        m, n-m, (unsigned)(a->array[m] + a->array[n-m]), n, (unsigned)a->array[n]);
      a->array[n] = a->array[m] + a->array[n-m];
    }
    else if (a->array[n-1]+1 > a->array[m] + a->array[n-m])
      printf("Counterexample to [2]: A(%u) + A(%u) = %u < A(%u-1)+1 = %u\n",
        m, n-m, (unsigned)(a->array[m] + a->array[n-m]), n, (unsigned)a->array[n-1]+1);
  }
}

void
A005245_multiplications_from_n(A005245_array_t *a, unsigned n)
{
  unsigned m, mn;
  for (m = 2, mn = 2*n; (m <= n) && (mn < a->size); m++, mn += n)
    if (a->array[mn] > a->array[m] + a->array[n])
      a->array[mn] = a->array[m] + a->array[n];
}


int main()
{
  int t1=clock();
  unsigned n = 0;
  A005245_array_t A005245;
  A005245_value_t A005520_n, A000792_n;
  unsigned A000792_value;
  unsigned A005421[MAX_A005245_VALUE+1];

  //if (argc > 1)
  //  sscanf(argv[1], "%u", &n);
  n = 100000000;
  if (n < 2) 
  {
    printf("Syntax: A005245 limit\n");
    return 1;
  }

  if (!A005245_init(&A005245, n))
  {
    A005245_free(&A005245);
    printf("Not enough memory\n");
    return 2;
  }

  A005520_n = 1;

  A000792_n = A000792_value = 2;
  for (n = 0; n <= MAX_A005245_VALUE; n++)
    A005421[n] = 0;

  for (n = 2; n < A005245.size; n++)
  {
    A005245_additions_to_n(&A005245, n);
    A005245_multiplications_from_n(&A005245, n);

    if (A005520_n < A005245.array[n])
    {
      A005520_n = A005245.array[n];
      //printf("A005520(%u) = %u\n", (unsigned)A005520_n, n);
    }

    A005421[A005245.array[n]]++;
    if (n >= A000792_value)
    {
      //printf("A005421(%u) = %u\n", (unsigned)A000792_n, A005421[A000792_n]);
      A000792_n++;
      A000792_value = A000792(A000792_n);
    }
  }

  A005245_free(&A005245);
  printf("%I64d\n",S);
  printf("%d\n",clock()-t1);
  return 0;
}

