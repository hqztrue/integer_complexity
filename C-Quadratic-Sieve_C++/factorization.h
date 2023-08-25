typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned __int128 u128;

#pragma region template

#pragma GCC optimize("O3")
#pragma GCC target("avx2")
#pragma GCC optimize("fast-math")
#pragma GCC optimize("unroll-loops")

#define _GNU_SOURCE
#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef   int8_t      i8;
typedef   int16_t     i16;
typedef   int32_t     i32;
typedef   int64_t     i64;
typedef __int128_t    i128;
typedef   uint8_t     u8;
typedef   uint16_t    u16;
typedef   uint32_t    u32;
typedef   uint64_t    u64;
typedef __uint128_t   u128;
typedef   float       f32;
typedef   double      f64;
typedef   long double f80;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SWAP(a, b)          \
            do {            \
                (a) ^= (b); \
                (b) ^= (a); \
                (a) ^= (b); \
            }               \
            while (0);
#define CTZ32(a)       ((a) ? __builtin_ctz((a)) : (32))
#define CTZ64(a)       ((a) ? __builtin_ctzll((a)) : (64))
#define CLZ32(a)       ((a) ? __builtin_clz((a)) : (32))
#define CLZ64(a)       ((a) ? __builtin_clzll((a)) : (64))
#define POPCNT32(a)    ((a) ? __builtin_popcount((a)) : (0))
#define POPCNT64(a)    ((a) ? __builtin_popcountll((a)) : (0))
#define BIT_WIDTH32(a) ((32) - CLZ32((a)))
#define BIT_WIDTH64(a) ((64) - CLZ64((a)))
#define BIT_FLOOR32(a) ((a) ? ((1u) << (BIT_WIDTH32((a)) - (1))) : (0))
#define BIT_FLOOR64(a) ((a) ? ((1ul) << (BIT_WIDTH64((a)) - (1))) : (0))
#define BIT_CEIL32(a)  (((a) <= 1) ? (1u) : ((1u) << BIT_WIDTH32((a) - (1))))
#define BIT_CEIL64(a)  (((a) <= 1) ? (1ul) : ((1ul) << BIT_WIDTH64((a) - (1))))
#define LSBit(a)       ((a) & (-(a)))
#define CLSBit(a)      ((a) & ((a) - (1)))
#define HAS_SINGLE_BIT32(a) (((a) != (0)) && (CLSBit((a)) == (0)))
#define HAS_SINGLE_BIT64(a) (((a) != (0)) && (CLSBit((a)) == (0)))
    #define _ROTL32_INNER(x, l) (((x) << (l)) | ((x) >> ((-l) & (31))))
    #define _ROTR32_INNER(x, r) (((x) >> (r)) | ((x) << ((-r) & (31))))
    #define _ROTL64_INNER(x, l) (((x) << (l)) | ((x) >> ((-l) & (63))))
    #define _ROTR64_INNER(x, r) (((x) >> (r)) | ((x) << ((-r) & (63))))
#define ROTR32(x, r) (((r) < (0)) ? (_ROTL32_INNER((x), ((u64)(-r) % (32)))) : (_ROTR32_INNER((x), ((r) % (32)))))
#define ROTL32(x, l) ROTR32((x), (-l))
#define ROTR64(x, r) (((r) < (0)) ? (_ROTL64_INNER((x), ((u64)(-r) % (64)))) : (_ROTR64_INNER((x), ((r) % (64)))))
#define ROTL64(x, l) ROTR64((x), (-l))

#pragma region gcd

u32 gcd32(u32 a, u32 b)
{
    if (!a || !b)
        return a | b;
    u32 sh = CTZ32(a | b);
    a >>= CTZ32(a);
    do
    {
        b >>= CTZ32(b);
        if (a > b)
            SWAP(a, b);
        b -= a;
    } while (b);
    return a << sh;
}

u64 gcd64(u64 a, u64 b)
{
    if (!a || !b)
        return a | b;
    u64 sh = CTZ64(a | b);
    a >>= CTZ64(a);
    do
    {
        b >>= CTZ64(b);
        if (a > b)
            SWAP(a, b);
        b -= a;
    } while (b);
    return a << sh;
}

#pragma endregion gcd

#pragma region integer square root

int isqrt(u64 n)
{
    if (n == 0)
        return 0;

    u64 a = n;
    u64 b = 1;
    u64 c;

    if ((c = (a >> 32)) != 0)
    {
        a = c;
        b <<= 16;
    }
    if ((c = (a >> 16)) != 0)
    {
        a = c;
        b <<= 8;
    }
    if ((c = (a >> 8)) != 0)
    {
        a = c;
        b <<= 4;
    }
    if ((c = (a >> 4)) != 0)
    {
        a = c;
        b <<= 2;
    }
    if ((c = (a >> 2)) != 0)
    {
        a = c;
        b <<= 1;
    }

    if (a <= 1)
        b += b >> 1;
    else
        b <<= 1;

    do
    {
        a = b;
        b = (b + n / b) >> 1;
    } while (b < a);

    return (int)a;
}

#pragma endregion integer square root

#pragma region mod inverse

typedef struct
{
    i32 f;
    i32 s;
    u32 t;
} Bezout32;

Bezout32 bezout32(u32 x, u32 y)
{
    bool swap = x < y;
    if (swap)
        SWAP(x, y);
    if (y == 0)
    {
        if (x == 0)
            return (Bezout32){0, 0, 0};
        else if (swap)
            return (Bezout32){0, 1, x};
        else
            return (Bezout32){1, 0, x};
    }
    i32 s0 = 1, s1 = 0, t0 = 0, t1 = 1;
    while (true)
    {
        u32 q = x / y, r = x % y;
        if (r == 0)
        {
            if (swap)
                return (Bezout32){t1, s1, y};
            else
                return (Bezout32){s1, t1, y};
        }
        i32 s2 = s0 - (i32)(q)*s1, t2 = t0 - (i32)(q)*t1;
        x = y, y = r;
        s0 = s1, s1 = s2, t0 = t1, t1 = t2;
    }
}

u32 mod_inverse32(u32 x, u32 mod)
{
    assert(gcd32(x, mod) == 1);
    Bezout32 b = bezout32(x, mod);
    assert(b.t == 1);
    return b.f < 0 ? mod + b.f : (u32)b.f;
}

typedef struct
{
    i64 f;
    i64 s;
    u64 t;
} Bezout64;

Bezout64 bezout64(u64 x, u64 y)
{
    bool swap = x < y;
    if (swap)
        SWAP(x, y);
    if (y == 0)
    {
        if (x == 0)
            return (Bezout64){0, 0, 0};
        else if (swap)
            return (Bezout64){0, 1, x};
        else
            return (Bezout64){1, 0, x};
    }
    i64 s0 = 1, s1 = 0, t0 = 0, t1 = 1;
    while (true)
    {
        u64 q = x / y, r = x % y;
        if (r == 0)
        {
            if (swap)
                return (Bezout64){t1, s1, y};
            else
                return (Bezout64){s1, t1, y};
        }
        i64 s2 = s0 - (i64)(q)*s1, t2 = t0 - (i64)(q)*t1;
        x = y, y = r;
        s0 = s1, s1 = s2, t0 = t1, t1 = t2;
    }
}

u64 mod_inverse64(u64 x, u64 mod)
{
    assert(gcd64(x, mod) == 1);
    Bezout64 b = bezout64(x, mod);
    assert(b.t == 1);
    return b.f < 0 ? mod + b.f : (u64)b.f;
}

#pragma endregion mod inverse

#pragma region quadratic residue

int jacobi_symbol32(i32 a, i32 n)
{
    int j = 1;
    while (a)
    {
        if (a < 0)
        {
            a = -a;
            if ((n & 3) == 3)
                j = -j;
        }
        int s = CTZ32(a);
        a >>= s;
        if (((n & 7) == 3 || (n & 7) == 5) && (s & 1))
            j = -j;
        if ((a & n & 3) == 3)
            j = -j;
        SWAP(a, n);
        a %= n;
        if (a > n / 2)
            a -= n;
    }
    return n == 1 ? j : 0;
}

int jacobi_symbol64(i64 a, i64 n)
{
    int j = 1;
    while (a)
    {
        if (a < 0)
        {
            a = -a;
            if ((n & 3) == 3)
                j = -j;
        }
        int s = CTZ64(a);
        a >>= s;
        if (((n & 7) == 3 || (n & 7) == 5) && (s & 1))
            j = -j;
        if ((a & n & 3) == 3)
            j = -j;
        SWAP(a, n);
        a %= n;
        if (a > n / 2)
            a -= n;
    }
    return n == 1 ? j : 0;
}

#pragma endregion quadratic residue

#pragma region RNGs

u32 lcg_rand(void)
{
    static u64 lcg_state = 14534622846793005ull;
    lcg_state = 6364136223846793005ull * lcg_state + 1442695040888963407ull;
    return (u32)lcg_state;
}
u32 lcg_range(u32 l, u32 r)
{
    return l + lcg_rand() % (r - l + 1);
}
f32 lcg_randf(void)
{
    u32 a = 0x3F800000u | (lcg_rand() >> 9);
    return (*((f32 *)(&a))) - 1;
}

u32 pcg_rand(void)
{
    static u64 pcg_state = 0x853c49e6748fea9bull;
    u64 t = pcg_state;
    pcg_state = t * 0x5851f42d4c957f2dull + 0xda3e39cb94b95bdbull;
    u32 sh = ((t >> 18u) ^ t) >> 27u;
    u32 ro = t >> 59u;
    return (sh >> ro) | (sh << ((-ro) & 31));
}
u32 pcg_range(u32 l, u32 r)
{
    return l + pcg_rand() % (r - l + 1);
}
f32 pcg_randf(void)
{
    u32 a = 0x3F800000u | (pcg_rand() >> 9);
    return (*((f32 *)(&a))) - 1;
}

u64 msws_rand(void)
{
    static u64 msws_state1 = 0;
    static u64 msws_state2 = 0;
    static u64 msws_state3 = 0xb5ad4eceda1ce2a9ul;
    static u64 msws_state4 = 0;
    static u64 msws_state5 = 0;
    static u64 msws_state6 = 0x278c5a4d8419fe6bul;
    u64 ret;
    msws_state1 *= msws_state1;
    ret = msws_state1 += (msws_state2 += msws_state3);
    msws_state1 = (msws_state1 >> 32) | (msws_state1 << 32);
    msws_state4 *= msws_state4;
    msws_state4 += (msws_state5 += msws_state6);
    msws_state4 = (msws_state4 >> 32) | (msws_state4 << 32);
    return ret ^ msws_state4;
}
u64 msws_range(u64 l, u64 r)
{
    return l + msws_rand() % (r - l + 1);
}
f64 msws_randf(void)
{
    u64 a = 0x3FF0000000000000ull | (msws_rand() >> 12);
    return (*((f64 *)(&a))) - 1;
}

u64 xrsr128p_rand(void)
{
    static u64 xrsr128p_state1 = 0x1ull;
    static u64 xrsr128p_state2 = 0x2ull;

    const u64 s0 = xrsr128p_state1;
    u64 s1 = xrsr128p_state2;
    const u64 ret = s0 + s1;

    s1 ^= s0;
    xrsr128p_state1 = ROTL64(s0, 24) ^ s1 ^ (s1 << 16);
    xrsr128p_state2 = ROTL64(s1, 37);

    return ret;
}
u64 xrsr128p_range(u64 l, u64 r)
{
    return l + xrsr128p_rand() % (r - l + 1);
}
f64 xrsr128p_randf(void)
{
    u64 a = 0x3FF0000000000000ull | (xrsr128p_rand() >> 12);
    return (*((f64 *)(&a))) - 1;
}

#pragma endregion RNGs

#pragma region sort

void comb_sort11(size_t a_len, u64 *a)
{
    size_t g = a_len;
    while (true)
    {
        bool no_swap = true;
        g = (g * 10) / 13 > 1 ? (g * 10) / 13 : 1;
        if (g == 9 || g == 10)
            g = 11;
        for (size_t i = 0; i + g < a_len; ++i)
        {
            if (a[i] > a[i + g])
            {
                SWAP(a[i + g], a[i]);
                no_swap = false;
            }
        }
        if (g == 1 && no_swap)
            break;
    }
}

#pragma endregion sort

// clang-format off

#pragma region Montgomery ModInt

static u32 N_32, N2_32, NI_32, R1_32, R2_32, R3_32;
void Montgomery32(u32 mod)
{
    assert(mod < 1073741824u);
    N_32 = mod;
    N2_32 = mod << 1;
    NI_32 = mod;
    NI_32 *= 2 - NI_32 * mod;
    NI_32 *= 2 - NI_32 * mod;
    NI_32 *= 2 - NI_32 * mod;
    NI_32 *= 2 - NI_32 * mod;
    R1_32 = (u32)(i32)-1 % mod + 1;
    R2_32 = (u64)(i64)-1 % mod + 1;
    R3_32 = (u32)(((u64)R1_32 * (u64)R2_32) % mod);
}
u32 mr32(u64 A) { u32 y = (u32)(A >> 32) - (u32)(((u64)((u32)A * NI_32) * N_32) >> 32); return (i32)y < 0 ? y + N_32 : y;}
u32 To32(u32 a) { return mr32((u64)a * R2_32); }
u32 From32(u32 A) { return mr32((u64)A); }
u32 Add32(u32 A, u32 B) { A += B - N2_32; A += N2_32 & -(A >> 31); return A; }
u32 Sub32(u32 A, u32 B) { A -= B; A += N2_32 & -(A >> 31); return A; }
u32 SAdd32(u32 A, u32 B) { A += B; A -= (A >= N_32 ? N_32 : 0); return A; }
u32 SSub32(u32 A, u32 B) { A += (A < B ? N_32 : 0); A -= B; return A; }
u32 Min32(u32 A) { return SSub32(0, A); }
u32 Mul32(u32 A, u32 B) { return mr32((u64)A * B); }
u32 Square32(u32 A) { return mr32((u64)A * A); }
u32 Twice32(u32 A) { return (A <<= 1) >= N_32 ? A - N_32 : A; }
u32 Power32(u32 A, size_t k) { return k ? Mul32(Power32(Square32(A), k >> 1), k & 1 ? A : R1_32) : R1_32; }
u32 Inverse32(u32 A) { return mr32((u64)R3_32 * mod_inverse32(A, N_32)); }
u32 Div32(u32 A, u32 B) { return Mul32(A, Inverse32(B)); }
u32 Half32(u32 A) { return (A & 1) ? ((A >> 1) + (N_32 >> 1) + 1) : (A >> 1); }
int Equal32(u32 A, u32 B) { return (((A >= N_32) ? (A - N_32) : A) == ((B >= N_32) ? (B - N_32) : B)) ? 1 : 0; }
int NotEqual32(u32 A, u32 B) { return (((A >= N_32) ? (A - N_32) : A) != ((B >= N_32) ? (B - N_32) : B)) ? 1 : 0; }

static u64 N_64, N2_64, NI_64, R1_64, R2_64, R3_64;
void Montgomery64(u64 mod)
{
    assert(mod < 4611686018427387904ull);
    N_64 = mod;
    N2_64 = mod << 1;
    NI_64 = mod;
    NI_64 *= 2 - NI_64 * mod;
    NI_64 *= 2 - NI_64 * mod;
    NI_64 *= 2 - NI_64 * mod;
    NI_64 *= 2 - NI_64 * mod;
    NI_64 *= 2 - NI_64 * mod;
    R1_64 = (u64)(i64)-1 % mod + 1;
    R2_64 = (u128)(i128)-1 % mod + 1;
    R3_64 = (u64)(((u128)R1_64 * (u128)R2_64) % mod);
}
u64 mr64(u128 A) { u64 y = (u64)(A >> 64) - (u64)(((u128)((u64)A * NI_64) * N_64) >> 64); return (i64)y < 0 ? y + N_64 : y; }
u64 To64(u64 a) { return mr64((u128)a * R2_64); }
u64 From64(u64 A) { return mr64((u128)A); }
u64 Add64(u64 A, u64 B) { A += B - N2_64; A += N2_64 & -(A >> 63); return A; }
u64 Sub64(u64 A, u64 B) { A -= B; A += N2_64 & -(A >> 63); return A; }
u64 SAdd64(u64 A, u64 B) { A += B; A -= (A >= N_64 ? N_64 : 0); return A; }
u64 SSub64(u64 A, u64 B) { A += (A < B ? N_64 : 0); A -= B; return A; }
u64 Min64(u64 A) { return SSub64(0, A); }
u64 Mul64(u64 A, u64 B) { return mr64((u128)A * B); }
u64 Square64(u64 A) { return mr64((u128)A * A); }
u64 Twice64(u64 A) { return (A <<= 1) >= N_64 ? A - N_64 : A; }
u64 Power64(u64 A, size_t k) { return k ? Mul64(Power64(Square64(A), k >> 1), k & 1 ? A : R1_64) : R1_64; }
u64 Inverse64(u64 A) { return mr64((u128)R3_64 * mod_inverse64(A, N_64)); }
u64 Div64(u64 A, u64 B) { return Mul64(A, Inverse64(B)); }
u64 Half64(u64 A) { return (A & 1) ? ((A >> 1) + (N_64 >> 1) + 1) : (A >> 1); }
int Equal64(u64 A, u64 B) { return (((A >= N_64) ? (A - N_64) : A) == ((B >= N_64) ? (B - N_64) : B)) ? 1 : 0; }
int NotEqual64(u64 A, u64 B) { return (((A >= N_64) ? (A - N_64) : A) != ((B >= N_64) ? (B - N_64) : B)) ? 1 : 0; }


#pragma endregion Montgomery ModInt

#pragma region Barrett ModInt

u64 m_b64;
u64 im_b64;
u64 divrem64[2] = {0};
void new_br64(u64 m) { m_b64 = m; im_b64 = (~((u64)0ul)) / m; }
void div_rem_br64(u64 lhs) { if (m_b64 == 1) { divrem64[0] = lhs; divrem64[1] = 0; return; } u64 q = (u64)(((u128)lhs * (u128)im_b64) >> 64); u64 r = lhs - q * m_b64; if (m_b64 <= r) { r -= m_b64; q += 1ul; } divrem64[0] = q; divrem64[1] = r; }
u32 add_br32(u32 a, u32 b) { a += b; a -= (a >= (u32)m_b64 ? (u32)m_b64 : 0); return a; }
u32 sub_br32(u32 a, u32 b) { a += (a < b ? (u32)m_b64 : 0); a -= b; return a; }
u32 mul_br32(u32 a, u32 b) { div_rem_br64((u64)a * b); return (u32)divrem64[1]; }
u32 sqr_br32(u32 a) { div_rem_br64((u64)a * a); return (u32)divrem64[1]; }
u32 pow_br32(u32 a, u32 k) { return k ? mul_br32(pow_br32(sqr_br32(a), k >> 1), k & 1 ? a : 1) : 1; }

u128 m_b128;
u128 im_b128;
u128 divrem128[2] = {0};
void new_br128(u128 m) { m_b128 = m; im_b128 = (~((u128)0ull)) / m; }
void div_rem_br128(u128 lhs) { if (m_b128 == 1) { divrem128[0] = lhs; divrem128[1] = 0; return; } u128 t = (lhs >> 64) * (im_b128 >> 64); u128 x = ((lhs & 0xffffffffffffffffull) * (im_b128 & 0xffffffffffffffffull)) >> 64; u8 flag; u128 auil = (lhs >> 64) * (im_b128 & 0xffffffffffffffffull); if (auil <= (u128)((i128)(-1L)) - x) flag = 0; else flag = 1; x += auil; t += flag; u128 aliu = (lhs & 0xffffffffffffffffull) * (im_b128 >> 64); if (aliu <= (u128)((i128)(-1L)) - x) flag = 0; else flag = 1; x += aliu; t += flag; u128 q = t + (x >> 64); u128 r = lhs - q * m_b128; if (m_b128 <= r) { r -= m_b128; q += 1; } divrem128[0] = q; divrem128[1] = r; }
u64 add_br64(u64 a, u64 b) { a += b; a -= (a >= (u64)m_b128 ? (u64)m_b128 : 0); return a; }
u64 sub_br64(u64 a, u64 b) { a += (a < b ? (u64)m_b128 : 0); a -= b; return a; }
u64 mul_br64(u64 a, u64 b) { div_rem_br128((u128)a * b); return (u64)divrem128[1]; }
u64 sqr_br64(u64 a) { div_rem_br128((u128)a * a); return (u64)divrem128[1]; }
u64 pow_br64(u64 a, u64 k) { return k ? mul_br64(pow_br64(sqr_br64(a), k >> 1), k & 1 ? a : 1) : 1; }

#pragma endregion Barrett ModInt

// clang-format on

#pragma region prime number

#pragma region primality test

bool miller_rabin(u64 n, size_t base_len, u64 *bases)
{
    u64 s = CTZ64(n - 1);
    u64 d = (n - 1) >> s;
    Montgomery64(n);
    for (int i = 0; i < base_len; ++i)
    {
        if (n <= bases[i])
            return true;
        u64 a = Power64(To64(bases[i]), d);
        if (a == R1_64)
            continue;
        u64 r = 1;
        while (a != n - R1_64)
        {
            if (r == s)
                return false;
            a = Square64(a);
            r++;
        }
    }
    return true;
}

bool miller_rabin_br(u64 n)
{
    u64 s = CTZ64(n - 1);
    u64 d = (n - 1) >> s;
    new_br128((u128)n);
    u64 bases[7] = {2ul, 325ul, 9375ul, 28178ul, 450775ul, 9780504ul, 1795265022ul};
    for (int i = 0; i < 7; ++i)
    {
        if (n <= bases[i])
            return true;
        u64 a = pow_br64(bases[i], d);
        if (a == 1)
            continue;
        u64 r = 1;
        while (a != n - 1)
        {
            if (r == s)
                return false;
            a = sqr_br64(a);
            r++;
        }
    }
    return true;
}

bool baillie_psw(u64 n)
{
    assert(n < 4611686018427387904ull);
    if (n < 64ul)
        return (1ull << n) & 2891462833508853932ull;
    if (!(n & 1))
        return false;
    Montgomery64(n);
    {
        u64 d = (n - 1) << CLZ64(n - 1);
        u64 t = Twice64(R1_64);
        for (d <<= 1; d; d <<= 1)
        {
            t = Square64(t);
            if (d >> 63)
                t = Twice64(t);
        }
        if (t != R1_64)
        {
            u64 x = LSBit(n - 1);
            u64 rev = Min64(R1_64);
            for (x >>= 1; t != rev; x >>= 1)
            {
                if (x == 0)
                    return false;
                t = Square64(t);
            }
        }
    }
    {
        i64 D = 5;
        for (int i = 0; jacobi_symbol64(D, n) != -1 && i < 64; ++i)
        {
            if (i == 32)
            {
                u64 sqrt_n = (u64)sqrtl((f80)n);
                return sqrt_n * sqrt_n == n;
            }
            if (i & 1)
                D -= 2;
            else
                D += 2;
            D = -D;
        }
        u64 Q = To64((D < 0) ? ((1 - D) / 4 % n) : (n - (D - 1) / 4 % n));
        u64 u = R1_64, v = R1_64, Qn = Q;
        u64 k = (n + 1) << CLZ64(n + 1);
        D %= (i64)n;
        D = To64((D < 0) ? (D + n) : D);
        for (k <<= 1; k; k <<= 1)
        {
            u = Mul64(u, v);
            v = SSub64(Square64(v), SAdd64(Qn, Qn));
            Qn = Square64(Qn);
            if (k >> 63)
            {
                u64 uu = SAdd64(u, v);
                uu = Half64(uu);
                v = Half64(SAdd64(Mul64(D, u), v));
                u = uu;
                Qn = Mul64(Qn, Q);
            }
        }
        if (u == 0 || v == 0)
            return true;
        u64 x = (n + 1) & ~n;
        for (x >>= 1; x; x >>= 1)
        {
            u = Mul64(u, v);
            v = SSub64(Square64(v), SAdd64(Qn, Qn));
            if (v == 0)
                return true;
            Qn = Square64(Qn);
        }
    }
    return false;
}

bool is_prime(u64 n)
{
    if (n < 64ul)
        return (1ull << n) & 2891462833508853932ull;
    if (!(n & 1))
        return false;
    u64 bases[12] = {2,3,5,7,11,13,17,19,23,29,31,37};
    return miller_rabin(n, 12, bases);
    if (n < 4759123141ul)
    {
        u64 bases[3] = {2ul, 7ul, 61ul};
        return miller_rabin(n, 3, bases);
    }
    else if (n < 4611686018427387904ull)
    {
        return baillie_psw(n);  //err?
    }
    else
    {
        return miller_rabin_br(n);
    }
}

#pragma endregion primality test

#pragma region prime factorization

u32 pollard_brent_rho_32(u32 n, u32 c)
{
    Montgomery32(n);
    const u32 one = R1_32;
    const u32 two = To32(2u);
    const u32 cc = To32(c);
    const u32 m = 1ull << ((31 - CLZ32(n)) / 5);

    u32 x = one, y = two, z = one, q = one;
    u32 g = 1u;

    for (u32 r = 1; g == 1; r <<= 1)
    {
        x = y;
        for (int i = 0; i < r; ++i)
            y = Add32(Square32(y), cc);
        for (u32 k = 0; k < r && g == 1u; k += m)
        {
            z = y;
            for (int _ = 0; _ < MIN(m, r - k); ++_)
            {
                y = Add32(Square32(y), cc);
                q = Mul32(q, Sub32(x, y));
            }
            g = gcd32(From32(q), n);
        }
    }
    if (g == n)
    {
        do
        {
            z = Add32(Square32(z), cc);
            g = gcd32(From32(Sub32(x, z)), n);
        } while (g == 1);
    }
    return g;
}

u64 pollard_brent_rho_64(u64 n, u64 c)
{
    assert(n < 4611686018427387904ull);
    Montgomery64(n);
    const u64 one = R1_64;
    const u64 two = To64(2ull);
    const u64 cc = To64(c);
    const u64 m = 1ull << ((63 - CLZ64(n)) / 5);

    u64 x = one, y = two, z = one, q = one;
    u64 g = 1ull;

    for (u64 r = 1; g == 1; r <<= 1)
    {
        x = y;
        for (int i = 0; i < r; ++i)
            y = Add64(Square64(y), cc);
        for (u64 k = 0; k < r && g == 1ull; k += m)
        {
            z = y;
            for (int _ = 0; _ < MIN(m, r - k); ++_)
            {
                y = Add64(Square64(y), cc);
                q = Mul64(q, Sub64(x, y));
            }
            g = gcd64(From64(q), n);
        }
    }
    if (g == n)
    {
        do
        {
            z = Add64(Square64(z), cc);
            g = gcd64(From64(Sub64(x, z)), n);
        } while (g == 1);
    }
    return g;
}

typedef struct
{
    u64 x;
    u64 z;
} MontgomeryCoordinates;

u64 check(MontgomeryCoordinates p)
{
    return gcd64(From64(p.z), N_64);
}

static u64 a24;

MontgomeryCoordinates create_curve_and_point(void)
{
    while (true)
    {
        u64 a = msws_range(0, N_64 - 1ull);
        u64 x = msws_range(0, N_64 - 1ull);
        u64 m1 = R1_64;
        u64 y2 = Mul64(x, Add64(Mul64(x, Add64(x, a)), m1));
        if (jacobi_symbol64(From64(y2), N_64) == 1)
        {
            a24 = Div64(Add64(a, To64(2ull)), To64(4ull));
            return (MontgomeryCoordinates){x, m1};
        }
    }
}

MontgomeryCoordinates double_EC(MontgomeryCoordinates p)
{
    u64 s = Add64(p.x, p.z);
    u64 d = Sub64(p.x, p.z);
    u64 s2 = Square64(s);
    u64 d2 = Square64(d);
    u64 t = Sub64(s2, d2);
    u64 new_x = Mul64(s2, d2);
    u64 new_z = Mul64(t, Add64(d2, Mul64(a24, t)));
    return (MontgomeryCoordinates){new_x, new_z};
}

MontgomeryCoordinates add_EC(MontgomeryCoordinates p, MontgomeryCoordinates q, MontgomeryCoordinates diff)
{
    u64 u = Mul64(Sub64(p.x, p.z), Add64(q.x, q.z));
    u64 v = Mul64(Add64(p.x, p.z), Sub64(q.x, q.z));
    u64 upv = Add64(u, v);
    u64 umv = Sub64(u, v);
    u64 new_x = Mul64(Mul64(diff.z, upv), upv);
    u64 new_z = Mul64(Mul64(diff.x, umv), umv);
    return (MontgomeryCoordinates){new_x, new_z};
}

MontgomeryCoordinates power_EC(MontgomeryCoordinates p, size_t k)
{
    MontgomeryCoordinates p0 = p;
    MontgomeryCoordinates p1 = double_EC(p);
    for (int b = BIT_WIDTH64(k) - 2; b >= 0; --b)
    {
        MontgomeryCoordinates tmp = add_EC(p1, p0, p);
        if ((k >> b) & 1)
        {
            p1 = double_EC(p1);
            p0 = tmp;
        }
        else
        {
            p0 = double_EC(p0);
            p1 = tmp;
        }
    }
    return p0;
}

u64 ecm(u64 n)
{
    for (int k = 2; k < BIT_WIDTH64(n); ++k)
    {
        u64 r = roundl(powl(n, 1.0l / k));
        u64 pw = r;
        for (int i = 1; i < k; ++i)
            pw *= r;
        if (pw == n)
            return r;
    }

    Montgomery64(n);
    u64 ecm_blocks[10] = {5690199479092128000ull, 810162134158954261ull, 326580695497527083ull, 13784092967194631821ull, 1107997261359193637ull, 6532397423431938467ull, 96265407405451883ull, 260006624961107813ull, 707992818804600227ull, 22417030981ull};

    while (true)
    {
        MontgomeryCoordinates point = create_curve_and_point();
        u64 f = 1ull;

        for (size_t block = 0; block < 10; ++block)
        {
            MontgomeryCoordinates new_point = power_EC(point, ecm_blocks[block]);
            f = check(new_point);
            if (f != 1)
            {
                if (f != N_64)
                    return f;
                else
                    break;
            }
            point = new_point;
        }
        if (f == N_64)
            continue;
        MontgomeryCoordinates six = double_EC(add_EC(double_EC(point), point, point));
        MontgomeryCoordinates q0 = six;
        MontgomeryCoordinates q1 = double_EC(six);
        for (int _ = 6; _ < 400; _ += 6)
        {
            q0 = add_EC(q1, six, q0);
            SWAP(q0.x, q1.x);
            SWAP(q0.z, q1.z);
        }
        u64 xprod = R1_64;
        u64 x_norm = Div64(point.x, point.z);
        for (int i = 396; i < 3000; i += 6)
        {
            xprod = Mul64(xprod, Sub64(q0.x, Mul64(q0.z, x_norm)));
            if (i % 300 == 0)
            {
                f = gcd64(From64(xprod), N_64);
                if (f != 1)
                {
                    if (f != N_64)
                        return f;
                    else
                        break;
                }
            }
            q0 = add_EC(q1, six, q0);
            SWAP(q0.x, q1.x);
            SWAP(q0.z, q1.z);
        }
        if (f == 1)
        {
            f = gcd64(From64(xprod), N_64);
            if (f != 1 && f != N_64)
                return f;
        }
    }
}

unordered_map<u64,u64> P; 
u64 find_prime_factor(u64 n)
{
	u64 n0=n;
	auto it=P.find(n);
	if (it!=P.end())return it->second;
    if (is_prime(n)){
    	P[n0]=n;
        return n;
    }
    for (int _ = 0; _ < 200; ++_)
    {
        u64 m;
        if (n < 1073741824ull)
        {
            m = pollard_brent_rho_32((u32)n, pcg_range(1u, (u32)n - 1));
            if (is_prime(m)){
            	P[n0]=m;
                return m;
            }
        }
        else if (n < 288230376151711744ull)
        {
            m = pollard_brent_rho_64(n, xrsr128p_range(1ull, n - 1));
            if (is_prime(m)){
            	P[n0]=m;
                return m;
            }
        }
        else
        {
            m = ecm(n);
            if (is_prime(m)){
                P[n0]=m;
                return m;
            }
        }
        n = m;
    }
    return -1;
}

u64 *factorize(u64 n)
{
    u64 *ret = (u64 *)calloc(65, sizeof(u64));
    if (ret == NULL)
    {
        exit(EXIT_FAILURE);
    }
    ret[64] = 0;
    for (int i = 0; i < 64; ++i)
        ret[i] = (u64)(i64)-1;
    u64 s = CTZ64(n);
    n >>= s;
    ret[64] += s;
    for (u64 i = 0; i < s; ++i)
        ret[i] = 2;
    for (u64 i = 3ull; i <= 100ull && i * i <= n; i += 2ull)
    {
        if (n % i == 0)
        {
            do
            {
                n /= i;
                ret[ret[64]++] = i;
            } while (n % i == 0);
        }
    }
    while (n > 1)
    {
        u64 p = find_prime_factor(n);
        do
        {
            n /= p;
            ret[ret[64]++] = p;
        } while (n % p == 0);
    }
    comb_sort11(64, ret);
    return ret;
}

#pragma endregion prime factorization
// https://judge.yosupo.jp/submission/107098

