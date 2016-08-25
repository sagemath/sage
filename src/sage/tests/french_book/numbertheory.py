"""
Test file for Chapter Number Theory.
"""

r"""
sage: a=IntegerModRing(15)(3); b=IntegerModRing(17)(3); a, b
(3, 3)
sage: a == b
False
sage: R=a.base_ring(); R
Ring of integers modulo 15
sage: R.characteristic()
15
sage: a+a, a-17, a*a+1, a^3
(6, 1, 10, 12)
sage: 1/(a+1)
4
sage: 1/a
Traceback (most recent call last):
...
ZeroDivisionError: Inverse does not exist.
sage: z=lift(a); y=ZZ(a); y, type(y), y==z
(3, <type 'sage.rings.integer.Integer'>, True)
sage: [Mod(x,15).additive_order() for x in range(0,15)]
[1, 15, 15, 5, 15, 3, 5, 15, 15, 5, 3, 15, 5, 15, 15]
sage: [[x,Mod(x,15).multiplicative_order()] for x in range(1,15) if gcd(x,15)==1]
[[1, 1], [2, 4], [4, 2], [7, 4], [8, 4], [11, 2], [13, 4], [14, 2]]
sage: p=10^20+39; mod(2,p).multiplicative_order()
50000000000000000019
sage: mod(3,p).multiplicative_order()
100000000000000000038
sage: n=3^100000; a=n-1; e=100
sage: timeit('(a^e) % n') # random long time
5 loops, best of 3: 387 ms per loop
sage: timeit('power_mod(a,e,n)') # random
125 loops, best of 3: 3.46 ms per loop
sage: R = GF(17); [1/R(x) for x in range(1,17)]
[1, 9, 6, 13, 7, 3, 5, 15, 2, 12, 14, 10, 4, 11, 8, 16]
sage: R = GF(9,name='x'); R
Finite Field in x of size 3^2
sage: R.polynomial()
x^2 + 2*x + 2
sage: Set([r for r in R])
{0, 1, 2, x, x + 1, x + 2, 2*x, 2*x + 1, 2*x + 2}
sage: Q.<x> = PolynomialRing(GF(3))
sage: R2 = GF(9,name='x',modulus=x^2+1); R2
Finite Field in x of size 3^2
sage: p = R(x+1); R2(p)
Traceback (most recent call last):
...
TypeError: unable to coerce from a finite field other than the prime subfield
sage: rational_reconstruction(411,1000)
-13/17
sage: rational_reconstruction(409,1000)
Traceback (most recent call last):
...
ArithmeticError: rational reconstruction of 409 (mod 1000) does not exist
sage: def harmonic(n):
...    return sum([1/x for x in range(1,n+1)])
sage: def harmonic_mod(n,m):
...    return add([1/x % m for x in range(1,n+1)])
sage: def harmonic2(n):
...    q = lcm(range(1,n+1))
...    pmax = RR(q*(log(n)+1))
...    m = ZZ(2*pmax^2)
...    m = ceil(m/q)*q + 1
...    a = harmonic_mod(n,m)
...    return rational_reconstruction(a,m)
sage: harmonic(100) == harmonic2(100)
True
sage: a=2; b=3; m=5; n=7; lambda0=(b-a)/m % n; a+lambda0*m
17
sage: crt(2,3,5,7)
17
sage: def harmonic3(n):
...    q = lcm(range(1,n+1))
...    pmax = RR(q*(log(n)+1))
...    B = ZZ(2*pmax^2)
...    m = 1; a = 0; p = 2^63
...    while m<B:
...        p = next_prime(p)
...        b = harmonic_mod(n,p)
...        a = crt(a,b,m,p)
...        m = m*p
...    return rational_reconstruction(a,m)
sage: harmonic(100) == harmonic3(100)
True
sage: crt(15,1,30,4)
45
sage: crt(15,2,30,4)
Traceback (most recent call last):
...
ValueError: No solution to crt problem since gcd(30,4) does not divide 15-2
sage: p=previous_prime(2^400)
sage: timeit('is_pseudoprime(p)') # random
625 loops, best of 3: 1.07 ms per loop
sage: timeit('is_prime(p)') # random long time
5 loops, best of 3: 485 ms per loop
sage: [560 % (x-1) for x in [3,11,17]]
[0, 0, 0]
sage: def count_primes1(n):
...    return add([1 for p in range(n+1) if is_prime(p)])
sage: def count_primes2(n):
...    return add([1 for p in range(n+1) if is_pseudoprime(p)])
sage: def count_primes3(n):
...    s=0; p=2
...    while p <= n: s+=1; p=next_prime(p)
...    return s
sage: def count_primes4(n):
...    s=0; p=2
...    while p <= n: s+=1; p=next_probable_prime(p)
...    return s
sage: def count_primes5(n):
...    s=0
...    for p in prime_range(n): s+=1
...    return s
sage: timeit('count_primes1(10^5)') # random, not tested
5 loops, best of 3: 674 ms per loop
sage: timeit('count_primes2(10^5)') # random, not tested
5 loops, best of 3: 256 ms per loop
sage: timeit('count_primes3(10^5)') # random
5 loops, best of 3: 49.2 ms per loop
sage: timeit('count_primes4(10^5)') # random
5 loops, best of 3: 48.6 ms per loop
sage: timeit('count_primes5(10^5)') # random
125 loops, best of 3: 2.67 ms per loop
sage: p=(2^42737+1)//3; a=3^42737
sage: timeit('a.gcd(p)') # random
125 loops, best of 3: 4.3 ms per loop
sage: timeit('a.jacobi(p)') # random
25 loops, best of 3: 26.1 ms per loop
sage: p=10^10+19; a=mod(17,p); a.log(2)
6954104378
sage: mod(2,p)^6954104378
17
sage: p=10^20+39; a=mod(17,p)
sage: time r=a.log(3) # not tested
CPU times: user 89.63 s, sys: 1.70 s, total: 91.33 s
"""
