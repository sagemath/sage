# sage_setup: distribution = sagemath-repl
# sage.doctest: needs sage.libs.pari
"""
This file contains all the example code from the published book
'Elementary Number Theory: Primes, Congruences, and Secrets' by
William Stein, Springer-Verlag, 2009.
"""

"""
sage: prime_range(10,50)
[11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
sage: [n for n in range(10,30) if not is_prime(n)]
[10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28]
sage: gcd(97,100)
1
sage: gcd(97 * 10^15, 19^20 * 97^2)
97
sage: factor(1275)
3 * 5^2 * 17
sage: factor(2007)
3^2 * 223
sage: factor(31415926535898)
2 * 3 * 53 * 73 * 2531 * 534697
sage: n = 7403756347956171282804679609742957314259318888\
....: 9231289084936232638972765034028266276891996419625117\
....: 8439958943305021275853701189680982867331732731089309\
....: 0055250511687706329907239638078671008609696253793465\
....: 0563796359
sage: len(n.str(2))
704
sage: len(n.str(10))
212
sage: n.is_prime()              # this is instant
False
sage: p = 2^32582657 - 1
sage: p.ndigits()                                                                       # needs sage.rings.real_interval_field
9808358
sage: s = p.str(10) # this takes a long time
sage: len(s)        # s is a very long string  (long time)
9808358
sage: s[:20]        # the first 20 digits of p (long time)
'12457502601536945540'
sage: s[-20:]       # the last 20 digits       (long time)
'11752880154053967871'

sage: # needs sage.symbolic
sage: prime_pi(6)
3
sage: prime_pi(100)
25
sage: prime_pi(3000000)
216816

sage: # needs sage.plot sage.symbolic
sage: plot(prime_pi, 1,1000, rgbcolor=(0,0,1))
Graphics object consisting of 1 graphics primitive
sage: P = plot(Li, 2,10000, rgbcolor='purple')
sage: Q = plot(prime_pi, 2,10000, rgbcolor='black')
sage: R = plot(sqrt(x)*log(x), 2,10000, rgbcolor='red')
sage: show(P + Q + R, xmin=0, figsize=[8,3])

sage: R = Integers(3)
sage: list(R)
[0, 1, 2]
sage: R = Integers(10)
sage: a = R(3)                 # create an element of Z/10Z
sage: a.multiplicative_order()
4
sage: [a^i for i in range(15)]
[1, 3, 9, 7, 1, 3, 9, 7, 1, 3, 9, 7, 1, 3, 9]
sage: euler_phi(2007)
1332
sage: n = 20
sage: k = euler_phi(n); k
8
sage: [Mod(x,n)^k for x in range(n) if gcd(x,n) == 1]
[1, 1, 1, 1, 1, 1, 1, 1]
sage: for n in range(1,10):
....:     print("{} {} {}".format(n, factorial(n-1) % n, -1 % n))
1 0 0
2 1 1
3 2 2
4 2 3
5 4 4
6 0 5
7 6 6
8 0 7
9 0 8
sage: CRT(2,3, 3, 5)
8
sage: CRT_list([2,3,2], [3,5,7])
23
sage: xgcd(5,7)
(1, 3, -2)
sage: xgcd(130,61)
(1, 23, -49)
sage: a = Mod(17, 61)
sage: a^(-1)
18
sage: 100.str(2)
'1100100'
sage: 0*2^0 + 0*2^1 + 1*2^2 + 0*2^3 + 0*2^4 + 1*2^5 + 1*2^6
100
sage: Mod(7,100)^91
43
sage: 7^91
80153343160247310515380886994816022539378033762994852007501964604841680190743
sage: n = 95468093486093450983409583409850934850938459083
sage: Mod(2,n)^(n-1)
34173444139265553870830266378598407069248687241
sage: factor(n)                # takes up to a few seconds.
1610302526747 * 59285812386415488446397191791023889
sage: n = 95468093486093450983409583409850934850938459083
sage: is_prime(n)
False
sage: for p in primes(100):
....:     if is_prime(2^p - 1):
....:         print("{} {}".format(p, 2^p - 1))
2 3
3 7
5 31
7 127
13 8191
17 131071
19 524287
31 2147483647
61 2305843009213693951
89 618970019642690137449562111
sage: def is_prime_lucas_lehmer(p):
....:     s = Mod(4, 2^p - 1)
....:     for i in range(3, p+1):
....:         s = s^2 - 2
....:     return s == 0
sage: # Check primality of 2^9941 - 1
sage: is_prime_lucas_lehmer(9941)
True
sage: # Check primality of 2^next_prime(1000)-1
sage: is_prime_lucas_lehmer(next_prime(1000))
False
sage: for p in primes(20):
....:     print("{} {}".format(p, primitive_root(p)))
2 1
3 2
5 2
7 3
11 2
13 2
17 3
19 2
sage: R.<x> = PolynomialRing(Integers(13))
sage: f = x^15 + 1
sage: f.roots()
[(12, 1), (10, 1), (4, 1)]
sage: f(12)
0
sage: R.<x> = PolynomialRing(Integers(13))
sage: f = x^6 + 1
sage: f.roots()
[(11, 1), (8, 1), (7, 1), (6, 1), (5, 1), (2, 1)]
sage: log(19683.0)
9.88751059801299
sage: log(3.0)
1.09861228866811
sage: log(19683.0) / log(3.0)
9.00000000000000
sage: plot(log, 0.1, 10, rgbcolor=(0,0,1))                                              # needs sage.plot
Graphics object consisting of 1 graphics primitive
sage: p = 53
sage: R = Integers(p)
sage: a = R.multiplicative_generator()
sage: v = sorted([(a^n, n) for n in range(p-1)])
sage: G = plot(point(v,pointsize=50,rgbcolor=(0,0,1)))                                  # needs sage.plot
sage: H = plot(line(v,rgbcolor=(0.5,0.5,0.5)))                                          # needs sage.plot
sage: G + H                                                                             # needs sage.plot
Graphics object consisting of 2 graphics primitives
sage: q = 93450983094850938450983409623
sage: q.is_prime()
True
sage: is_prime((q-1)//2)
True
sage: g = Mod(-2, q)
sage: g.multiplicative_order()
93450983094850938450983409622
sage: n = 18319922375531859171613379181
sage: m = 82335836243866695680141440300
sage: g^n
45416776270485369791375944998
sage: g^m
15048074151770884271824225393
sage: (g^n)^m
85771409470770521212346739540
sage: (g^m)^n
85771409470770521212346739540
sage: def rsa(bits):
....:     # only prove correctness up to 1024 bits
....:     proof = (bits <= 1024)
....:     while True:
....:         p = next_prime(ZZ.random_element(2**(bits//2 +1)),
....:                        proof=proof)
....:         q = next_prime(ZZ.random_element(2**(bits//2 +1)),
....:                        proof=proof)
....:         if (p != q and p > 2 and q > 2): break
....:     n = p * q
....:     phi_n = (p-1) * (q-1)
....:     while True:
....:         e = ZZ.random_element(1,phi_n)
....:         if gcd(e,phi_n) == 1: break
....:     d = lift(Mod(e,phi_n)^(-1))
....:     return e, d, n
sage: def encrypt(m,e,n):
....:     return lift(Mod(m,n)^e)
sage: def decrypt(c,d,n):
....:     return lift(Mod(c,n)^d)
sage: e,d,n = rsa(20)
sage: c = encrypt(123, e, n)
sage: decrypt(c, d, n)
123
sage: def encode(s):
....:     s = str(s)      # make input a string
....:     return sum(ord(s[i])*256^i for i in range(len(s)))
sage: def decode(n):
....:     n = Integer(n)  # make input an integer
....:     v = []
....:     while n != 0:
....:         v.append(chr(n % 256))
....:         n //= 256    # this replaces n by floor(n/256).
....:     return ''.join(v)
sage: m = encode('Run Nikita!'); m
40354769014714649421968722
sage: decode(m)
'Run Nikita!'
sage: def crack_rsa(n, phi_n):
....:     R.<x> = PolynomialRing(QQ)
....:     f = x^2 - (n+1 -phi_n)*x + n
....:     return [b for b, _ in f.roots()]
sage: crack_rsa(31615577110997599711, 31615577098574867424)
[8850588049, 3572144239]
sage: def crack_when_pq_close(n):
....:     t = Integer(ceil(sqrt(n)))
....:     while True:
....:        k = t^2 - n
....:        if k > 0:
....:           s = Integer(int(round(sqrt(t^2 - n))))
....:           if s^2 + n == t^2:
....:              return t+s, t-s
....:        t += 1
sage: crack_when_pq_close(23360947609)
(153649, 152041)
sage: p = next_prime(2^128); p
340282366920938463463374607431768211507
sage: q = next_prime(p)
sage: crack_when_pq_close(p*q)
(340282366920938463463374607431768211537,
      340282366920938463463374607431768211507)
sage: def crack_given_decrypt(n, m):
....:     n = Integer(n); m = Integer(m);  # some type checking
....:     # Step 1: divide out powers of 2
....:     while True:
....:         if is_odd(m): break
....:         divide_out = True
....:         for i in range(5):
....:            a = randrange(1,n)
....:            if gcd(a,n) == 1:
....:               if Mod(a,n)^(m//2) != 1:
....:                   divide_out = False
....:                   break
....:         if divide_out:
....:             m = m//2
....:         else:
....:             break
....:     # Step 2: Compute GCD
....:     while True:
....:         a = randrange(1,n)
....:         g = gcd(lift(Mod(a, n)^(m//2)) - 1, n)
....:         if g != 1 and g != n:
....:             return g
sage: n=32295194023343; e=29468811804857; d=11127763319273
sage: p = crack_given_decrypt(n, e*d - 1)
sage: p in (737531, n/737531)  # could be other prime divisor
True
sage: factor(n)
737531 * 43788253
sage: e = 22601762315966221465875845336488389513
sage: d = 31940292321834506197902778067109010093
sage: n = 268494924039590992469444675130990465673
sage: p = crack_given_decrypt(n, e*d - 1)  # not tested, known bug (see :issue:`32097`)
sage: p   # random output (could be other prime divisor)  # not tested
13432418150982799907
sage: n % p  # not tested
0
sage: set_random_seed(0)
sage: p = next_prime(randrange(2^96))
sage: q = next_prime(randrange(2^97))
sage: n = p * q
sage: qsieve(n)  # long time (8s on sage.math, 2011)
[(6340271405786663791648052309, 1), (46102313108592180286398757159, 1)]
sage: legendre_symbol(2,3)
-1
sage: legendre_symbol(1,3)
1
sage: legendre_symbol(3,5)
-1
sage: legendre_symbol(Mod(3,5), 5)
-1
sage: legendre_symbol(69,389)
1
sage: def kr(a, p):
....:  if Mod(a,p)^((p-1)//2) == 1:
....:     return 1
....:  else:
....:     return -1
sage: for a in range(1,5):
....:     print("{} {}".format(a, kr(a,5)))
1 1
2 -1
3 -1
4 1
sage: p = 726377359
sage: Mod(3, p)^((p-1)//2)
726377358
sage: def gauss(a, p):
....:     # make the list of numbers reduced modulo p
....:     v = [(n*a)%p for n in range(1, (p-1)//2 + 1)]
....:     # normalize them to be in the range -p/2 to p/2
....:     v = [(x if (x < p/2) else x - p) for x in v]
....:     # sort and print the resulting numbers
....:     v.sort()
....:     print(v)
....:     # count the number that are negative
....:     num_neg = len([x for x in v if x < 0])
....:     return (-1)^num_neg
sage: gauss(2, 13)
[-5, -3, -1, 2, 4, 6]
-1
sage: legendre_symbol(2,13)
-1
sage: gauss(4, 13)
[-6, -5, -2, -1, 3, 4]
1
sage: legendre_symbol(4,13)
1
sage: gauss(2,31)
[-15, -13, -11, -9, -7, -5, -3, -1, 2, 4, 6, 8, 10, 12, 14]
1
sage: legendre_symbol(2,31)
1
sage: K.<zeta> = CyclotomicField(5)
sage: zeta^5
1
sage: 1/zeta
-zeta^3 - zeta^2 - zeta - 1
sage: def gauss_sum(a,p):
....:     K.<zeta> = CyclotomicField(p)
....:     return sum(legendre_symbol(n,p) * zeta^(a*n) for n in range(1,p))
sage: g2 = gauss_sum(2,5); g2
2*zeta^3 + 2*zeta^2 + 1
sage: g2.complex_embedding()
-2.236067977... + ...e-16*I
sage: g2^2
5
sage: [gauss_sum(a, 7)^2 for a in range(1,7)]
[-7, -7, -7, -7, -7, -7]
sage: [gauss_sum(a, 13)^2 for a in range(1,13)]
[13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13]
sage: S.<x> = PolynomialRing(GF(13))
sage: R.<alpha> = S.quotient(x^2 - 3)
sage: (2+3*alpha)*(1+2*alpha)
7*alpha + 7
sage: def find_sqrt(a, p):
....:     assert (p-1)%4 == 0
....:     assert legendre_symbol(a,p) == 1
....:     S.<x> = PolynomialRing(GF(p))
....:     R.<alpha> = S.quotient(x^2 - a)
....:     while True:
....:         z = GF(p).random_element()
....:         w = (1 + z*alpha)^((p-1)//2)
....:         (u, v) = (w[0], w[1])
....:         if v != 0: break
....:     if (-u/v)^2 == a: return -u/v
....:     if ((1-u)/v)^2 == a: return (1-u)/v
....:     if ((-1-u)/v)^2 == a: return (-1-u)/v
sage: b = find_sqrt(3,13)
sage: b                        # random: either 9 or 3
9
sage: b^2
3
sage: b = find_sqrt(3,13)
sage: b                        # see, it's random
4
sage: find_sqrt(5,389)         # random: either 303 or 86
303
sage: find_sqrt(5,389)         # see, it's random
86

# Several of the examples below had to be changed due to improved
# behavior of the continued_fraction function #8017 and #14567.

sage: continued_fraction(17/23)
[0; 1, 2, 1, 5]
sage: reset('e')
sage: continued_fraction(e)
[2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, ...]
sage: continued_fraction_list(e, bits=21)
[2, 1, 2, 1, 1, 4, 1, 1, 6]
sage: continued_fraction_list(e, bits=30)
[2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8]
sage: a = continued_fraction(17/23); a
[0; 1, 2, 1, 5]
sage: a.value()
17/23
sage: b = continued_fraction(6/23); b
[0; 3, 1, 5]
sage: c = continued_fraction(pi); c
[3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
sage: [c.convergent(i) for i in range(5)]
[3, 22/7, 333/106, 355/113, 103993/33102]
sage: [c.p(n)*c.q(n-1) - c.q(n)*c.p(n-1) for n in range(-1, 13)]
[1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1]
sage: [c.p(n)*c.q(n-2) - c.q(n)*c.p(n-2) for n in range(13)]
[3, -7, 15, -1, 292, -1, 1, -1, 2, -1, 3, -1, 14]
sage: c = continued_fraction([1,2,3,4,5])
sage: c.convergents()
[1, 3/2, 10/7, 43/30, 225/157]
sage: [c.p(n) for n in range(len(c))]
[1, 3, 10, 43, 225]
sage: [c.q(n) for n in range(len(c))]
[1, 2, 7, 30, 157]
sage: c = continued_fraction([1,1,1,1,1,1,1,1])
sage: v = [(i, c.p(i)/c.q(i)) for i in range(len(c))]

sage: # needs sage.plot
sage: P = point(v, rgbcolor=(0,0,1), pointsize=40)
sage: L = line(v, rgbcolor=(0.5,0.5,0.5))
sage: L2 = line([(0,c.value()), (len(c)-1,c.value())],
....:           thickness=0.5, rgbcolor=(0.7,0,0))
sage: (L + L2 + P).show(xmin=0, ymin=1)

sage: def cf(bits):
....:     x = (1 + sqrt(RealField(bits)(5))) / 2
....:     return continued_fraction(x)
sage: cf(10)
[1; 1, 1, 1, 1, 1, 1, 2]
sage: cf(30)
[1; 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
sage: cf(50)
 [1; 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
sage: def cf_sqrt_d(d, bits):
....:     x = sqrt(RealField(bits)(d))
....:     return continued_fraction(x)
sage: cf_sqrt_d(389,50)
[19; 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38]
sage: cf_sqrt_d(389,100)
[19; 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 1, 1, 1, 1, 2, 1, 38, 1, 2, 2]
sage: def newton_root(f, iterates=2, x0=0, prec=53):
....:     x = RealField(prec)(x0)
....:     R = PolynomialRing(ZZ,'x')
....:     f = R(f)
....:     g = f.derivative()
....:     for i in range(iterates):
....:         x = x - f(x)/g(x)
....:     return x
sage: reset('x')
sage: a = newton_root(3847*x^2 - 14808904*x + 36527265); a
2.46815700480740
sage: cf = continued_fraction(a); cf
[2; 2, 7, 2, 1, 5, 1, 1, 1, 1, 1, 1, 103, 8, 1, 2, 3, 2]
sage: c = cf[:12]; c
[2; 2, 7, 2, 1, 5, 1, 1, 1, 1, 2]
sage: c.value()
9495/3847
sage: def sum_of_two_squares_naive(n):
....:     for i in range(int(sqrt(n))):
....:         if is_square(n - i^2):
....:             return i, (Integer(n-i^2)).sqrt()
....:     return "%s is not a sum of two squares"%n
sage: sum_of_two_squares_naive(23)
'23 is not a sum of two squares'
sage: sum_of_two_squares_naive(389)
(10, 17)
sage: sum_of_two_squares_naive(2007)
'2007 is not a sum of two squares'
sage: sum_of_two_squares_naive(2008)
'2008 is not a sum of two squares'
sage: sum_of_two_squares_naive(2009)
(28, 35)
sage: 28^2 + 35^2
2009
sage: sum_of_two_squares_naive(2*3^4*5*7^2*13)
(189, 693)
sage: def sum_of_two_squares(p):
....:     p = Integer(p)
....:     assert p%4 == 1, "p must be 1 modulo 4"
....:     r = Mod(-1,p).sqrt().lift()
....:     v = continued_fraction(-r/p)
....:     n = floor(sqrt(p))
....:     for x in v.convergents():
....:         c = r*x.denominator() + p*x.numerator()
....:         if -n <= c and c <= n:
....:             return (abs(x.denominator()),abs(c))
sage: p = next_prime(next_prime(10^10))
sage: sum_of_two_squares(p)
(55913, 82908)
sage: sum_of_two_squares_naive(p)
(55913, 82908)
sage: E = EllipticCurve([-5, 4])
sage: E
Elliptic Curve defined by y^2  = x^3 - 5*x + 4
over Rational Field
sage: P = E.plot(thickness=4,rgbcolor=(0.1,0.7,0.1))                                    # needs sage.plot
sage: P.show(figsize=[4,6])                                                             # needs sage.plot
sage: E = EllipticCurve(GF(37), [1,0])
sage: E
Elliptic Curve defined by y^2  = x^3 + x over
Finite Field of size 37
sage: E.plot(pointsize=45)                                                              # needs sage.plot
Graphics object consisting of 1 graphics primitive
sage: E = EllipticCurve([-5,4])
sage: P = E([1,0]); Q = E([0,2])
sage: P + Q
(3 : 4 : 1)
sage: P + P
(0 : 1 : 0)
sage: P + Q + Q + Q + Q
(350497/351649 : 16920528/208527857 : 1)
sage: R.<x1,y1,x2,y2,x3,y3,a,b> = QQ[]
sage: rels = [y1^2 - (x1^3 + a*x1 + b),
....:         y2^2 - (x2^3 + a*x2 + b),
....:         y3^2 - (x3^3 + a*x3 + b)]
...
sage: Q = R.quotient(rels)
sage: def op(P1,P2):
....:     x1,y1 = P1;  x2,y2 = P2
....:     lam = (y1 - y2)/(x1 - x2); nu  = y1 - lam*x1
....:     x3 = lam^2 - x1 - x2; y3 = -lam*x3 - nu
....:     return (x3, y3)
sage: P1 = (x1,y1); P2 = (x2,y2); P3 = (x3,y3)
sage: Z = op(P1, op(P2,P3)); W = op(op(P1,P2),P3)
sage: (Q(Z[0].numerator()*W[0].denominator() -
....:        Z[0].denominator()*W[0].numerator())) == 0
True
sage: (Q(Z[1].numerator()*W[1].denominator() -
....:        Z[1].denominator()*W[1].numerator())) == 0
True
sage: def lcm_upto(B):
....:     return prod([p^int(math.log(B)/math.log(p))
....:                  for p in prime_range(B+1)])
sage: lcm_upto(10^2)
69720375229712477164533808935312303556800
sage: LCM([1..10^2])
69720375229712477164533808935312303556800
sage: def pollard(N, B=10^5, stop=10):
....:     m = prod([p^int(math.log(B)/math.log(p))
....:               for p in prime_range(B+1)])
....:     for a in [2..stop]:
....:         x = (Mod(a,N)^m - 1).lift()
....:         if x == 0: continue
....:         g = gcd(x, N)
....:         if g != 1 or g != N: return g
....:     return 1
sage: pollard(5917,5)
61
sage: pollard(779167,5)
1
sage: pollard(779167,15)
2003
sage: pollard(4331,7)
1
sage: pollard(4331,5)
61
sage: pollard(187, 15, 2)
1
sage: pollard(187, 15)
11
sage: def ecm(N, B=10^3, trials=10):
....:     m = prod([p^int(math.log(B)/math.log(p))
....:               for p in prime_range(B+1)])
....:     R = Integers(N)
....:     # Make Sage think that R is a field:
....:     R.is_field = lambda : True
....:     for _ in range(trials):
....:         while True:
....:             a = R.random_element()
....:             if gcd(4*a.lift()^3 + 27, N) == 1: break
....:         try:
....:             m * EllipticCurve([a, 1])([0,1])
....:         except ZeroDivisionError as msg:
....:             # msg: "Inverse of <int> does not exist"
....:             return gcd(Integer(str(msg).split()[2]), N)
....:     return 1
sage: set_random_seed(2)
sage: ecm(5959, B=20)
101
sage: ecm(next_prime(10^20)*next_prime(10^7), B=10^3)
10000019
sage: p = 785963102379428822376694789446897396207498568951
sage: E = EllipticCurve(GF(p), \
....:  [317689081251325503476317476413827693272746955927,
....:   79052896607878758718120572025718535432100651934])
sage: E.cardinality()
785963102379428822376693024881714957612686157429
sage: E.cardinality().is_prime()
True
sage: B = E([
....:  771507216262649826170648268565579889907769254176,
....:  390157510246556628525279459266514995562533196655])
sage: n=670805031139910513517527207693060456300217054473
sage: r=70674630913457179596452846564371866229568459543
sage: P = E([14489646124220757767,
....:  669337780373284096274895136618194604469696830074])
sage: encrypt = (r*B, P + r*(n*B))
sage: encrypt[1] - n*encrypt[0] == P   # decrypting works
True
sage: T = lambda v: EllipticCurve(v
....:        ).torsion_subgroup().invariants()
sage: T([-5,4])
(2,)
sage: T([-43,166])
(7,)
sage: T([-4,0])
(2, 2)
sage: T([-1386747, 368636886])
(2, 8)
sage: r = lambda v: EllipticCurve(v).rank()
sage: r([-5,4])
1
sage: r([0,1])
0
sage: r([-3024, 46224])
2
sage: r([-112, 400])
3
sage: r([-102627, 12560670])
4
sage: def cong(n):
....:     G = EllipticCurve([-n^2,0]).gens()
....:     if len(G) == 0: return False
....:     x,y,_ = G[0]
....:     return ((n^2-x^2)/y,-2*x*n/y,(n^2+x^2)/y)
sage: cong(6)
(3, 4, 5)
sage: cong(5)
(3/2, 20/3, 41/6)
sage: cong(1)
False
sage: cong(13)
(323/30, 780/323, 106921/9690)
sage: (323/30 * 780/323)/2
13
sage: (323/30)^2 + (780/323)^2 == (106921/9690)^2
True
"""
