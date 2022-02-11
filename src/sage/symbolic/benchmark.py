r"""
Benchmarks

Tests that will take a long time if something is wrong, but be very
quick otherwise.  See https://wiki.sagemath.org/symbench.  The
parameters chosen below are such that with pynac most of these take
well less than a second, but would not even be feasible using Sage's
Maxima-based symbolics.

Problem R1

Important note. Below we do s.expand().real() because s.real() takes forever (TODO?). ::

    sage: f(z) = sqrt(1/3)*z^2 + i/3
    sage: s = f(f(f(f(f(f(f(f(f(f(i/2))))))))))
    sage: s.expand().real()
    -15323490199844318074242473679071410934833494247466385771803570370858961112774390851798166656796902695599442662754502211584226105508648298600018090510170430216881977761279503642801008178271982531042720727178135881702924595044672634313417239327304576652633321095875724771887486594852083526001648217317718794685379391946143663292907934545842931411982264788766619812559999515408813796287448784343854980686798782575952258163992236113752353237705088451481168691158059505161807961082162315225057299394348203539002582692884735745377391416638540520323363224931163680324690025802009761307137504963304640835891588925883135078996398616361571065941964628043214890356454145039464055430143/160959987592246947739944859375773744043416001841910423046466880402863187009126824419781711398533250016237703449459397319370100476216445123130147322940019839927628599479294678599689928643570237983736966305423831947366332466878486992676823215303312139985015592974537721140932243906832125049776934072927576666849331956351862828567668505777388133331284248870175178634054430823171923639987569211668426477739974572402853248951261366399284257908177157179099041115431335587887276292978004143353025122721401971549897673882099546646236790739903146970578001092018346524464799146331225822142880459202800229013082033028722077703362360159827236163041299500992177627657014103138377287073792*sqrt(1/3)


Problem R2::

    sage: def hermite(n,y):
    ....:     if n == 1: return 2*y
    ....:     if n == 0: return 1
    ....:     return expand(2*y*hermite(n-1,y) - 2*(n-1)*hermite(n-2,y))
    sage: hermite(15,var('y'))
    32768*y^15 - 1720320*y^13 + 33546240*y^11 - 307507200*y^9 + 1383782400*y^7 - 2905943040*y^5 + 2421619200*y^3 - 518918400*y

Problem R3::

    sage: f = sum(var('x,y,z')); a = [bool(f==f) for _ in range(100000)]

Problem R4::

    sage: u=[e,pi,sqrt(2)]; Tuples(u,3).cardinality()
    27

Problem R5::

    sage: def blowup(L,n):
    ....:    for i in [0..n]:
    ....:        L.append( (L[i] + L[i+1]) * L[i+2] )
    sage: L = list(var('x,y,z'))
    sage: blowup(L,15)
    sage: len(set(L))
    19

Problem R6::

    sage: sum(((x+sin(i))/x+(x-sin(i))/x) for i in range(100)).expand()
    200

Problem R7::

    sage: f = x^24+34*x^12+45*x^3+9*x^18 +34*x^10+ 32*x^21
    sage: a = [f(x=random()) for _ in range(10^4)]

Problem R10::

    sage: v = [float(z) for z in [-pi,-pi+1/100..,pi]]

Problem R11::

    sage: a = [random() + random()*I for w in [0..100]]
    sage: a.sort()


Problem W3::

    sage: acos(cos(x))
    arccos(cos(x))

PROBLEM S1::

    sage: _=var('x,y,z')
    sage: f = (x+y+z+1)^10
    sage: g = expand(f*(f+1))


PROBLEM S2::

    sage: _=var('x,y')
    sage: a = expand((x^sin(x) + y^cos(y) - z^(x+y))^100)

PROBLEM S3::

    sage: _=var('x,y,z')
    sage: f = expand((x^y + y^z + z^x)^50)
    sage: g = f.diff(x)

PROBLEM S4::

    w = (sin(x)*cos(x)).series(x,400)



"""
