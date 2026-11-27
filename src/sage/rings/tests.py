"""
Tests for rings

TESTS::

    sage: K.<x> = FractionField(QQ['x'])
    sage: V.<z> = K[]
    sage: x+z
    z + x
"""

import sage.misc.prandom as random

from sage.misc.random_testing import random_testing


def prime_finite_field():
    """
    Create a random prime finite field with cardinality at most 10^20.

    OUTPUT: a prime finite field

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.prime_finite_field(); K                              # needs sage.rings.finite_rings
        Finite Field of size ...
        sage: K.cardinality().is_prime()                                                # needs sage.rings.finite_rings
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    return GF(ZZ.random_element(x=2, y=10**20 - 12).next_prime())


def finite_field():
    """
    Create a random finite field with degree at most 20 and prime at most 10^6.

    OUTPUT: a finite field

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.finite_field(); K                                    # needs sage.rings.finite_rings
        Finite Field...of size ...
        sage: K.cardinality().is_prime_power()                                          # needs sage.rings.finite_rings
        True
        sage: while K.cardinality().is_prime():                                         # needs sage.rings.finite_rings
        ....:     K = sage.rings.tests.finite_field()
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    p = ZZ.random_element(x=2, y=10**6 - 18).next_prime()
    d = ZZ.random_element(x=1, y=20)
    return GF(p**d, 'a')


def small_finite_field():
    """
    Create a random finite field with cardinality at most 2^16.

    OUTPUT: a finite field

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.small_finite_field(); K
        Finite Field...of size ...
        sage: q = K.cardinality()
        sage: q.is_prime_power()
        True
        sage: q <= 2^16
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.finite_field_constructor import GF
    while True:
        q = ZZ.random_element(x=2, y=2**16)
        if q.is_prime_power():
            return GF(q, 'a')


def integer_mod_ring():
    """
    Return a random ring of integers modulo n with n at most 50000.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: R = sage.rings.tests.integer_mod_ring(); R
        Ring of integers modulo ...
        sage: R.cardinality() <= 50000
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    n = ZZ.random_element(x=2, y=50000)
    return IntegerModRing(n)


def padic_field():
    """
    Return a random `p`-adic field modulo n with p at most 10000
    and precision between 10 and 100.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: sage.rings.tests.padic_field()                                            # needs sage.rings.padics
        ...-adic Field with capped relative precision ...
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.padics.factory import Qp
    prec = ZZ.random_element(x=10, y=100)
    p = ZZ.random_element(x=2, y=10**4 - 30).next_prime()
    return Qp(p, prec)


def quadratic_number_field():
    """
    Return a quadratic extension of QQ.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.quadratic_number_field(); K                          # needs sage.rings.number_field
        Number Field in a with defining polynomial x^2 ... with a = ...
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.number_field.number_field import QuadraticField
    while True:
        d = ZZ.random_element(x=-10**5, y=10**5)
        if not d.is_square():
            return QuadraticField(d, 'a')


def absolute_number_field(maxdeg=10):
    """
    Return an absolute extension of QQ of degree at most 10.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.absolute_number_field(); K                           # needs sage.rings.number_field
        Number Field in a with defining polynomial ...
        sage: K.degree() <= 10                                                          # needs sage.rings.number_field
        True
    """
    from sage.rings.integer_ring import ZZ
    from sage.rings.number_field.number_field import NumberField
    R = ZZ['x']
    while True:
        f = R.random_element(degree=ZZ.random_element(x=1, y=maxdeg),
                             x=-100, y=100)
        if f.degree() <= 0:
            continue
        f = f + R.gen()**(f.degree() + 1)  # make monic
        if f.is_irreducible():
            return NumberField(f, 'a')


def relative_number_field(n=2, maxdeg=2):
    """
    Return a tower of at most n extensions each of degree at most maxdeg.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: K = sage.rings.tests.relative_number_field(3); K
        Number Field in aaa with defining polynomial x^2 ... over its base field
        sage: K.relative_degree()
        2
        sage: L = K.base_ring()
        sage: L.relative_degree()
        2
        sage: M = L.base_ring()
        sage: M.relative_degree()
        2
        sage: M.base_ring() is QQ
        True

    TESTS:

    Check that :issue:`32117` is fixed::

        sage: set_random_seed(3030)
        sage: from sage.rings.tests import relative_number_field
        sage: _ = relative_number_field(3)                                              # needs sage.rings.number_field
    """
    from sage.rings.integer_ring import ZZ
    K = absolute_number_field(maxdeg)
    n -= 1
    var = 'aa'
    R = ZZ['x']
    R1 = K['x']
    while n >= 1:
        while True:
            f = R.random_element(degree=ZZ.random_element(x=1, y=maxdeg),
                                 x=-100, y=100)
            if f.degree() <= 0:
                continue
            f = f * f.denominator()  # bug trac #4781
            f = f + R.gen()**maxdeg  # make monic
            if R1(f).is_irreducible():
                break
        K = K.extension(f, var)
        R1 = K['x']
        var += 'a'
        n -= 1
    return K


def rings0():
    """
    Return a list of pairs (f, desc), where f is a function that when
    called creates a random ring of a certain representative type
    described by desc.

    RINGS:

    - ZZ
    - QQ
    - ZZ/nZZ
    - GF(p)
    - GF(q)
    - p-adic fields
    - quadratic number fields
    - absolute number fields
    - relative number fields

    EXAMPLES::

        sage: import sage.rings.tests
        sage: type(sage.rings.tests.rings0())
        <... 'list'>
    """
    from sage.rings.integer_ring import IntegerRing
    from sage.rings.rational_field import RationalField

    v = [(IntegerRing, 'ring of integers'),
         (RationalField, 'field of rational numbers'),
         (integer_mod_ring, 'integers modulo n for n at most 50000')]
    try:
        v += [(prime_finite_field, 'a prime finite field with cardinality at most 10^20'),
              (finite_field, 'finite field with degree at most 20 and prime at most 10^6'),
              (small_finite_field, 'finite field with cardinality at most 2^16')]
    except ImportError:
        pass

    try:
        v += [(padic_field, 'a p-adic field')]
    except ImportError:
        pass

    try:
        v += [(quadratic_number_field, 'a quadratic number field'),
              (absolute_number_field, 'an absolute number field of degree at most 10'),
              (relative_number_field, 'a tower of at most 2 extensions each of degree at most 2')]
    except ImportError:
        pass

    return v


def rings1():
    """
    Return an iterator over random rings.

    Return a list of pairs (f, desc), where f is a function that
    outputs a random ring that takes a ring and possibly
    some other data as constructor.

    RINGS:

    - polynomial ring in one variable over a rings0() ring.
    - polynomial ring over a rings1() ring.
    - multivariate polynomials
    - power series rings in one variable over a rings0() ring.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: type(sage.rings.tests.rings0())
        <... 'list'>
    """
    v = rings0()
    X = random_rings(level=0)
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.rings.power_series_ring import PowerSeriesRing
    from sage.rings.integer_ring import ZZ

    v = [(lambda: PolynomialRing(next(X), names='x'),
          'univariate polynomial ring over level 0 ring'),
         (lambda: PowerSeriesRing(next(X), names='x'),
          'univariate power series ring over level 0 ring')]

    try:
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
    except ImportError:
        pass
    else:
        v += [(lambda: LaurentPolynomialRing(next(X), names='x'),
               'univariate Laurent polynomial ring over level 0 ring')]

    v += [(lambda: PolynomialRing(next(X), abs(ZZ.random_element(x=2, y=10)),
                                  names='x'),
           'multivariate polynomial ring in between 2 and 10 variables over a level 0 ring')]

    return v


MAX_LEVEL = 99999


def random_rings(level=MAX_LEVEL):
    """
    Return an iterator over random rings up to the given "level" of complexity.

    EXAMPLES::

        sage: import sage.rings.tests
        sage: type(sage.rings.tests.random_rings())
        <... 'generator'>
    """
    v = rings0()
    if level >= 1:
        v += rings1()
    while True:
        yield random.choice(v)[0]()


@random_testing
def check_random_elements(level=MAX_LEVEL, trials=1):
    """
    Create random elements of random rings until a crash occurs, in
    which case an exception is raised.  Defaults to running a single
    trial, but more can be specified.  To run tests in an infinite
    loop, you could use::

        while True: check_random_elements(trials=100, print_seed=True)

    INPUT:

    - ``level`` -- (default: ``MAX_LEVEL``) controls the types of rings to use
    - ``trials`` -- a positive integer (default: 1); the number of trials to run
    - ``seed`` -- the random seed to use; if not specified, uses a truly random seed
    - ``print_seed`` -- if ``True`` (default: ``False``), prints the random seed chosen

    EXAMPLES::

        sage: import sage.rings.tests
        sage: sage.rings.tests.check_random_elements(trials=2, seed=0)                   # needs sage.rings.number_field
        survived 0 tests
        Rational Field
        -1/2
        ----
        survived 1 tests
        Multivariate Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7 over Number Field in a with defining polynomial x^5 + 39*x^4 - 46*x^3 + 82*x^2 + 17*x - 15
        (-4*a^4 - 1/5*a^3 - 1/14*a^2 - a - 1/4)*x1*x2 + (-a^4 + 3/2*a^3 + a^2 + 1/9*a + 5/2)*x1*x3 + (a^4 + 1/3*a^3 - a^2 + a - 1)*x0*x5 + (a^4 - a^3 - a)*x4*x6 + (-3*a^4 - 1/4*a^3 - 1/2*a^2 + 4*a)*x1*x7
        ----

        sage: sage.rings.tests.check_random_elements(trials=10)
        survived 0 tests...
        sage: sage.rings.tests.check_random_elements(trials=1000)  # long time (5 seconds)
        survived 0 tests...
    """
    r = random_rings(level)
    i = 0
    for R in r:
        print("survived %s tests" % i)
        i += 1
        print(R)
        print(R.random_element())
        print("----")
        if i >= trials:
            return


@random_testing
def check_random_arith(level=MAX_LEVEL, trials=1):
    """
    Create random elements of random rings and do some arithmetic with them.

    Repeats until a crash occurs, in which case an exception is
    raised.  Defaults to running a single trial, but more can be
    specified.  To run tests in an infinite loop, you could use::

        while True: check_random_arith(trials=100, print_seed=True)

    INPUT:

    - ``level`` -- (default: ``MAX_LEVEL``) controls the types of rings to use
    - ``trials`` -- positive integer (default: 1); the number of trials to run
    - ``seed`` -- the random seed to use; if not specified, uses a truly random seed
    - ``print_seed`` -- if ``True`` (default: ``False``), prints the random seed chosen

    EXAMPLES::

        sage: import sage.rings.tests
        sage: sage.rings.tests.check_random_arith(trials=2, seed=0)
        survived 0 tests
        Rational Field
        -1/2 -1/95
        49/95
        survived 1 tests
        Multivariate Polynomial Ring in x0, x1, x2, x3, x4, x5 over Number Field in a with defining polynomial x^8 + 82*x^7 + 17*x^6 - 15*x^5 - 79*x^4 + 87*x^3 + 95*x^2 + 5*x - 1
        (1/4*a^7 - a^6 - 1/3*a^5 + 1/2*a^4 - a^3 + 1/2*a^2 - a - 2/27)*x3^2 + (2*a^7 + 2*a^5 + a^4 - 1/357*a^3 - 1/14*a^2 + 2)*x0*x4 + (22*a^7 - 3/1955*a^6 + 1/2*a^5 - 2/7*a^4 - a^3 - 1/2*a^2 - a + 3/2)*x1*x4 + (-4*a^7 - 1/5*a^6 - 1/14*a^5 - a^4 - 1/4*a^3 + a^2 - a)*x4*x5 + (a^7 + 1/9*a^6 + 5/2*a^5 - 3*a^4 - 1/4*a^3 - 1/2*a^2 + 4*a)*x1 (-2*a^7 - 1/10*a^6 + 2/5*a^5 - 1/27*a^4 - 2*a^3 + 2*a^2 - a + 1)*x0*x3 + (-1/2*a^5 - 1/3*a^2 + a + 1/6)*x3^2 + (a^7 - 1/2*a^6 + 94*a^5 - 1/2*a^4 - a^2 - 1/3*a + 1)*x0*x5 + (a^7 - 3*a^6 + a^5 + 1/4*a^4 - a^3 - 4/5*a^2 + 2/3)*x1*x5 + (-13*a^7 + 2*a^5 - 4/3*a^4 - a^3 - 1/2*a + 1/2)
        (8335830563728609/648*a^7 + 4375169223294811/1620*a^6 - 3572531535458341/1620*a^5 - 367129730038025827/29160*a^4 + 43742564579945269/3240*a^3 + 4031628460406779/270*a^2 + 851388560942053/1080*a - 101917410191525/648)*x0*x3^3 + (8658913703/18*a^7 + 807954417/8*a^6 - 17812785073/216*a^5 - 11299528315/24*a^4 + 1514598493/3*a^3 + 45229188080/81*a^2 + 6367586285/216*a - 1905616129/324)*x3^4 + (1892073453153160475/19278*a^7 + 1324105854787133021/64260*a^6 - 1081194733520142317/64260*a^5 - 6172696163274863233/64260*a^4 + 2206382690878441393/21420*a^3 + 2440272627775145449/21420*a^2 + 386497705816640749/64260*a - 77110982828047163/64260)*x0^2*x3*x4 + (159541839636254386607/147798*a^7 + 1674752242788361695463/7389900*a^6 - 683757004135171199879/3694950*a^5 - 3903667031650238859241/3694950*a^4 + 4186007123190116706733/3694950*a^3 + 264557099846797814957/211140*a^2 + 162949685960759150951/2463300*a - 696652283545469812/52785)*x0*x1*x3*x4 + (1747028552159/476*a^7 + 16181482733/21*a^6 - 673860035066/1071*a^5 - 1282387952246/357*a^4 + 4125417243814/1071*a^3 + 6083648636327/1428*a^2 + 80295545534/357*a - 32039871065/714)*x0*x3^2*x4 + (331450840052539/8211*a^7 + 23195505266229/2737*a^6 - 284103319866802/41055*a^5 - 3243973419377599/82110*a^4 + 1739299959336677/41055*a^3 + 7694705191339393/164220*a^2 + 406236837637037/164220*a - 40524589539473/82110)*x1*x3^2*x4 + (-1418230900929283/216*a^7 - 297750782807311/216*a^6 + 243127524163409/216*a^5 + 694024993587073/108*a^4 - 13781882631763/2*a^3 - 823113668644273/108*a^2 - 260734428555329/648*a + 8669947126141/108)*x0*x3^2*x5 + (-3603487699045651/540*a^7 - 1008713303801291/720*a^6 + 494197125293329/432*a^5 + 14107212163813613/2160*a^4 - 15127542930084671/2160*a^3 - 1045697135021857/135*a^2 - 49072785777299/120*a + 528693280606921/6480)*x1*x3^2*x5 + (-35767908689480873/714*a^7 - 7509301062897307/714*a^6 + 721376194378387/84*a^5 + 210040460249221963/4284*a^4 - 75077344424128273/1428*a^3 - 83035997935705561/1428*a^2 - 4383830224147987/1428*a + 291542443891301/476)*x0^2*x4*x5 + (-4942036926751918354/8211*a^7 - 2766818600332080223/21896*a^6 + 2824046967567076043/27370*a^5 + 32245780231960140987/54740*a^4 - 15560107737622455316/24633*a^3 - 114730489727287249373/164220*a^2 - 1211423963103680813/32844*a + 402822970704135781/54740)*x0*x1*x4*x5 + (-919575723486092563307/1642200*a^7 - 16088376449609089341/136850*a^6 + 6305719860814926017/65688*a^5 + 450003264731284780307/821100*a^4 - 137871600241618847881/234600*a^3 - 213481555602902273941/328440*a^2 - 11270616588162289157/328440*a + 468463497709203519/68425)*x1^2*x4*x5 + (-296570993234981827/1512*a^7 - 518863826148202811/12600*a^6 + 1271030184371047943/37800*a^5 + 7256494042407260587/37800*a^4 - 7781333680517159641/37800*a^3 - 63749645992630777/280*a^2 - 454358715458875733/37800*a + 755417315081218/315)*x0*x3*x4*x5 + (-14669790437/2*a^7 - 258707585821/168*a^6 + 1056234575803/840*a^5 + 6030195094033/840*a^4 - 1077723389123/140*a^3 - 7151811607883/840*a^2 - 125858377409/280*a + 75330821071/840)*x3^2*x4*x5 + (1501714795107104/15*a^7 + 5885186092945767/280*a^6 - 14416594080171413/840*a^5 - 82306408098976987/840*a^4 + 88259374529418377/840*a^3 + 97615403121218651/840*a^2 + 5153540213684083/840*a - 514096907873707/420)*x0*x4*x5^2 + (854695816128330157/8400*a^7 + 12817092263543953/600*a^6 - 14652062481291307/840*a^5 - 418253655231120007/4200*a^4 + 448504640867919821/4200*a^3 + 99209770187625719/840*a^2 + 748244823055711/120*a - 870822883368243/700)*x1*x4*x5^2 + (158818466671888433/3240*a^7 + 33343176167568581/3240*a^6 - 81678816706902197/9720*a^5 - 155438702168107007/3240*a^4 + 166681099904177029/3240*a^3 + 553050919967402617/9720*a^2 + 29197955088382571/9720*a - 2912673192557701/4860)*x0*x1*x3 + (197968812379/108*a^7 + 83125208285/216*a^6 - 22625214433/72*a^5 - 387511804901/216*a^4 + 69256554457/36*a^3 + 459589015369/216*a^2 + 24263696167/216*a - 4840901819/216)*x1*x3^2 + (-675521339326351/27*a^7 - 378193249936939/72*a^6 + 308812583669081/72*a^5 + 1763055434274023/72*a^4 - 1890571748692697/72*a^3 - 232331534640807/8*a^2 - 110392098127027/72*a + 305896712356)*x0*x1*x5 + (-18308126159975299/720*a^7 - 2882777528847821/540*a^6 + 34872940673869/8*a^5 + 8959258419605051/360*a^4 - 640483580504903/24*a^3 - 708378721654847/24*a^2 - 186991914791501/120*a + 4663394540741/15)*x1^2*x5 + (2*a^7 + 1/10*a^6 - 2/5*a^5 + 1/27*a^4 + 2*a^3 - 2*a^2 + a - 1)*x0*x3 + (18072135893255449/216*a^7 + 1264719448645567/72*a^6 - 3098108815811495/216*a^5 - 53062645815387341/648*a^4 + 18966833967284401/216*a^3 + 2330825923921351/24*a^2 + 1107489624724565/216*a - 73652540191333/72)*x3^2 + (1367342594277744215/2142*a^7 + 31896336292953897/238*a^6 - 11162082599960971/102*a^5 - 63725934193782685/102*a^4 + 318896819516971189/476*a^3 + 352701815042788259/476*a^2 + 27930986234873423/714*a - 2786285361486205/357)*x0*x4 + (288239804512522153973/41055*a^7 + 242058263923052911661/164220*a^6 - 11626586002859102347/9660*a^5 - 282105866310724657066/41055*a^4 + 1210038823799157812587/164220*a^3 + 11246304810899299411/1380*a^2 + 3364533322153197667/7820*a - 14096569821296035687/164220)*x1*x4 + (-a^7 + 1/2*a^6 - 94*a^5 + 1/2*a^4 + a^2 + 1/3*a - 1)*x0*x5 + (-a^7 + 3*a^6 - a^5 - 1/4*a^4 + a^3 + 4/5*a^2 - 2/3)*x1*x5 + (-535806565281138601/420*a^7 - 2678333840353439/10*a^6 + 91853395374613427/420*a^5 + 1048808477574711929/840*a^4 - 1124665531766446127/840*a^3 - 103657241175713587/70*a^2 - 32835090185999321/420*a + 2183666325731897/140)*x4*x5 + (1275257363094063/4*a^7 + 1204803780912301/18*a^6 - 5902673860120745/108*a^5 - 22466142413184839/72*a^4 + 24091048598585143/72*a^3 + 39967268628811033/108*a^2 + 2110045336326155/108*a - 420979652006525/108)*x1 + (13*a^7 - 2*a^5 + 4/3*a^4 + a^3 + 1/2*a + 1/2)
        sage: sage.rings.tests.check_random_arith(trials=10)
        survived 0 tests...
        sage: sage.rings.tests.check_random_arith(trials=1000)   # long time (5 seconds?)
        survived 0 tests...
    """
    i = 0
    for x in random_rings(level):
        print("survived %s tests" % i)
        i += 1
        print(x)
        a = x.random_element()
        b = x.random_element()
        print(a, b)
        print(a * b + a - b + 1)
        if i >= trials:
            return


@random_testing
def check_karatsuba_multiplication(base_ring, maxdeg1, maxdeg2,
                                  ref_mul=lambda f, g: f._mul_generic(g),
                                  base_ring_random_elt_args=[],
                                  numtests=10, verbose=False):
    """
    Test univariate Karatsuba multiplication against other multiplication algorithms.

    EXAMPLES:

    First check that random tests are reproducible::

        sage: from sage.rings.tests import check_karatsuba_multiplication
        sage: check_karatsuba_multiplication(ZZ, 6, 5, verbose=True, seed=42)
        check_karatsuba_multiplication: ring=Univariate Polynomial Ring in x over Integer Ring, threshold=2
          (x + 4)*(4*x - 3)
          (-x^6 - x^5 + 2*x^2 + 3*x - 20)*(-2*x^3 + x^2 + 4*x)
          (16*x^5 + x^4 - 41*x^3 + x + 1)*(-x + 8)
          (2*x^2 + x + 3)*(-4)
          (-x^3 - 8*x^2 - x)*(-x)
          (-1)*(-x^2 + x + 1)
          (2*x)*(-x^3 + 3*x^2 - x + 1)
          (x^2 + x + 1)*(-x^4 - x^3 + 76*x^2 + 4*x + 1)
          (6*x^4 - x^3 - 5*x^2 + 1)*(4*x^3 - x^2 + 5*x - 2)
          (3*x - 1)*(x^5 + 21*x^4 + x^3 + 4*x^2 - x)

    Test Karatsuba multiplication of polynomials of small degree over some common rings::

        sage: rings = [QQ]
        sage: rings += [ZZ[I], ZZ[I, sqrt(2)]]                                          # needs sage.rings.number_field sage.symbolic
        sage: rings += [GF(49, 'a')]                                                    # needs sage.rings.finite_rings
        sage: rings += [MatrixSpace(GF(17), 3)]                                         # needs sage.modules
        sage: for C in rings:                                                           # needs sage.modules
        ....:     check_karatsuba_multiplication(C, 10, 10)

    Zero-tests over ``QQbar`` are currently very slow, so we test only very small examples::

        sage: check_karatsuba_multiplication(QQbar, 3, 3, numtests=2)    # long time, needs sage.rings.number_field

    Larger degrees (over ``ZZ``, using FLINT)::

        sage: check_karatsuba_multiplication(ZZ, 1000, 1000,
        ....:                               ref_mul=lambda f,g: f*g,
        ....:                               base_ring_random_elt_args=[1000])

    Some more aggressive tests::

        sage: testrings = [ZZ[I, sqrt(2)], ZZ[I, sqrt(2), sqrt(3)]]     # long time
        sage: for C in testrings:                                       # long time
        ....:     check_karatsuba_multiplication(C, 100, 100)
        sage: check_karatsuba_multiplication(ZZ, 10000, 10000,           # long time
        ....:                               ref_mul=lambda f,g: f*g,
        ....:                               base_ring_random_elt_args=[100000])
    """
    from sage.misc.prandom import randint
    from sage.misc.sage_input import sage_input
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    threshold = randint(0, min(maxdeg1, maxdeg2))
    R = PolynomialRing(base_ring, 'x')
    if verbose:
        print(f"check_karatsuba_multiplication: ring={R}, threshold={threshold}")
    for _ in range(numtests):
        f = R.random_element(randint(0, maxdeg1), False, *base_ring_random_elt_args)
        g = R.random_element(randint(0, maxdeg2), False, *base_ring_random_elt_args)
        if verbose:
            print("  ({})*({})".format(f, g))
        if ref_mul(f, g) - f._mul_karatsuba(g, threshold) != 0:
            msg = "Multiplication failed for elements defined by\n"
            msg += f"{sage_input(f)}\n"
            msg += "and\n"
            msg += f"{sage_input(g)}"
            raise ValueError(msg)
