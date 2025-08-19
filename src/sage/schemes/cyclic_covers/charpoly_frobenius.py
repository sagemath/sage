# sage.doctest: needs sage.rings.padics
r"""
Computation of the Frobenius polynomial using Newton's identities
"""
# *****************************************************************************
#  Copyright (C) 2018 Edgar Costa <edgarc@mit.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.rings.integer_ring import ZZ
from sage.misc.lazy_import import lazy_import
lazy_import("sage.functions.log", "log")


def charpoly_frobenius(frob_matrix, charpoly_prec, p, weight, a=1, known_factor=[1]):
    """
    Return the characteristic polynomial of the given Frobenius matrix.

    INPUT:

    - ``frob_matrix`` -- a matrix representing the Frobenius matrix up to some precision

    - ``charpoly_prec`` -- a vector ai, such that, ``frob_matrix.change_ring(ZZ).charpoly()[i]``
      will be correct mod `p^ai`, this can be easily deduced from the
      Hodge numbers and knowing the q-adic precision of ``frob_matrix``

    - ``p`` -- prime `p`

    - ``weight`` -- weight of the motive

    - ``a`` -- `p = q^a`

    - ``known_factor`` -- the list of coefficients of the known factor

    OUTPUT: list of integers corresponding to the characteristic polynomial
    of the Frobenius action

    EXAMPLES::

        sage: from sage.schemes.cyclic_covers.charpoly_frobenius import charpoly_frobenius
        sage: M = Matrix([[O(17), 8 + O(17)], [O(17), 15 + O(17)]])
        sage: charpoly_frobenius(M, [2, 1, 1], 17, 1, 1)
        [17, 2, 1]

        sage: R = Zq(17**2, names=('a',))
        sage: M = Matrix(R, [[8*17 + 16*17**2 + O(17**3), 8 + 11*17 + O(17**2)],
        ....:                [7*17**2 + O(17**3), 15 + 8*17 + O(17**2)]])
        sage: charpoly_frobenius(M*M, [3, 2, 2], 17, 1, 2)
        [289, 30, 1]

        sage: M = Matrix([[8*31 + 8*31**2 + O(31**3), O(31**3), O(31**3), O(31**3)],
        ....:             [O(31**3), 23*31 + 22*31**2 + O(31**3), O(31**3), O(31**3)],
        ....:             [O(31**3), O(31**3), 27 + 7*31 + O(31**3), O(31**3)],
        ....:             [O(31**3), O(31**3), O(31**3), 4 + 23*31 + O(31**3)]])
        sage: charpoly_frobenius(M, [4, 3, 2, 2, 2], 31, 1, 1)
        [961, 0, 46, 0, 1]

        sage: M = Matrix([(4*43^2 + O(43^3), 17*43 + 11*43^2 + O(43^3), O(43^3), O(43^3), 17 + 37*43 + O(43^3), O(43^3)),
        ....:             (30*43 + 23*43^2 + O(43^3), 5*43 + O(43^3), O(43^3), O(43^3), 3 + 38*43 + O(43^3), O(43^3)),
        ....:             (O(43^3), O(43^3), 9*43 + 32*43^2 + O(43^3), 13 + 25*43 + O(43^3), O(43^3), 17 + 18*43 + O(43^3)),
        ....:             (O(43^3), O(43^3), 22*43 + 25*43^2 + O(43^3), 11 + 24*43 + O(43^3), O(43^3), 36 + 5*43 + O(43^3)),
        ....:             (42*43 + 15*43^2 + O(43^3), 22*43 + 8*43^2 + O(43^3), O(43^3), O(43^3), 29 + 4*43 + O(43^3), O(43^3)),
        ....:             (O(43^3), O(43^3), 6*43 + 19*43^2 + O(43^3), 8 + 24*43 + O(43^3), O(43^3), 31 + 42*43 + O(43^3))])
        sage: charpoly_frobenius(M, [5, 4, 3, 2, 2, 2, 2], 43, 1, 1)
            [79507, 27735, 6579, 1258, 153, 15, 1]

        sage: M = Matrix([(1 + O(4999), O(4999), 0, 0),
        ....:             (O(4999), 4860 + O(4999), 0, 0),
        ....:             (0, 0, O(4999), O(4999)),
        ....:             (0, 0, O(4999), 1 + O(4999))])
        sage: charpoly_frobenius(M, [2, 1, 1], 4999, 1, 1, [1, -2, 1])
        [4999, 139, 1]

    TESTS::

        sage: M = Matrix([[-149196156000219, 0, 0, 0, 0, 0, 0, 0],
        ....:             [0, 76324364094257, 0, 0, 0, 0, 0, 0],
        ....:             [0, 0, 76324364094257, 0, 0, 0, 0, 0],
        ....:             [0, 0, 0, -149196156000219, 0, 0, 0, 0],
        ....:             [0, 0, 0, 0, 281855171388275, 0, 0, 0],
        ....:             [0, 0, 0, 0, 0, -208983379482579, 0, 0],
        ....:             [0, 0, 0, 0, 0, 0, -208983379482579, 0],
        ....:             [0, 0, 0, 0, 0, 0, 0, 281855171388275]])
        sage: charpoly_frobenius(M, [9, 8, 7, 6, 5, 5, 5, 5, 5], 1009, 1, 2)
        [1074309286591662654798721,
         561382189105547134612,
         -2982540407204025062,
         -247015136050256,
         4390163797795,
         -242628176,
         -2877542,
         532,
         1]
        sage: M = Matrix([[0, 0, 0, -338082603, 0, 0, 0, 0],
        ....:             [0, 0, -317436968, 0, 0, 0, 0, 0],
        ....:             [0, -120741807, 0, 0, 0, 0, 0, 0],
        ....:             [200618482, 0, 0, 0, 0, 0, 0, 0],
        ....:             [0, 0, 0, 0, 0, 0, 0, 123492519],
        ....:             [0, 0, 0, 0, 0, 0, 426826171, 0],
        ....:             [0, 0, 0, 0, 0, 157417117, 0, 0],
        ....:             [0, 0, 0, 0, 373415235, 0, 0, 0]])
        sage: charpoly_frobenius(M, [7, 6, 5, 4, 3, 3, 3, 3, 3], 1009, 1, 1)
        [1036488922561, 0, 270809546, 0, -1474149, 0, 266, 0, 1]

        sage: M = Matrix({(0, 31): 1814236329200021268558465351501717,
        ....: (1, 30): 3268331092352160631300311212049390,
        ....: (2, 29): 1002349136486054751305109007707560,
        ....: (3, 28): 1789497403160078628636360424523308,
        ....: (4, 19): 919866278512654133838788268427125,
        ....: (5, 18): 2918980842679879118243999587726673,
        ....: (6, 17): 2062741569795231121341967954037400,
        ....: (7, 16): 3562554496811633214919332352788305,
        ....: (8, 7): 287823825201170974551150606916601,
        ....: (9, 6): 2657175570144838727074228404244845,
        ....: (10, 5): 3200631048273888400670606576807785,
        ....: (11, 4): 707085630754978281870563133348521,
        ....: (12, 39): 679572779843478608532167180287595,
        ....: (13, 38): 510867456922807824071915371084390,
        ....: (14, 37): 3300741705093235469798877501619286,
        ....: (15, 36): 1374430202827161695034370373469332,
        ....: (16, 27): 1897240889699239396313755822318254,
        ....: (17, 26): 3171751877741319729745976757727266,
        ....: (18, 25): 1151779650995750952707414056498421,
        ....: (19, 24): 1309748952162524211332312241346156,
        ....: (20, 15): 2914640274871541651939754878647777,
        ....: (21, 14): 2524322227034087814555116576604052,
        ....: (22, 13): 693999428630644346611319813759997,
        ....: (23, 12): 2093267437436875555592094407087011,
        ....: (24, 3): 101158112439244133585487537448909,
        ....: (25, 2): 638873050956374173808321501215560,
        ....: (26, 1): 3529335795023815426485172749287314,
        ....: (27, 0): 618726320422582798159865537548600,
        ....: (28, 35): 2510605595766272594980682702750921,
        ....: (29, 34): 2978146199632282120435531158312695,
        ....: (30, 33): 1724161588290366191539756998844438,
        ....: (31, 32): 516507426627993787229114955328811,
        ....: (32, 23): 1716672265998537901154333190869011,
        ....: (33, 22): 3787144776814278856737374038432424,
        ....: (34, 21): 3765560528316833596614887925578722,
        ....: (35, 20): 1628311006615824767735977131865996,
        ....: (36, 11): 3638935478569769465046956942756848,
        ....: (37, 10): 1878821491042105813643148323053706,
        ....: (38, 9): 1187568624951630613061547491748348,
        ....: (39, 8): 2538351040819233009959661983810741}
        ....: )
        sage: charpoly_frobenius(M,
        ....: [31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16,
        ....:  15, 14, 13, 12] + [11]*21, 1129, 1, 1)
        [11320844849639649951608809973589776933203136765026963553258401,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         24687045654725446027864774006541463602997309796,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         20187877911930897108199045855206,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         7337188909826596,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         0,
         1]

        sage: F = Matrix(Qp(17),
        ....: [(28442601332527957763, 729848492961404015, 70994086070709920),
        ....: (24928804992606688137, 1345506389644311177, 147442915782003034),
        ....: (7562462964206075698, 1262441299395996535, 92309755559576133)])
        sage: F+= F.base_ring()(0).add_bigoh(6)*ones_matrix(*F.dimensions())
        sage: charpoly_frobenius(F, [6, 5, 4, 4], 17, 2)
        [-4913, -221, 13, 1]
    """
    assert known_factor[-1] == 1
    try:
        cp = frob_matrix.change_ring(ZZ).charpoly().list()
    except ValueError:
        # the given matrix was not integral
        cp = frob_matrix.charpoly().change_ring(ZZ).list()
    assert len(charpoly_prec) == len(cp) - (len(known_factor) - 1)
    assert cp[-1] == 1

    # reduce cp mod prec
    degree = len(charpoly_prec) - 1
    mod = [0] * (degree + 1)
    for i in range(len(charpoly_prec)):
        mod[-i] = p**charpoly_prec[-i]
        cp[-i] = cp[-i] % mod[-i]

    # figure out the sign
    # i.e., if it is a reciprocal or an antireciprocal polynomial
    if weight % 2:
        # for odd weight the sign is always 1
        # it's the charpoly of a USp matrix
        # and charpoly of a symplectic matrix is reciprocal
        sign = 1
    else:
        # For the moment I will not worry about this case
        if known_factor != [1]:
            raise NotImplementedError()
        # we compare ith coefficient and  (degree - i)th coefficient to deduce the sign
        # note, if degree is even, the middle coefficient will not help us determine the sign
        for i in range((degree + 1) // 2):
            # Note: degree*weight is even
            p_power = p**min(
                charpoly_prec[i],
                charpoly_prec[degree - i] + ((a * (degree - 2 * i) * weight) // 2),
            )
            if cp[i] % p_power != 0 and cp[degree - i] % p_power != 0:
                other = cp[degree - i] * p**((a * (degree - 2 * i) * weight) // 2)
                if (cp[i] + other) % p_power == 0:
                    sign = -1
                else:
                    sign = 1
                assert (-sign * cp[i] + other) % p_power == 0
                break
    # halfdegree is the number of coefficients that we will compute
    # the rest will be deduced using the functional equation
    # as up to scaling of the variable
    # the polynomial is either reciprocal or antireciprocal polynomial
    # note, this includes the middle coefficient if degree is even
    halfdegree = degree // 2 + 1

    cp[0] = sign * p**((a * degree * weight) // 2)
    # Note: degree*weight is even

    # calculate the i-th power sum of the roots and correct cp along the way
    e = cp[-halfdegree:]
    e.reverse()
    for k in range(halfdegree):
        if k % 2:
            e[k] = -e[k] % mod[degree - k]
        # e[k] = cp[degree - k] if (k%2 ==0) else -cp[degree - k]
        if k > 0:
            # verify if p^charpoly_prec[degree - k] > 2*degree/k * q^(w*k/2)
            assert (
                log(k, p) + charpoly_prec[degree - k]
                > log(2 * degree, p) + a * 0.5 * weight * k
            ), (
                "log(k)/log(p) + charpoly_prec[degree - k] <= log(2*degree)/log(p) + a*0.5*weight*k, k = %d"
                % k
            )

    fix_e = known_factor[:]
    fix_e.reverse()
    if len(fix_e) < halfdegree:
        fix_e.extend([0] * (halfdegree - len(fix_e)))
    for i in range(halfdegree):
        if i % 2:
            fix_e[i] *= -1

    # e[k] = \sum x_{i_1} x_{i_2} ... x_{i_k} # where x_* are eigenvalues
    # and i_1 < i_2 ... < i_k

    # s[k] = \sum x_i ^k for k>0
    s = [None] * (halfdegree)
    res = [None] * len(charpoly_prec)
    res[0] = sign * p**((a * degree * weight) // 2)
    # Note: degree*weight is even

    res[-1] = 1
    e[1] -= fix_e[1]
    e[1] = e[1] % mod[degree - 1]
    for k in range(1, halfdegree):
        # assume that s[i] and e[i] are correct for i < k
        # e[k] correct modulo mod[degree - k]
        # S = sum (-1)^i e[k-i] * s[i]
        # s[k] = (-1)^(k-1) (k*e[k] + S) ==> (-1)^(k-1) s[k] - S = k*e[k]
        S = sum((-1)**i * e[k - i] * s[i] for i in range(1, k))
        s[k] = (-1)**(k - 1) * (S + k * e[k])
        # hence s[k] is correct modulo k*mod[degree - k]
        localmod = k * mod[degree - k]
        # s[k] +=   (-1)**k * fix_power_sum[k]
        s[k] = s[k] % localmod

        # |x_i| = p^(w*0.5)
        # => s[k] <= degree*p^(a*w*k*0.5)
        # recall, 2*degree*p^(a*w*k*0.5) /k < mod[degree - k]
        if s[k]**2 > degree**2 * p**(a * weight * k):
            s[k] = -(-s[k] % localmod)

        # now correct e[k] with:
        # (-1)^(k-1) s[k] - S = k*e[k]
        e[k] = (-S + (-1)**(k - 1) * s[k]) // k
        assert (-S + (-1)**(k - 1) * s[k]) % k == 0
        res[degree - k] = e[k] if not k % 2 else -e[k]
        # Note: degree*weight is even

        res[k] = sign * res[degree - k] * p**((a * (degree - 2 * k) * weight) // 2)
        # fix e[k + 1]
        if k + 1 < halfdegree:
            e[k + 1] -= sum([fix_e[k + 1 - i] * e[i] for i in range(k + 1)])
            e[k + 1] = e[k + 1] % mod[degree - (k + 1)]
    return res
