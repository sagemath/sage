"""
Cython interface to Cremona's ``eclib`` library (also known as ``mwrank``)

EXAMPLES::

    sage: from sage.libs.eclib.mwrank import _Curvedata, _mw
    sage: c = _Curvedata(1,2,3,4,5)

    sage: print(c)
    [1,2,3,4,5]
    b2 = 9       b4 = 11         b6 = 29         b8 = 35
    c4 = -183           c6 = -3429
    disc = -10351       (# real components = 1)
    #torsion not yet computed

    sage: t= _mw(c)
    sage: t.search(10)
    sage: t
    [[1:2:1]]
"""

import os
import sys

from cysignals.memory cimport sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.cpython.string import FS_ENCODING
from sage.libs.eclib cimport bigint, Curvedata, mw, two_descent

cdef extern from "wrap.cpp":
    ### misc functions ###
    long mwrank_get_precision()
    void mwrank_set_precision(long n)
    void mwrank_initprimes(char* pfilename, int verb)

    ### bigint ###
    bigint* str_to_bigint(char* s)
    char* bigint_to_str(bigint* x)

    ### Curvedata ###
    char* Curvedata_repr(Curvedata* curve)
    double Curvedata_silverman_bound(Curvedata* curve)
    double Curvedata_cps_bound(Curvedata* curve)
    double Curvedata_height_constant(Curvedata* curve)
    char* Curvedata_getdiscr(Curvedata* curve)
    char* Curvedata_conductor(Curvedata* m)
    char* Curvedata_isogeny_class(Curvedata* E, int verbose)

    ## mw ##
    int mw_process(Curvedata* curve, mw* m,
                   bigint* x, bigint* y,
                   bigint* z, int sat)
    char* mw_getbasis(mw* m)
    double mw_regulator(mw* m)
    int mw_rank(mw* m)
    int mw_saturate(mw* m, bigint* index, char** unsat,
                    long sat_bd, int odd_primes_only)
    void mw_search(mw* m, char* h_lim, int moduli_option, int verb)

    ### two_descent ###
    int two_descent_ok(two_descent* t)
    long two_descent_get_certain(two_descent* t)
    char* two_descent_get_basis(two_descent* t)
    double two_descent_regulator(two_descent* t)
    long two_descent_get_rank(two_descent* t)
    long two_descent_get_rank_bound(two_descent* t)
    long two_descent_get_selmer_rank(two_descent* t)
    void two_descent_saturate(two_descent* t, long sat_bd)


cdef object string_sigoff(char* s):
    sig_off()
    # Makes a python string and deletes what is pointed to by s.
    t = char_to_str(s)
    sig_free(s)
    return t


# set the default bit precision
mwrank_set_precision(150)

def get_precision():
    """
    Returns the working floating point bit precision of mwrank, which is
    equal to the global NTL real number precision.

    OUTPUT:

    (int) The current precision in bits.

    See also :meth:`set_precision`.

    EXAMPLES::

        sage: mwrank_get_precision()
        150
    """
    return mwrank_get_precision()


def set_precision(n):
    """
    Sets the working floating point bit precision of mwrank, which is
    equal to the global NTL real number precision.

    NTL real number bit precision.  This has a massive effect on the
    speed of mwrank calculations.  The default (used if this function is
    not called) is ``n=150``, but it might have to be increased if a
    computation fails.

    INPUT:

    - ``n`` -- a positive integer: the number of bits of precision.

    .. warning::

       This change is global and affects *all* future calls of eclib
       functions by Sage.

    .. note::

        The minimal value to which the precision may be set is 53.
        Lower values will be increased to 53.

    See also :meth:`get_precision`.

    EXAMPLES::

        sage: from sage.libs.eclib.mwrank import set_precision, get_precision
        sage: old_prec = get_precision(); old_prec
        150
        sage: set_precision(50)
        sage: get_precision()
        53
        sage: set_precision(old_prec)
        sage: get_precision()
        150
    """
    mwrank_set_precision(n)


def initprimes(filename, verb=False):
    """
    Initialises mwrank/eclib's internal prime list.

    INPUT:

    - ``filename`` (string) -- the name of a file of primes.

    - ``verb`` (bool: default ``False``) -- verbose or not?

    EXAMPLES::

        sage: file = os.path.join(SAGE_TMP, 'PRIMES')
        sage: with open(file, 'w') as fobj:
        ....:     _ = fobj.write(' '.join([str(p) for p in prime_range(10^7,10^7+20)]))
        sage: mwrank_initprimes(file, verb=True)
        Computed 78519 primes, largest is 1000253
        reading primes from file ...
        read extra prime 10000019
        finished reading primes from file ...
        Extra primes in list: 10000019

        sage: mwrank_initprimes("x" + file, True)
        Traceback (most recent call last):
        ...
        IOError: No such file or directory: ...
    """
    if not os.path.exists(filename):
        raise IOError('No such file or directory: %s' % filename)
    filename = str_to_bytes(filename, FS_ENCODING, 'surrogateescape')
    mwrank_initprimes(filename, verb)

############# bigint ###########################################
#
# In mwrank (and eclib) bigint is synonymous with NTL's ZZ class.
#
################################################################

cdef class _bigint:
    """
    Cython class wrapping eclib's bigint class.
    """
    cdef bigint* x

    def __init__(self, x="0"):
        """
        Constructor for bigint class.

        INPUT:

        - ``x`` -- string or int: a string representing a decimal
          integer, or a Sage integer

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _bigint
            sage: _bigint(123)
            123
            sage: _bigint('123')
            123
            sage: type(_bigint(123))
            <type 'sage.libs.eclib.mwrank._bigint'>
        """
        s = str(x)
        if s.isdigit() or s[0] == "-" and s[1:].isdigit():
            self.x = str_to_bigint(str_to_bytes(s))
        else:
            raise ValueError("invalid _bigint: %r"%x)

    def __dealloc__(self):
        """
        Destructor for bigint class (releases memory).
        """
        del self.x

    def __repr__(self):
        """
        String representation of bigint.

        OUTPUT:

        (string) the bigint as a string (decimal)

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _bigint
            sage: a = _bigint('123')
            sage: a.__repr__()
            '123'
            sage: a = _bigint('-456')
            sage: a.__repr__()
            '-456'
        """
        sig_on()
        return string_sigoff(bigint_to_str(self.x))


cdef make_bigint(bigint* x):
    cdef _bigint y
    sig_off()
    y = _bigint.__new__(_bigint)
    y.x = x
    return y

############# Curvedata #################

cdef class _Curvedata:   # cython class wrapping eclib's Curvedata class
    cdef Curvedata* x

    def __init__(self, a1, a2, a3,
                 a4, a6, min_on_init=0):
        """
        Constructor for Curvedata class.

        INPUT:

        - ``a1``, ``a2``, ``a3``, ``a4``, ``a6`` (int) -- integer
          coefficients of a Weierstrass equation (must be
          nonsingular).

        - ``min_on_init`` (int, default 0) -- flag controlling whether
          the constructed curve is replaced by a global minimal model.
          If nonzero then this minimisation does take place.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: _Curvedata(1,2,3,4,5)
            [1,2,3,4,5]
            b2 = 9       b4 = 11         b6 = 29         b8 = 35
            c4 = -183           c6 = -3429
            disc = -10351       (# real components = 1)
            #torsion not yet computed

        A non-minimal example::

            sage: _Curvedata(0,0,0,0,64)
            [0,0,0,0,64]
            b2 = 0       b4 = 0  b6 = 256        b8 = 0
            c4 = 0              c6 = -55296
            disc = -1769472     (# real components = 1)
            #torsion not yet computed

            sage: _Curvedata(0,0,0,0,64,min_on_init=1)
            [0,0,0,0,1] (reduced minimal model)
            b2 = 0       b4 = 0  b6 = 4  b8 = 0
            c4 = 0              c6 = -864
            disc = -432 (# real components = 1)
            #torsion not yet computed
        """
        cdef _bigint _a1, _a2, _a3, _a4, _a6
        _a1 = _bigint(a1)
        _a2 = _bigint(a2)
        _a3 = _bigint(a3)
        _a4 = _bigint(a4)
        _a6 = _bigint(a6)
        self.x = new Curvedata(_a1.x[0], _a2.x[0], _a3.x[0], _a4.x[0], _a6.x[0], min_on_init)
        if self.discriminant() == 0:
            msg = "Invariants (= {},{},{},{},{}) do not describe an elliptic curve."
            raise ArithmeticError(msg.format(a1, a2, a3, a4, a6))

    def __dealloc__(self):
        """
        Destructor for Curvedata class.
        """
        del self.x

    def __repr__(self):
        r"""
        String representation of Curvedata

        OUTPUT:

        (string) the Curvedata as a string

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,2,3,4,5)
            sage: E.__repr__()
            '[1,2,3,4,5]\nb2 = 9\t b4 = 11\t b6 = 29\t b8 = 35\nc4 = -183\t\tc6 = -3429\ndisc = -10351\t(# real components = 1)\n#torsion not yet computed'
            sage: E
            [1,2,3,4,5]
            b2 = 9       b4 = 11         b6 = 29         b8 = 35
            c4 = -183           c6 = -3429
            disc = -10351       (# real components = 1)
            #torsion not yet computed
            """
        sig_on()
        return string_sigoff(Curvedata_repr(self.x))[:-1]

    def silverman_bound(self):
        """
        The Silverman height bound for this elliptic curve.

        OUTPUT:

        (float) A non-negative real number `B` such that for every
        rational point on this elliptic curve `E`, `h(P)\le\hat{h}(P)
        + B`, where `h(P)` is the naive height and `\hat{h}(P)` the
        canonical height.

        .. note::

            Since eclib can compute this to arbitrary precision, we
            could return a Sage real, but this is only a bound and in
            the contexts in which it is used extra precision is
            irrelevant.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,2,3,4,5)
            sage: E.silverman_bound()
            6.52226179519101...
            sage: type(E.silverman_bound())
            <type 'float'>
        """
        return Curvedata_silverman_bound(self.x)

    def cps_bound(self):
        """
        The Cremona-Prickett-Siksek height bound for this elliptic curve.

        OUTPUT:

        (float) A non-negative real number `B` such that for every
        rational point on this elliptic curve `E`, `h(P)\le\hat{h}(P)
        + B`, where `h(P)` is the naive height and `\hat{h}(P)` the
        canonical height.

        .. note::

            Since eclib can compute this to arbitrary precision, we
            could return a Sage real, but this is only a bound and in
            the contexts in which it is used extra precision is
            irrelevant.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,2,3,4,5)
            sage: E.cps_bound()
            0.11912451909250982...

        Note that this is a better bound than Silverman's in this case::

            sage: E.silverman_bound()
            6.52226179519101...
        """
        sig_on()
        x = Curvedata_cps_bound(self.x)
        sig_off()
        return x

    def height_constant(self):
        """
        A height bound for this elliptic curve.

        OUTPUT:

        (float) A non-negative real number `B` such that for every
        rational point on this elliptic curve `E`, `h(P)\le\hat{h}(P)
        + B`, where `h(P)` is the naive height and `\hat{h}(P)` the
        canonical height.  This is the minimum of the Silverman and
        Cremona_Prickett-Siksek height bounds.

        .. note::

            Since eclib can compute this to arbitrary precision, we
            could return a Sage real, but this is only a bound and in
            the contexts in which it is used extra precision is
            irrelevant.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,2,3,4,5)
            sage: E.height_constant()
            0.119124519092509...
        """
        return Curvedata_height_constant(self.x)

    def discriminant(self):
        """
        The discriminant of this elliptic curve.

        OUTPUT:

        (Integer) The discriminant.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,2,3,4,5)
            sage: E.discriminant()
            -10351
            sage: E = _Curvedata(100,200,300,400,500)
            sage: E.discriminant()
            -1269581104000000
            sage: ZZ(E.discriminant())
            -1269581104000000
        """
        sig_on()
        from sage.rings.all import Integer
        return Integer(string_sigoff(Curvedata_getdiscr(self.x)))

    def conductor(self):
        """
        The conductor of this elliptic curve.

        OUTPUT:

        (Integer) The conductor.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,2,3,4,5)
            sage: E.discriminant()
            -10351
            sage: E = _Curvedata(100,200,300,400,500)
            sage: E.conductor()
            126958110400
        """
        sig_on()
        from sage.rings.all import Integer
        return Integer(string_sigoff(Curvedata_conductor(self.x)))

    def isogeny_class(self, verbose=False):
        """
        The isogeny class of this elliptic curve.

        OUTPUT:

        (tuple) A tuple consisting of (1) a list of the curves in the
        isogeny class, each as a list of its Weierstrass coefficients;
        (2) a matrix of the degrees of the isogenies between the
        curves in the class (prime degrees only).

        .. warning::

           The list may not be complete, if the precision is too low.
           Use ``mwrank_set_precision()`` to increase the precision.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,0,1,4,-6)
            sage: E.conductor()
            14
            sage: E.isogeny_class()
            ([[1, 0, 1, 4, -6], [1, 0, 1, -36, -70], [1, 0, 1, -1, 0], [1, 0, 1, -171, -874], [1, 0, 1, -11, 12], [1, 0, 1, -2731, -55146]], [[0, 2, 3, 3, 0, 0], [2, 0, 0, 0, 3, 3], [3, 0, 0, 0, 2, 0], [3, 0, 0, 0, 0, 2], [0, 3, 2, 0, 0, 0], [0, 3, 0, 2, 0, 0]])

        """
        sig_on()
        s = string_sigoff(Curvedata_isogeny_class(self.x, verbose))
        return eval(s)


############# _mw #################

cdef class _mw:
    """
    Cython class wrapping eclib's mw class.
    """
    cdef mw* x
    cdef Curvedata* curve
    cdef int verb

    def __init__(self, _Curvedata curve, verb=False, pp=1, maxr=999):
        """
        Constructor for mw class.

        INPUT:

        - ``curve`` (_Curvedata) -- an elliptic curve

        - ``verb`` (bool, default ``False``) -- verbosity flag (controls
          amount of output produced in point searches)

        - ``pp`` (int, default 1) -- process points flag (if nonzero,
          the points found are processed, so that at all times only a
          `\ZZ`-basis for the subgroup generated by the points found
          so far is stored; if zero, no processing is done and all
          points found are stored).

        - ``maxr`` (int, default 999) -- maximum rank (quit point
          searching once the points found generate a subgroup of this
          rank; useful if an upper bound for the rank is already
          known).

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _mw
            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(1,0,1,4,-6)
            sage: EQ = _mw(E)
            sage: EQ
            []
            sage: type(EQ)
            <type 'sage.libs.eclib.mwrank._mw'>

            sage: E = _Curvedata(0,0,1,-7,6)
            sage: EQ = _mw(E)
            sage: EQ.search(2)
            sage: EQ
            [[1:-1:1], [-2:3:1], [-14:25:8]]

        Example to illustrate the verbose parameter::

            sage: from sage.libs.eclib.mwrank import _mw
            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(0,0,1,-7,6)
            sage: EQ = _mw(E, verb=False)
            sage: EQ.search(1)
            sage: EQ = _mw(E, verb=True)
            sage: EQ.search(1)
            P1 = [0:1:0]         is torsion point, order 1
            P1 = [-3:0:1]         is generator number 1
            ...
            P4 = [12:35:27]      = 1*P1 + -1*P2 + -1*P3 (mod torsion)

        The previous command produces the following output::

            P1 = [0:1:0]         is torsion point, order 1
            P1 = [-3:0:1]         is generator number 1
            saturating up to 20...Checking 2-saturation
            Points have successfully been 2-saturated (max q used = 7)
            Checking 3-saturation
            Points have successfully been 3-saturated (max q used = 7)
            Checking 5-saturation
            Points have successfully been 5-saturated (max q used = 23)
            Checking 7-saturation
            Points have successfully been 7-saturated (max q used = 41)
            Checking 11-saturation
            Points have successfully been 11-saturated (max q used = 17)
            Checking 13-saturation
            Points have successfully been 13-saturated (max q used = 43)
            Checking 17-saturation
            Points have successfully been 17-saturated (max q used = 31)
            Checking 19-saturation
            Points have successfully been 19-saturated (max q used = 37)
            done
            P2 = [-2:3:1]         is generator number 2
            saturating up to 20...Checking 2-saturation
            possible kernel vector = [1,1]
            This point may be in 2E(Q): [14:-52:1]
            ...and it is!
            Replacing old generator #1 with new generator [1:-1:1]
            Points have successfully been 2-saturated (max q used = 7)
            Index gain = 2^1
            Checking 3-saturation
            Points have successfully been 3-saturated (max q used = 13)
            Checking 5-saturation
            Points have successfully been 5-saturated (max q used = 67)
            Checking 7-saturation
            Points have successfully been 7-saturated (max q used = 53)
            Checking 11-saturation
            Points have successfully been 11-saturated (max q used = 73)
            Checking 13-saturation
            Points have successfully been 13-saturated (max q used = 103)
            Checking 17-saturation
            Points have successfully been 17-saturated (max q used = 113)
            Checking 19-saturation
            Points have successfully been 19-saturated (max q used = 47)
            done (index = 2).
            Gained index 2, new generators = [ [1:-1:1] [-2:3:1] ]
            P3 = [-14:25:8]       is generator number 3
            saturating up to 20...Checking 2-saturation
            Points have successfully been 2-saturated (max q used = 11)
            Checking 3-saturation
            Points have successfully been 3-saturated (max q used = 13)
            Checking 5-saturation
            Points have successfully been 5-saturated (max q used = 71)
            Checking 7-saturation
            Points have successfully been 7-saturated (max q used = 101)
            Checking 11-saturation
            Points have successfully been 11-saturated (max q used = 127)
            Checking 13-saturation
            Points have successfully been 13-saturated (max q used = 151)
            Checking 17-saturation
            Points have successfully been 17-saturated (max q used = 139)
            Checking 19-saturation
            Points have successfully been 19-saturated (max q used = 179)
            done (index = 1).
            P4 = [-1:3:1]        = -1*P1 + -1*P2 + -1*P3 (mod torsion)
            P4 = [0:2:1]         = 2*P1 + 0*P2 + 1*P3 (mod torsion)
            P4 = [2:13:8]        = -3*P1 + 1*P2 + -1*P3 (mod torsion)
            P4 = [1:0:1]         = -1*P1 + 0*P2 + 0*P3 (mod torsion)
            P4 = [2:0:1]         = -1*P1 + 1*P2 + 0*P3 (mod torsion)
            P4 = [18:7:8]        = -2*P1 + -1*P2 + -1*P3 (mod torsion)
            P4 = [3:3:1]         = 1*P1 + 0*P2 + 1*P3 (mod torsion)
            P4 = [4:6:1]         = 0*P1 + -1*P2 + -1*P3 (mod torsion)
            P4 = [36:69:64]      = 1*P1 + -2*P2 + 0*P3 (mod torsion)
            P4 = [68:-25:64]     = -2*P1 + -1*P2 + -2*P3 (mod torsion)
            P4 = [12:35:27]      = 1*P1 + -1*P2 + -1*P3 (mod torsion)
            sage: EQ
            [[1:-1:1], [-2:3:1], [-14:25:8]]

        Example to illustrate the process points ``pp`` parameter::

            sage: from sage.libs.eclib.mwrank import _mw
            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: E = _Curvedata(0,0,1,-7,6)
            sage: EQ = _mw(E, pp=1)
            sage: EQ.search(1); EQ
            [[1:-1:1], [-2:3:1], [-14:25:8]]
            sage: EQ = _mw(E, pp=0)
            sage: EQ.search(1); EQ
            [[-3:0:1], [-2:3:1], [-14:25:8], [-1:3:1], [0:2:1], [2:13:8], [1:0:1], [2:0:1], [18:7:8], [3:3:1], [4:6:1], [36:69:64], [68:-25:64], [12:35:27]]

        """
        self.curve = curve.x
        self.x = new mw(curve.x, verb, pp, maxr)
        self.verb = verb

    def __dealloc__(self):
        """
        Destructor for mw class.
        """
        del self.x

    def __repr__(self):
        """
        String representation of the current basis of this mw group.

        OUTPUT:

        (string) the current basis of this mw as a string

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,0,1,-7,6)
            sage: EQ = _mw(E)
            sage: EQ # indirect doctest
            []
            sage: EQ.search(2)
            sage: EQ
            [[1:-1:1], [-2:3:1], [-14:25:8]]
        """
        sig_on()
        return string_sigoff(mw_getbasis(self.x))

    def process(self, point, sat=0):
        """
        Processes the given point, adding it to the mw group.

        INPUT:

        - ``point`` (tuple or list) -- tuple or list of 3 integers.
          An ``ArithmeticError`` is raised if the point is not on the
          curve.

        - ``sat`` (int, default 0) --saturate at primes up to ``sat``.
          No saturation is done if ``sat=0``.  (Note that it is more
          efficient to add several points at once and then saturate
          just once at the end).

        .. NOTE::

           The eclib function which implements this only carries out
           any saturation if the rank of the points increases upon
           adding the new point.  This is because it is assumed that
           one saturates as ones goes along.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,1,1,-2,0)
            sage: EQ = _mw(E)

        Initially the list of gens is empty::

            sage: EQ
            []

        We process a point of infinite order::

            sage: EQ.process([-1,1,1])
            sage: EQ
            [[-1:1:1]]

        And another independent one::

            sage: EQ.process([0,-1,1])
            sage: EQ
            [[-1:1:1], [0:-1:1]]

        Processing a point dependent on the current basis will not
        change the basis::

            sage: EQ.process([4,8,1])
            sage: EQ
            [[-1:1:1], [0:-1:1]]

        """
        if not isinstance(point, (tuple, list)) and len(point) == 3:
            raise TypeError("point must be a list or tuple of length 3.")
        cdef _bigint x,y,z
        sig_on()
        x,y,z = _bigint(point[0]), _bigint(point[1]), _bigint(point[2])
        r = mw_process(self.curve, self.x, x.x, y.x, z.x, sat)
        sig_off()
        if r != 0:
            raise ArithmeticError("point (=%s) not on curve." % point)

    def getbasis(self):
        """
        Returns the current basis of the mw structure.

        OUTPUT:

        (string) String representation of the points in the basis of
        the mw group.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,1,1,-2,0)
            sage: EQ = _mw(E)
            sage: EQ.search(3)
            sage: EQ.getbasis()
            '[[0:-1:1], [-1:1:1]]'
            sage: EQ.rank()
            2
        """
        sig_on()
        s = string_sigoff(mw_getbasis(self.x))
        return s

    def regulator(self):
        """
        Returns the regulator of the current basis of the mw group.

        OUTPUT:

        (double) The current regulator.

        .. TODO::

            ``eclib`` computes the regulator to arbitrary precision, and
            the full precision value should be returned.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,1,1,-2,0)
            sage: EQ = _mw(E)
            sage: EQ.search(3)
            sage: EQ.getbasis()
            '[[0:-1:1], [-1:1:1]]'
            sage: EQ.rank()
            2
            sage: EQ.regulator()
            0.15246017794314376
        """
        sig_on()
        f = mw_regulator(self.x)
        sig_off()
        return f

    def rank(self):
        """
        Returns the rank of the current basis of the mw group.

        OUTPUT:

        (Integer) The current rank.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,1,1,-2,0)
            sage: EQ = _mw(E)
            sage: EQ.search(3)
            sage: EQ.getbasis()
            '[[0:-1:1], [-1:1:1]]'
            sage: EQ.rank()
            2
        """
        sig_on()
        r = mw_rank(self.x)
        sig_off()
        from sage.rings.all import Integer
        return Integer(r)

    def saturate(self, int sat_bd=-1, int odd_primes_only=0):
        """
        Saturates the current subgroup of the mw group.

        INPUT:

        - ``sat_bnd`` (int, default -1) -- bound on primes at which to
          saturate.  If -1 (default), compute a bound for the primes
          which may not be saturated, and use that.

        - ``odd_primes_only`` (bool, default ``False``) -- only do
          saturation at odd primes.  (If the points have been found
          via 2-descent they should already be 2-saturated.)

        OUTPUT:

        (tuple) (success flag, index, list) The success flag will be 1
        unless something failed (usually an indication that the points
        were not saturated but the precision is not high enough to
        divide out successfully).  The index is the index of the mw
        group before saturation in the mw group after.  The list is a
        string representation of the primes at which saturation was
        not proved or achieved.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,1,1,-2,0)
            sage: EQ = _mw(E)
            sage: EQ.process([494, -5720, 6859]) # 3 times another point
            sage: EQ
            [[494:-5720:6859]]
            sage: EQ.saturate()
            (1, 3, '[ ]')
            sage: EQ
            [[-1:1:1]]

        If we set the saturation bound at 2, then saturation will fail::

            sage: EQ = _mw(E)
            sage: EQ.process([494, -5720, 6859]) # 3 times another point
            sage: EQ.saturate(sat_bd=2)
            Saturation index bound = 10
            WARNING: saturation at primes p > 2 will not be done;
            points may be unsaturated at primes between 2 and index bound
            Failed to saturate MW basis at primes [ ]
            (0, 1, '[ ]')
            sage: EQ
            [[494:-5720:6859]]

        The following output is also seen in the preceding example::

            Saturation index bound = 10
            WARNING: saturation at primes p > 2 will not be done;
            points may be unsaturated at primes between 2 and index bound
            Failed to saturate MW basis at primes [ ]


        """
        cdef _bigint index
        cdef char* s
        cdef int ok
        sig_on()
        index = _bigint()
        ok = mw_saturate(self.x, index.x, &s, sat_bd, odd_primes_only)
        unsat = string_sigoff(s)
        return ok, index, unsat

    def search(self, h_lim, int moduli_option=0, int verb=0):
        """
        Search for points in the mw group.

        INPUT:

        - ``h_lim`` (int) -- bound on logarithmic naive height of points

        - ``moduli_option`` (int, default 0) -- option for sieving
          strategy.  The default (0) uses an adapted version of
          Stoll's ratpoints code and is recommended.

        - ``verb`` (int, default 0) -- level of verbosity.  If 0, no
          output.  If positive, the points are output as found and
          some details of the processing, finding linear relations,
          and partial saturation are output.

        .. NOTE::

           The effect of the search is also governed by the class
           options, notably whether the points found are processed:
           meaning that linear relations are found and saturation is
           carried out, with the result that the list of generators
           will always contain a `\ZZ`-span of the saturation of the
           points found, modulo torsion.

        OUTPUT:

        None.  The effect of the search is to update the list of
        generators.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: from sage.libs.eclib.mwrank import _mw
            sage: E = _Curvedata(0,0,1,-19569,-4064513) # 873c1
            sage: EQ = _mw(E)
            sage: EQ = _mw(E)
            sage: for i in [1..11]:
            ....:     print("{} {} {}".format(i, EQ.search(i), EQ))
            1 None []
            2 None []
            3 None []
            4 None []
            5 None []
            6 None []
            7 None []
            8 None []
            9 None []
            10 None []
            11 None [[3639568:106817593:4096]]
        """

        h_lim = str_to_bytes(str(h_lim))

        sig_on()
        mw_search(self.x, h_lim, moduli_option, verb)
        sig_off()


############# two_descent #################
cdef class _two_descent:
    """
    Cython class wrapping eclib's two_descent class.
    """
    cdef two_descent* x

    def __init__(self):
        """
        Constructor for two_descent class.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
        """
        self.x = <two_descent*> 0

    def __dealloc__(self):
        """
        Destructor for two_descent class.
        """
        del self.x

    def do_descent(self, _Curvedata curve,
                 int verb = 1,
                 int sel = 0,
                 int firstlim = 20,
                 int secondlim = 8,
                 int n_aux = -1,
                 int second_descent = 1):
        """
        Carry out a 2-descent.

        INPUT:

        - ``curvedata`` (_Curvedata) -- the curve on which to do descent.

        - ``verb`` (int, default 1) -- verbosity level.

        - ``sel`` (int, default 0) -- Selmer-only flag.  If 1, only
          the 2-Selmer group will be computed, with no rational
          points.  Useful as a faster way of getting an upper bound on
          the rank.

        - ``firstlim`` (int, default 20) -- naive height bound on
          first point search on quartic homogeneous spaces (before
          testing local solubility; very simple search with no
          overheads).

        - ``secondlim`` (int, default 8) -- naive height bound on
          second point search on quartic homogeneous spaces (after
          testing local solubility; sieve-assisted search)

        - ``n_aux`` (int, default -1) -- If positive, the number of
          auxiliary primes used in sieve-assisted search for quartics.
          If -1 (the default) use a default value (set in the eclib
          code in ``src/qrank/mrank1.cc`` in DEFAULT_NAUX: currently 8).
          Only relevant for curves with no 2-torsion, where full
          2-descent is carried out.  Worth increasing for curves
          expected to be of rank > 6 to one or two more than the
          expected rank.

        - ``second_descent`` (int, default 1) -- flag specifying
          whether or not a second descent will be carried out (yes if
          1, the default; no if 0).  Only relevant for curves with
          2-torsion.  Recommended left as the default except for
          experts interested in details of Selmer groups.

        OUTPUT:

        None

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.getrank()
            3
            sage: D2.getcertain()
            1
            sage: D2.ok()
            1
        """
        sig_on()
        self.x = new two_descent(curve.x, verb, sel, firstlim, secondlim, n_aux, second_descent)
        sig_off()

    def getrank(self):
        """
        Returns the rank (after doing a 2-descent).

        OUTPUT:

        (Integer) the rank (or an upper bound).

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.getrank()
            3
        """
        cdef int r
        sig_on()
        r = two_descent_get_rank(self.x)
        sig_off()
        from sage.rings.all import Integer
        return Integer(r)

    def getrankbound(self):
        """
        Returns the rank upper bound (after doing a 2-descent).

        OUTPUT:

        (Integer) an upper bound on the rank.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.getrankbound()
            3
        """
        cdef int r
        sig_on()
        r = two_descent_get_rank_bound(self.x)
        sig_off()
        from sage.rings.all import Integer
        return Integer(r)

    def getselmer(self):
        """
        Returns the 2-Selmer rank (after doing a 2-descent).

        OUTPUT:

        (Integer) The 2-Selmer rank.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.getselmer()
            3
        """
        sig_on()
        r = two_descent_get_selmer_rank(self.x)
        sig_off()
        from sage.rings.all import Integer
        return Integer(r)

    def ok(self):
        """
        Returns the success flag (after doing a 2-descent).

        OUTPUT:

        (bool) Flag indicating whether or not 2-descent was successful.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.ok()
            1
        """
        return two_descent_ok(self.x)

    def getcertain(self):
        """
        Returns the certainty flag (after doing a 2-descent).

        OUTPUT:

        (bool) ``True`` if the rank upper and lower bounds are equal.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.getcertain()
            1
        """
        return two_descent_get_certain(self.x)

    def saturate(self, saturation_bound=0):
        """
        Carries out saturation of the points found by a 2-descent.

        OUTPUT:

        None.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.saturate()
            Searching for points (bound = 8)...done:
              found points which generate a subgroup of rank 3
              and regulator 0.417...
            Processing points found during 2-descent...done:
              now regulator = 0.417...
            No saturation being done
            sage: D2.getbasis()
            '[[1:-1:1], [-2:3:1], [-14:25:8]]'
        """
        sig_on()
        two_descent_saturate(self.x, saturation_bound)
        sig_off()

    def getbasis(self):
        """
        Returns the basis of points found by doing a 2-descent.

        If the success and certain flags are 1, this will be a
        `\ZZ/2\ZZ`-basis for `E(\QQ)/2E(\QQ)` (modulo torsion),
        otherwise possibly only for a proper subgroup.

        .. NOTE::

           You must call ``saturate()`` first, or a RunTimeError will be raised.

        OUTPUT:

        (string) String representation of the list of points after
        saturation.

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0
            sage: D2.saturate()
            Searching for points (bound = 8)...done:
              found points which generate a subgroup of rank 3
              and regulator 0.417...
            Processing points found during 2-descent...done:
              now regulator = 0.417...
            No saturation being done
            sage: D2.getbasis()
            '[[1:-1:1], [-2:3:1], [-14:25:8]]'
        """
        sig_on()
        return string_sigoff(two_descent_get_basis(self.x))

    def regulator(self):
        """
        Returns the regulator of the points found by doing a 2-descent.

        OUTPUT:

        (double) The regulator (of the subgroup found by 2-descent).

        EXAMPLES::

            sage: from sage.libs.eclib.mwrank import _Curvedata
            sage: CD = _Curvedata(0,0,1,-7,6)
            sage: from sage.libs.eclib.mwrank import _two_descent
            sage: D2 = _two_descent()
            sage: D2.do_descent(CD)
            Basic pair: I=336, J=-10800
            disc=35092224
            ...
            Mordell rank contribution from B=im(eps) = 3
            Selmer  rank contribution from B=im(eps) = 3
            Sha     rank contribution from B=im(eps) = 0
            Mordell rank contribution from A=ker(eps) = 0
            Selmer  rank contribution from A=ker(eps) = 0
            Sha     rank contribution from A=ker(eps) = 0

        If called before calling ``saturate()``, a bogus value of 1.0
        is returned::

            sage: D2.regulator()
            1.0

        After saturation, both ``getbasis()`` and ``regulator()``
        return the basis and regulator of the subgroup found by
        2-descent::

            sage: D2.saturate()
            Searching for points (bound = 8)...done:
              found points which generate a subgroup of rank 3
              and regulator 0.417...
            Processing points found during 2-descent...done:
              now regulator = 0.417...
            No saturation being done
            sage: D2.getbasis()
            '[[1:-1:1], [-2:3:1], [-14:25:8]]'
            sage: D2.regulator()
            0.417143558758384
        """
        sig_on()
        reg = two_descent_regulator(self.x)
        sig_off()
        return reg
