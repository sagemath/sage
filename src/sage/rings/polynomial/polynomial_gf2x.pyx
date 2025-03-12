# distutils: libraries = gmp NTL_LIBRARIES
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: extra_compile_args = NTL_CFLAGS M4RI_CFLAGS
# distutils: include_dirs = NTL_INCDIR M4RI_INCDIR
# distutils: language = c++
"""
Univariate Polynomials over GF(2) via NTL's GF2X

AUTHOR:
- Martin Albrecht (2008-10) initial implementation
"""


# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

cdef inline cparent get_cparent(parent) noexcept:
    return 0

# first we include the definitions
include "sage/libs/ntl/decl.pxi"
include "sage/libs/ntl/ntl_GF2X_linkage.pxi"

# and then the interface
include "polynomial_template.pxi"

from sage.libs.pari import pari

from sage.libs.m4ri cimport mzd_write_bit, mzd_read_bit
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense

from sage.misc.cachefunc import cached_method

cdef class Polynomial_GF2X(Polynomial_template):
    r"""
    Univariate Polynomials over `\GF{2}` via NTL's GF2X.

    EXAMPLES::

        sage: P.<x> = GF(2)[]
        sage: x^3 + x^2 + 1
        x^3 + x^2 + 1
    """
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        r"""
        Create a new univariate polynomials over `\GF{2}`.

        EXAMPLES::

            sage: P.<x> = GF(2)[]
            sage: x^3 + x^2 + 1
            x^3 + x^2 + 1

        We check that the bug noted at :issue:`12724` is fixed::

            sage: R.<x> = Zmod(2)[]
            sage: R([2^80])
            0
        """
        try:
            if (isinstance(x, int)
                or isinstance(x, Integer)):
                x = int(x % 2)
            elif (x.parent() is parent.base_ring()
                or x.parent() == parent.base_ring()):
                x = int(x)
        except AttributeError:
            pass
        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    cdef get_unsafe(self, Py_ssize_t i):
        """
        Return the `i`-th coefficient of ``self``.

        EXAMPLES::

            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1; f
            x^3 + x^2 + 1
            sage: f[0]
            1
            sage: f[1]
            0
            sage: f[:50] == f
            True
            sage: f[:3]
            x^2 + 1
        """
        cdef long c = GF2_conv_to_long(GF2X_coeff(self.x, i))
        return self._parent._base(c)

    def __pari__(self, variable=None):
        """
        EXAMPLES::

            sage: P.<x> = GF(2)[]
            sage: f = x^3 + x^2 + 1
            sage: pari(f)
            Mod(1, 2)*x^3 + Mod(1, 2)*x^2 + Mod(1, 2)
        """
        # TODO: put this in a superclass
        parent = self._parent
        if variable is None:
            variable = parent.variable_name()
        return pari(self.list()).Polrev(variable) * pari(1).Mod(2)

    def modular_composition(Polynomial_GF2X self, Polynomial_GF2X g, Polynomial_GF2X h, algorithm=None):
        r"""
        Compute `f(g) \pmod h`.

        Both implementations use Brent-Kung's Algorithm 2.1 (*Fast Algorithms
        for Manipulation of Formal Power Series*, JACM 1978).

        INPUT:

        - ``g`` -- a polynomial
        - ``h`` -- a polynomial
        - ``algorithm`` -- either ``'native'`` or ``'ntl'`` (default: ``'native'``)

        EXAMPLES::

            sage: P.<x> = GF(2)[]
            sage: r = 279
            sage: f = x^r + x +1
            sage: g = x^r
            sage: g.modular_composition(g, f) == g(g) % f
            True

            sage: P.<x> = GF(2)[]
            sage: f = x^29 + x^24 + x^22 + x^21 + x^20 + x^16 + x^15 + x^14 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^2
            sage: g = x^31 + x^30 + x^28 + x^26 + x^24 + x^21 + x^19 + x^18 + x^11 + x^10 + x^9 + x^8 + x^5 + x^2 + 1
            sage: h = x^30 + x^28 + x^26 + x^25 + x^24 + x^22 + x^21 + x^18 + x^17 + x^15 + x^13 + x^12 + x^11 + x^10 + x^9 + x^4
            sage: f.modular_composition(g, h) == f(g) % h
            True

        AUTHORS:

        - Paul Zimmermann (2008-10) initial implementation
        - Martin Albrecht (2008-10) performance improvements
        """
        if g.parent() is not self.parent() or h.parent() is not self.parent():
            raise TypeError("Parents of the first three parameters must match.")

        from sage.misc.timing import cputime
        from sage.misc.verbose import verbose
        from sage.arith.misc import integer_ceil as ceil
        from sage.matrix.constructor import Matrix
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

        cdef Polynomial_GF2X res
        cdef GF2XModulus_c modulus
        GF2XModulus_build(modulus, (<Polynomial_GF2X>h).x)

        res = <Polynomial_GF2X>Polynomial_GF2X.__new__(Polynomial_GF2X)
        res._parent = self._parent
        res._cparent = self._cparent

        if algorithm == "ntl":
            t = cputime()
            sig_on()
            GF2X_CompMod(res.x, self.x, g.x, modulus)
            sig_off()
            verbose("NTL %5.3f s"%cputime(t),level=1)
            return res

        cdef Py_ssize_t i, j, k, l, n, maxlength
        cdef Matrix_mod2_dense F, G, H

        if g.degree() >= h.degree():
            g = g % h

        cdef GF2X_c _f = (<Polynomial_GF2X>self).x
        cdef GF2X_c _g = (<Polynomial_GF2X>g).x
        cdef GF2X_c gpow, g2, tt
        GF2X_conv_long(gpow, 1)

        maxlength = GF2X_NumBits(_f)

        t = cputime()

        n = h.degree()

        k = ceil(Integer(n+1).sqrt(prec=Integer(n).log(2,prec=30)+1))
        l = ceil((self.degree() + 1) / k)

        # we store all matrices transposed for performance reasons
        G = <Matrix_mod2_dense>Matrix(GF(2), k, n)

        # first compute g^j mod h, 2 <= j < k
        # first deal with j=0
        for i from 0 <= i < GF2X_NumBits(gpow):
            mzd_write_bit(G._entries, 0, i, GF2_conv_to_long(GF2X_coeff(gpow, i)))
        # precompute g^2
        GF2X_SqrMod_pre(g2, _g, modulus)
        gpow = _g
        for j in range(1, k, 2):
            if j > 1:
                GF2X_MulMod_pre(gpow, gpow, g2, modulus) # gpow = g^j
            for i from 0 <= i < GF2X_NumBits(gpow):
                mzd_write_bit(G._entries, j, i, GF2_conv_to_long(GF2X_coeff(gpow, i)))
            # we now process 2j, 4j, 8j, ... by squaring each time
            if 2*j < k:
                tt = gpow
                jj = j
                while 2*jj < k:
                    GF2X_SqrMod_pre(tt, tt, modulus)
                    jj = 2*jj
                    for i from 0 <= i < GF2X_NumBits(tt):
                        mzd_write_bit(G._entries, jj, i, GF2_conv_to_long(GF2X_coeff(tt, i)))
        # we need that gpow = g^k at the end
        if k % 2 == 1: # k is odd, last j is k-2
            GF2X_MulMod_pre(gpow, gpow, g2, modulus)
        else:          # k is even, last j is k-1
            GF2X_MulMod_pre(gpow, gpow, _g, modulus)
        verbose("G %d x %d %5.3f s"%(G.nrows(), G.ncols(),cputime(t)),level=1)

        # split f in chunks of degree < k
        t = cputime()
        F = <Matrix_mod2_dense>Matrix(GF(2), l, k)
        for j in range(0, l):
            if j*k+k <= maxlength:
                for i from j*k <= i < j*k+k:
                    mzd_write_bit(F._entries, j, i-j*k, GF2_conv_to_long(GF2X_coeff(_f, i)))
            else:
                for i from j*k <= i < maxlength:
                    mzd_write_bit(F._entries, j, i-j*k, GF2_conv_to_long(GF2X_coeff(_f, i)))

        verbose("F %d x %d %5.3f s"%(F.nrows(), F.ncols(), cputime(t)),level=1)

        t = cputime()
        H = <Matrix_mod2_dense>(F * G)
        verbose("H %d x %d %5.3f s"%(H.nrows(), H.ncols(), cputime(t)),level=1)

        t = cputime()
        # H is a n x l matrix now H[i,j] = sum(G[i,m]*F[m,j],
        # m=0..k-1) = sum(g^m[i] * f[j*k+m], m=0..k-1) where g^m[i] is
        # the coefficient of degree i in g^m and f[j*k+m] is the
        # coefficient of degree j*k+m in f thus f[j*k+m]*g^m[i] should
        # be multiplied by g^(j*k) gpow = (g^k) % h

        GF2X_conv_long(res.x, 0)
        j = l - 1
        while j >= 0:
            #res = (res * gpow) % h
            GF2X_MulMod_pre(res.x, res.x, gpow, modulus)

            # res = res + parent([H[j,i] for i in range(0,n)])
            GF2X_conv_long(tt, 0)
            for i from 0<= i < n:
                GF2X_SetCoeff_long(tt, i, mzd_read_bit(H._entries, j, i))
            GF2X_add(res.x, res.x, tt)
            j = j - 1

        verbose("Res %5.3f s"%cputime(t),level=1)
        return res

    # Other polynomials have compose_mod as methods following the naming of
    # NTL/Flint bindings but the above method predates these. We expose
    # compose_mod here so all polynomial ring elements which support this can
    # use either name
    compose_mod = modular_composition

    @cached_method
    def is_irreducible(self):
        r"""
        Return whether this polynomial is irreducible over `\GF{2}`.

        EXAMPLES::

            sage: R.<x> = GF(2)[]
            sage: (x^2 + 1).is_irreducible()
            False
            sage: (x^3 + x + 1).is_irreducible()
            True

        Test that caching works::

            sage: R.<x> = GF(2)[]
            sage: f = x^2 + 1
            sage: f.is_irreducible()
            False
            sage: f.is_irreducible.cache
            False
        """
        return 0 != GF2X_IterIrredTest(self.x)


# The three functions below are used in polynomial_ring.py, but are in
# this Cython file since they call C++ functions.  They return
# polynomials as lists so that no variable has to be specified.
# AUTHOR: Peter Bruin (June 2013)

def GF2X_BuildIrred_list(n):
    """
    Return the list of coefficients of the lexicographically smallest
    irreducible polynomial of degree `n` over the field of 2 elements.

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list
        sage: GF2X_BuildIrred_list(2)
        [1, 1, 1]
        sage: GF2X_BuildIrred_list(3)
        [1, 1, 0, 1]
        sage: GF2X_BuildIrred_list(4)
        [1, 1, 0, 0, 1]
        sage: GF(2)['x'](GF2X_BuildIrred_list(33))
        x^33 + x^6 + x^3 + x + 1
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    cdef GF2X_c f
    GF2 = FiniteField(2)
    GF2X_BuildIrred(f, int(n))
    return [GF2(not GF2_IsZero(GF2X_coeff(f, i))) for i in range(n + 1)]


def GF2X_BuildSparseIrred_list(n):
    """
    Return the list of coefficients of an irreducible polynomial of
    degree `n` of minimal weight over the field of 2 elements.

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildIrred_list, GF2X_BuildSparseIrred_list
        sage: all([GF2X_BuildSparseIrred_list(n) == GF2X_BuildIrred_list(n)
        ....:      for n in range(1,33)])
        True
        sage: GF(2)['x'](GF2X_BuildSparseIrred_list(33))
        x^33 + x^10 + 1
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    cdef GF2X_c f
    GF2 = FiniteField(2)
    GF2X_BuildSparseIrred(f, int(n))
    return [GF2(not GF2_IsZero(GF2X_coeff(f, i))) for i in range(n + 1)]


def GF2X_BuildRandomIrred_list(n):
    """
    Return the list of coefficients of an irreducible polynomial of
    degree `n` of minimal weight over the field of 2 elements.

    EXAMPLES::

        sage: from sage.rings.polynomial.polynomial_gf2x import GF2X_BuildRandomIrred_list
        sage: GF2X_BuildRandomIrred_list(2)
        [1, 1, 1]
        sage: GF2X_BuildRandomIrred_list(3) in [[1, 1, 0, 1], [1, 0, 1, 1]]
        True
    """
    from sage.misc.randstate import current_randstate
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    cdef GF2X_c tmp, f
    GF2 = FiniteField(2)
    current_randstate().set_seed_ntl(False)
    GF2X_BuildSparseIrred(tmp, int(n))
    GF2X_BuildRandomIrred(f, tmp)
    return [GF2(not GF2_IsZero(GF2X_coeff(f, i))) for i in range(n + 1)]
