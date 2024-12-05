"""
libSingular: Functions

Sage implements a C wrapper around the Singular interpreter which
allows to call any function directly from Sage without string parsing
or interprocess communication overhead. Users who do not want to call
Singular functions directly, usually do not have to worry about this
interface, since it is handled by higher level functions in Sage.

EXAMPLES:

The direct approach for loading a Singular function is to call the
function :func:`singular_function` with the function name as
parameter::

    sage: from sage.libs.singular.function import singular_function
    sage: P.<a,b,c,d> = PolynomialRing(GF(7))
    sage: std = singular_function('std')
    sage: I = sage.rings.ideal.Cyclic(P)
    sage: std(I)
    [a + b + c + d,
     b^2 + 2*b*d + d^2,
     b*c^2 + c^2*d - b*d^2 - d^3,
     b*c*d^2 + c^2*d^2 - b*d^3 + c*d^3 - d^4 - 1,
     b*d^4 + d^5 - b - d,
     c^3*d^2 + c^2*d^3 - c - d,
     c^2*d^4 + b*c - b*d + c*d - 2*d^2]

If a Singular library needs to be loaded before a certain function is
available, use the :func:`lib` function as shown below::

    sage: from sage.libs.singular.function import singular_function, lib as singular_lib
    sage: primdecSY = singular_function('primdecSY')
    Traceback (most recent call last):
    ...
    NameError: Singular library function 'primdecSY' is not defined

    sage: singular_lib('primdec.lib')
    sage: primdecSY = singular_function('primdecSY')

There is also a short-hand notation for the above::

    sage: import sage.libs.singular.function_factory
    sage: primdecSY = sage.libs.singular.function_factory.ff.primdec__lib.primdecSY

The above line will load "primdec.lib" first and then load the
function ``primdecSY``.

TESTS::

    sage: from sage.libs.singular.function import singular_function
    sage: std = singular_function('std')
    sage: loads(dumps(std)) == std
    True

AUTHORS:

- Michael Brickenstein (2009-07): initial implementation, overall design
- Martin Albrecht (2009-07): clean up, enhancements, etc.
- Michael Brickenstein (2009-10): extension to more Singular types
- Martin Albrecht (2010-01): clean up, support for attributes
- Simon King (2011-04): include the documentation provided by Singular as a code block
- Burcin Erocal, Michael Brickenstein, Oleksandr Motsak, Alexander Dreyer, Simon King (2011-09): plural support
"""

#*****************************************************************************
#       Copyright (C) 2009 Michael Brickenstein <brickenstein@mfo.de>
#       Copyright (C) 2009,2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.string cimport memcpy
from cysignals.signals cimport sig_on, sig_off

from sage.cpython.string cimport str_to_bytes, char_to_str

from sage.structure.sage_object cimport SageObject
from sage.structure.richcmp cimport richcmp
from sage.structure.sequence import Sequence, Sequence_generic

from sage.modules.free_module_element cimport FreeModuleElement_generic_dense

from sage.rings.integer cimport Integer
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, new_MP
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomialRing_libsingular
from sage.rings.polynomial.plural cimport NCPolynomialRing_plural, NCPolynomial_plural, new_NCP
from sage.rings.polynomial.multi_polynomial_ideal import NCPolynomialIdeal
from sage.rings.polynomial.multi_polynomial_ideal import MPolynomialIdeal
from sage.rings.polynomial.multi_polynomial_ideal_libsingular cimport sage_ideal_to_singular_ideal, singular_ideal_to_sage_sequence
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence_generic

from sage.libs.singular.decl cimport *
from sage.libs.singular.option import opt_ctx
from sage.libs.singular.polynomial cimport singular_vector_maximal_component
from sage.libs.singular.singular cimport sa2si, si2sa, si2sa_intvec, si2sa_bigintvec, start_catch_error, check_error
from sage.libs.singular.singular import error_messages

from sage.interfaces.singular import get_docstring

from sage.misc.verbose import get_verbose


cdef poly* sage_vector_to_poly(v, ring *r) except <poly*> -1:
    """
    Convert a vector or list of multivariate polynomials to a
    polynomial by adding them all up.
    """
    cdef poly *res = NULL
    cdef poly *poly_component
    cdef poly *p_iter
    cdef int component

    for (i, p) in enumerate(v):
        component = <int>i+1
        poly_component = copy_sage_polynomial_into_singular_poly(p)
        p_iter = poly_component
        while p_iter!=NULL:
            p_SetComp(p_iter, component, r)
            p_Setm(p_iter, r)
            p_iter=pNext(p_iter)
        res=p_Add_q(res, poly_component, r)
    return res


cdef class RingWrap:
    """
    A simple wrapper around Singular's rings.
    """
    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P)
            <RingWrap>
        """
        if not self.is_commutative():
            return "<noncommutative RingWrap>"
        return "<RingWrap>"

    def __dealloc__(self):
        if self._ring != NULL:
            self._ring.ref -= 1

    def ngens(self):
        """
        Get number of generators.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).ngens()
            3
        """
        return self._ring.N

    def var_names(self):
        """
        Get names of variables.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).var_names()
            ['x', 'y', 'z']
        """
        return [char_to_str(self._ring.names[i]) for i in range(self.ngens())]

    def npars(self):
        """
        Get number of parameters.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).npars()
            0
        """
        return n_NumberOfParameters(self._ring.cf)

    def ordering_string(self):
        """
        Get Singular string defining monomial ordering.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).ordering_string()
            'dp(3),C'
        """
        return char_to_str(rOrderingString(self._ring))

    def par_names(self):
        """
        Get parameter names.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).par_names()
            []
        """
        return [char_to_str(n_ParameterNames(self._ring.cf)[i])
                for i in range(self.npars())]

    def characteristic(self):
        """
        Get characteristic.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).characteristic()
            0
        """
        return self._ring.cf.ch

    def is_commutative(self):
        """
        Determine whether a given ring is commutative.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P).is_commutative()
            True
        """
        return not rIsPluralRing(self._ring)

    def _output(self):
        """
        Use Singular output.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: ringlist = singular_function("ringlist")
            sage: l = ringlist(P)
            sage: ring = singular_function("ring")
            sage: ring(l, ring=P)._output()
            //   coefficients: QQ
            //   number of vars : 3
            //        block   1 : ordering dp
            //                  : names    x y z
            //        block   2 : ordering C
        """
        rPrint(self._ring)


cdef class Resolution:
    """
    A simple wrapper around Singular's resolutions.
    """
    def __init__(self, base_ring):
        """
        EXAMPLES::

           sage: from sage.libs.singular.function import singular_function
           sage: mres = singular_function("mres")
           sage: syz = singular_function("syz")
           sage: P.<x,y,z> = PolynomialRing(QQ)
           sage: I = P.ideal([x+y,x*y-y, y*2,x**2+1])
           sage: M = syz(I)
           sage: resolution = mres(M, 0)
        """
        # FIXME: still not working noncommutative
        assert is_sage_wrapper_for_singular_ring(base_ring)
        self.base_ring = base_ring

    def __repr__(self):
        """
        EXAMPLES::

           sage: from sage.libs.singular.function import singular_function
           sage: mres = singular_function("mres")
           sage: syz = singular_function("syz")
           sage: P.<x,y,z> = PolynomialRing(QQ)
           sage: I = P.ideal([x+y,x*y-y, y*2,x**2+1])
           sage: M = syz(I)
           sage: resolution = mres(M, 0)
           sage: resolution
           <Resolution>
        """
        return "<Resolution>"

    def __dealloc__(self):
        """
        EXAMPLES::

           sage: from sage.libs.singular.function import singular_function
           sage: mres = singular_function("mres")
           sage: syz = singular_function("syz")
           sage: P.<x,y,z> = PolynomialRing(QQ)
           sage: I = P.ideal([x+y,x*y-y, y*2,x**2+1])
           sage: M = syz(I)
           sage: resolution = mres(M, 0)
           sage: del resolution
        """
        if self._resolution != NULL:
            self._resolution.references -= 1


cdef leftv* new_leftv(void *data, res_type) noexcept:
    """
    INPUT:

    - ``data`` -- some Singular data this interpreter object points to
    - ``res_type`` -- the type of that data
    """
    cdef leftv* res
    res = <leftv*>omAllocBin(sleftv_bin)
    res.Init()
    res.data = data
    res.rtyp = res_type
    return res

cdef free_leftv(leftv *args, ring *r=NULL):
    """
    Kills this ``leftv`` and all ``leftv``s in the tail.

    INPUT:

    - ``args`` -- list of Singular arguments
    """
    args.CleanUp(r)
    omFreeBin(args, sleftv_bin)


# =====================================
# = Singular/Plural Abstraction Layer =
# =====================================

def is_sage_wrapper_for_singular_ring(ring):
    """
    Check whether wrapped ring arises from Singular or Singular/Plural.

    EXAMPLES::

        sage: from sage.libs.singular.function import is_sage_wrapper_for_singular_ring
        sage: P.<x,y,z> = QQ[]
        sage: is_sage_wrapper_for_singular_ring(P)
        True

    ::

        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: P = A.g_algebra(relations={y*x:-x*y}, order = 'lex')
        sage: is_sage_wrapper_for_singular_ring(P)
        True
    """
    if isinstance(ring, MPolynomialRing_libsingular):
        return True
    if isinstance(ring, NCPolynomialRing_plural):
        return True
    return False


cdef new_sage_polynomial(ring,  poly *p):
    if isinstance(ring, MPolynomialRing_libsingular):
        return new_MP(ring, p)
    else:
        if isinstance(ring, NCPolynomialRing_plural):
            return new_NCP(ring, p)
    raise ValueError("not a singular or plural ring")


def is_singular_poly_wrapper(p):
    """
    Check if ``p`` is some data type corresponding to some singular ``poly``.

    EXAMPLES::

        sage: from sage.libs.singular.function import is_singular_poly_wrapper
        sage: A.<x,y,z> = FreeAlgebra(QQ, 3)
        sage: H.<x,y,z> = A.g_algebra({z*x:x*z+2*x, z*y:y*z-2*y})
        sage: is_singular_poly_wrapper(x+y)
        True
    """
    return isinstance(p, (MPolynomial_libsingular, NCPolynomial_plural))


def all_singular_poly_wrapper(s):
    """
    Test for a sequence ``s``, whether it consists of
    singular polynomials.

    EXAMPLES::

        sage: from sage.libs.singular.function import all_singular_poly_wrapper
        sage: P.<x,y,z> = QQ[]
        sage: all_singular_poly_wrapper([x+1, y])
        True
        sage: all_singular_poly_wrapper([x+1, y, 1])
        False
    """
    return all(is_singular_poly_wrapper(p) for p in s)


cdef poly* access_singular_poly(p) except <poly*> -1:
    """
    Get the raw ``poly`` pointer from a wrapper object.
    """
    if isinstance(p, MPolynomial_libsingular):
        return (<MPolynomial_libsingular> p)._poly
    else:
        if isinstance(p, NCPolynomial_plural):
            return (<NCPolynomial_plural> p)._poly
    raise ValueError("not a singular polynomial wrapper")

cdef ring* access_singular_ring(r) except <ring*> -1:
    """
    Get the singular ``ring`` pointer from a wrapper object.
    """
    if isinstance(r, MPolynomialRing_libsingular):
        return (<MPolynomialRing_libsingular> r )._ring
    if isinstance(r, NCPolynomialRing_plural):
        return (<NCPolynomialRing_plural> r )._ring
    raise ValueError("not a singular polynomial ring wrapper")

cdef poly* copy_sage_polynomial_into_singular_poly(p) noexcept:
    return p_Copy(access_singular_poly(p), access_singular_ring(p.parent()))


def all_vectors(s):
    """
    Check if a sequence ``s`` consists of free module elements
    over a singular ring.

    EXAMPLES::

        sage: from sage.libs.singular.function import all_vectors
        sage: P.<x,y,z> = QQ[]
        sage: M = P**2
        sage: all_vectors([x])
        False
        sage: all_vectors([(x,y)])
        False
        sage: all_vectors([M(0), M((x,y))])
        True
        sage: all_vectors([M(0), M((x,y)),(0,0)])
        False
    """
    return all(isinstance(p, FreeModuleElement_generic_dense)
               and is_sage_wrapper_for_singular_ring(p.parent().base_ring())
               for p in s)


cdef class Converter(SageObject):
    """
    A :class:`Converter` interfaces between Sage objects and Singular
    interpreter objects.
    """

    def __init__(self, args, ring, attributes=None):
        """
        Create a new argument list.

        INPUT:

        - ``args`` -- list of Python objects
        - ``ring`` -- a multivariate polynomial ring
        - ``attributes`` -- an optional dictionary of Singular
          attributes (default: ``None``)

        EXAMPLES::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: Converter([a,b,c],ring=P)
            Singular Converter in Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
        """
        cdef leftv *v
        self.args = NULL
        self._sage_ring = ring
        if ring is not None:
            self._singular_ring = access_singular_ring(ring)

        from sage.matrix.matrix_mpolynomial_dense import Matrix_mpolynomial_dense
        from sage.matrix.matrix_integer_dense import Matrix_integer_dense
        from sage.matrix.matrix_generic_dense import Matrix_generic_dense
        for a in args:
            if is_singular_poly_wrapper(a):
                v = self.append_polynomial(a)

            elif is_sage_wrapper_for_singular_ring(a):
                v = self.append_ring(a)

            elif isinstance(a, MPolynomialIdeal) or \
                    isinstance(a, NCPolynomialIdeal):
                v = self.append_ideal(a)

            elif isinstance(a, int):
                v = self.append_int(a)

            elif isinstance(a, str):
                v = self.append_str(a)

            elif isinstance(a, Matrix_mpolynomial_dense):
                v = self.append_matrix(a)

            elif isinstance(a, Matrix_integer_dense):
                v = self.append_intmat(a)

            elif isinstance(a, Matrix_generic_dense) and\
                is_sage_wrapper_for_singular_ring(a.parent().base_ring()):
                v = self.append_matrix(a)

            elif isinstance(a, Resolution):
                v = self.append_resolution(a)

            elif isinstance(a, FreeModuleElement_generic_dense)\
                and is_sage_wrapper_for_singular_ring(
                    a.parent().base_ring()):
                v = self.append_vector(a)

            # as output ideals get converted to sequences
            # sequences of polynomials should get converted to ideals
            # this means, that Singular lists should not be converted to Sequences,
            # as we do not want ambiguities
            elif isinstance(a, Sequence_generic)\
                and all_singular_poly_wrapper(a):
                v = self.append_ideal(ring.ideal(a))
            elif isinstance(a, PolynomialSequence_generic):
                v = self.append_ideal(ring.ideal(a))
            elif isinstance(a, Sequence_generic)\
                and all_vectors(a):
                v = self.append_module(a)
            elif isinstance(a, list):
                v = self.append_list(a)

            elif isinstance(a, tuple):
                is_intvec = True
                for i in a:
                    if not (isinstance(i, int)
                        or isinstance(i, Integer)):
                        is_intvec = False
                        break
                if is_intvec:
                    v = self.append_intvec(a)
                else:
                    v = self.append_list(a)
            elif a.parent() is self._sage_ring.base_ring():
                v = self.append_number(a)

            elif isinstance(a, Integer):
                v = self.append_int(a)

            else:
                raise TypeError("unknown argument type '%s'"%(type(a),))

            if attributes and a in attributes:
                for attrib in attributes[a]:
                    if attrib == "isSB":
                        val = <long>(attributes[a][attrib])
                        atSet(v, omStrDup("isSB"), <void*>val, INT_CMD)
                        setFlag(v, FLAG_STD)
                    else:
                        raise NotImplementedError("Support for attribute '%s' not implemented yet."%attrib)

    def ring(self):
        """
        Return the ring in which the arguments of this list live.

        EXAMPLES::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: Converter([a,b,c],ring=P).ring()
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
        """
        return self._sage_ring

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: Converter([a,b,c],ring=P) # indirect doctest
            Singular Converter in Multivariate Polynomial Ring in a, b, c over Finite Field of size 127
        """
        return "Singular Converter in %s"%(self._sage_ring)

    def __dealloc__(self):
        cdef ring *r = access_singular_ring(self._sage_ring)
        if self.args:
            free_leftv(self.args, r)

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import Converter
            sage: P.<a,b,c> = PolynomialRing(GF(127))
            sage: len(Converter([a,b,c],ring=P))
            3
        """
        cdef leftv * v
        v=self.args
        cdef int l
        l=0
        while v != NULL:
            l=l+1
            v=v.next
        return l

    cdef leftv* pop_front(self) except NULL:
        """
        Pop a Singular element from the front of the list.
        """
        assert(self.args != NULL)
        cdef leftv *res = self.args
        self.args = self.args.next
        res.next = NULL
        return res

    cdef leftv *_append_leftv(self, leftv *v) noexcept:
        """
        Append a new Singular element to the list.
        """
        cdef leftv* last
        if not self.args == NULL:
            last = self.args
            while not last.next == NULL:
                last=last.next
            last.next=v
        else:
            self.args = v
        return v

    cdef leftv *_append(self, void* data, int res_type) noexcept:
        """
        Create a new ``leftv`` and append it to the list.

        INPUT:

        - ``data`` -- the raw data
        - ``res_type`` -- the type of the data
        """
        return self._append_leftv( new_leftv(data, res_type) )

    cdef to_sage_matrix(self, matrix* mat):
        """
        Convert singular matrix to matrix over the polynomial ring.
        """
        from sage.matrix.constructor import Matrix
        ncols = mat.ncols
        nrows = mat.nrows
        result = Matrix(self._sage_ring, nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                p = new_sage_polynomial(self._sage_ring, mat.m[i*ncols+j])
                mat.m[i*ncols+j]=NULL
                result[i,j] = p
        return result

    cdef to_sage_vector_destructive(self, poly *p, free_module=None):
        cdef int rank
        if free_module:
            rank = free_module.rank()
        else:
            rank = singular_vector_maximal_component(p, self._singular_ring)
            free_module = self._sage_ring**rank
        cdef poly *acc
        cdef poly *p_iter
        cdef poly *first
        cdef poly *previous
        cdef int i
        result = []
        for i from 1 <= i <= rank:
            previous = NULL
            acc = NULL
            first = NULL
            p_iter = p
            while p_iter != NULL:
                if p_GetComp(p_iter, self._singular_ring) == i:
                    p_SetComp(p_iter, 0, self._singular_ring)
                    p_Setm(p_iter, self._singular_ring)
                    if acc == NULL:
                        first = p_iter
                    else:
                        acc.next = p_iter
                    acc = p_iter
                    if p_iter == p:
                        p = pNext(p_iter)
                    if previous != NULL:
                        previous.next = pNext(p_iter)
                    p_iter = pNext(p_iter)
                    acc.next = NULL
                else:
                    previous = p_iter
                    p_iter = pNext(p_iter)

            result.append(new_sage_polynomial(self._sage_ring, first))
        return free_module(result)

    cdef object to_sage_module_element_sequence_destructive( self, ideal *i):
        """
        Convert a SINGULAR module to a Sage Sequence (the format Sage
        stores a Groebner basis in).

        INPUT:

        - ``i`` -- a SINGULAR ideal
        - ``r`` -- a SINGULAR ring
        - ``sage_ring`` -- a Sage ring matching r
        """
        cdef int j
        cdef int rank=i.rank
        free_module = self._sage_ring ** rank
        l = []

        for j from 0 <= j < IDELEMS(i):
            p = self.to_sage_vector_destructive(i.m[j], free_module)
            i.m[j] = NULL  # save it from getting freed
            l.append( p )

        return Sequence(l, check=False, immutable=True)

    cdef to_sage_integer_matrix(self, intvec* mat):
        """
        Convert Singular matrix to matrix over the polynomial ring.
        """
        from sage.matrix.constructor import Matrix
        from sage.rings.integer_ring import ZZ

        ncols = mat.cols()
        nrows = mat.rows()

        result = Matrix(ZZ, nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                result[i,j] = mat.get(i*ncols+j)
        return result

    cdef leftv *append_polynomial(self, p) except NULL:
        """
        Append the polynomial ``p`` to the list.
        """
        cdef poly* _p
        _p = copy_sage_polynomial_into_singular_poly(p)

        return self._append(_p, POLY_CMD)

    cdef leftv *append_ideal(self,  i) except NULL:
        """
        Append the ideal ``i`` to the list.
        """
        cdef ideal* singular_ideal = sage_ideal_to_singular_ideal(i)
        return self._append(singular_ideal, IDEAL_CMD)

    cdef leftv *append_module(self, m) except NULL:
        """
        Append sequence ``m`` of vectors over the polynomial ring to
        the list
        """
        rank = max([v.parent().rank() for v in m])
        cdef ideal *result
        cdef ring *r = self._singular_ring
        cdef ideal *i
        cdef int j = 0

        i = idInit(len(m),rank)
        for f in m:
            i.m[j] = sage_vector_to_poly(f, r)
            j+=1
        return self._append(<void*> i, MODUL_CMD)

    cdef leftv *append_number(self, n) except NULL:
        """
        Append the number ``n`` to the list.
        """
        cdef number *_n =  sa2si(n, self._singular_ring)
        return self._append(<void *>_n, NUMBER_CMD)

    cdef leftv *append_ring(self, r) except NULL:
        """
        Append the ring ``r`` to the list.
        """
        cdef ring *_r =  access_singular_ring(r)
        _r.ref+=1
        return self._append(<void *>_r, RING_CMD)

    cdef leftv *append_matrix(self, mat) except NULL:
        sage_ring = mat.base_ring()
        cdef ring *r=<ring*> access_singular_ring(sage_ring)

        cdef poly *p
        ncols = mat.ncols()
        nrows = mat.nrows()
        cdef matrix* _m=mpNew(nrows, ncols)
        for i in range(nrows):
            for j in range(ncols):
                #FIXME
                p = copy_sage_polynomial_into_singular_poly(mat[i,j])
                _m.m[ncols*i+j]=p
        return self._append(_m, MATRIX_CMD)

    cdef leftv *append_int(self, n) except NULL:
        """
        Append the integer ``n`` to the list.
        """
        cdef long _n =  n
        return self._append(<void*>_n, INT_CMD)

    cdef leftv *append_list(self, l) except NULL:
        """
        Append the list ``l`` to the list.
        """

        cdef Converter c = Converter(l, self._sage_ring)
        n = len(c)

        cdef lists *singular_list=<lists*>omAlloc0Bin(slists_bin)
        singular_list.Init(n)
        cdef leftv* iv
        for i in range(n):
            iv=c.pop_front()
            memcpy(&singular_list.m[i],iv,sizeof(leftv))
            omFreeBin(iv, sleftv_bin)

        return self._append(<void*>singular_list, LIST_CMD)

    cdef leftv *append_intvec(self, a) except NULL:
        """
        Append ``a`` to the list as intvec.
        """
        s = len(a)
        cdef intvec *iv = new intvec()
        iv.resize(s)

        for i in range(s):
            iv.ivGetVec()[i]=<int>a[i]
        return self._append(<void*>iv, INTVEC_CMD)

    cdef leftv *append_vector(self, v) except NULL:
        """
        Append vector ``v`` from free
        module over polynomial ring.
        """
        cdef ring *r = self._singular_ring
        cdef poly *p = sage_vector_to_poly(v, r)
        return self._append(<void*> p, VECTOR_CMD)

    cdef leftv *append_resolution(self, Resolution resolution) except NULL:
        """
        Append free resolution ``r`` to the list.
        """
        resolution._resolution.references += 1
        return self._append(<void*> resolution._resolution, RESOLUTION_CMD)

    cdef leftv *append_intmat(self, a) except NULL:
        """
        Append ``a`` to the list as intvec.
        """
        cdef int nrows = <int> a.nrows()
        cdef int ncols = <int> a.ncols()
        cdef intvec *iv = new intvec(nrows, ncols, 0)

        for i in range(nrows):
            for j in range(ncols):
                iv.ivGetVec()[i*ncols+j]=<int>a[i,j]
        return self._append(<void*>iv, INTMAT_CMD)

    cdef leftv *append_str(self, n) except NULL:
        """
        Append the string ``n`` to the list.
        """
        b = str_to_bytes(n)
        return self._append(omStrDup(b), STRING_CMD)

    cdef to_python(self, leftv* to_convert):
        """
        Convert the ``leftv`` to a Python object.

        INPUT:

        - ``to_convert`` -- a Singular ``leftv``

        TESTS:

        Check that negative integers come through unscathed::

            sage: P.<x,y,z> = QQ[]
            sage: C = Curve((x-y)*(y-z)*(z-x))
            sage: I = C.defining_ideal()
            sage: import sage.libs.singular.function_factory
            sage: freerank = sage.libs.singular.function_factory.ff.polylib__lib.freerank
            sage: freerank(I, true)
            [-1, [x^2*y - x*y^2 - x^2*z + y^2*z + x*z^2 - y*z^2]]

        Singular's genus function is prone to crashing, see :issue:`12851` and :issue:`19750` ::

            sage: sing_genus = sage.libs.singular.function_factory.ff.normal__lib.genus  # known bug
            sage: sing_genus(I)  # known bug
            -2
        """
        # FIXME
        cdef MPolynomial_libsingular res_poly
        cdef int rtyp = to_convert.rtyp
        cdef lists *singular_list
        cdef Resolution res_resolution

        if rtyp == IDEAL_CMD:
            return singular_ideal_to_sage_sequence(<ideal*>to_convert.data, self._singular_ring, self._sage_ring)
        elif rtyp == POLY_CMD:
            # FIXME
            res_poly = MPolynomial_libsingular(self._sage_ring)
            res_poly._poly = <poly*>to_convert.data
            to_convert.data = NULL  # prevent it getting free, when cleaning the leftv
            return res_poly
        elif rtyp == INT_CMD:
            return <int><long>to_convert.data
        elif rtyp == NUMBER_CMD:
            return si2sa(<number *>to_convert.data, self._singular_ring, self._sage_ring.base_ring())
        elif rtyp == INTVEC_CMD:
            return si2sa_intvec(<intvec *> to_convert.data)
        elif rtyp == BIGINTVEC_CMD:
            return si2sa_bigintvec(<bigintmat *> to_convert.data)
        elif rtyp == STRING_CMD:
            # TODO: Need to determine what kind of data can be returned by a
            # STRING_CMD--is it just ASCII strings or can it be an arbitrary
            # binary?
            ret = char_to_str(<char *>to_convert.data)
            return ret
        elif rtyp == VECTOR_CMD:
            result = self.to_sage_vector_destructive(<poly *> to_convert.data)
            to_convert.data = NULL
            return result
        elif rtyp == RING_CMD or rtyp==QRING_CMD:
            return new_RingWrap( <ring*> to_convert.data )
        elif rtyp == MATRIX_CMD:
            return self.to_sage_matrix(<matrix*> to_convert.data )
        elif rtyp == LIST_CMD:
            singular_list = <lists*> to_convert.data
            ret = []
            for i in range(singular_list.nr+1):
                ret.append(self.to_python(&(singular_list.m[i])))
            return ret
        elif rtyp == MODUL_CMD:
            return self.to_sage_module_element_sequence_destructive(<ideal*> to_convert.data)
        elif rtyp == INTMAT_CMD:
            return self.to_sage_integer_matrix(<intvec*> to_convert.data)
        elif rtyp == RESOLUTION_CMD:
            res_resolution = Resolution(self._sage_ring)
            res_resolution._resolution = <syStrategy *> to_convert.data
            res_resolution._resolution.references += 1
            return res_resolution
        elif rtyp == NONE:
            return None
        else:
            raise NotImplementedError("rtyp %d not implemented."%(rtyp))


cdef class BaseCallHandler:
    """
    A call handler is an abstraction which hides the details of the
    implementation differences between kernel and library functions.
    """
    cdef leftv* handle_call(self, Converter argument_list, ring *_ring=NULL) noexcept:
        """
        Actual function call.
        """
        return NULL

    cdef bint free_res(self) noexcept:
        """
        Do we need to free the result object.
        """
        return False


cdef class LibraryCallHandler(BaseCallHandler):
    """
    A call handler is an abstraction which hides the details of the
    implementation differences between kernel and library functions.

    This class implements calling a library function.

    .. NOTE::

        Do not construct this class directly, use
        :func:`singular_function` instead.
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import LibraryCallHandler
            sage: LibraryCallHandler()
            <sage.libs.singular.function.LibraryCallHandler object at 0x...>
        """
        super().__init__()

    cdef leftv* handle_call(self, Converter argument_list, ring *_ring=NULL) noexcept:
        if _ring != currRing: rChangeCurrRing(_ring)
        cdef bint error = iiMake_proc(self.proc_idhdl, NULL, argument_list.args)
        cdef leftv * res
        if not error:
            res = <leftv*> omAllocBin(sleftv_bin)
            res.Init()
            res.Copy(&iiRETURNEXPR)
            iiRETURNEXPR.Init()
            return res
        raise RuntimeError("Error raised calling singular function")

    cdef bint free_res(self) noexcept:
        """
        We do not need to free the result object for library
        functions.
        """
        return False

# mapping int --> string for function arity
arity_dict = {
        CMD_1: "CMD_1",
        CMD_2: "CMD_2",
        CMD_3: "CMD_3",
        CMD_12: "CMD_12",
        CMD_13: "CMD_13",
        CMD_23: "CMD_23",
        CMD_123: "CMD_123",
        CMD_M: "CMD_M"
}

cdef class KernelCallHandler(BaseCallHandler):
    """
    A call handler is an abstraction which hides the details of the
    implementation differences between kernel and library functions.

    This class implements calling a kernel function.

    .. NOTE::

        Do not construct this class directly, use
        :func:`singular_function` instead.
    """
    def __init__(self, cmd_n, arity):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import KernelCallHandler
            sage: KernelCallHandler(0,0)
            <sage.libs.singular.function.KernelCallHandler object at 0x...>
        """
        super().__init__()
        self.cmd_n = cmd_n
        self.arity = arity

    cdef leftv* handle_call(self, Converter argument_list, ring *_ring=NULL) noexcept:
        cdef leftv * res
        res = <leftv*> omAllocBin(sleftv_bin)
        res.Init()
        cdef leftv *arg1
        cdef leftv *arg2
        cdef leftv *arg3

        cdef Py_ssize_t number_of_arguments = len(argument_list)

        # Handle functions with an arbitrary number of arguments, sent
        # by an argument list.
        if self.arity in [CMD_M, ROOT_DECL_LIST, RING_DECL_LIST]:
            if _ring != currRing: rChangeCurrRing(_ring)
            iiExprArithM(res, argument_list.args, self.cmd_n)
            return res

        if number_of_arguments == 1:
            if self.arity in [CMD_1, CMD_12, CMD_13, CMD_123, RING_CMD]:
                arg1 = argument_list.pop_front()
                if _ring != currRing: rChangeCurrRing(_ring)
                iiExprArith1(res, arg1, self.cmd_n)
                free_leftv(arg1, _ring)
                return res

        elif number_of_arguments == 2:
            if self.arity in [CMD_2, CMD_12, CMD_23, CMD_123]:
                arg1 = argument_list.pop_front()
                arg2 = argument_list.pop_front()
                if _ring != currRing: rChangeCurrRing(_ring)
                iiExprArith2(res, arg1, self.cmd_n, arg2, True)
                free_leftv(arg1, _ring)
                free_leftv(arg2, _ring)
                return res

        elif number_of_arguments == 3:
            if self.arity in [CMD_3, CMD_13, CMD_23, CMD_123, RING_CMD]:
                arg1 = argument_list.pop_front()
                arg2 = argument_list.pop_front()
                arg3 = argument_list.pop_front()
                if _ring != currRing: rChangeCurrRing(_ring)
                iiExprArith3(res, self.cmd_n, arg1, arg2, arg3)
                free_leftv(arg1, _ring)
                free_leftv(arg2, _ring)
                free_leftv(arg3, _ring)
                return res

        global errorreported
        global error_messages

        errorreported += 1
        error_messages.append(
                "Wrong number of arguments (got {} arguments, arity is {})"
                .format(number_of_arguments,
                        arity_dict.get(self.arity) or self.arity))
        return NULL

    cdef bint free_res(self) noexcept:
        """
        We need to free the result object for kernel functions.
        """
        return True


# The Sage ring used as a dummy for singular function calls.
cdef object dummy_ring


cdef class SingularFunction(SageObject):
    """
    The base class for Singular functions either from the kernel or
    from the library.
    """
    def __init__(self, name):
        """
        INPUT:

        - ``name`` -- the name of the function

        EXAMPLES::

            sage: from sage.libs.singular.function import SingularFunction
            sage: SingularFunction('foobar')
            foobar (singular function)
        """
        self._name = name
        global currRingHdl
        if currRingHdl == NULL:
            currRingHdl = ggetid("my_awesome_sage_ring")
            if currRingHdl == NULL:
                currRingHdl = enterid("my_awesome_sage_ring", 0, RING_CMD, &IDROOT, 1)
                currRingHdl.data.uring = <ring *>omAlloc0Bin(sip_sring_bin)
            currRingHdl.data.uring.ref += 1

    cdef BaseCallHandler get_call_handler(self):
        """
        Return a call handler which does the actual work.
        """
        raise NotImplementedError

    cdef bint function_exists(self) noexcept:
        """
        Return ``True`` if the function exists in this interface.
        """
        raise NotImplementedError

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import SingularFunction
            sage: SingularFunction('foobar') # indirect doctest
            foobar (singular function)
        """
        return "%s (singular function)" %(self._name)

    def __call__(self, *args, ring=None, bint interruptible=True, attributes=None):
        """
        Call this function with the provided arguments ``args`` in the
        ``ring``.

        INPUT:

        - ``args`` -- list of arguments
        - ``ring`` -- a multivariate polynomial ring
        - ``interruptible`` -- if ``True`` pressing :kbd:`Ctrl` + :kbd:`C`
          during the execution of this function will interrupt the computation
          (default: ``True``)

        - ``attributes`` -- dictionary of optional Singular
          attributes assigned to Singular objects (default: ``None``)

        If ``ring`` is not specified, it is guessed from the given arguments.
        If this is not possible, then a dummy ring, univariate polynomial ring
        over ``QQ``, is used.

        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: size = singular_function('size')
            sage: P.<a,b,c> = PolynomialRing(QQ)
            sage: size(a, ring=P)
            1
            sage: size(2r,ring=P)
            1
            sage: size(2, ring=P)
            1
            sage: size(2)
            1
            sage: size(Ideal([a*b + c, a + 1]))
            2
            sage: size(Ideal([a*b + c, a + 1]))
            2
            sage: size(1,2)
            Traceback (most recent call last):
            ...
            RuntimeError: error in Singular function call 'size':
            Wrong number of arguments (got 2 arguments, arity is CMD_1)
            sage: size('foobar', ring=P)
            6

        Show the usage of the optional ``attributes`` parameter::

            sage: P.<x,y,z> = PolynomialRing(QQ)
            sage: I = Ideal([x^3*y^2 + 3*x^2*y^2*z + y^3*z^2 + z^5])
            sage: I = Ideal(I.groebner_basis())
            sage: hilb = sage.libs.singular.function_factory.ff.hilb
            sage: from sage.misc.sage_ostools import redirection
            sage: out = tmp_filename()
            sage: with redirection(sys.stdout,  open(out, 'w')):
            ....:     hilb(I) # Singular will print // ** _ is no standard basis
            sage: with open(out) as f:
            ....:     'is no standard basis' in f.read()
            True

        So we tell Singular that ``I`` is indeed a Groebner basis::

            sage: out = tmp_filename()
            sage: with redirection(sys.stdout,  open(out, 'w')):
            ....:     hilb(I,attributes={I:{'isSB':1}}) # no complaint from Singular
            sage: with open(out) as f:
            ....:     'is no standard basis' in f.read()
            False


        TESTS:

        We show that the interface recovers gracefully from errors::

            sage: P.<e,d,c,b,a> = PolynomialRing(QQ,5,order='lex')
            sage: I = sage.rings.ideal.Cyclic(P)

            sage: triangL = sage.libs.singular.function_factory.ff.triang__lib.triangL
            sage: _ = triangL(I)
            Traceback (most recent call last):
            ...
            RuntimeError: error in Singular function call 'triangL':
            The input is no groebner basis.
            leaving triang.lib::triangL (0)

        Flush any stray output -- see :issue:`28622`::

            sage: sys.stdout.flush()
            ...

            sage: G= Ideal(I.groebner_basis())
            sage: triangL(G,attributes={G:{'isSB':1}})
            [[e + d + c + b + a, ...]]
        """
        global dummy_ring

        if ring is None:
            ring = self.common_ring(args, ring)
            if ring is None:
                if dummy_ring is None:
                    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                    from sage.rings.rational_field import QQ
                    dummy_ring = PolynomialRing(QQ, "dummy", implementation='singular') # seems a reasonable default
                ring = dummy_ring
        if not (isinstance(ring, MPolynomialRing_libsingular) or isinstance(ring, NCPolynomialRing_plural)):
            raise TypeError("cannot call Singular function '%s' with ring parameter of type '%s'" % (self._name,type(ring)))
        return call_function(self, args, ring, interruptible, attributes)

    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: groebner = singular_function('groebner')
            sage: 'groebner' in groebner.__doc__  # needs info
            True
        """

        prefix = \
"""
This function is an automatically generated C wrapper around the Singular
function '%s'.

This wrapper takes care of converting Sage datatypes to Singular
datatypes and vice versa. In addition to whatever parameters the
underlying Singular function accepts when called, this function also
accepts the following keyword parameters:

INPUT:

- ``args`` -- list of arguments
- ``ring`` -- a multivariate polynomial ring
- ``interruptible`` -- if ``True`` pressing :kbd:`Ctrl` + :kbd:`C` during the
  execution of this function will interrupt the computation
  (default: ``True``)
- ``attributes`` -- dictionary of optional Singular attributes
  assigned to Singular objects (default: ``None``)

If ``ring`` is not specified, it is guessed from the given arguments.
If this is not possible, then a dummy ring, univariate polynomial ring
over ``QQ``, is used.

EXAMPLES::

    sage: groebner = sage.libs.singular.function_factory.ff.groebner
    sage: P.<x, y> = PolynomialRing(QQ)
    sage: I = P.ideal(x^2-y, y+x)
    sage: groebner(I)
    [x + y, y^2 - y]
    sage: triangL = sage.libs.singular.function_factory.ff.triang__lib.triangL
    sage: P.<x1, x2> = PolynomialRing(QQ, order='lex')
    sage: f1 = 1/2*((x1^2 + 2*x1 - 4)*x2^2 + 2*(x1^2 + x1)*x2 + x1^2)
    sage: f2 = 1/2*((x1^2 + 2*x1 + 1)*x2^2 + 2*(x1^2 + x1)*x2 - 4*x1^2)
    sage: I = Ideal(Ideal(f1,f2).groebner_basis()[::-1])
    sage: triangL(I, attributes={I:{'isSB':1}})
    [[x2^4 + 4*x2^3 - 6*x2^2 - 20*x2 + 5, 8*x1 - x2^3 + x2^2 + 13*x2 - 5],
     [x2, x1^2],
     [x2, x1^2],
     [x2, x1^2]]

"""%(self._name)
        from sage.interfaces.singular import get_docstring
        return prefix + get_docstring(self._name, prefix=True, code=True)

    cdef common_ring(self, tuple args, ring=None):
        """
        Return the common ring for the argument list ``args``.

        If ``ring`` is not ``None`` this routine checks whether it is
        the parent/ring of all members of ``args`` instead.

        If no common ring was found, None is returned.

        INPUT:

        - ``args`` -- list of Python objects
        - ``ring`` -- an optional ring to check
        """
        from sage.matrix.matrix_mpolynomial_dense import Matrix_mpolynomial_dense
        from sage.matrix.matrix_integer_dense import Matrix_integer_dense
        ring2 = None
        for a in args:
            if isinstance(a, MPolynomialIdeal) or \
                    isinstance(a, NCPolynomialIdeal):
                ring2 = a.ring()
            elif is_singular_poly_wrapper(a):
                ring2 = a.parent()
            elif is_sage_wrapper_for_singular_ring(a):
                ring2 = a
            elif isinstance(a, (int, str)):
                continue
            elif isinstance(a, Matrix_integer_dense):
                continue
            elif isinstance(a, Matrix_mpolynomial_dense):
                ring2 = a.base_ring()
            elif isinstance(a, (list, tuple, Sequence_generic)):
                #TODO: catch exception, if recursion finds no ring
                ring2 = self.common_ring(tuple(a), ring)
            elif isinstance(a, Resolution):
                ring2 = (<Resolution> a).base_ring
            elif isinstance(a, FreeModuleElement_generic_dense)\
                and is_sage_wrapper_for_singular_ring(
                    a.parent().base_ring()):
                ring2 = a.parent().base_ring()
            elif ring is not None:
                a.parent() is ring
                continue

            if ring is None:
                ring = ring2
            elif ring is not ring2:
                raise ValueError("Rings do not match up.")
        return ring

    def __reduce__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: groebner = singular_function('groebner')
            sage: groebner == loads(dumps(groebner))
            True
        """
        return singular_function, (self._name,)

    def __richcmp__(self, other, op):
        """
        EXAMPLES::

            sage: from sage.libs.singular.function import singular_function
            sage: groebner = singular_function('groebner')
            sage: groebner == singular_function('groebner')
            True
            sage: groebner == singular_function('std')
            False
            sage: groebner == 1
            False
            sage: groebner == None
            False
        """
        try:
            lx = <SingularFunction?>self
            rx = <SingularFunction?>other
        except TypeError:
            return NotImplemented
        return richcmp(lx._name, rx._name, op)


cdef inline call_function(SingularFunction self, tuple args, object R, bint signal_handler=True, attributes=None):
    global currRingHdl
    global errorreported
    global currentVoice
    global myynest

    cdef ring *si_ring
    if isinstance(R, MPolynomialRing_libsingular):
        si_ring = (<MPolynomialRing_libsingular>R)._ring
    else:
        si_ring = (<NCPolynomialRing_plural>R)._ring

    if si_ring != currRing: rChangeCurrRing(si_ring)

    if currRingHdl.data.uring!= currRing:
        currRingHdl.data.uring.ref -= 1
        currRingHdl.data.uring = currRing # ref counting?
        currRingHdl.data.uring.ref += 1

    cdef Converter argument_list = Converter(args, R, attributes)

    cdef leftv * _res

    currentVoice = NULL
    myynest = 0
    start_catch_error()

    with opt_ctx: # we are preserving the global options state here
        if signal_handler:
            try:
                sig_on()
                _res = self.call_handler.handle_call(argument_list, si_ring)
                sig_off()
            finally:
                s = check_error()
        else:
            _res = self.call_handler.handle_call(argument_list, si_ring)
            s = check_error()

    if myynest:
        myynest = 0

    if currentVoice:
        currentVoice = NULL

    if s:
        raise RuntimeError("error in Singular function call %r:\n%s" %
            (self._name, "\n".join(s)))

    res = argument_list.to_python(_res)

    if self.call_handler.free_res():
        free_leftv(_res, si_ring)
    else:
        _res.CleanUp(si_ring)

    return res


cdef class SingularLibraryFunction(SingularFunction):
    """
    EXAMPLES::

        sage: from sage.libs.singular.function import SingularLibraryFunction
        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: I = R.ideal(x, x+1)
        sage: f = SingularLibraryFunction("groebner")
        sage: f(I)
        [1]
    """
    def __init__(self, name):
        """
        Construct a new Singular library function.

        EXAMPLES::

            sage: from sage.libs.singular.function import SingularLibraryFunction
            sage: R.<x,y> = PolynomialRing(QQ, order='lex')
            sage: I = R.ideal(x + 1, x*y + 1)
            sage: f = SingularLibraryFunction("groebner")
            sage: f(I)
            [y - 1, x + 1]
        """
        super().__init__(name)
        self.call_handler = self.get_call_handler()

    cdef BaseCallHandler get_call_handler(self):
        cdef idhdl* singular_idhdl = ggetid(str_to_bytes(self._name))
        if singular_idhdl==NULL:
            raise NameError("Singular library function {!r} is not defined".format(self._name))
        if singular_idhdl.typ!=PROC_CMD:
            raise ValueError("Not a procedure")

        cdef LibraryCallHandler res = LibraryCallHandler()
        res.proc_idhdl = singular_idhdl
        return res

    cdef bint function_exists(self) noexcept:
        cdef idhdl* singular_idhdl = ggetid(str_to_bytes(self._name))
        return singular_idhdl!=NULL


cdef class SingularKernelFunction(SingularFunction):
    """
    EXAMPLES::

        sage: from sage.libs.singular.function import SingularKernelFunction
        sage: R.<x,y> = PolynomialRing(QQ, order='lex')
        sage: I = R.ideal(x, x+1)
        sage: f = SingularKernelFunction("std")
        sage: f(I)
        [1]
    """
    def __init__(self, name):
        """
        Construct a new Singular kernel function.

        EXAMPLES::

            sage: from sage.libs.singular.function import SingularKernelFunction
            sage: R.<x,y> = PolynomialRing(QQ, order='lex')
            sage: I = R.ideal(x + 1, x*y + 1)
            sage: f = SingularKernelFunction("std")
            sage: f(I)
            [y - 1, x + 1]
            sage: SingularKernelFunction("no_such_function")
            Traceback (most recent call last):
            ...
            NameError: Singular kernel function 'no_such_function' is not defined
        """
        super().__init__(name)
        self.call_handler = self.get_call_handler()

    cdef BaseCallHandler get_call_handler(self):
        cdef int cmd_n = 0
        arity = IsCmd(str_to_bytes(self._name), cmd_n) # call by reverence for CMD_n
        if not cmd_n:
            raise NameError("Singular kernel function {!r} is not defined".format(self._name))

        return KernelCallHandler(cmd_n, arity)

    cdef bint function_exists(self) noexcept:
        cdef int cmd_n = -1
        arity = IsCmd(str_to_bytes(self._name), cmd_n) # call by reverence for CMD_n
        return cmd_n != -1


def singular_function(name):
    """
    Construct a new libSingular function object for the given
    ``name``.

    This function works both for interpreter and built-in functions.

    INPUT:

    - ``name`` -- the name of the function

    EXAMPLES::

        sage: P.<x,y,z> = PolynomialRing(QQ)
        sage: f = 3*x*y + 2*z + 1
        sage: g = 2*x + 1/2
        sage: I = Ideal([f,g])

    ::

        sage: from sage.libs.singular.function import singular_function
        sage: std = singular_function("std")
        sage: std(I)
        [3*y - 8*z - 4, 4*x + 1]
        sage: size = singular_function("size")
        sage: size([2, 3, 3])
        3
        sage: size("sage")
        4
        sage: size(["hello", "sage"])
        2
        sage: factorize = singular_function("factorize")
        sage: factorize(f)
        [[1, 3*x*y + 2*z + 1], (1, 1)]
        sage: factorize(f, 1)
        [3*x*y + 2*z + 1]

    We give a wrong number of arguments::

        sage: factorize()
        Traceback (most recent call last):
        ...
        RuntimeError: error in Singular function call 'factorize':
        Wrong number of arguments (got 0 arguments, arity is CMD_12)
        sage: factorize(f, 1, 2)
        Traceback (most recent call last):
        ...
        RuntimeError: error in Singular function call 'factorize':
        Wrong number of arguments (got 3 arguments, arity is CMD_12)
        sage: factorize(f, 1, 2, 3)
        Traceback (most recent call last):
        ...
        RuntimeError: error in Singular function call 'factorize':
        Wrong number of arguments (got 4 arguments, arity is CMD_12)

    The Singular function ``list`` can be called with any number of
    arguments::

        sage: singular_list = singular_function("list")
        sage: singular_list(2, 3, 6)
        [2, 3, 6]
        sage: singular_list()
        []
        sage: singular_list(1)
        [1]
        sage: singular_list(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    We try to define a non-existing function::

        sage: number_foobar = singular_function('number_foobar')
        Traceback (most recent call last):
        ...
        NameError: Singular library function 'number_foobar' is not defined

    ::

        sage: from sage.libs.singular.function import lib as singular_lib
        sage: singular_lib('general.lib')
        sage: number_e = singular_function('number_e')
        sage: number_e(10r)
        67957045707/25000000000
        sage: RR(number_e(10r))
        2.71828182828000

    ::

        sage: singular_lib('primdec.lib')
        sage: primdecGTZ = singular_function("primdecGTZ")
        sage: primdecGTZ(I)
        [[[y - 8/3*z - 4/3, x + 1/4], [y - 8/3*z - 4/3, x + 1/4]]]
        sage: singular_list((1,2,3),3,[1,2,3], ring=P)
        [(1, 2, 3), 3, [1, 2, 3]]
        sage: ringlist=singular_function("ringlist")
        sage: l = ringlist(P)
        sage: l[3].__class__
        <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>
        sage: l
        [0, ['x', 'y', 'z'], [['dp', (1, 1, 1)], ['C', (0,)]], [0]]
        sage: ring=singular_function("ring")
        sage: ring(l)
        <RingWrap>
        sage: matrix = Matrix(P,2,2)
        sage: matrix.randomize(terms=1)
        sage: det = singular_function("det")
        sage: det(matrix) == matrix[0, 0] * matrix[1, 1] - matrix[0, 1] * matrix[1, 0]
        True
        sage: coeffs = singular_function("coeffs")
        sage: coeffs(x*y+y+1,y)
        [    1]
        [x + 1]
        sage: intmat = Matrix(ZZ, 2,2, [100,2,3,4])
        sage: det(intmat)
        394
        sage: random = singular_function("random")
        sage: A = random(10,2,3); A.nrows(), max(A.list()) <= 10
        (2, True)
        sage: P.<x,y,z> = PolynomialRing(QQ)
        sage: M=P**3
        sage: leadcoef = singular_function("leadcoef")
        sage: v=M((100*x,5*y,10*z*x*y))
        sage: leadcoef(v)
        10
        sage: v = M([x+y,x*y+y**3,z])
        sage: lead = singular_function("lead")
        sage: lead(v)
        (0, y^3)
        sage: jet = singular_function("jet")
        sage: jet(v, 2)
        (x + y, x*y, z)
        sage: syz = singular_function("syz")
        sage: I = P.ideal([x+y,x*y-y, y*2,x**2+1])
        sage: M = syz(I)
        sage: M
        [(-2*y, 2, y + 1, 0), (0, -2, x - 1, 0), (x*y - y, -y + 1, 1, -y), (x^2 + 1, -x - 1, -1, -x)]
        sage: singular_lib("mprimdec.lib")
        sage: syz(M)
        [(-x - 1, y - 1, 2*x, -2*y)]
        sage: GTZmod = singular_function("GTZmod")
        sage: GTZmod(M)
        [[[(-2*y, 2, y + 1, 0), (0, x + 1, 1, -y), (0, -2, x - 1, 0), (x*y - y, -y + 1, 1, -y), (x^2 + 1, 0, 0, -x - y)], [0]]]
        sage: mres = singular_function("mres")
        sage: resolution = mres(M, 0)
        sage: resolution
        <Resolution>
        sage: singular_list(resolution)
        [[(-2*y, 2, y + 1, 0), (0, -2, x - 1, 0), (x*y - y, -y + 1, 1, -y), (x^2 + 1, -x - 1, -1, -x)], [(-x - 1, y - 1, 2*x, -2*y)], [(0)]]

        sage: A.<x,y> = FreeAlgebra(QQ, 2)
        sage: P.<x,y> = A.g_algebra({y*x:-x*y})
        sage: I= Sequence([x*y,x+y], check=False, immutable=True)
        sage: twostd = singular_function("twostd")
        sage: twostd(I)
        [x + y, y^2]
        sage: M=syz(I)
        doctest...
        sage: M
        [(x + y, x*y)]
        sage: syz(M)
        [(0)]
        sage: mres(I, 0)
        <Resolution>
        sage: M=P**3
        sage: v=M((100*x,5*y,10*y*x*y))
        sage: leadcoef(v)
        -10
        sage: v = M([x+y,x*y+y**3,x])
        sage: lead(v)
        (0, y^3)
        sage: jet(v, 2)
        (x + y, x*y, x)
        sage: l = ringlist(P)
        sage: len(l)
        6
        sage: ring(l)
        <noncommutative RingWrap>
        sage: I=twostd(I)
        sage: l[3]=I
        sage: ring(l)
        <noncommutative RingWrap>
    """
    try:
        return SingularKernelFunction(name)
    except NameError:
        return SingularLibraryFunction(name)


def lib(name):
    """
    Load the Singular library ``name``.

    INPUT:

    - ``name`` -- a Singular library name

    EXAMPLES::

        sage: from sage.libs.singular.function import singular_function
        sage: from sage.libs.singular.function import lib as singular_lib
        sage: singular_lib('general.lib')
        sage: primes = singular_function('primes')
        sage: primes(2,10, ring=GF(127)['x,y,z'])
        (2, 3, 5, 7)
    """
    global si_opt_2

    cdef int vv = si_opt_2

    if get_verbose() <= 0:
        si_opt_2 &= ~Sy_bit(V_LOAD_LIB)
        si_opt_2 &= ~Sy_bit(V_REDEFINE)

    cdef char* cname = omStrDup(str_to_bytes(name))
    sig_on()
    cdef bint failure = iiLibCmd(cname, 1, 1, 1)
    sig_off()
    si_opt_2 = vv

    if failure:
        raise NameError("Singular library {!r} not found".format(name))


def get_printlevel():
    """
    Return the value of the variable ``printlevel``.

    This is useful to switch off and back the comments.

    EXAMPLES::

        sage: from sage.libs.singular.function import get_printlevel, set_printlevel
        sage: l = get_printlevel()
        sage: set_printlevel(-1)
        sage: get_printlevel()
        -1
        sage: set_printlevel(l)
    """
    global printlevel
    cdef int pl = printlevel
    return pl


def set_printlevel(l):
    """
    Set the  value of the variable ``printlevel``.

    This is useful to switch off and back the comments.

    EXAMPLES::

        sage: from sage.libs.singular.function import get_printlevel, set_printlevel
        sage: l = get_printlevel()
        sage: set_printlevel(2)
        sage: get_printlevel()
        2
        sage: set_printlevel(l)
    """
    global printlevel
    printlevel = <int>l


def list_of_functions(packages=False):
    """
    Return a list of all function names currently available.

    INPUT:

    - ``packages`` -- include local functions in packages

    EXAMPLES::

        sage: from sage.libs.singular.function import list_of_functions
        sage: 'groebner' in list_of_functions()
        True
    """
    cdef list l = []
    cdef idhdl *h=IDROOT
    cdef idhdl *ph = NULL
    while h!=NULL:
        if PROC_CMD == IDTYP(h):
            l.append(char_to_str(h.id))
        if PACKAGE_CMD == IDTYP(h):
            if packages:
                ph = IDPACKAGE(h).idroot
                while ph != NULL:
                    if PROC_CMD == IDTYP(ph):
                        l.append(char_to_str(ph.id))
                    ph = IDNEXT(ph)
        h = IDNEXT(h)
    return l


cdef inline RingWrap new_RingWrap(ring* r):
    cdef RingWrap ring_wrap_result = RingWrap.__new__(RingWrap)
    ring_wrap_result._ring = r
    ring_wrap_result._ring.ref += 1

    return ring_wrap_result

# Add support for _instancedoc_
from sage.misc.instancedoc import instancedoc
instancedoc(SingularFunction)
instancedoc(SingularLibraryFunction)
instancedoc(SingularKernelFunction)
