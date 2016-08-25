r"""
Combinatorial Algebras

A combinatorial algebra is an algebra whose basis elements are
indexed by a combinatorial class. Some examples of combinatorial
algebras are the symmetric group algebra of order n (indexed by
permutations of size n) and the algebra of symmetric functions
(indexed by integer partitions).

The CombinatorialAlgebra base class makes it easy to define and
work with new combinatorial algebras in Sage. For example, the
following code constructs an algebra which models the power-sum
symmetric functions.

::

    sage: class PowerSums(CombinatorialAlgebra):
    ...     def __init__(self, R):
    ...         self._one = Partition([])
    ...         self._name = 'Power-sum symmetric functions'
    ...         CombinatorialAlgebra.__init__(self, R, Partitions())
    ...         self.print_options(prefix='p')
    ...
    ...     def _multiply_basis(self, a, b):
    ...         l = list(a)+list(b)
    ...         l.sort(reverse=True)
    ...         return Partition(l)
    ...

::

    sage: ps = PowerSums(QQ); ps
    Power-sum symmetric functions over Rational Field
    sage: ps([2,1])^2
    p[2, 2, 1, 1]
    sage: ps([2,1])+2*ps([1,1,1])
    2*p[1, 1, 1] + p[2, 1]
    sage: ps(2)
    2*p[]

The important things to define are ._indices which
specifies the combinatorial class that indexes the basis elements,
._one which specifies the identity element in the algebra, ._name
which specifies the name of the algebra, .print_options is used to set
the print options for the elements, and finally a _multiply
or _multiply_basis method that defines the multiplication in the
algebra.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.misc import repr_lincomb
from sage.misc.cachefunc import cached_method
from sage.categories.all import AlgebrasWithBasis
from sage.structure.element import Element

# TODO: migrate this completely to the combinatorial free module + categories framework

# for backward compatibility
CombinatorialAlgebraElement = CombinatorialFreeModule.Element

# Not used anymore
class CombinatorialAlgebraElementOld(CombinatorialFreeModule.Element):
#     def __init__(self, A, x):
#         """
#         Create a combinatorial algebra element x.  This should never
#         be called directly, but only through the parent combinatorial
#         algebra's __call__ method.

#         TESTS::
#
#             sage: s = SFASchur(QQ)
#             sage: a = s._element_class(s, {Partition([2,1]):QQ(2)}); a
#             2*s[2, 1]
#             sage: a == loads(dumps(a))
#             True
#         """
#         AlgebraElement.__init__(self, A)
#         self._monomial_coefficients = x


    def _mul_(self, y):
        """
        EXAMPLES::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: a = s([2])
            sage: a._mul_(a) #indirect doctest
            s[2, 2] + s[3, 1] + s[4]
        """
        return self.parent().product(self, y)


    def __invert__(self):
        """
        EXAMPLES::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: ~s(2)
            1/2*s[]
            sage: ~s([2,1])
            Traceback (most recent call last):
            ...
            ValueError: cannot invert self (= s[2, 1])
        """
        mcs = self._monomial_coefficients
        one = self.parent()._one
        if len(mcs) == 1 and one in mcs:
            return self.parent()( ~mcs[ one ] )
        else:
            raise ValueError("cannot invert self (= %s)"%self)





    def __repr__(self):
        """
        EXAMPLES::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: a = 2 + s([3,2,1])
            sage: print(a.__repr__())
            2*s[] + s[3, 2, 1]
        """
        v = sorted(self._monomial_coefficients.items())
        prefix = self.parent().prefix()
        retur = repr_lincomb( [(prefix + repr(m), c) for m,c in v ], strip_one = True)

class CombinatorialAlgebra(CombinatorialFreeModule):
    """

    Deprecated! Don't use!

    """

    # For backward compatibility
    @cached_method
    def one_basis(self):
        """
        TESTS::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: s.one_basis()
            []
        """
        return self._one

    def __init__(self, R, cc = None, element_class = None, category = None):
        """
        TESTS::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: TestSuite(s).run()
        """
        #Check to make sure that the user defines the necessary
        #attributes / methods to make the combinatorial module
        #work
        required = ['_one',]
        for r in required:
            if not hasattr(self, r):
                raise ValueError("%s is required"%r)
        if not hasattr(self, '_multiply') and not hasattr(self, '_multiply_basis'):
            raise ValueError("either _multiply or _multiply_basis is required")

        #Create a class for the elements of this combinatorial algebra
        #We need to do this so to distinguish between element of different
        #combinatorial algebras
#         if element_class is None:
#             if not hasattr(self, '_element_class'):
#                 class CAElement(CombinatorialAlgebraElement):
#                     pass
#                 self._element_class = CAElement
#         else:
#             self._element_class = element_class

        if category is None:
            category = AlgebrasWithBasis(R)

        # for backward compatibility
        if cc is None:
            if hasattr(self, "_indices"):
                cc = self._indices
        assert(cc is not None)

        CombinatorialFreeModule.__init__(self, R, cc, element_class, category = category)

    # see sage.combinat.free_module._repr_term
    # this emulates the _repr_ of CombinatorialAlgebraElement did not add brackets around the basis indices
    _repr_option_bracket = False

    def __call__(self, x):
        """
        TESTS::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: s([2])
            s[2]
        """
        try:
            return CombinatorialFreeModule.__call__(self, x)
        except TypeError:
            pass

        R = self.base_ring()
        eclass = self._element_class
        #Coerce elements of the base ring
        if not isinstance(x, Element):
            x = R(x)
        if x.parent() == R:
            if x == R(0):
                return eclass(self, {})
            else:
                return eclass(self, {self._one:x})
        #Coerce things that coerce into the base ring
        elif R.has_coerce_map_from(x.parent()):
            rx = R(x)
            if rx == R(0):
                return eclass(self, {})
            else:
                return eclass(self, {self._one:R(x)})

        raise TypeError("do not know how to make x (= %s) an element of self (=%s)"%(x,self))

    def _an_element_impl(self):
        """
        Returns an element of self, namely the unit element.

        EXAMPLES::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: s._an_element_impl()
            s[]
            sage: _.parent() is s
            True
        """
        return self._element_class(self, {self._one:self.base_ring()(1)})

    def _coerce_impl(self, x):
        """
        EXAMPLES::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: s(2)          # indirect doctest
            2*s[]
        """
        try:
            R = x.parent()
            if R.__class__ is self.__class__:
                #Only perform the coercion if we can go from the base
                #ring of x to the base ring of self
                if self.base_ring().has_coerce_map_from( R.base_ring() ):
                    return self(x)
        except AttributeError:
            pass

        # any ring that coerces to the base ring
        return self._coerce_try(x, [self.base_ring()])

    def product(self, left, right):
        """
        Returns left\*right where left and right are elements of self.
        product() uses either _multiply or _multiply basis to carry out
        the actual multiplication.

        EXAMPLES::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: a = s([2])
            sage: s.product(a,a)
            s[2, 2] + s[3, 1] + s[4]
            sage: ZS3 = SymmetricGroupAlgebra(ZZ, 3)
            sage: a = 2 + ZS3([2,1,3])
            sage: a*a
            5*[1, 2, 3] + 4*[2, 1, 3]
            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ,3)
            sage: j2 = H3.jucys_murphy(2)
            sage: j2*j2
            (q^3-q^2+q)*T[1, 2, 3] + (q^3-q^2+q-1)*T[2, 1, 3]
            sage: X = SchubertPolynomialRing(ZZ)
            sage: X([1,2,3])*X([2,1,3])
            X[2, 1]
        """
        A = left.parent()
        ABR = A.base_ring()
        ABRzero = ABR(0)
        z_elt = {}

        #Do the case where the user specifies how to multiply basis elements
        if hasattr(self, '_multiply_basis'):
            for (left_m, left_c) in left._monomial_coefficients.iteritems():
                for (right_m, right_c) in right._monomial_coefficients.iteritems():
                    res = self._multiply_basis(left_m, right_m)
                    coeffprod = left_c * right_c
                    #Handle the case where the user returns a dictionary
                    #where the keys are the monomials and the values are
                    #the coefficients.  If res is not a dictionary, then
                    #it is assumed to be an element of self
                    if not isinstance(res, dict):
                        if isinstance(res, self._element_class):
                            res = res._monomial_coefficients
                        else:
                            z_elt[res] = z_elt.get(res, ABRzero) + coeffprod
                            continue
                    for m, c in res.iteritems():
                        z_elt[m] = z_elt.get(m, ABRzero) + coeffprod * c

        #We assume that the user handles the multiplication correctly on
        #his or her own, and returns a dict with monomials as keys and
        #coefficients as values
        else:
            m = self._multiply(left, right)
            if isinstance(m, self._element_class):
                return m
            if not isinstance(m, dict):
                z_elt = m.monomial_coefficients()
            else:
                z_elt = m

        #Remove all entries that are equal to 0
        BR = self.base_ring()
        zero = BR(0)
        del_list = []
        for m, c in z_elt.iteritems():
            if c == zero:
                del_list.append(m)
        for m in del_list:
            del z_elt[m]

        return self._from_dict(z_elt)

class TestAlgebra(CombinatorialAlgebra):

    _name = "TestAlgebra"

    def __init__(self, R):
        """
        TESTS::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: TestSuite(s).run()
        """
        from sage.combinat.partition import Partition, Partitions
        self._one = Partition([])
        CombinatorialAlgebra.__init__(self, R, cc=Partitions())

    def _multiply_basis(self, part1, part2):
        """
        TESTS::

            sage: s = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: a = s([2])
            sage: a._mul_(a) #indirect doctest
            s[2, 2] + s[3, 1] + s[4]
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        S = SymmetricFunctions(self.base_ring()).schur()
        return self.sum_of_terms(S(part1) * S(part2))

    def prefix(self):
        """
        TESTS::

            sage: sa = sage.combinat.combinatorial_algebra.TestAlgebra(QQ)
            sage: x = sa([2]); x  # indirect doctest
            s[2]
        """
        return "s"
