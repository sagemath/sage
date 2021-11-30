r"""
Relative finite field extensions

Considering a *absolute field* `F_{q^m}` and a *relative_field* `F_q`, with
`q = p^s`, `p` being a prime and `s, m` being integers, this file
contains a class to take care of the representation of `F_{q^m}`-elements
as `F_q`-elements.

.. WARNING::

    As this code is experimental, a warning is thrown when a
    relative finite field extension is created for the first time
    in a session (see :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: from sage.coding.relative_finite_field_extension import *
        sage: Fqm.<aa> = GF(16)
        sage: Fq.<a> = GF(4)
        sage: RelativeFiniteFieldExtension(Fqm, Fq)
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/20284 for details.
        Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
"""

# ****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.categories.homset import Hom
from sage.matrix.constructor import column_matrix
from sage.modules.free_module_element import vector
from sage.misc.superseded import experimental


class RelativeFiniteFieldExtension(SageObject):
    r"""
    Considering `p` a prime number, n an integer and three finite fields
    `F_p`, `F_q` and `F_{q^m}`, this class contains a set of methods
    to manage the representation of elements of the relative extension
    `F_{q^m}` over `F_q`.

    INPUT:

    - ``absolute_field``, ``relative_field`` -- two finite fields, ``relative_field``
      being a subfield of ``absolute_field``

    - ``embedding`` -- (default: ``None``) an homomorphism from ``relative_field`` to
      ``absolute_field``. If ``None`` is provided, it will default to the first
      homomorphism of the list of homomorphisms Sage can build.

    EXAMPLES::

        sage: from sage.coding.relative_finite_field_extension import *
        sage: Fqm.<aa> = GF(16)
        sage: Fq.<a> = GF(4)
        sage: RelativeFiniteFieldExtension(Fqm, Fq)
        Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2

    It is possible to specify the embedding to use
    from ``relative_field`` to ``absolute_field``::

        sage: Fqm.<aa> = GF(16)
        sage: Fq.<a> = GF(4)
        sage: FE = RelativeFiniteFieldExtension(Fqm, Fq, embedding=Hom(Fq, Fqm)[1])
        sage: FE.embedding() == Hom(Fq, Fqm)[1]
        True
    """

    @experimental(trac_number=20284)
    def __init__(self, absolute_field, relative_field, embedding=None):
        r"""
        TESTS:

        If ``absolute_field`` is not a finite field, an error is raised::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm = RR
            sage: Fq.<a> = GF(4)
            sage: RelativeFiniteFieldExtension(Fqm, Fq)
            Traceback (most recent call last):
            ...
            ValueError: absolute_field has to be a finite field

        Same for ``relative_field``::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq = RR
            sage: RelativeFiniteFieldExtension(Fqm, Fq)
            Traceback (most recent call last):
            ...
            ValueError: relative_field has to be a finite field

        If ``relative_field`` is not a subfield of ``absolute_field``, an exception
        is raised::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(8)
            sage: RelativeFiniteFieldExtension(Fqm, Fq)
            Traceback (most recent call last):
            ...
            ValueError: relative_field has to be a subfield of absolute_field
        """
        if not absolute_field.is_finite():
            raise ValueError("absolute_field has to be a finite field")
        if not relative_field.is_finite():
            raise ValueError("relative_field has to be a finite field")
        s = relative_field.degree()
        sm = absolute_field.degree()
        if not s.divides(sm):
            raise ValueError("relative_field has to be a subfield of absolute_field")
        H = Hom(relative_field, absolute_field)
        if embedding is not None and embedding not in H:
            raise ValueError("embedding has to be an embedding from relative_field to absolute_field")
        elif embedding is not None:
            self._phi = embedding
        else:
            self._phi = H[0]
        self._prime_field = relative_field.base_ring()
        self._relative_field = relative_field
        self._absolute_field = absolute_field
        alpha = relative_field.gen()
        beta = absolute_field.gen()
        self._alphas = [alpha ** i for i in range(s)]
        self._betas = [beta ** i for i in range(sm)]
        self._relative_field_degree = s
        self._absolute_field_degree = sm

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: RelativeFiniteFieldExtension(Fqm, Fq)
            Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
        """
        return "Relative field extension between %s and %s" % (self.absolute_field(), self.relative_field())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: latex(RelativeFiniteFieldExtension(Fqm, Fq))
            \textnormal{Relative field extension between \Bold{F}_{2^{4}} and \Bold{F}_{2^{2}}}
        """
        return "\\textnormal{Relative field extension between %s and %s}" % (self.absolute_field()._latex_(),
                self.relative_field()._latex_())

    def __eq__(self, other):
        r"""
        Tests equality between embeddings.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fq = GF(4)
            sage: FQ = GF(4**3)
            sage: H = Hom(Fq, FQ)
            sage: E1 = RelativeFiniteFieldExtension(FQ, Fq)
            sage: E2 = RelativeFiniteFieldExtension(FQ, Fq, H[0])
            sage: E3 = RelativeFiniteFieldExtension(FQ, Fq, H[1])
            sage: E1 == E2
            True
            sage: E1 == E3
            False
        """
        return isinstance(other, RelativeFiniteFieldExtension) \
                and self.embedding() == other.embedding()

    @cached_method
    def _representation_matrix(self):
        r"""
        Returns the matrix used to represents elements of the absolute field
        as vectors in the basis of the relative field over the prime field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE._representation_matrix()
            [1 0 0 0]
            [0 0 1 1]
            [0 1 1 1]
            [0 0 0 1]
        """
        s = self.relative_field_degree()
        m = self.extension_degree()
        betas = self.absolute_field_basis()
        phi_alphas = [ self._phi(self._alphas[i]) for i in range(s) ]
        A = column_matrix([vector(betas[i] * phi_alphas[j])
            for i in range(m) for j in range(s)])
        return A.inverse()

    def _flattened_relative_field_representation(self, b):
        r"""
        Returns a vector representation of ``b`` in the basis of
        the relative field over the prime field.

        INPUT:

        - ``b`` -- an element of the absolute field

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: b = aa^3 + aa^2 + aa + 1
            sage: FE._flattened_relative_field_representation(b)
            (1, 0, 1, 1)
        """
        if b not in self.absolute_field():
            raise ValueError("The input has to be an element of the absolute field")
        return self._representation_matrix() * vector(b)

    def relative_field_representation(self, b):
        r"""
        Returns a vector representation of the field element ``b`` in the basis
        of the absolute field over the relative field.

        INPUT:

        - ``b`` -- an element of the absolute field

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: b = aa^3 + aa^2 + aa + 1
            sage: FE.relative_field_representation(b)
            (1, a + 1)
        """
        if b not in self.absolute_field():
            raise ValueError("The input has to be an element of the absolute field")
        s = self.relative_field_degree()
        if s == 1:
            return vector(b)
        Fq = self.relative_field()
        vect = self._flattened_relative_field_representation(b)
        sm = self.absolute_field_degree()
        list_elts = []
        for i in range(0, sm, s):
            list_elts.append(Fq(vect[i:i + s]))
        return vector(Fq, list_elts)

    def absolute_field_representation(self, a):
        r"""
        Returns an absolute field representation of the relative field
        vector ``a``.

        INPUT:

        - ``a`` -- a vector in the relative extension field

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: b = aa^3 + aa^2 + aa + 1
            sage: rel = FE.relative_field_representation(b)
            sage: FE.absolute_field_representation(rel) == b
            True
        """
        s = self.relative_field_degree()
        m = self.extension_degree()
        if len(a) != m:
            raise ValueError("The input has to be a vector with length equal to the order of the absolute field")
        if not a.base_ring() == self.relative_field():
            raise ValueError("The input has to be over the prime field")
        alphas = self.relative_field_basis()
        betas = self.absolute_field_basis()
        phi = self.embedding()
        b = self.absolute_field().zero()
        flattened_relative_field_rep_list = []
        for i in a:
            tmp = vector(i).list()
            for j in tmp:
                flattened_relative_field_rep_list.append(j)

        flattened_relative_field_rep = vector(flattened_relative_field_rep_list)
        for i in range(m):
            b += betas[i] * phi(sum([flattened_relative_field_rep[j] * alphas[j%s] for j in range(i*s, i*s + s)]))
        return b

    def is_in_relative_field(self, b):
        r"""
        Returns ``True`` if ``b`` is in the relative field.

        INPUT:

        - ``b`` -- an element of the absolute field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.is_in_relative_field(aa^2 + aa)
            True
            sage: FE.is_in_relative_field(aa^3)
            False
        """
        vect = self.relative_field_representation(b)
        return vect[1:vect.length()].is_zero()

    def cast_into_relative_field(self, b, check=True):
        r"""
        Casts an absolute field element into the relative field (if possible).
        This is the inverse function of the field embedding.

        INPUT:

        - ``b`` -- an element of the absolute field which also lies in the
          relative field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: phi = FE.embedding()
            sage: b = aa^2 + aa
            sage: FE.is_in_relative_field(b)
            True
            sage: FE.cast_into_relative_field(b)
            a
            sage: phi(FE.cast_into_relative_field(b)) == b
            True
        """
        if check:
            if not self.is_in_relative_field(b):
                raise ValueError("%s does not belong to the relative field" % b)
        return self.relative_field_representation(b)[0]

    def embedding(self):
        r"""
        Returns the embedding which is used to go from the
        relative field to the absolute field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.embedding()
            Ring morphism:
             From: Finite Field in a of size 2^2
             To:   Finite Field in aa of size 2^4
             Defn: a |--> aa^2 + aa
        """
        return self._phi

    def relative_field_basis(self):
        r"""
        Returns a basis of the relative field over the prime field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.relative_field_basis()
            [1, a]
        """
        return self._alphas

    def absolute_field_basis(self):
        r"""
        Returns a basis of the absolute field over the prime field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.absolute_field_basis()
            [1, aa, aa^2, aa^3]
        """
        return self._betas

    def relative_field_degree(self):
        r"""
        Let `F_p` be the base field of our relative field `F_q`.
        Returns `s` where `p^s = q`

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.relative_field_degree()
            2
        """
        return self._relative_field_degree

    def absolute_field_degree(self):
        r"""
        Let `F_p` be the base field of our absolute field `F_{q^m}`.
        Returns `sm` where `p^{sm} = q^{m}`

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.absolute_field_degree()
            4
        """
        return self._absolute_field_degree


    def extension_degree(self):
        r"""
        Return `m`, the extension degree of the absolute field over
        the relative field.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(64)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.extension_degree()
            3
        """
        return self.absolute_field_degree() // self.relative_field_degree()

    def prime_field(self):
        r"""
        Returns the base field of our absolute and relative fields.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.prime_field()
            Finite Field of size 2
        """
        return self._prime_field

    def relative_field(self):
        r"""
        Returns the relative field of ``self``.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.relative_field()
            Finite Field in a of size 2^2
        """
        return self._relative_field

    def absolute_field(self):
        r"""
        Returns the absolute field of ``self``.

        EXAMPLES::

            sage: from sage.coding.relative_finite_field_extension import *
            sage: Fqm.<aa> = GF(16)
            sage: Fq.<a> = GF(4)
            sage: FE = RelativeFiniteFieldExtension(Fqm, Fq)
            sage: FE.absolute_field()
            Finite Field in aa of size 2^4
        """
        return self._absolute_field
