# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.combinat sage.graphs
r"""
Modular forms for Hecke triangle groups

AUTHORS:

- Jonas Jermann (2013): initial version
"""
# ****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import infinity

from sage.modules.module import Module
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from .hecke_triangle_groups import HeckeTriangleGroup
from .abstract_space import FormsSpace_abstract


def canonical_parameters(group, base_ring, k, ep, n=None):
    r"""
    Return a canonical version of the parameters.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.space import canonical_parameters
        sage: canonical_parameters(5, ZZ, 20/3, int(1))
        (Hecke triangle group for n = 5, Integer Ring, 20/3, 1, 5)

        sage: canonical_parameters(infinity, ZZ, 2, int(-1))
        (Hecke triangle group for n = +Infinity, Integer Ring, 2, -1, +Infinity)
    """
    if n is not None:
        group = n

    if group == infinity:
        group = HeckeTriangleGroup(infinity)
    else:
        try:
            group = HeckeTriangleGroup(ZZ(group))
        except TypeError:
            group = HeckeTriangleGroup(group.n())

    n = group.n()
    k = QQ(k)
    if ep is None:
        if n == infinity:
            ep = (-1)**(k / 2)
        elif ZZ(2).divides(n):
            ep = (-1)**(k*ZZ(n-2)/ZZ(4))
        else:
            ep = (-1)**(k*ZZ(n-2)/ZZ(2))
    ep = ZZ(ep)

    if n == infinity:
        num = (k-(1-ep)) / ZZ(4)
    else:
        num = (k-(1-ep)*n/(n-2)) * (n-2) / ZZ(4)

    try:
        num = ZZ(num)
    except TypeError:
        raise ValueError(f"Invalid or non-occurring weight k={k}, ep={ep}!")

    return (group, base_ring, k, ep, n)


class QuasiMeromorphicModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi meromorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiMeromorphicModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, ZZ, 20/3, int(1))
            sage: QuasiMeromorphicModularForms(5, ZZ, 20/3, int(1)) == QuasiMeromorphicModularForms(group, base_ring, k, ep, n)
            True
        """
        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi meromorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MF = QuasiMeromorphicModularForms(5, ZZ, 20/3, 1)
            sage: MF
            QuasiMeromorphicModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi meromorphic modular
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.ambient_space() == MF
            True
        """
        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["quasi", "mero"])


class QuasiWeakModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi weakly holomorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiWeakModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(4, ZZ, 8, -1)
            sage: QuasiWeakModularForms(4, ZZ, 8, -1) == QuasiWeakModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi weakly holomorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiWeakModularForms
            sage: MF = QuasiWeakModularForms(4, ZZ, 8, 1)
            sage: MF
            QuasiWeakModularForms(n=4, k=8, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi weakly holomorphic modular
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.is_ambient()
            True
        """
        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["quasi", "weak"])


class QuasiModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, ZZ, 10/3, -1)
            sage: QuasiModularForms(5, ZZ, 10/3) == QuasiModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(5, ZZ, 20/3, 1)
            sage: MF
            QuasiModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi modular
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["quasi", "holo"])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(n=5, k=6, ep=-1)
            sage: MF.default_prec(2)
            sage: MF.gens()
            [1 - 37/(200*d)*q + O(q^2),
             1 + 33/(200*d)*q + O(q^2),
             1 - 27/(200*d)*q + O(q^2)]

            sage: MF = QuasiModularForms(n=infinity, k=2, ep=-1)
            sage: MF.default_prec(2)
            sage: MF.gens()
            [1 - 24*q + O(q^2), 1 - 8*q + O(q^2)]
        """

        return self.quasi_part_gens()

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(n=5, k=6, ep=-1)
            sage: MF.dimension()
            3
            sage: len(MF.gens()) == MF.dimension()
            True
        """

        return self.quasi_part_dimension()

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.

        INPUT:

        - ``v`` -- an element of ``self``

        OUTPUT:

        An element of ``self.module()``, namely the
        corresponding coordinate vector of ``v`` with respect
        to the basis ``self.gens()``.

        The module is the free module over the coefficient
        ring of ``self`` with the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            22
            sage: el = MF(MF.E4()^2*MF.E6()^2 + MF.E4()*MF.E2()^2*MF.Delta() + MF.E2()^3*MF.E4()^2*MF.E6())
            sage: el
            2 + 25*q - 2478*q^2 - 82731*q^3 - 448484*q^4 + O(q^5)
            sage: vec = el.coordinate_vector()    # long time
            sage: vec    # long time
            (1, 1/(9*d), -11/(81*d^2), -4499/(104976*d^3), 0, 0, 0, 0, 1, 1/(2*d), 1, 5/(18*d), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: vec.parent()    # long time
            Vector space of dimension 22 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec.parent() == MF.module()    # long time
            True
            sage: el == MF(sum([vec[l]*MF.gen(l) for l in range(0,22)]))    # long time
            True
            sage: el == MF.element_from_coordinates(vec)    # long time
            True
            sage: MF.gen(1).coordinate_vector() == vector([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])    # long time
            True

            sage: MF = QuasiModularForms(n=infinity, k=4, ep=1)
            sage: el2 = MF.E4() + MF.E2()^2
            sage: el2
            2 + 160*q^2 + 512*q^3 + 1632*q^4 + O(q^5)
            sage: el2.coordinate_vector()
            (1, 1/(4*d), 0, 1)
            sage: el2 == MF.element_from_coordinates(el2.coordinate_vector())
            True
        """
        x, y, z, d = self.pol_ring().gens()
        k = self._weight
        rmax = (QQ(k) / 2).floor()
        partlist = v.rat().numerator().polynomial(z).list()
        denom = self.coeff_ring()(v.rat().denominator())
        partlist = [part / denom for part in partlist]
        parts = partlist + [0] * (rmax + 1 - len(partlist))
        E2 = self.E2()
        coord_vector = []

        for r in range(rmax + 1):
            gens = [v / E2**r for v in self.quasi_part_gens(r)]

            if gens:
                ambient_space = self.graded_ring().reduce_type("holo", degree=(gens[0].weight(), gens[0].ep()))
                subspace = ambient_space.subspace(gens)
                vector_part_in_subspace = subspace(parts[r])
                coord_part = list(vector_part_in_subspace.coordinate_vector())
                coord_vector += coord_part

        return self._module(vector(self.coeff_ring(), coord_vector))


class QuasiCuspForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi cusp forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiCuspForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(8, ZZ, 16/3, None)
            sage: QuasiCuspForms(8, ZZ, 16/3) == QuasiCuspForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi cusp forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(8, ZZ, 16/3)
            sage: MF
            QuasiCuspForms(n=8, k=16/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi cuspidal
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.is_ambient()
            True

            sage: QuasiCuspForms(n=infinity)
            QuasiCuspForms(n=+Infinity, k=0, ep=1) over Integer Ring
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["quasi", "cusp"])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(n=8, k=46/3, ep=-1)
            sage: MF.default_prec(4)
            sage: MF.dimension()
            7
            sage: MF.gens()
            [q - 17535/(262144*d^2)*q^3 + O(q^4),
             q^2 - 47/(128*d)*q^3 + O(q^4),
             q - 9/(128*d)*q^2 + 15633/(262144*d^2)*q^3 + O(q^4),
             q^2 - 7/(128*d)*q^3 + O(q^4),
             q - 23/(64*d)*q^2 - 3103/(262144*d^2)*q^3 + O(q^4),
             q - 3/(64*d)*q^2 - 4863/(262144*d^2)*q^3 + O(q^4),
             q - 27/(64*d)*q^2 + 17217/(262144*d^2)*q^3 + O(q^4)]

            sage: MF = QuasiCuspForms(n=infinity, k=10, ep=-1)
            sage: MF.gens()
            [q - 16*q^2 - 156*q^3 - 256*q^4 + O(q^5), q - 60*q^3 - 256*q^4 + O(q^5)]
        """

        return self.quasi_part_gens()

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(n=8, k=46/3, ep=-1)
            sage: MF.default_prec(3)
            sage: MF.dimension()
            7
            sage: len(MF.gens()) == MF.dimension()
            True

            sage: QuasiCuspForms(n=infinity, k=10, ep=-1).dimension()
            2
        """

        return self.quasi_part_dimension()

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.

        INPUT:

        - ``v`` -- an element of ``self``

        OUTPUT:

        An element of ``self.module()``, namely the
        corresponding coordinate vector of ``v`` with respect
        to the basis ``self.gens()``.

        The module is the free module over the coefficient
        ring of ``self`` with the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            12
            sage: el = MF(MF.E4()^2*MF.Delta() + MF.E4()*MF.E2()^2*MF.Delta())
            sage: el
            2*q + 120*q^2 + 3402*q^3 + 61520*q^4 + O(q^5)
            sage: vec = el.coordinate_vector()    # long time
            sage: vec    # long time
            (1, 13/(18*d), 103/(432*d^2), 0, 0, 1, 1/(2*d), 0, 0, 0, 0, 0)
            sage: vec.parent()    # long time
            Vector space of dimension 12 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec.parent() == MF.module()    # long time
            True
            sage: el == MF(sum([vec[l]*MF.gen(l) for l in range(0,12)]))    # long time
            True
            sage: el == MF.element_from_coordinates(vec)    # long time
            True
            sage: MF.gen(1).coordinate_vector() == vector([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])    # long time
            True

            sage: MF = QuasiCuspForms(n=infinity, k=10, ep=-1)
            sage: el2 = MF(MF.E4()*MF.f_inf()*(MF.f_i() - MF.E2()))
            sage: el2.coordinate_vector()
            (1, -1)
            sage: el2 == MF.element_from_coordinates(el2.coordinate_vector())
            True
        """
        x, y, z, d = self.pol_ring().gens()
        k = self._weight
        rmax = (QQ(k) / 2).floor()
        partlist = v.rat().numerator().polynomial(z).list()
        denom = self.coeff_ring()(v.rat().denominator())
        partlist = [part / denom for part in partlist]
        parts = partlist + [0] * (rmax + 1 - len(partlist))
        E2 = self.E2()
        coord_vector = []

        for r in range(rmax + 1):
            gens = [v / E2**r for v in self.quasi_part_gens(r)]

            if gens:
                ambient_space = self.graded_ring().reduce_type("cusp", degree=(gens[0].weight(), gens[0].ep()))
                subspace = ambient_space.subspace(gens)
                vector_part_in_subspace = subspace(parts[r])
                coord_part = list(vector_part_in_subspace.coordinate_vector())
                coord_vector += coord_part

        return self._module(vector(self.coeff_ring(), coord_vector))


class MeromorphicModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) meromorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, MeromorphicModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(3, ZZ, 0, 1)
            sage: MeromorphicModularForms() == MeromorphicModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) meromorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import MeromorphicModularForms
            sage: MF = MeromorphicModularForms()
            sage: MF
            MeromorphicModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF.analytic_type()
            meromorphic modular
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["mero"])


class WeakModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) weakly holomorphic modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, WeakModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, CC, 20/3, None)
            sage: WeakModularForms(5, CC, 20/3) == WeakModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) weakly holomorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(5, CC, 20/3)
            sage: MF
            WeakModularForms(n=5, k=20/3, ep=1) over Complex Field with 53 bits of precision
            sage: MF.analytic_type()
            weakly holomorphic modular
            sage: MF.category()
            Category of vector spaces over Complex Field with 53 bits of precision
            sage: MF in MF.category()
            True
        """
        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["weak"])


class ModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) modular forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, ModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(3, ZZ, 0, None)
            sage: ModularForms() == ModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms()
            sage: MF
            ModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF.analytic_type()
            modular
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.module()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True

            sage: MF = ModularForms(n=infinity, k=8)
            sage: MF
            ModularForms(n=+Infinity, k=8, ep=1) over Integer Ring
            sage: MF.analytic_type()
            modular
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["holo"])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            4
            sage: MF.gens()
            [1 + 360360*q^4 + O(q^5),
             q + 21742*q^4 + O(q^5),
             q^2 + 702*q^4 + O(q^5),
             q^3 - 6*q^4 + O(q^5)]

            sage: ModularForms(n=infinity, k=4).gens()
            [1 + 240*q^2 + 2160*q^4 + O(q^5), q - 8*q^2 + 28*q^3 - 64*q^4 + O(q^5)]
        """
        return [self.F_basis(m) for m in range(self.dimension())]

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            4
            sage: len(MF.gens()) == MF.dimension()
            True

            sage: ModularForms(n=infinity, k=8).dimension()
            3
        """
        return max(self._l1 + 1, ZZ.zero())

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.

        INPUT:

        - ``v`` -- an element of ``self``

        OUTPUT:

        An element of ``self.module()``, namely the
        corresponding coordinate vector of ``v`` with respect
        to the basis ``self.gens()``.

        The module is the free module over the coefficient
        ring of ``self`` with the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            4
            sage: el = MF.E4()^2*MF.Delta()
            sage: el
            q + 78*q^2 + 2781*q^3 + 59812*q^4 + O(q^5)
            sage: vec = el.coordinate_vector()
            sage: vec
            (0, 1, 13/(18*d), 103/(432*d^2))
            sage: vec.parent()
            Vector space of dimension 4 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec.parent() == MF.module()
            True
            sage: el == vec[0]*MF.gen(0) + vec[1]*MF.gen(1) + vec[2]*MF.gen(2) + vec[3]*MF.gen(3)
            True
            sage: el == MF.element_from_coordinates(vec)
            True

            sage: MF = ModularForms(n=infinity, k=8, ep=1)
            sage: (MF.E4()^2).coordinate_vector()
            (1, 1/(2*d), 15/(128*d^2))
        """

        vec = v.q_expansion_vector(min_exp=0, max_exp=self.degree() - 1)
        return self._module(vec)


class CuspForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) cusp forms
    for the given group, base ring, weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, CuspForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(6, ZZ, 6, 1)
            sage: CuspForms(6, ZZ, 6, 1) == CuspForms(group, base_ring, k, ep, n)
            True
        """
        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) cusp forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(6, ZZ, 6, 1)
            sage: MF
            CuspForms(n=6, k=6, ep=1) over Integer Ring
            sage: MF.analytic_type()
            cuspidal
            sage: MF.category()
            Category of modules over Integer Ring
            sage: MF in MF.category()
            True
            sage: MF.module()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True
        """
        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT(["cusp"])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF=CuspForms(n=12, k=72/5, ep=1)
            sage: MF
            CuspForms(n=12, k=72/5, ep=1) over Integer Ring
            sage: MF.dimension()
            3
            sage: MF.gens()
            [q + 296888795/(10319560704*d^3)*q^4 + O(q^5),
             q^2 + 6629/(221184*d^2)*q^4 + O(q^5),
             q^3 - 25/(96*d)*q^4 + O(q^5)]

            sage: MF = CuspForms(n=infinity, k=8, ep=1)
            sage: MF.gen(0) == MF.E4()*MF.f_inf()
            True
        """
        return [self.F_basis(m, order_1=ZZ.one())
                for m in range(1, self.dimension() + 1)]

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(n=12, k=72/5, ep=1)
            sage: MF.dimension()
            3
            sage: len(MF.gens()) == MF.dimension()
            True

            sage: CuspForms(n=infinity, k=8).dimension()
            1
        """
        if self.hecke_n() == infinity:
            return max(self._l1 - 1, ZZ.zero())
        return max(self._l1, ZZ.zero())

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.

        INPUT:

        - ``v`` -- an element of ``self``

        OUTPUT:

        An element of ``self.module()``, namely the
        corresponding coordinate vector of ``v`` with respect
        to the basis ``self.gens()``.

        The module is the free module over the coefficient
        ring of ``self`` with the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(n=12, k=72/5, ep=-1)
            sage: MF.default_prec(4)
            sage: MF.dimension()
            2
            sage: el = MF(MF.f_i()*MF.Delta())
            sage: el
            q - 1/(288*d)*q^2 - 96605/(1327104*d^2)*q^3 + O(q^4)
            sage: vec = el.coordinate_vector()
            sage: vec
            (1, -1/(288*d))
            sage: vec.parent()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec.parent() == MF.module()
            True
            sage: el == vec[0]*MF.gen(0) + vec[1]*MF.gen(1)
            True
            sage: el == MF.element_from_coordinates(vec)
            True

            sage: MF = CuspForms(n=infinity, k=16)
            sage: el2 = MF(MF.Delta()*MF.E4())
            sage: vec2 = el2.coordinate_vector()
            sage: vec2
            (1, 5/(8*d), 187/(1024*d^2))
            sage: el2 == MF.element_from_coordinates(vec2)
            True
        """

        vec = v.q_expansion_vector(min_exp=1, max_exp=self.degree())
        return self._module(vec)


class ZeroForm(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Zero Module for the zero form for the given group, base ring
    weight and multiplier
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, ZeroForm)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(6, CC, 3, -1)
            sage: ZeroForm(6, CC, 3, -1) == ZeroForm(group, base_ring, k, ep, n)
            True
        """
        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the zero Module for the zero form of weight ``k`` with multiplier ``ep``
        for the given ``group`` and ``base_ring``.

        The ZeroForm space coerces into any forms space or ring with a compatible group.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: MF = ZeroForm(6, CC, 3, -1)
            sage: MF
            ZeroForms(n=6, k=3, ep=-1) over Complex Field with 53 bits of precision
            sage: MF.analytic_type()
            zero
            sage: MF.category()
            Category of vector spaces over Complex Field with 53 bits of precision
            sage: MF in MF.category()
            True
            sage: MF.module()
            Vector space of dimension 0 over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=base_ring)
        self._analytic_type = self.AT([])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    def _change_degree(self, k, ep):
        r"""
        Return the Zeroform space with a different weight ``k`` and multiplier ``ep``
        for the same group and base_ring.

        All those spaces coerce into each other.

        INPUT:

        - ``k`` -- a rational number; the weight

        - ``ep`` -- ``1`` or ``-1``; the multiplier

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: MF = ZeroForm()
            sage: MF
            ZeroForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF._change_degree(14, -1)
            ZeroForms(n=3, k=14, ep=-1) over Integer Ring
        """
        return ZeroForm(group=self.group(), base_ring=self.base_ring(),
                        k=k, ep=ep)

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.
        Since this is the zero module an empty list is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: ZeroForm(6, CC, 3, -1).gens()
            []
        """
        return []

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.
        Since this is the zero module ``0`` is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: ZeroForm(6, CC, 3, -1).dimension()
            0
        """
        return 0

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.

        Since this is the zero module which only contains
        the zero form the trivial vector in the trivial
        module of dimension ``0`` is returned.

        INPUT:

        - ``v`` -- an element of ``self``, i.e. in this case the zero vector

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: MF = ZeroForm(6, QQ, 3, -1)
            sage: el = MF(0)
            sage: el
            O(q^5)
            sage: vec = el.coordinate_vector()
            sage: vec
            ()
            sage: vec.parent()
            Vector space of dimension 0 over Fraction Field of Univariate Polynomial Ring in d over Rational Field
            sage: vec.parent() == MF.module()
            True
        """
        return self._module(vector(self.coeff_ring(), []))
