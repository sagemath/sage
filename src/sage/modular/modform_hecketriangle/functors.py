r"""
Functor construction for all spaces

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

from sage.categories.functor                     import Functor
from sage.categories.pushout                     import ConstructionFunctor
from sage.categories.sets_cat                    import Sets
from sage.structure.parent                       import Parent
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.rings                       import Rings

from .constructor                                 import FormsSpace, FormsRing
from .abstract_space                              import FormsSpace_abstract
from .subspace                                    import SubSpaceForms


def _get_base_ring(ring, var_name="d"):
    r"""
    Return the base ring of the given ``ring``:

    If ``ring`` is of the form ``FractionField(PolynomialRing(R,'d'))``:
    Return ``R``.

    If ``ring`` is of the form ``FractionField(R)``:
    Return ``R``.

    If ``ring`` is of the form ``PolynomialRing(R,'d')``:
    Return ``R``.

    Otherwise return ``ring``.

    The base ring is used in the construction of the corresponding
    ``FormsRing`` or ``FormsSpace``. In particular in the construction
    of holomorphic forms of degree (0, 1). For (binary)
    operations a general ring element is considered (coerced to)
    a (constant) holomorphic form of degree (0, 1)
    whose construction should be based on the returned base ring
    (and not on ``ring``!).

    If ``var_name`` (default: "d") is specified then this variable
    name is used for the polynomial ring.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.functors import _get_base_ring
        sage: _get_base_ring(ZZ) == ZZ
        True
        sage: _get_base_ring(QQ) == ZZ
        True
        sage: _get_base_ring(PolynomialRing(CC, 'd')) == CC
        True
        sage: _get_base_ring(PolynomialRing(QQ, 'd')) == ZZ
        True
        sage: _get_base_ring(FractionField(PolynomialRing(CC, 'd'))) == CC
        True
        sage: _get_base_ring(FractionField(PolynomialRing(QQ, 'd'))) == ZZ
        True
        sage: _get_base_ring(PolynomialRing(QQ, 'x')) == PolynomialRing(QQ, 'x')
        True
    """

    #from sage.rings.fraction_field import is_FractionField
    from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
    from sage.categories.pushout import FractionField as FractionFieldFunctor

    base_ring = ring
    #if (is_FractionField(base_ring)):
    #    base_ring = base_ring.base()
    if (base_ring.construction() and base_ring.construction()[0] == FractionFieldFunctor()):
        base_ring = base_ring.construction()[1]
    if (is_PolynomialRing(base_ring) and base_ring.ngens()==1 and base_ring.variable_name()==var_name):
        base_ring = base_ring.base()
    if (base_ring.construction() and base_ring.construction()[0] == FractionFieldFunctor()):
        base_ring = base_ring.construction()[1]

    return base_ring


def _common_subgroup(group1, group2):
    r"""
    Return a common (Hecke triangle) subgroup of both given groups
    ``group1`` and ``group2`` if it exists. Otherwise return ``None``.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.functors import _common_subgroup
        sage: from sage.modular.modform_hecketriangle.hecke_triangle_groups import HeckeTriangleGroup
        sage: _common_subgroup(HeckeTriangleGroup(n=3), HeckeTriangleGroup(n=infinity))
        Hecke triangle group for n = +Infinity
        sage: _common_subgroup(HeckeTriangleGroup(n=infinity), HeckeTriangleGroup(n=3))
        Hecke triangle group for n = +Infinity
        sage: _common_subgroup(HeckeTriangleGroup(n=4), HeckeTriangleGroup(n=infinity)) is None
        True
        sage: _common_subgroup(HeckeTriangleGroup(n=4), HeckeTriangleGroup(n=4))
        Hecke triangle group for n = 4
    """

    if group1 == group2:
        return group1
    elif (group1.n() == 3) and (group2.n() == infinity):
        return group2
    elif (group1.n() == infinity) and (group2.n() == 3):
        return group1
    else:
        return None


def ConstantFormsSpaceFunctor(group):
    r"""
    Construction functor for the space of constant forms.

    When determining a common parent between a ring
    and a forms ring or space this functor is first
    applied to the ring.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.functors import (ConstantFormsSpaceFunctor, FormsSpaceFunctor)
        sage: ConstantFormsSpaceFunctor(4) == FormsSpaceFunctor("holo", 4, 0, 1)
        True
        sage: ConstantFormsSpaceFunctor(4)
        ModularFormsFunctor(n=4, k=0, ep=1)
    """
    return FormsSpaceFunctor("holo", group, QQ.zero(), ZZ.one())


class FormsSubSpaceFunctor(ConstructionFunctor):
    r"""
    Construction functor for forms sub spaces.
    """

    rank = 10

    def __init__(self, ambient_space_functor, generators):
        r"""
        Construction functor for the forms sub space
        for the given ``generators`` inside the ambient space
        which is constructed by the ``ambient_space_functor``.

        The functor can only be applied to rings for which the generators
        can be converted into the corresponding forms space
        given by the ``ambient_space_functor`` applied to the ring.

        See :meth:`__call__` for a description of the functor.

        INPUT:

        - ``ambient_space_functor`` -- A FormsSpaceFunctor

        - ``generators``            -- A list of elements of some ambient space
                                       over some base ring.

        OUTPUT:

        The construction functor for the corresponding forms sub space.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsSubSpaceFunctor)
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ambient_space = ModularForms(n=4, k=12, ep=1)
            sage: ambient_space_functor = FormsSpaceFunctor("holo", group=4, k=12, ep=1)
            sage: ambient_space_functor
            ModularFormsFunctor(n=4, k=12, ep=1)
            sage: el = ambient_space.gen(0).full_reduce()
            sage: FormsSubSpaceFunctor(ambient_space_functor, [el])
            FormsSubSpaceFunctor with 1 generator for the ModularFormsFunctor(n=4, k=12, ep=1)
        """

        Functor.__init__(self, Rings(), CommutativeAdditiveGroups())
        if not isinstance(ambient_space_functor, FormsSpaceFunctor):
            raise ValueError("{} is not a FormsSpaceFunctor!".format(ambient_space_functor))
        # TODO: canonical parameters? Some checks?
        # The generators should have an associated base ring
        # self._generators_ring = ...
        # on call check if there is a coercion from self._generators_ring to R

        self._ambient_space_functor = ambient_space_functor
        self._generators = generators

    def __call__(self, R):
        r"""
        Return the corresponding subspace of the ambient space
        constructed by ``self._ambient_space`` with the generators ``self._generators``.
        If the ambient space is not a forms space the ambient space is returned.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsSubSpaceFunctor, BaseFacade)
            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: ambient_space = CuspForms(n=4, k=12, ep=1)
            sage: ambient_space_functor = FormsSpaceFunctor("cusp", group=4, k=12, ep=1)
            sage: el = ambient_space.gen(0)
            sage: F = FormsSubSpaceFunctor(ambient_space_functor, [el])
            sage: F
            FormsSubSpaceFunctor with 1 generator for the CuspFormsFunctor(n=4, k=12, ep=1)

            sage: F(BaseFacade(ZZ))
            Subspace of dimension 1 of CuspForms(n=4, k=12, ep=1) over Integer Ring
            sage: F(BaseFacade(QQ))
            Subspace of dimension 1 of CuspForms(n=4, k=12, ep=1) over Integer Ring
            sage: F(QQ)
            ModularFormsRing(n=4) over Integer Ring

            sage: ambient_space_functor = FormsSpaceFunctor("holo", group=4, k=0, ep=1)
            sage: F = FormsSubSpaceFunctor(ambient_space_functor, [1])
            sage: F
            FormsSubSpaceFunctor with 1 generator for the ModularFormsFunctor(n=4, k=0, ep=1)
            sage: F(BaseFacade(ZZ))
            Subspace of dimension 1 of ModularForms(n=4, k=0, ep=1) over Integer Ring
        """

        ambient_space = self._ambient_space_functor(R)
        if isinstance(ambient_space, FormsSpace_abstract):
            return SubSpaceForms(ambient_space, self._generators)
        else:
            return ambient_space

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsSubSpaceFunctor)
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ambient_space = ModularForms(n=4, k=12, ep=1)
            sage: ambient_space_functor = FormsSpaceFunctor("holo", group=4, k=12, ep=1)
            sage: FormsSubSpaceFunctor(ambient_space_functor, ambient_space.gens())
            FormsSubSpaceFunctor with 2 generators for the ModularFormsFunctor(n=4, k=12, ep=1)
            sage: FormsSubSpaceFunctor(ambient_space_functor, [ambient_space.gen(0)])
            FormsSubSpaceFunctor with 1 generator for the ModularFormsFunctor(n=4, k=12, ep=1)
        """

        return "FormsSubSpaceFunctor with {} generator{} for the {}".format(len(self._generators), 's' if len(self._generators) != 1 else '', self._ambient_space_functor)

    def merge(self, other):
        r"""
        Return the merged functor of ``self`` and ``other``.

        If ``other`` is a ``FormsSubSpaceFunctor`` then
        first the common ambient space functor is constructed by
        merging the two corresponding functors.

        If that ambient space functor is a FormsSpaceFunctor
        and the generators agree the corresponding ``FormsSubSpaceFunctor``
        is returned.

        If ``other`` is not a ``FormsSubSpaceFunctor`` then ``self``
        is merged as if it was its ambient space functor.


        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsSubSpaceFunctor)
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ambient_space = ModularForms(n=4, k=12, ep=1)
            sage: ambient_space_functor1 = FormsSpaceFunctor("holo", group=4, k=12, ep=1)
            sage: ambient_space_functor2 = FormsSpaceFunctor("cusp", group=4, k=12, ep=1)
            sage: ss_functor1 = FormsSubSpaceFunctor(ambient_space_functor1, [ambient_space.gen(0)])
            sage: ss_functor2 = FormsSubSpaceFunctor(ambient_space_functor2, [ambient_space.gen(0)])
            sage: ss_functor3 = FormsSubSpaceFunctor(ambient_space_functor2, [2*ambient_space.gen(0)])
            sage: merged_ambient = ambient_space_functor1.merge(ambient_space_functor2)
            sage: merged_ambient
            ModularFormsFunctor(n=4, k=12, ep=1)
            sage: functor4 = FormsSpaceFunctor(["quasi", "cusp"], group=4, k=10, ep=-1)

            sage: ss_functor1.merge(ss_functor1) is ss_functor1
            True
            sage: ss_functor1.merge(ss_functor2)
            FormsSubSpaceFunctor with 2 generators for the ModularFormsFunctor(n=4, k=12, ep=1)
            sage: ss_functor1.merge(ss_functor2) == FormsSubSpaceFunctor(merged_ambient, [ambient_space.gen(0), ambient_space.gen(0)])
            True
            sage: ss_functor1.merge(ss_functor3) == FormsSubSpaceFunctor(merged_ambient, [ambient_space.gen(0), 2*ambient_space.gen(0)])
            True
            sage: ss_functor1.merge(ambient_space_functor2) == merged_ambient
            True
            sage: ss_functor1.merge(functor4)
            QuasiModularFormsRingFunctor(n=4, red_hom=True)
        """

        if (self == other):
            return self
        elif isinstance(other, FormsSubSpaceFunctor):
            merged_ambient_space_functor = self._ambient_space_functor.merge(other._ambient_space_functor)
            if isinstance(merged_ambient_space_functor, FormsSpaceFunctor):
                generators = self._generators + other._generators
                return FormsSubSpaceFunctor(merged_ambient_space_functor, generators)
            # This includes the case when None is returned
            else:
                return merged_ambient_space_functor
        else:
            return self._ambient_space_functor.merge(other)

    def __eq__(self, other):
        r"""
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsSubSpaceFunctor)
            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: ambient_space = ModularForms(n=4, k=12, ep=1)
            sage: ambient_space_functor1 = FormsSpaceFunctor("holo", group=4, k=12, ep=1)
            sage: ss_functor1 = FormsSubSpaceFunctor(ambient_space_functor1, [ambient_space.gen(0)])
            sage: ss_functor2 = FormsSubSpaceFunctor(ambient_space_functor1, [ambient_space.gen(1)])
            sage: ss_functor1 == ss_functor2
            False
        """

        if (type(self) is type(other) and
            self._ambient_space_functor == other._ambient_space_functor and
            self._generators == other._generators):
                return True
        else:
            return False


class FormsSpaceFunctor(ConstructionFunctor):
    r"""
    Construction functor for forms spaces.

    NOTE:

    When the base ring is not a ``BaseFacade`` the functor is first
    merged with the ConstantFormsSpaceFunctor.  This case occurs in
    the pushout constructions (when trying to find a common parent
    between a forms space and a ring which is not a ``BaseFacade``).
    """

    from .analytic_type import AnalyticType
    AT = AnalyticType()

    rank = 10

    def __init__(self, analytic_type, group, k, ep):
        r"""
        Construction functor for the forms space
        (or forms ring, see above) with
        the given ``analytic_type``, ``group``,
        weight ``k`` and multiplier ``ep``.

        See :meth:`__call__` for a description of the functor.

        INPUT:

        - ``analytic_type``  -- An element of ``AnalyticType()``.

        - ``group``          -- The index of a Hecke Triangle group.

        - ``k``              -- A rational number, the weight of the space.

        - ``ep``             -- `1` or `-1`, the multiplier of the space.

        OUTPUT:

        The construction functor for the corresponding forms space/ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import FormsSpaceFunctor
            sage: FormsSpaceFunctor(["holo", "weak"], group=4, k=0, ep=-1)
            WeakModularFormsFunctor(n=4, k=0, ep=-1)
        """

        Functor.__init__(self, Rings(), CommutativeAdditiveGroups())
        from .space import canonical_parameters
        (self._group, R, self._k, self._ep, n) = canonical_parameters(group, ZZ, k, ep)

        self._analytic_type = self.AT(analytic_type)

    def __call__(self, R):
        r"""
        If ``R`` is a ``BaseFacade(S)`` then return the corresponding
        forms space with base ring ``_get_base_ring(S)``.

        If not then we first merge the functor with the ConstantFormsSpaceFunctor.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, BaseFacade)
            sage: F = FormsSpaceFunctor(["holo", "weak"], group=4, k=0, ep=-1)
            sage: F(BaseFacade(ZZ))
            WeakModularForms(n=4, k=0, ep=-1) over Integer Ring
            sage: F(BaseFacade(CC))
            WeakModularForms(n=4, k=0, ep=-1) over Complex Field with 53 bits of precision
            sage: F(CC)
            WeakModularFormsRing(n=4) over Complex Field with 53 bits of precision
            sage: F(CC).has_reduce_hom()
            True
        """

        if (isinstance(R, BaseFacade)):
            R = _get_base_ring(R._ring)
            return FormsSpace(self._analytic_type, self._group, R, self._k, self._ep)
        else:
            R = BaseFacade(_get_base_ring(R))
            merged_functor = self.merge(ConstantFormsSpaceFunctor(self._group))
            return merged_functor(R)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import FormsSpaceFunctor
            sage: F = FormsSpaceFunctor(["cusp", "quasi"], group=5, k=10/3, ep=-1)
            sage: str(F)
            'QuasiCuspFormsFunctor(n=5, k=10/3, ep=-1)'
            sage: F
            QuasiCuspFormsFunctor(n=5, k=10/3, ep=-1)
        """

        return "{}FormsFunctor(n={}, k={}, ep={})".format(self._analytic_type.analytic_space_name(), self._group.n(), self._k, self._ep)

    def merge(self, other):
        r"""
        Return the merged functor of ``self`` and ``other``.

        It is only possible to merge instances of ``FormsSpaceFunctor``
        and ``FormsRingFunctor``. Also only if they share the same group.
        An ``FormsSubSpaceFunctors`` is replaced by its ambient space functor.

        The analytic type of the merged functor is the extension
        of the two analytic types of the functors.
        The ``red_hom`` parameter of the merged functor
        is the logical ``and`` of the two corresponding ``red_hom``
        parameters (where a forms space is assumed to have it
        set to ``True``).

        Two ``FormsSpaceFunctor`` with different (k,ep) are merged to a
        corresponding ``FormsRingFunctor``. Otherwise the corresponding
        (extended) ``FormsSpaceFunctor`` is returned.

        A ``FormsSpaceFunctor`` and ``FormsRingFunctor``
        are merged to a corresponding (extended) ``FormsRingFunctor``.

        Two ``FormsRingFunctors`` are merged to the corresponding
        (extended) ``FormsRingFunctor``.


        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsRingFunctor)
            sage: functor1 = FormsSpaceFunctor("holo", group=5, k=0, ep=1)
            sage: functor2 = FormsSpaceFunctor(["quasi", "cusp"], group=5, k=10/3, ep=-1)
            sage: functor3 = FormsSpaceFunctor(["quasi", "mero"], group=5, k=0, ep=1)
            sage: functor4 = FormsRingFunctor("cusp", group=5, red_hom=False)
            sage: functor5 = FormsSpaceFunctor("holo", group=4, k=0, ep=1)

            sage: functor1.merge(functor1) is functor1
            True
            sage: functor1.merge(functor5) is None
            True
            sage: functor1.merge(functor2)
            QuasiModularFormsRingFunctor(n=5, red_hom=True)
            sage: functor1.merge(functor3)
            QuasiMeromorphicModularFormsFunctor(n=5, k=0, ep=1)
            sage: functor1.merge(functor4)
            ModularFormsRingFunctor(n=5)
        """

        if (self == other):
            return self

        if isinstance(other, FormsSubSpaceFunctor):
            other = other._ambient_space_functor

        if isinstance(other, FormsSpaceFunctor):
            group = _common_subgroup(self._group, other._group)
            if group is None:
                return None
            analytic_type = self._analytic_type + other._analytic_type
            if (self._k == other._k) and (self._ep == other._ep):
                return FormsSpaceFunctor(analytic_type, group, self._k, self._ep)
            else:
                return FormsRingFunctor(analytic_type, group, True)
        elif isinstance(other, FormsRingFunctor):
            group = _common_subgroup(self._group, other._group)
            if group is None:
                return None
            red_hom = other._red_hom
            analytic_type = self._analytic_type + other._analytic_type
            return FormsRingFunctor(analytic_type, group, red_hom)

    def __eq__(self, other):
        r"""
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import FormsSpaceFunctor
            sage: functor1 = FormsSpaceFunctor("holo", group=4, k=12, ep=1)
            sage: functor2 = FormsSpaceFunctor("holo", group=4, k=12, ep=-1)
            sage: functor1 == functor2
            False
        """

        if (type(self) is type(other) and
            self._group == other._group and
            self._analytic_type == other._analytic_type and
            self._k == other._k and
            self._ep == other._ep):
                return True
        else:
            return False


class FormsRingFunctor(ConstructionFunctor):
    r"""
    Construction functor for forms rings.

    NOTE:

    When the base ring is not a ``BaseFacade`` the functor is first
    merged with the ConstantFormsSpaceFunctor.  This case occurs in
    the pushout constructions.  (when trying to find a common parent
    between a forms ring and a ring which is not a ``BaseFacade``).
    """

    from .analytic_type import AnalyticType
    AT = AnalyticType()

    rank = 10

    def __init__(self, analytic_type, group, red_hom):
        r"""
        Construction functor for the forms ring
        with the given ``analytic_type``, ``group``
        and variable ``red_hom``

        See :meth:`__call__` for a description of the functor.

        INPUT:

        - ``analytic_type``  -- An element of ``AnalyticType()``.

        - ``group``          -- The index of a Hecke Triangle group.

        - ``red_hom``        -- A boolean variable for the parameter ``red_hom``
                                (also see ``FormsRing_abstract``).

        OUTPUT:

        The construction functor for the corresponding forms ring.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import FormsRingFunctor
            sage: FormsRingFunctor(["quasi", "mero"], group=6, red_hom=False)
            QuasiMeromorphicModularFormsRingFunctor(n=6)
            sage: FormsRingFunctor(["quasi", "mero"], group=6, red_hom=True)
            QuasiMeromorphicModularFormsRingFunctor(n=6, red_hom=True)
        """

        Functor.__init__(self, Rings(), Rings())
        from .graded_ring import canonical_parameters
        (self._group, R, red_hom, n) = canonical_parameters(group, ZZ, red_hom)
        self._red_hom = bool(red_hom)
        self._analytic_type = self.AT(analytic_type)

    def __call__(self, R):
        r"""
        If ``R`` is a ``BaseFacade(S)`` then return the corresponding
        forms ring with base ring ``_get_base_ring(S)``.

        If not then we first merge the functor with the ConstantFormsSpaceFunctor.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsRingFunctor, BaseFacade)
            sage: F = FormsRingFunctor(["quasi", "mero"], group=6, red_hom=False)
            sage: F(BaseFacade(ZZ))
            QuasiMeromorphicModularFormsRing(n=6) over Integer Ring
            sage: F(BaseFacade(CC))
            QuasiMeromorphicModularFormsRing(n=6) over Complex Field with 53 bits of precision
            sage: F(CC)
            QuasiMeromorphicModularFormsRing(n=6) over Complex Field with 53 bits of precision
            sage: F(CC).has_reduce_hom()
            False
        """

        if (isinstance(R, BaseFacade)):
            R = _get_base_ring(R._ring)
            return FormsRing(self._analytic_type, self._group, R, self._red_hom)
        else:
            R = BaseFacade(_get_base_ring(R))
            merged_functor = self.merge(ConstantFormsSpaceFunctor(self._group))
            return merged_functor(R)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import FormsRingFunctor
            sage: str(FormsRingFunctor(["quasi", "mero"], group=6, red_hom=True))
            'QuasiMeromorphicModularFormsRingFunctor(n=6, red_hom=True)'
            sage: FormsRingFunctor(["quasi", "mero"], group=6, red_hom=False)
            QuasiMeromorphicModularFormsRingFunctor(n=6)
        """

        if (self._red_hom):
            red_arg = ", red_hom=True"
        else:
            red_arg = ""
        return "{}FormsRingFunctor(n={}{})".format(self._analytic_type.analytic_space_name(), self._group.n(), red_arg)

    def merge(self, other):
        r"""
        Return the merged functor of ``self`` and ``other``.

        It is only possible to merge instances of ``FormsSpaceFunctor``
        and ``FormsRingFunctor``. Also only if they share the same group.
        An ``FormsSubSpaceFunctors`` is replaced by its ambient space functor.

        The analytic type of the merged functor is the extension
        of the two analytic types of the functors.
        The ``red_hom`` parameter of the merged functor
        is the logical ``and`` of the two corresponding ``red_hom``
        parameters (where a forms space is assumed to have it
        set to ``True``).

        Two ``FormsSpaceFunctor`` with different (k,ep) are merged to a
        corresponding ``FormsRingFunctor``. Otherwise the corresponding
        (extended) ``FormsSpaceFunctor`` is returned.

        A ``FormsSpaceFunctor`` and ``FormsRingFunctor``
        are merged to a corresponding (extended) ``FormsRingFunctor``.

        Two ``FormsRingFunctors`` are merged to the corresponding
        (extended) ``FormsRingFunctor``.


        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import (FormsSpaceFunctor, FormsRingFunctor)
            sage: functor1 = FormsRingFunctor("mero", group=6, red_hom=True)
            sage: functor2 = FormsRingFunctor(["quasi", "cusp"], group=6, red_hom=False)
            sage: functor3 = FormsSpaceFunctor("weak", group=6, k=0, ep=1)
            sage: functor4 = FormsRingFunctor("mero", group=5, red_hom=True)

            sage: functor1.merge(functor1) is functor1
            True
            sage: functor1.merge(functor4) is None
            True
            sage: functor1.merge(functor2)
            QuasiMeromorphicModularFormsRingFunctor(n=6)
            sage: functor1.merge(functor3)
            MeromorphicModularFormsRingFunctor(n=6, red_hom=True)
        """

        if (self == other):
            return self

        if isinstance(other, FormsSubSpaceFunctor):
            other = other._ambient_space_functor

        if isinstance(other, FormsSpaceFunctor):
            group = _common_subgroup(self._group, other._group)
            if group is None:
                return None
            red_hom = self._red_hom
            analytic_type = self._analytic_type + other._analytic_type
            return FormsRingFunctor(analytic_type, group, red_hom)
        elif isinstance(other, FormsRingFunctor):
            group = _common_subgroup(self._group, other._group)
            if group is None:
                return None
            red_hom = self._red_hom & other._red_hom
            analytic_type = self._analytic_type + other._analytic_type
            return FormsRingFunctor(analytic_type, group, red_hom)

    def __eq__(self, other):
        r"""
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import FormsRingFunctor
            sage: functor1 = FormsRingFunctor("holo", group=4, red_hom=True)
            sage: functor2 = FormsRingFunctor("holo", group=4, red_hom=False)
            sage: functor1 == functor2
            False
        """

        if (type(self) is type(other) and
            self._group == other._group and
            self._analytic_type == other._analytic_type and
            self._red_hom == other._red_hom):
                return True
        else:
            return False


from sage.structure.unique_representation import UniqueRepresentation
class BaseFacade(Parent, UniqueRepresentation):
    r"""
    BaseFacade of a ring.

    This class is used to distinguish the construction of
    constant elements (modular forms of weight 0) over the given ring
    and the construction of ``FormsRing`` or ``FormsSpace``
    based on the BaseFacade of the given ring.

    If that distinction was not made then ring elements
    couldn't be considered as constant modular forms
    in e.g. binary operations. Instead the coercion model would
    assume that the ring element lies in the common parent
    of the ring element and e.g. a ``FormsSpace`` which
    would give the ``FormsSpace`` over the ring. However
    this is not correct, the ``FormsSpace`` might
    (and probably will) not even contain the (constant)
    ring element. Hence we use the ``BaseFacade`` to
    distinguish the two cases.

    Since the ``BaseFacade`` of a ring embeds into that ring,
    a common base (resp. a coercion) between the two (or even a
    more general ring) can be found, namely the ring
    (not the ``BaseFacade`` of it).
    """

    def __init__(self, ring):
        r"""
        BaseFacade of ``ring`` (see above).

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import BaseFacade
            sage: BaseFacade(ZZ)
            BaseFacade(Integer Ring)
            sage: ZZ.has_coerce_map_from(BaseFacade(ZZ))
            True
            sage: CC.has_coerce_map_from(BaseFacade(ZZ))
            True
        """

        Parent.__init__(self, facade=ring, category=Rings())
        self._ring = _get_base_ring(ring)
        # The BaseFacade(R) coerces/embeds into R, used in pushout
        self.register_embedding(self.Hom(self._ring,Sets())(lambda x: x))

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.functors import BaseFacade
            sage: BaseFacade(ZZ)
            BaseFacade(Integer Ring)
        """

        return "BaseFacade({})".format(self._ring)
