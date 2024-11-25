# sage.doctest: needs sage.combinat sage.graphs
r"""
Graded rings of modular forms for Hecke triangle groups

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
from sage.rings.infinity import infinity

from sage.structure.parent import Parent
from sage.categories.algebras import Algebras
from sage.structure.unique_representation import UniqueRepresentation

from .hecke_triangle_groups import HeckeTriangleGroup
from .abstract_ring import FormsRing_abstract


def canonical_parameters(group, base_ring, red_hom, n=None):
    r"""
    Return a canonical version of the parameters.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import canonical_parameters
        sage: canonical_parameters(4, ZZ, 1)
        (Hecke triangle group for n = 4, Integer Ring, True, 4)
        sage: canonical_parameters(infinity, RR, 0)
        (Hecke triangle group for n = +Infinity, Real Field with 53 bits of precision, False, +Infinity)
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

    red_hom = bool(red_hom)
    n = group.n()

    return (group, base_ring, red_hom, n)


class QuasiMeromorphicModularFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi meromorphic modular forms
    for the given group and base ring.
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiMeromorphicModularFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(4, ZZ, 1)
            sage: QuasiMeromorphicModularFormsRing(4, ZZ, 1) == QuasiMeromorphicModularFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) quasi meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(4, ZZ, 1)
            sage: MR
            QuasiMeromorphicModularFormsRing(n=4) over Integer Ring
            sage: MR.analytic_type()
            quasi meromorphic modular
            sage: MR.category()
            Category of commutative algebras over Integer Ring
            sage: MR in MR.category()
            True

            sage: QuasiMeromorphicModularFormsRing(n=infinity)
            QuasiMeromorphicModularFormsRing(n=+Infinity) over Integer Ring
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["quasi", "mero"])


class QuasiWeakModularFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi weakly holomorphic modular forms
    for the given group and base ring.
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiWeakModularFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(5, CC, 0)
            sage: QuasiWeakModularFormsRing(5, CC, 0) == QuasiWeakModularFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) quasi weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiWeakModularFormsRing
            sage: MR = QuasiWeakModularFormsRing(5, CC, 0)
            sage: MR
            QuasiWeakModularFormsRing(n=5) over Complex Field with 53 bits of precision
            sage: MR.analytic_type()
            quasi weakly holomorphic modular
            sage: MR.category()
            Category of commutative algebras over Complex Field with 53 bits of precision
            sage: MR in MR.category()
            True
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["quasi", "weak"])


class QuasiModularFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiModularFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(6, ZZ, True)
            sage: QuasiModularFormsRing(6, ZZ, True) == QuasiModularFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) quasi modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: MR = QuasiModularFormsRing(6, ZZ, True)
            sage: MR
            QuasiModularFormsRing(n=6) over Integer Ring
            sage: MR.analytic_type()
            quasi modular
            sage: MR.category()
            Category of commutative algebras over Integer Ring
            sage: MR in MR.category()
            True
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["quasi", "holo"])


class QuasiCuspFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi cusp forms
    for the given group and base ring.
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiCuspFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(7, ZZ, 1)
            sage: QuasiCuspFormsRing(7, ZZ, 1) == QuasiCuspFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) quasi cusp forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi cusp forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiCuspFormsRing
            sage: MR = QuasiCuspFormsRing(7, ZZ, 1)
            sage: MR
            QuasiCuspFormsRing(n=7) over Integer Ring
            sage: MR.analytic_type()
            quasi cuspidal
            sage: MR.category()
            Category of commutative algebras over Integer Ring
            sage: MR in MR.category()
            True
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["quasi", "cusp"])


class MeromorphicModularFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) meromorphic modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, MeromorphicModularFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(4, ZZ, 1)
            sage: MeromorphicModularFormsRing(4, ZZ, 1) == MeromorphicModularFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import MeromorphicModularFormsRing
            sage: MR = MeromorphicModularFormsRing(4, ZZ, 1)
            sage: MR
            MeromorphicModularFormsRing(n=4) over Integer Ring
            sage: MR.analytic_type()
            meromorphic modular
            sage: MR.category()
            Category of commutative algebras over Integer Ring
            sage: MR in MR.category()
            True
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["mero"])


class WeakModularFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) weakly holomorphic modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, WeakModularFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(5, ZZ, 0)
            sage: WeakModularFormsRing(5, ZZ, 0) == WeakModularFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import WeakModularFormsRing
            sage: MR = WeakModularFormsRing(5, ZZ, 0)
            sage: MR
            WeakModularFormsRing(n=5) over Integer Ring
            sage: MR.analytic_type()
            weakly holomorphic modular
            sage: MR.category()
            Category of commutative algebras over Integer Ring
            sage: MR in MR.category()
            True
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["weak"])


class ModularFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing(3, ZZ, 0) == ModularFormsRing()
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: MR = ModularFormsRing()
            sage: MR
            ModularFormsRing(n=3) over Integer Ring
            sage: MR.analytic_type()
            modular
            sage: MR.category()
            Category of commutative algebras over Integer Ring
            sage: MR in MR.category()
            True
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["holo"])


class CuspFormsRing(FormsRing_abstract, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) cusp forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group=HeckeTriangleGroup(3), base_ring=ZZ, red_hom=False, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, CuspFormsRing)
            sage: (group, base_ring, red_hom, n) = canonical_parameters(5, CC, True)
            sage: CuspFormsRing(5, CC, True) == CuspFormsRing(group, base_ring, red_hom, n)
            True
        """
        (group, base_ring, red_hom, n) = canonical_parameters(group, base_ring, red_hom, n)
        return super().__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n)

    def __init__(self, group, base_ring, red_hom, n):
        r"""
        Return the graded ring of (Hecke) cusp forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group`` -- the Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring`` -- the base_ring (default: ``ZZ``)

        - ``red_hom`` -- if ``True`` then results of binary operations are
          considered homogeneous whenever it makes sense (default: ``False``).
          This is mainly used by the spaces of homogeneous elements.

        OUTPUT:

        The corresponding graded ring of (Hecke) cusp forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import CuspFormsRing
            sage: MR = CuspFormsRing(5, CC, True)
            sage: MR
            CuspFormsRing(n=5) over Complex Field with 53 bits of precision
            sage: MR.analytic_type()
            cuspidal
            sage: MR.category()
            Category of commutative algebras over Complex Field with 53 bits of precision
            sage: MR in MR.category()
            True

            sage: CuspFormsRing(n=infinity, base_ring=CC, red_hom=True)
            CuspFormsRing(n=+Infinity) over Complex Field with 53 bits of precision
        """
        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n)
        cat = Algebras(base_ring).Commutative()
        Parent.__init__(self, category=cat)
        self._analytic_type = self.AT(["cusp"])
