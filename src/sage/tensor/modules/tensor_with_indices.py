# -*- coding: utf-8 -*-
r"""
Index notation for tensors

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version
- Léo Brunswic (2019): add multiple symmetries and multiple contractions

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************
from __future__ import annotations

from typing import Optional, TYPE_CHECKING, Union, List, Tuple
from sage.structure.sage_object import SageObject
from sage.groups.perm_gps.permgroup import PermutationGroup
import re
from itertools import combinations

if TYPE_CHECKING:
    from sage.tensor.modules.free_module_tensor import (
        IndexConfigurationNormalized,
        IndexCharacterNormalized,
        FreeModuleTensor,
    )

    IndicesWithCharacter = List[Tuple[str, IndexCharacterNormalized]]

# Regular expression for the allowed characters in index notation.
# This includes Unicode word constituents but excludes digits and underscores.
# Compare with https://docs.python.org/3/reference/lexical_analysis.html#identifiers
# The dot is special syntax for unnamed index positions.
_alph_or_dot_pattern = r"([.]|[^\d\W_])"

class TensorWithIndices(SageObject):
    r"""
    Index notation for tensors.

    This is a technical class to allow one to write some tensor operations
    (contractions and symmetrizations) in index notation.

    INPUT:

    - ``tensor`` -- a tensor (or a tensor field)
    - ``indices`` -- string containing the indices, as single letters; the
      contravariant indices must be stated first and separated from the
      covariant indices by the character ``_``

    EXAMPLES:

    Index representation of tensors on a rank-3 free module::

        sage: M = FiniteRankFreeModule(QQ, 3, name='M')
        sage: e = M.basis('e')
        sage: a = M.tensor((2,0), name='a')
        sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
        sage: b = M.tensor((0,2), name='b')
        sage: b[:] = [[-1,2,-3], [-4,5,6], [7,-8,9]]
        sage: t = a*b ; t.set_name('t') ; t
        Type-(2,2) tensor t on the 3-dimensional vector space M over the
         Rational Field
        sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
        sage: T = TensorWithIndices(t, '^ij_kl') ; T
        t^i^j_k_l

    The :class:`TensorWithIndices` object is returned by the square
    bracket operator acting on the tensor and fed with the string specifying
    the indices::

        sage: a['^ij']
        a^i^j
        sage: type(a['^ij'])
        <class 'sage.tensor.modules.tensor_with_indices.TensorWithIndices'>
        sage: b['_ef']
        b_e_f
        sage: t['^ij_kl']
        t^i^j_k_l

    The symbol '^' may be omitted, since the distinction between covariant
    and contravariant indices is performed by the index position relative to
    the symbol '_'::

        sage: t['ij_kl']
        t^i^j_k_l

    Also, LaTeX notation may be used::

        sage: t['^{ij}_{kl}']
        t^i^j_k_l

    If some operation is asked in the index notation, the resulting tensor
    is returned, not a :class:`TensorWithIndices` object; for instance, for
    a symmetrization::

        sage: s = t['^(ij)_kl'] ; s  # the symmetrization on i,j is indicated by parentheses
        Type-(2,2) tensor on the 3-dimensional vector space M over the
         Rational Field
        sage: s.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: s == t.symmetrize(0,1)
        True

    The letters denoting the indices can be chosen freely; since they carry no
    information, they can even be replaced by dots::

        sage: t['^(..)_..'] == t.symmetrize(0,1)
        True

    Similarly, for an antisymmetrization::

        sage: s = t['^ij_[kl]'] ; s # the symmetrization on k,l is indicated by square brackets
        Type-(2,2) tensor on the 3-dimensional vector space M over the Rational
         Field
        sage: s.symmetries()
        no symmetry;  antisymmetry: (2, 3)
        sage: s == t.antisymmetrize(2,3)
        True

    One can also perform multiple symmetrization-antisymmetrizations::

        sage: aa = a*a
        sage: aa['(..)(..)'] == aa.symmetrize(0,1).symmetrize(2,3)
        True
        sage: aa == aa['(..)(..)'] + aa['[..][..]'] + aa['(..)[..]'] + aa['[..](..)']
        True

    Another example of an operation indicated by indices is a contraction::

        sage: s = t['^ki_kj'] ; s  # contraction on the repeated index k
        Type-(1,1) tensor on the 3-dimensional vector space M over the Rational
         Field
        sage: s == t.trace(0,2)
        True

    Indices not involved in the contraction may be replaced by dots::

        sage: s == t['^k._k.']
        True

    The contraction of two tensors is indicated by repeated indices and
    the ``*`` operator::

        sage: s = a['^ik'] * b['_kj'] ; s
        Type-(1,1) tensor on the 3-dimensional vector space M over the Rational
         Field
        sage: s == a.contract(1, b, 0, swap_indices=False)
        True
        sage: s = t['^.k_..'] * b['_.k'] ; s
        Type-(1,3) tensor on the 3-dimensional vector space M over the Rational
         Field
        sage: s == t.contract(1, b, 1, swap_indices=False)
        True
        sage: t['^{ik}_{jl}']*b['_{mk}'] == s # LaTeX notation
        True

    Contraction on two indices::

        sage: s = a['^kl'] * b['_kl'] ; s
        105
        sage: s == (a*b)['^kl_kl']
        True
        sage: s == a.contract(0,1, b, 0,1, swap_indices=False)
        True

    The square bracket operator acts in a similar way on :class:`TensorWithIndices`::

        sage: b = +a["ij"] ; b._tensor.set_name("b") # create a copy of a["ij"]
        sage: b
        b^i^j
        sage: b[:]
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: b[0,0] == 1
        True
        sage: b["ji"]
        b^j^i
        sage: b["(ij)"][:]
        [1 3 5]
        [3 5 7]
        [5 7 9]
        sage: b["(ij)"] == b["(ij)"]["ij"]
        True

    However, it keeps track of indices::

        sage: b["ij"] = a["ji"]
        sage: b[:] == a[:]
        False
        sage: b[:] == a[:].transpose()
        True


    Arithmetics::

        sage: 2*a['^ij']
        X^i^j
        sage: (2*a['^ij'])._tensor == 2*a
        True
        sage: 2*t['ij_kl']
        X^i^j_k_l
        sage: +a['^ij']
        +a^i^j
        sage: +t['ij_kl']
        +t^i^j_k_l
        sage: -a['^ij']
        -a^i^j
        sage: -t['ij_kl']
        -t^i^j_k_l
        sage: a["^(..)"]["ij"] == 1/2*(a["^ij"] + a["^ji"])
        True

    The output indices are the ones of the left term of the addition::

        sage: a["^(..)"]["ji"] == 1/2*(a["^ij"] + a["^ji"])
        False
        sage: (a*a)["^..(ij)"]["abij"] == 1/2*((a*a)["^abij"] + (a*a)["^abji"])
        True
        sage: c = 1/2*((a*a)["^abij"] + (a*a)["^ijab"])
        sage: from itertools import product
        sage: all(c[i,j,k,l] == c[k,l,i,j] for i,j,k,l in product(range(3),repeat=4))
        True

    Non-digit unicode identifier characters are allowed::

        sage: a['^μξ']
        a^μ^ξ
    """

    _tensor: FreeModuleTensor
    _indices: IndicesWithCharacter

    @staticmethod
    def _parse_indices(
        indices: str,
        index_configuration: Optional[IndexConfigurationNormalized] = None,
        allow_contraction: bool = True,
        allow_symmetries: bool = True,
    ):
        r"""
        Parse index notation for tensors, enforces conventions and return
        indices.

        Parse ``indices`` checking usual conventions on repeating indices,
        wildcard, balanced parentheses/brackets and raises a ValueError if not.
        Also check whether the usual convention for indices, symmetries and
        contractions are respected. This includes restrictions on the
        indices symbols used, non nested (anti)symmetries,
        (co/contra)variant  identification of repeated indices, as well
        as checking the number of covariant and contravariant indices.
        Latex notations '{' and '}' are totally ignored.

        A valid input is for example: "^{ijkl}_{ib(cd)}"

        INPUT:

        - ``indices`` -- a string of index notation
        - ``index_configuration`` -- (default : ``None``) a valid index configuration.
            If not ``None``, the indices are checked to have the correct type.
        - ``allow_contraction`` -- (default : ``True``) Determines if
          repeated indices are allowed in the index notation.
        - ``allow_symmetries`` -- (default : ``True``) Determines if
          symmetries ()/[] are allowed in the index notation.

        OUTPUT:

        - A list of indices and their configuration.
        - A list of symmetries.
        - A list of antisymmetries.
        - A list of contractions.
        """
        # Suppress all '{' and '}' coming from LaTeX notations:
        indices = indices.replace("{", "").replace("}", "")

        output: IndicesWithCharacter = []
        last_type: IndexCharacterNormalized = "UP"
        in_symmetrization: Optional[list[int]] = None
        in_antisymmetrization: Optional[list[int]] = None
        symmetries: list[tuple[int, ...]] = []
        antisymmetries: list[tuple[int, ...]] = []
        position = 0
        for char in indices:
            if char == "^":
                last_type = "UP"
            elif char == "_":
                last_type = "DOWN"
            elif char == "(":
                if not allow_symmetries:
                    raise ValueError("symmetries not allowed")
                if in_symmetrization or in_antisymmetrization:
                    raise ValueError("nested symmetries not allowed (yet)")
                in_symmetrization = []
            elif char == "[":
                if not allow_symmetries:
                    raise ValueError("symmetries not allowed")
                if in_symmetrization or in_antisymmetrization:
                    raise ValueError("nested symmetries not allowed (yet)")
                in_antisymmetrization = []
            elif char == ")":
                if not in_symmetrization:
                    raise ValueError("unbalanced symmetrization")
                symmetries.append(tuple(in_symmetrization))
                in_symmetrization = None
            elif char == "]":
                if not in_antisymmetrization:
                    print(output, in_antisymmetrization)
                    raise ValueError("unbalanced antisymmetrization")
                antisymmetries.append(tuple(in_antisymmetrization))
                in_antisymmetrization = None
            elif re.match(_alph_or_dot_pattern, char):
                output.append((char, last_type))
                if in_symmetrization is not None:
                    in_symmetrization.append(position)
                elif in_antisymmetrization is not None:
                    in_antisymmetrization.append(position)
                position += 1
            else:
                raise ValueError(f"{char} is not a valid character")

        if in_symmetrization is not None or in_antisymmetrization is not None:
            raise ValueError("unbalanced symmetrization or antisymmetrization")

        contraction_pairs: list[tuple[int, int]] = []
        for i, (ind1, type1) in enumerate(output):
            for j, (ind2, type2) in enumerate(output):
                if i < j and ind1 == ind2 and ind1 != ".":
                    if not allow_contraction:
                        raise ValueError(
                            f"contraction of '{ind1}' and '{ind2}' not allowed"
                        )
                    if type1 == type2:
                        raise ValueError(
                            f"contraction of '{ind1}' with same index character not allowed"
                        )
                    contraction_pairs.append((i, j))

        if index_configuration is not None:
            if len(index_configuration) != len(output):
                raise IndexError(
                    "number of indices not compatible with the index configuration"
                )
            for i, (character, index_type) in enumerate(output):
                if index_type != index_configuration[i] and character != ".":
                    raise IndexError(
                        f"index configuration not satisfied for index '{character}'"
                    )
        return output, symmetries, antisymmetries, contraction_pairs

    def __init__(
        self, tensor: FreeModuleTensor, indices: Union[str, IndicesWithCharacter]
    ):
        r"""
        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')

        We need to skip the pickling test because we can't check equality
        unless the tensor was defined w.r.t. a basis::

            sage: TestSuite(ti).run(skip="_test_pickling")

        ::

            sage: e = M.basis('e')
            sage: t[:] = [[[1,2,3], [-4,5,6], [7,8,-9]],
            ....:         [[10,-11,12], [13,14,-15], [16,17,18]],
            ....:         [[19,-20,-21], [-22,23,24], [25,26,-27]]]
            sage: ti = TensorWithIndices(t, 'ab_c')
            sage: TestSuite(ti).run()

        """
        self._tensor = tensor # may be changed below
        self._changed = False # indicates whether self contains an altered
                              # version of the original tensor (True if
                              # symmetries or contractions are indicated in the
                              # indices)

        if not isinstance(indices, str):
            self._indices = indices
            return

        parsed_indices, symmetries, antisymmetries, contractions = self._parse_indices(
            indices, index_configuration=self._tensor.config()
        )

        # Apply (anti)symmetrizations
        for sym in symmetries:
            self._tensor = self._tensor.symmetrize(*sym)
            self._changed = True  # self does no longer contain the original tensor
        for sym in antisymmetries:
            self._tensor = self._tensor.antisymmetrize(*sym)
            self._changed = True  # self does no longer contain the original tensor

        # Treatment of possible self-contractions:
        if contractions:
            self._tensor = self._tensor.trace(contractions)
            self._changed = True  # self does no longer contain the original tensor
            for pos in sorted(
                [pos for contraction in contractions for pos in contraction],
                reverse=True,
            ):
                parsed_indices.pop(pos)

        self._indices = parsed_indices

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')
            sage: ti._repr_()
            't^a^b_c'
            sage: t = M.tensor((0,2), name='t')
            sage: ti = TensorWithIndices(t, '_{ij}')
            sage: ti._repr_()
            't_i_j'

        """
        name = 'X'
        if hasattr(self._tensor, '_name'):
            if self._tensor._name is not None:
                name = self._tensor._name
        if len(self._indices) == 0:
            return "scalar"

        return name + "".join(
            [
                (("^" if character == "UP" else "_") + index)
                for (index, character) in self._indices
            ]
        )

    def update(self):
        r"""
        Return the tensor contains in ``self`` if it differs from that used
        for creating ``self``, otherwise return ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((1,1),  name='a')
            sage: a[:] = [[1,-2,3], [-4,5,-6], [7,-8,9]]
            sage: a_ind = TensorWithIndices(a, 'i_j') ; a_ind
            a^i_j
            sage: a_ind.update()
            a^i_j
            sage: a_ind.update() is a_ind
            True
            sage: a_ind = TensorWithIndices(a, 'k_k') ; a_ind
            scalar
            sage: a_ind.update()
            15

        """
        if self._changed:
            return self._tensor
        else:
            return self

    def __eq__(self, other):
        r"""
        Check equality.

        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')
            sage: ti == TensorWithIndices(t, '^{ab}_c')
            True
            sage: ti == TensorWithIndices(t, 'ac_b')
            False
            sage: tp = M.tensor((2,1))
            sage: ti == TensorWithIndices(tp, 'ab_c')
            Traceback (most recent call last):
            ...
            ValueError: no common basis for the comparison

        """
        if not isinstance(other, TensorWithIndices):
            return False
        return self._tensor == other._tensor and self._indices == other._indices

    def __ne__(self, other):
        r"""
        Check not equals.

        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')
            sage: ti != TensorWithIndices(t, '^{ab}_c')
            False
            sage: ti != TensorWithIndices(t, 'ac_b')
            True

        """
        return not self == other

    def __mul__(self, other):
        r"""
        Tensor product or contraction on specified indices.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,-2,3], [-4,5,-6], [7,-8,9]]
            sage: b = M.linear_form(name='b')
            sage: b[:] = [4,2,1]
            sage: ai = TensorWithIndices(a, '^ij')
            sage: bi = TensorWithIndices(b, '_k')
            sage: s = ai.__mul__(bi) ; s  # no repeated indices ==> tensor product
            Type-(2,1) tensor a⊗b on the 3-dimensional vector space M over the
             Rational Field
            sage: s == a*b
            True
            sage: s[:]
            [[[4, 2, 1], [-8, -4, -2], [12, 6, 3]],
             [[-16, -8, -4], [20, 10, 5], [-24, -12, -6]],
             [[28, 14, 7], [-32, -16, -8], [36, 18, 9]]]
            sage: ai = TensorWithIndices(a, '^kj')
            sage: s = ai.__mul__(bi) ; s  # repeated index k ==> contraction
            Element of the 3-dimensional vector space M over the Rational Field
            sage: s == a.contract(0, b, swap_indices=False)
            True
            sage: s[:]
            [3, -6, 9]
        """
        if not isinstance(other, TensorWithIndices):
            raise TypeError("the second item of * must be a tensor with " +
                            "specified indices")
        contractions: list[tuple[int, int]] = []
        for position1, (index, character) in enumerate(self._indices):
            if index != ".":
                for position2, (other_index, other_character) in enumerate(
                    other._indices
                ):
                    if index == other_index:
                        if character == other_character:
                            raise IndexError(
                                f"the index {index} appears twice with the same character"
                            )
                        contractions.append((position1, position2))

        if not contractions:
            # No contraction is performed: the tensor product is returned
            return self._tensor * other._tensor
        ncontr = len(contractions)
        pos1 = [contractions[i][0] for i in range(ncontr)]
        pos2 = [contractions[i][1] for i in range(ncontr)]
        args = pos1 + [other._tensor] + pos2
        return self._tensor.contract(*args, swap_indices=False)

    def __rmul__(self, other):
        r"""
        Multiplication on the left by ``other``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,1), name='a')
            sage: a[0,2,1], a[1,2,0] = 7, -4
            sage: ai = TensorWithIndices(a, 'ij_k')
            sage: s = ai.__rmul__(3) ; s
            X^i^j_k
            sage: s._tensor == 3*a
            True

        """
        return TensorWithIndices(other * self._tensor, self._indices)

    def __add__(self, other):
        r"""
        Addition between tensors with indices.

        The underlying tensor of the output is the sum of the underlying tensor
        of ``self`` with the underlying tensor of ``other`` whose entries have
        be permuted to respect Einstein summation usual conventions. The
        indices names of the output are those of self.


        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b = M.tensor((0,2), name='b')
            sage: b[:] = [[-1,2,-3], [-4,5,6], [7,-8,9]]
            sage: T = a*a*b*b
            sage: 1/4*(T["ijkl_abcd"] + T["jikl_abcd"] + T["ijkl_abdc"]\
             + T["jikl_abdc"]) == T["(..).._..(..)"]["ijkl_abcd"]
            True

        """
        # Check tensor types are compatible
        if self._tensor.tensor_type() != other._tensor.tensor_type():
            raise ValueError("Tensors are not of the same type")

        permutation = _find_permutation(self._indices, other._indices)

        result = self.__pos__()
        result._tensor = result._tensor + other.permute_indices(permutation)._tensor
        return result

    def __sub__(self, other):
        r"""
        Subtraction between tensors with indices.

        The underlying tensor of the output is the underlying tensor of
        ``self`` minus the underlying tensor of ``other`` whose entries have
        be permuted to respect Einstein summation usual conventions. The
        indices names of the output are those of self.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b = M.tensor((0,2), name='b')
            sage: b[:] = [[-1,2,-3], [-4,5,6], [7,-8,9]]
            sage: a["^[..]"]["ij"] == 1/2*(a["^ij"]-a["^ji"])
            True
            sage: (a*a)["^..[ij]"]["abij"] == 1/2*((a*a)["^abij"]-(a*a)["^abji"])
            True
            sage: Riem = a*a
            sage: Riem = Riem["[ij][kl]"]
            sage: Riem = 1/2*(Riem["ijkl"]+Riem["klij"])
            sage: O = M.tensor((4,0), name='O')
            sage: O[0,0,0,0] = 0
            sage: (Riem["ijkl"]+Riem["iklj"]+Riem["iljk"]) == O["ijkl"]
            True

        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b = M.tensor((0,2), name='b')
            sage: b[:] = [[-1,2,-3], [-4,5,6], [7,-8,9]]
            sage: T = a*a*b*b
            sage: 1/4*(T["ijkl_abcd"]-T["jikl_abcd"] - T["ijkl_abdc"]\
                + T["jikl_abdc"] ) == T["[..].._..[..]"]["ijkl_abcd"]
            True

        """
        return self + (-other)

    def __getitem__(self, args):
        r"""
        Return a component of the underlying tensor w.r.t. some basis.

        NB: if ``args`` is a string, this method acts as a shortcut for
        tensor contractions and symmetrizations, the string containing
        abstract indices.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are returned. The basis can be passed
          as the first item of ``args``; if not, the free module's default
          basis is assumed.
          if ``args`` is a string, this method acts as a shortcut for
          tensor contractions and symmetrizations, the string containing
          abstract indices.

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b = a["ij"]
            sage: b
            a^i^j
            sage: b[:]
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: b[0,0] == 1
            True
            sage: b["ji"]
            a^j^i
            sage: b["(ij)"][:]
            [1 3 5]
            [3 5 7]
            [5 7 9]

        """
        if isinstance(args, str):
            result = +self
            result.__init__(self._tensor, args)
            return result
        else:
            return self._tensor[args]

    def __setitem__(self, args, value):
        r"""
         Set a component w.r.t. some basis.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are set. The basis can be passed
          as the first item of ``args``; if not, the free module's default
          basis is assumed  if ``args`` is a string and value is a tensor
          with indices, this method permutes the coefficients of ``value``
          before assigning the underlying tensor of ``value`` to ``self``.

        - ``value`` -- the value to be set or a list of values if
          ``args = [:]`` or a tensor with indices

        TESTS::

            sage: M = FiniteRankFreeModule(ZZ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')["ij"]
            sage: b = M.tensor((2,0), name='b')["ij"]
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b["ij"] = a["ji"]
            sage: b[:] == a[:].transpose()
            True

        """
        if isinstance(args, str):
            if not isinstance(value,TensorWithIndices):
                raise ValueError("The tensor provided should be with indices")
            elif self._tensor.tensor_type() != value._tensor.tensor_type():
                raise ValueError("The tensors are not of the same type")

            indices, *rest = self._parse_indices(
                args,
                index_configuration=self._tensor.config(),
                allow_symmetries=False,
                allow_contraction=False,
            )

            permutation = _find_permutation(value._indices, indices)
            self._tensor[:] = value.permute_indices(permutation)[:]

        else:
            self._tensor.__setitem__(args,value)

    def permute_indices(self, permutation):
        r"""
        Return a tensor with indices with permuted indices.

        INPUT:

        - ``permutation`` -- permutation that has to be applied to the indices
          the input should be a ``list`` containing the second line of the permutation
          in Cauchy notation.

        OUTPUT:

        - an instance of ``TensorWithIndices`` whose indices names and place
          are those of ``self`` but whose components have been permuted with
          ``permutation``.

        EXAMPLES::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: b = M.tensor((2,0), name='b')
            sage: b[:] = [[-1,2,-3], [-4,5,6], [7,-8,9]]
            sage: identity = [0,1]
            sage: transposition = [1,0]
            sage: a["ij"].permute_indices(identity) == a["ij"]
            True
            sage: a["ij"].permute_indices(transposition)[:] == a[:].transpose()
            True
            sage: cycle = [1,2,3,0] # the cyclic permutation sending 0 to 1
            sage: (a*b)[0,1,2,0] == (a*b)["ijkl"].permute_indices(cycle)[1,2,0,0]
            True

        TESTS::

            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,0), name='a')
            sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
            sage: identity = [0,1]
            sage: transposition = [1,0]
            sage: a["ij"].permute_indices(identity) == a["ij"]
            True
            sage: a["ij"].permute_indices(transposition)[:] == a[:].transpose()
            True
            sage: (a*a)["ijkl"].permute_indices([1,2,3,0])[0,1,2,1] == (a*a)[1,2,1,0]
            True
        """
        # Decomposition of the permutation of the components of self
        # into product of swaps given by the method
        # sage.tensor.modules.comp.Components.swap_adjacent_indices

        # A swap is determined by 3 distinct integers
        swap_params = list(combinations(range(self._tensor.tensor_rank()+1), 3))

        # The associated permutation is as follows
        def swap(param,N):
            i,j,k = param
            L = list(range(1,N+1))
            L = L[:i] + L[j:k] + L[i:j] + L[k:]
            return L

        # Construction of the permutation group generated by swaps
        perm_group = PermutationGroup(
            [swap(param, self._tensor.tensor_rank()) for param in swap_params],
            canonicalize = False
        )
        # Compute a decomposition of the permutation as a product of swaps
        decomposition_as_string = perm_group([x+1 for x in permutation]).word_problem(
            perm_group.gens(),
            display=False
        )[0]

        if decomposition_as_string != "<identity ...>":
            decomposition_as_string = [
                # Two cases whether the term appear with an exponent or not
                ("^" in term)*term.split("^") + ("^" not in term)*(term.split("^")+['1'])
                for term in decomposition_as_string.replace("x","").split("*")
            ]
            decomposition = [(swap_params[int(x)-1], int(y)) for x, y in decomposition_as_string]
            decomposition.reverse()  # /!\ The symmetric group acts on the right by default /!\.
        else:
            decomposition = []
        # Choice of a basis
        basis = self._tensor._fmodule._def_basis

        # Swap of components

        swaped_components = self._tensor.comp(basis)
        for swap_param,exponent in decomposition:
            if exponent > 0:
                for i in range(exponent):
                    # Apply the swap given by swap_param
                    swaped_components = swaped_components\
                        .swap_adjacent_indices(*swap_param)
            elif exponent < 0:
                for i in range(-exponent):
                    # Apply the opposite of the swap given by swap_param
                    swaped_components = swaped_components\
                        .swap_adjacent_indices(
                            swap_param[0],
                            swap_param[0] + swap_param[2] - swap_param[1],
                            swap_param[2]
                        )
            else:
                pass
        result = self.__pos__()
        result._tensor = self._tensor._fmodule.tensor_from_comp(
            self._tensor.tensor_type(),
            swaped_components
        )

        return result

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,1), name='a')
            sage: a[0,2,1], a[1,2,0] = 7, -4
            sage: ai = TensorWithIndices(a, 'ij_k')
            sage: s = ai.__pos__() ; s
            +a^i^j_k
            sage: s._tensor == a
            True

        """
        return TensorWithIndices(+self._tensor, self._indices)

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - negative of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,1), name='a')
            sage: a[0,2,1], a[1,2,0] = 7, -4
            sage: ai = TensorWithIndices(a, 'ij_k')
            sage: s = ai.__neg__() ; s
            -a^i^j_k
            sage: s._tensor == -a
            True

        """
        return TensorWithIndices(-self._tensor, self._indices)


def _find_permutation(first: IndicesWithCharacter, second: IndicesWithCharacter):
    permutation = list(range(len(first)))
    try:
        for i, (index, _) in enumerate(first):
            permutation[i] = next(
                j for j, (index2, _) in enumerate(second) if index2 == index
            )
    except StopIteration:
        raise ValueError(
            f"Indices '{first}' and '{second}' cannot be permutated to match each other"
        )
    return permutation
