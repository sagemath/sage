# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.groups
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

from sage.structure.sage_object import SageObject
import re
from itertools import combinations

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
        t^ij_kl

    The :class:`TensorWithIndices` object is returned by the square
    bracket operator acting on the tensor and fed with the string specifying
    the indices::

        sage: a['^ij']
        a^ij
        sage: type(a['^ij'])
        <class 'sage.tensor.modules.tensor_with_indices.TensorWithIndices'>
        sage: b['_ef']
        b_ef
        sage: t['^ij_kl']
        t^ij_kl

    The symbol '^' may be omitted, since the distinction between covariant
    and contravariant indices is performed by the index position relative to
    the symbol '_'::

        sage: t['ij_kl']
        t^ij_kl

    Also, LaTeX notation may be used::

        sage: t['^{ij}_{kl}']
        t^ij_kl

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
        sage: s == a.contract(1, b, 0)
        True
        sage: s = t['^.k_..'] * b['_.k'] ; s
        Type-(1,3) tensor on the 3-dimensional vector space M over the Rational
         Field
        sage: s == t.contract(1, b, 1)
        True
        sage: t['^{ik}_{jl}']*b['_{mk}'] == s # LaTeX notation
        True

    Contraction on two indices::

        sage: s = a['^kl'] * b['_kl'] ; s
        105
        sage: s == (a*b)['^kl_kl']
        True
        sage: s == (a*b)['_kl^kl']
        True
        sage: s == a.contract(0,1, b, 0,1)
        True

    The square bracket operator acts in a similar way on :class:`TensorWithIndices`::

        sage: b = +a["ij"] ; b._tensor.set_name("b") # create a copy of a["ij"]
        sage: b
        b^ij
        sage: b[:]
        [1 2 3]
        [4 5 6]
        [7 8 9]
        sage: b[0,0] == 1
        True
        sage: b["ji"]
        b^ji
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
        X^ij
        sage: (2*a['^ij'])._tensor == 2*a
        True
        sage: 2*t['ij_kl']
        X^ij_kl
        sage: +a['^ij']
        +a^ij
        sage: +t['ij_kl']
        +t^ij_kl
        sage: -a['^ij']
        -a^ij
        sage: -t['ij_kl']
        -t^ij_kl
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
        a^μξ

    Conventions are checked and non acceptable indices raise
    :exc:`ValueError`, for instance::

        sage: a['([..])']  # nested symmetries
        Traceback (most recent call last):
        ...
        ValueError: index conventions not satisfied
        sage: a['(..']  # unbalanced parenthis
        Traceback (most recent call last):
        ...
        ValueError: index conventions not satisfied
        sage: a['ii']  # repeated indices of the same type
        Traceback (most recent call last):
        ...
        ValueError: index conventions not satisfied: repeated indices of same type
        sage: (a*a)['^(ij)^(kl)']  # multiple indices group of the same type
        Traceback (most recent call last):
        ...
        ValueError: index conventions not satisfied
        sage: a["^\u2663\u2665"]  # non-word-constituent
        Traceback (most recent call last):
        ...
        ValueError: index conventions not satisfied
    """

    @staticmethod
    def _parse_indices(indices, tensor_type=None, allow_contraction=True,
                       allow_symmetries=True):
        r"""
        Parse index notation for tensors, enforces conventions and return
        indices.

        Parse ``indices`` checking usual conventions on repeating indices,
        wildcard, balanced parentheses/brackets and raises a :exc:`ValueError`
        if not. Return a couple contravariant/covariant indices.

        INPUT:

        - ``indices`` -- string of index notation
        - ``tensor_type`` -- (default: ``None``) a valid tensor type
          (a couple of nonnegative integers). If not ``None``, the indices
          are checked to have the correct type.
        - ``allow_contraction`` -- (default: ``True``) determines if
          repeated indices are allowed in the index notation
        - ``allow_symmetries`` -- (default: ``True``) determines if
          symmetries ()/[] are allowed in the index notation

        OUTPUT:

        - A couple of string corresponding to the contravariant and the
          covariant part

        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: TensorWithIndices._parse_indices('([..])')  # nested symmetries
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices('(..')  # unbalanced parenthis
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices('ii')  # repeated indices of the same type
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied: repeated indices of same type
            sage: TensorWithIndices._parse_indices('^(ij)^(kl)')  # multiple indices group of the same type
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices("^17")  # digits are not allowed as names
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices("^;")  # non-word-constituents are not allowed as names
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices("^\u00ae")  # non-word-constituents are not allowed as names
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices("^\u25e2")  # non-word-constituents are not allowed as names
            Traceback (most recent call last):
            ...
            ValueError: index conventions not satisfied
            sage: TensorWithIndices._parse_indices('^ij_kl')
            ('ij', 'kl')
            sage: TensorWithIndices._parse_indices('_kl^ij')
            ('ij', 'kl')
            sage: TensorWithIndices._parse_indices("(ij)_ik",tensor_type=(2,2))
            ('(ij)', 'ik')
            sage: TensorWithIndices._parse_indices("(ij)_ik",tensor_type=(2,0))
            Traceback (most recent call last):
            ...
            IndexError: number of covavariant indices not compatible with the tensor type
            sage: TensorWithIndices._parse_indices("(ij)_ik", allow_contraction=False)
            Traceback (most recent call last):
            ...
            IndexError: no contraction allowed
            sage: TensorWithIndices._parse_indices("(ij)_ik", allow_symmetries=False)
            Traceback (most recent call last):
            ...
            IndexError: no symmetry allowed
        """
        # Suppress all '{' and '}' coming from LaTeX notations:
        indices = indices.replace('{','').replace('}','')

        # Check index notation conventions and parse indices
        allowed_pattern = r"(\(" + _alph_or_dot_pattern + r"{2,}\)|\[" + _alph_or_dot_pattern + r"{2,}\]|" + _alph_or_dot_pattern + r"+)*"
        con_then_cov = r"^(\^|)" + allowed_pattern + r"(\_" + allowed_pattern + r"|)$"
        cov_then_con = r"^\_" + allowed_pattern + r"(\^" + allowed_pattern + r"|)$"
        if (re.match(con_then_cov,indices) is None
            and re.match(cov_then_con,indices) is None):
            raise ValueError("index conventions not satisfied")
        elif re.match(con_then_cov,indices):
            try:
                con,cov = indices.replace("^","").split("_")
            except ValueError:
                con = indices.replace("^","")
                cov = ""
        else:
            try:
                cov,con = indices.replace("_","").split("^")
            except ValueError:
                cov = indices.replace("_","")
                con = ""
        if not allow_contraction:
            for ind in con:
                if ind != '.' and ind in cov:
                    raise IndexError("no contraction allowed")
        con_without_sym = (con.replace("(","").replace(")","").replace("[","").replace("]",""))
        cov_without_sym = (cov.replace("(","").replace(")","").replace("[","").replace("]",""))
        if allow_symmetries:
            if len(con_without_sym) != len(set(con_without_sym)) \
                                       + max(con_without_sym.count(".")-1, 0):
                raise ValueError("index conventions not satisfied: "
                                 "repeated indices of same type")
            if len(cov_without_sym) != len(set(cov_without_sym)) \
                                       + max(cov_without_sym.count(".")-1, 0):
                raise ValueError("index conventions not satisfied: "
                                 "repeated indices of same type")
        else:
            if re.search(r"[()\[\]]",con) is not None:
                raise IndexError("no symmetry allowed")
            if re.search(r"[()\[\]]",cov) is not None:
                raise IndexError("no symmetry allowed")
        if tensor_type is not None:
            # Check number of (co/contra)variant indices
            if len(con_without_sym) != tensor_type[0]:
                raise IndexError("number of contravariant indices not compatible "
                                 "with the tensor type")
            if len(cov_without_sym) != tensor_type[1]:
                raise IndexError("number of covavariant indices not compatible "
                                 "with the tensor type")
        return con,cov

    def __init__(self, tensor, indices):
        r"""
        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')

        We need to skip the pickling test because we can't check equality
        unless the tensor was defined w.r.t. a basis::

            sage: TestSuite(ti).run(skip='_test_pickling')

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

        # Check whether the usual convention for indices, symmetries and
        # contractions are respected. This includes restrictions on the
        # indices symbols used, non nested (anti)symmetries,
        # (co/contra)variant  identification of repeated indices, as well
        # as checking the number of covariant and contravariant indices.
        # Latex notations '{' and '}' are totally ignored.
        # "^{ijkl}_{ib(cd)}"

        con,cov = self._parse_indices(
            indices,
            tensor_type=self._tensor.tensor_type()
        )

        # Apply (anti)symmetrizations on contravariant indices
        first_sym_regex = r"(\(|\[)" + _alph_or_dot_pattern + r"*[)\]]"
        while re.search(first_sym_regex,con):
            first_sym = re.search(first_sym_regex,con)
            sym1 = first_sym.span()[0]
            sym2 = first_sym.span()[1]-1
            if first_sym.groups()[0] == "(":
                self._tensor = self._tensor.symmetrize(*range(
                    sym1,
                    sym2-1
                ))
            else:
                self._tensor = self._tensor.antisymmetrize(*range(
                    sym1,
                    sym2-1
                ))
            self._changed = True # self does no longer contain the original tensor
            con = con[:sym1] + con[sym1+1:sym2] + con[sym2+1:]
        self._con = con

        # Apply (anti)symmetrizations on covariant indices
        while re.search(first_sym_regex,cov):
            first_sym = re.search(first_sym_regex,cov)
            sym1 = first_sym.span()[0]
            sym2 = first_sym.span()[1]-1
            if first_sym.groups()[0] == "(":
                self._tensor = self._tensor.symmetrize(*range(
                    self._tensor._tensor_type[0] + sym1,
                    self._tensor._tensor_type[0] + sym2-1
                ))
            else:
                self._tensor = self._tensor.antisymmetrize(*range(
                    self._tensor._tensor_type[0] + sym1,
                    self._tensor._tensor_type[0] + sym2-1
                ))
            self._changed = True # self does no longer contain the original tensor
            cov = cov[:sym1] + cov[sym1+1:sym2] + cov[sym2+1:]
        self._cov = cov

        # Treatment of possible self-contractions:
        # ---------------------------------------
        contraction_pair_list = []
        for ind in self._con:
            if ind != '.' and ind in self._cov:
                pos1 = self._con.index(ind)
                pos2 = self._tensor._tensor_type[0] + self._cov.index(ind)
                contraction_pair_list.append([pos1, pos2])
        while contraction_pair_list:
            pos1, pos2 = contraction_pair_list.pop()
            self._tensor = self._tensor.trace(pos1, pos2)
            for contraction_pair in contraction_pair_list:
                if contraction_pair[0] > pos1:
                    contraction_pair[0] = contraction_pair[0]-1
                if contraction_pair[1] > pos2:
                    contraction_pair[1] = contraction_pair[1]-1
                contraction_pair[1] = contraction_pair[1]-1
            self._changed = True # self does no longer contain the original
                                 # tensor
            ind = self._con[pos1]
            self._con = self._con.replace(ind, '')
            self._cov = self._cov.replace(ind, '')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')
            sage: ti._repr_()
            't^ab_c'
            sage: t = M.tensor((0,2), name='t')
            sage: ti = TensorWithIndices(t, '_{ij}')
            sage: ti._repr_()
            't_ij'
        """
        name = 'X'
        if hasattr(self._tensor, '_name'):
            if self._tensor._name is not None:
                name = self._tensor._name
        if self._con == '':
            if self._cov == '':
                return 'scalar'
            else:
                return name + '_' + self._cov
        elif self._cov == '':
            return name + '^' + self._con
        else:
            return name + '^' + self._con + '_' + self._cov

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
        return (self._tensor == other._tensor
                and self._con == other._con
                and self._cov == other._cov)

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
            sage: s == a.contract(0, b)
            True
            sage: s[:]
            [3, -6, 9]
        """
        if not isinstance(other, TensorWithIndices):
            raise TypeError("the second item of * must be a tensor with " +
                            "specified indices")
        contraction_pairs = []
        for ind in self._con:
            if ind != '.':
                if ind in other._cov:
                    pos1 = self._con.index(ind)
                    pos2 = other._tensor._tensor_type[0] + other._cov.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if ind in other._con:
                    raise IndexError(f"the index {ind} appears twice "
                                     + "in a contravariant position")
        for ind in self._cov:
            if ind != '.':
                if ind in other._con:
                    pos1 = self._tensor._tensor_type[0] + self._cov.index(ind)
                    pos2 = other._con.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if ind in other._cov:
                    raise IndexError(f"the index {ind} appears twice "
                                     + "in a covariant position")
        if not contraction_pairs:
            # No contraction is performed: the tensor product is returned
            return self._tensor * other._tensor
        ncontr = len(contraction_pairs)
        pos1 = [contraction_pairs[i][0] for i in range(ncontr)]
        pos2 = [contraction_pairs[i][1] for i in range(ncontr)]
        args = pos1 + [other._tensor] + pos2
        return self._tensor.contract(*args)

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
            X^ij_k
            sage: s._tensor == 3*a
            True
        """
        return TensorWithIndices(other*self._tensor,
                                 self._con + '_' + self._cov)

    def __add__(self, other):
        r"""
        Addition between tensors with indices.

        The underlying tensor of the output is the sum of the underlying tensor
        of ``self`` with the underlying tensor of ``other`` whose entries have
        be permuted to respect Einstein summation usual conventions. The
        indices names of the output are those of ``self``.


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
        # Check the set of indices are compatible
        if set(self._cov) != set(other._cov):
            raise ValueError("The covariant Indices sets are not identical")
        if set(self._con) != set(other._con):
            raise ValueError("The contravariant Indices sets are not identical")
        self_wild_card_indices = [match.span()[0] for match in re.finditer(r"\.", self._con)]
        other_wild_card_indices = [match.span()[0] for match in re.finditer(r"\.", self._cov)]
        if self_wild_card_indices != other_wild_card_indices:
            raise ValueError("Ambiguous wildcard notation")

        # Permutation of the components of self
        # -------------------------------------

        permutation = list(range(other._tensor.tensor_rank()))
        for other_index in range(other._tensor.tensor_type()[0]):
            if other._con[other_index] == self._con[other_index]:
                permutation[other_index] = other_index
            else:
                permutation[other_index] = self._con.index(other._con[other_index])
        for other_index in range(other._tensor.tensor_type()[1]):
            if other._cov[other_index] == self._cov[other_index]:
                permutation[other._tensor.tensor_type()[0] + other_index]\
                    = other._tensor.tensor_type()[0] + other_index
            else:
                permutation[other._tensor.tensor_type()[0] + other_index]\
                    = other._tensor.tensor_type()[0] + self._cov.index(other._cov[other_index])

        result = self.__pos__()
        result._tensor = result._tensor + other.permute_indices(permutation)._tensor
        return result

    def __sub__(self, other):
        r"""
        Subtraction between tensors with indices.

        The underlying tensor of the output is the underlying tensor of
        ``self`` minus the underlying tensor of ``other`` whose entries have
        be permuted to respect Einstein summation usual conventions. The
        indices names of the output are those of ``self``.

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
            a^ij
            sage: b[:]
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: b[0,0] == 1
            True
            sage: b["ji"]
            a^ji
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
            else:
                con,cov = self._parse_indices(
                    args,
                    tensor_type=self._tensor.tensor_type(),
                    allow_symmetries=False,
                    allow_contraction=False
                )

            permutation = list(range(value._tensor.tensor_rank()))
            for value_index in range(value._tensor.tensor_type()[0]):
                if value._con[value_index] == self._con[value_index]:
                    permutation[value_index] = value_index
                else:
                    permutation[value_index] = self._con.index(value._con[value_index])
            for value_index in range(value._tensor.tensor_type()[1]):
                if value._cov[value_index] == self._cov[value_index]:
                    permutation[value._tensor.tensor_type()[0] + value_index]\
                        = value._tensor.tensor_type()[0] + value_index
                else:
                    permutation[value._tensor.tensor_type()[0] + value_index]\
                        = value._tensor.tensor_type()[0] + self._cov.index(value._cov[value_index])
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
        from sage.groups.perm_gps.permgroup import PermutationGroup

        perm_group = PermutationGroup(
            [swap(param, self._tensor.tensor_rank()) for param in swap_params],
            canonicalize=False
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

        OUTPUT: an exact copy of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,1), name='a')
            sage: a[0,2,1], a[1,2,0] = 7, -4
            sage: ai = TensorWithIndices(a, 'ij_k')
            sage: s = ai.__pos__() ; s
            +a^ij_k
            sage: s._tensor == a
            True
        """
        return TensorWithIndices(+self._tensor,
                                 self._con + '_' + self._cov)

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT: negative of ``self``

        EXAMPLES::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: e = M.basis('e')
            sage: a = M.tensor((2,1), name='a')
            sage: a[0,2,1], a[1,2,0] = 7, -4
            sage: ai = TensorWithIndices(a, 'ij_k')
            sage: s = ai.__neg__() ; s
            -a^ij_k
            sage: s._tensor == -a
            True
        """
        return TensorWithIndices(-self._tensor,
                                 self._con + '_' + self._cov)
