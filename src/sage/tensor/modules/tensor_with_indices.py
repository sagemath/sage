r"""
Index notation for tensors

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

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
        sage: s == a.contract(0,1, b, 0,1)
        True

    Some minimal arithmetics::

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

    """
    def __init__(self, tensor, indices):
        r"""
        TESTS::

            sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
            sage: M = FiniteRankFreeModule(QQ, 3, name='M')
            sage: t = M.tensor((2,1), name='t')
            sage: ti = TensorWithIndices(t, 'ab_c')

        We need to skip the pickling test because we can't check equality
        unless the tensor was defined w.r.t. a basis::

            sage: TestSuite(ti).run(skip="_test_pickling")

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
        # Suppress all '{' and '}' coming from LaTeX notations:
        indices = indices.replace('{','').replace('}','')
        # Suppress the first '^':
        if indices[0] == '^':
            indices = indices[1:]
        if '^' in indices:
            raise IndexError("the contravariant indices must be placed first")
        con_cov = indices.split('_')
        con = con_cov[0]

        # Contravariant indices
        # ---------------------
        #  Search for (anti)symmetries:
        if '(' in con:
            sym1 = con.index('(')
            sym2 = con.index(')')-2
            if con.find('(', sym1+1) != -1 or '[' in con:
                raise NotImplementedError("Multiple symmetries are not " +
                                          "treated yet.")
            self._tensor = self._tensor.symmetrize(*range(sym1, sym2+1))
            self._changed = True # self does no longer contain the original tensor
            con = con.replace('(','').replace(')','')
        if '[' in con:
            sym1 = con.index('[')
            sym2 = con.index(']')-2
            if con.find('[', sym1+1) != -1 or '(' in con:
                raise NotImplementedError("multiple symmetries are not " +
                                          "treated yet")
            self._tensor = self._tensor.antisymmetrize(*range(sym1, sym2+1))
            self._changed = True # self does no longer contain the original tensor
            con = con.replace('[','').replace(']','')
        if len(con) != self._tensor._tensor_type[0]:
            raise IndexError("number of contravariant indices not compatible " +
                             "with the tensor type")
        self._con = con

        # Covariant indices
        # -----------------
        if len(con_cov) == 1:
            if tensor._tensor_type[1] != 0:
                raise IndexError("number of covariant indices not compatible " +
                                 "with the tensor type")
            self._cov = ''
        elif len(con_cov) == 2:
            cov = con_cov[1]
            #  Search for (anti)symmetries:
            if '(' in cov:
                sym1 = cov.index('(')
                sym2 = cov.index(')')-2
                if cov.find('(', sym1+1) != -1 or '[' in cov:
                    raise NotImplementedError("multiple symmetries are not " +
                                              "treated yet")
                csym1 = sym1 + self._tensor._tensor_type[0]
                csym2 = sym2 + self._tensor._tensor_type[0]
                self._tensor = self._tensor.symmetrize(*range(csym1, csym2+1))
                self._changed = True # self does no longer contain the original
                                     # tensor
                cov = cov.replace('(','').replace(')','')
            if '[' in cov:
                sym1 = cov.index('[')
                sym2 = cov.index(']')-2
                if cov.find('[', sym1+1) != -1 or '(' in cov:
                    raise NotImplementedError("multiple symmetries are not " +
                                              "treated yet")
                csym1 = sym1 + self._tensor._tensor_type[0]
                csym2 = sym2 + self._tensor._tensor_type[0]
                self._tensor = self._tensor.antisymmetrize(
                                                        *range(csym1, csym2+1))
                self._changed = True # self does no longer contain the original
                                     # tensor
                cov = cov.replace('[','').replace(']','')
            if len(cov) != tensor._tensor_type[1]:
                raise IndexError("number of covariant indices not " +
                                 "compatible with the tensor type")
            self._cov = cov
        else:
            raise IndexError("too many '_' in the list of indices")

        # Treatment of possible self-contractions:
        # ---------------------------------------
        contraction_pairs = []
        for ind in self._con:
            if ind != '.' and ind in self._cov:
                pos1 = self._con.index(ind)
                pos2 = self._tensor._tensor_type[0] + self._cov.index(ind)
                contraction_pairs.append((pos1, pos2))
        if len(contraction_pairs) > 1:
            raise NotImplementedError("multiple self-contractions are not " +
                                      "implemented yet")
        if len(contraction_pairs) == 1:
            pos1 = contraction_pairs[0][0]
            pos2 = contraction_pairs[0][1]
            self._tensor = self._tensor.trace(pos1, pos2)
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
        """
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
        """
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
            Type-(2,1) tensor a*b on the 3-dimensional vector space M over the
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
                if  ind in other._cov:
                    pos1 = self._con.index(ind)
                    pos2 = other._tensor._tensor_type[0] + other._cov.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if ind in other._con:
                    raise IndexError("the index {} appears twice ".format(ind)
                                     + "in a contravariant position")
        for ind in self._cov:
            if ind != '.':
                if ind in other._con:
                    pos1 = self._tensor._tensor_type[0] + self._cov.index(ind)
                    pos2 = other._con.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if ind in other._cov:
                    raise IndexError("the index {} appears twice ".format(ind)
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
            +a^ij_k
            sage: s._tensor == a
            True

        """
        return TensorWithIndices(+self._tensor,
                                 self._con + '_' + self._cov)

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
            -a^ij_k
            sage: s._tensor == -a
            True

        """
        return TensorWithIndices(-self._tensor,
                                 self._con + '_' + self._cov)
