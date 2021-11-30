r"""
Documentation for the matroids in the catalog

This module contains implementations for many of the functions accessible
through :mod:`matroids. <sage.matroids.matroids_catalog>` and
:mod:`matroids.named_matroids. <sage.matroids.matroids_catalog>`
(type those lines in Sage and hit ``tab`` for a list).

The docstrings include educational information about each named matroid with
the hopes that this class can be used as a reference. However, for a more
comprehensive list of properties we refer to the appendix of [Oxl2011]_.

.. TODO::

    Add optional argument ``groundset`` to each method so users can customize
    the groundset of the matroid. We probably want some means of relabeling to
    accomplish that.

    Add option to specify the field for represented matroids.

AUTHORS:

- Michael Welsh, Stefan van Zwam (2013-04-01): initial version

Functions
=========
"""
# ****************************************************************************
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz >
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import Matrix
from sage.graphs.all import graphs

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.schemes.all import ProjectiveSpace

import sage.matroids.matroid
import sage.matroids.basis_exchange_matroid
from sage.matroids.basis_matroid import BasisMatroid
from sage.matroids.circuit_closures_matroid import CircuitClosuresMatroid
from sage.matroids.linear_matroid import (LinearMatroid, RegularMatroid,
                                          BinaryMatroid, TernaryMatroid,
                                          QuaternaryMatroid)
from sage.matroids.rank_matroid import RankMatroid
from sage.matroids.constructor import Matroid

# The order is the same as in Oxley.


def Q6():
    """
    Return the matroid `Q_6`, represented over `GF(4)`.

    The matroid `Q_6` is a 6-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    four elements. It is the unique relaxation of the rank-3 whirl.
    See [Oxl2011]_, p. 641.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.Q6(); M
        Q6: Quaternary matroid of rank 3 on 6 elements
        sage: setprint(M.hyperplanes())
        [{'a', 'b', 'd'}, {'a', 'c'}, {'a', 'e'}, {'a', 'f'}, {'b', 'c', 'e'},
         {'b', 'f'}, {'c', 'd'}, {'c', 'f'}, {'d', 'e'}, {'d', 'f'},
         {'e', 'f'}]
        sage: M.nonspanning_circuits() == M.noncospanning_cocircuits()
        False
    """
    F = GF(4, 'x')
    x = F.gens()[0]
    A = Matrix(F, [
        [1, 0, 0, 1, 0, 1],
        [0, 1, 0, 1, 1, x],
        [0, 0, 1, 0, 1, 1]
    ])
    M = QuaternaryMatroid(A, 'abcdef')
    M.rename('Q6: ' + repr(M))
    return M


def P6():
    """
    Return the matroid `P_6`, represented as circuit closures.

    The matroid `P_6` is a 6-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    five elements.
    It is the unique relaxation of `Q_6`.
    It is an excluded minor for the class of quaternary matroids.
    See [Oxl2011]_, p. 641.

    EXAMPLES::

        sage: M = matroids.named_matroids.P6(); M
        P6: Matroid of rank 3 on 6 elements with circuit-closures
        {2: {{'a', 'b', 'c'}}, 3: {{'a', 'b', 'c', 'd', 'e', 'f'}}}
        sage: len(set(M.nonspanning_circuits()).difference(M.nonbases())) == 0
        True
        sage: Matroid(matrix=random_matrix(GF(4, 'a'), ncols=5,
        ....:                                          nrows=5)).has_minor(M)
        False
        sage: M.is_valid()
        True
    """
    E = 'abcdef'
    CC = {
        2: ['abc'],
        3: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('P6: ' + repr(M))
    return M


def R6():
    """
    Return the matroid `R_6`, represented over `GF(3)`.

    The matroid `R_6` is a 6-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    three elements.
    It is isomorphic to the 2-sum of two copies of `U_{2, 4}`.
    See [Oxl2011]_, p. 642.

    EXAMPLES::

        sage: M = matroids.named_matroids.R6(); M
        R6: Ternary matroid of rank 3 on 6 elements, type 2+
        sage: M.equals(M.dual())
        True
        sage: M.is_connected()
        True
        sage: M.is_3connected()
        False
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 2, 1],
        [0, 0, 1, 1, 0, 2]
    ])
    M = TernaryMatroid(A, 'abcdef')
    M.rename('R6: ' + repr(M))
    return M


def Fano():
    r"""
    Return the Fano matroid, represented over `GF(2)`.

    The Fano matroid, or Fano plane, or `F_7`, is a 7-element matroid of
    rank-3.
    It is representable over a field if and only if that field has
    characteristic two.
    It is also the projective plane of order two, i.e. `\mathrm{PG}(2, 2)`.
    See [Oxl2011]_, p. 643.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.Fano(); M
        Fano: Binary matroid of rank 3 on 7 elements, type (3, 0)
        sage: setprint(sorted(M.nonspanning_circuits()))
        [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'}, {'b', 'c', 'd'},
         {'b', 'e', 'g'}, {'c', 'f', 'g'}, {'d', 'e', 'f'}]
        sage: M.delete(M.groundset_list()[randrange(0,
        ....:                  7)]).is_isomorphic(matroids.CompleteGraphic(4))
        True
    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 1, 1, 0, 1]
    ])
    M = BinaryMatroid(A, 'abcdefg')
    M.rename('Fano: ' + repr(M))
    return M


def NonFano():
    """
    Return the non-Fano matroid, represented over `GF(3)`

    The non-Fano matroid, or `F_7^-`, is a 7-element matroid of rank-3.
    It is representable over a field if and only if that field has
    characteristic other than two.
    It is the unique relaxation of `F_7`. See [Oxl2011]_, p. 643.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.NonFano(); M
        NonFano: Ternary matroid of rank 3 on 7 elements, type 0-
        sage: setprint(M.nonbases())
        [{'a', 'b', 'f'}, {'a', 'c', 'e'}, {'a', 'd', 'g'}, {'b', 'c', 'd'},
         {'b', 'e', 'g'}, {'c', 'f', 'g'}]
        sage: M.delete('f').is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.delete('g').is_isomorphic(matroids.CompleteGraphic(4))
        False
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 0, 1, 1],
        [0, 0, 1, 1, 1, 0, 1]
    ])
    M = TernaryMatroid(A, 'abcdefg')
    M.rename('NonFano: ' + repr(M))
    return M


def O7():
    """
    Return the matroid `O_7`, represented over `GF(3)`.

    The matroid `O_7` is a 7-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    three elements.
    It is obtained by freely adding a point to any line of `M(K_4)`.
    See [Oxl2011]_, p. 644

    EXAMPLES::

        sage: M = matroids.named_matroids.O7(); M
        O7: Ternary matroid of rank 3 on 7 elements, type 0+
        sage: M.delete('e').is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.tutte_polynomial()
        y^4 + x^3 + x*y^2 + 3*y^3 + 4*x^2 + 5*x*y + 5*y^2 + 4*x + 4*y
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 1, 1, 1, 1],
        [0, 1, 0, 0, 1, 2, 2],
        [0, 0, 1, 1, 0, 1, 0]
    ])
    M = TernaryMatroid(A, 'abcdefg')
    M.rename('O7: ' + repr(M))
    return M


def P7():
    """
    Return the matroid `P_7`, represented over `GF(3)`.

    The matroid `P_7` is a 7-element matroid of rank-3.
    It is representable over a field if and only if that field has at least
    3 elements.
    It is one of two ternary 3-spikes, with the other being `F_7^-`.
    See [Oxl2011]_, p. 644.

    EXAMPLES::

        sage: M = matroids.named_matroids.P7(); M
        P7: Ternary matroid of rank 3 on 7 elements, type 1+
        sage: M.f_vector()
        [1, 7, 11, 1]
        sage: M.has_minor(matroids.CompleteGraphic(4))
        False
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 2, 1, 1, 0],
        [0, 1, 0, 1, 1, 0, 1],
        [0, 0, 1, 1, 0, 1, 1]
    ])
    M = TernaryMatroid(A, 'abcdefg')
    M.rename('P7: ' + repr(M))
    return M

# AG(3, 2) - use AG()


def AG32prime():
    """
    Return the matroid `AG(3, 2)'`, represented as circuit closures.

    The matroid `AG(3, 2)'` is a 8-element matroid of rank-4.
    It is a smallest non-representable matroid.
    It is the unique relaxation of `AG(3, 2)`. See [Oxl2011]_, p. 646.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.AG32prime(); M
        AG(3, 2)': Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'c', 'e', 'g'}, {'a', 'd', 'g', 'h'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'},
             {'b', 'e', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'},
             {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.contract('c').is_isomorphic(matroids.named_matroids.Fano())
        True
        sage: setprint(M.noncospanning_cocircuits())
        [{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
         {'a', 'c', 'd', 'f'}, {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'},
         {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'b', 'd', 'f', 'h'},
         {'b', 'e', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'},
         {'d', 'e', 'f', 'g'}]
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefgh'
    CC = {
        3: ['abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch', 'abed',
            'cfgh', 'bcef', 'adgh', 'acdf', 'begh', 'aceg'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('AG(3, 2)\': ' + repr(M))
    return M


def R8():
    """
    Return the matroid `R_8`, represented over `GF(3)`.

    The matroid `R_8` is a 8-element matroid of rank-4.
    It is representable over a field if and only if the characteristic of that
    field is not two.
    It is the real affine cube. See [Oxl2011]_, p. 646.

    EXAMPLES::

        sage: M = matroids.named_matroids.R8(); M
        R8: Ternary matroid of rank 4 on 8 elements, type 0+
        sage: M.contract(M.groundset_list()[randrange(0,
        ....:            8)]).is_isomorphic(matroids.named_matroids.NonFano())
        True
        sage: M.equals(M.dual())
        True
        sage: M.has_minor(matroids.named_matroids.Fano())
        False
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 2, 1, 1, 1],
        [0, 1, 0, 0, 1, 2, 1, 1],
        [0, 0, 1, 0, 1, 1, 2, 1],
        [0, 0, 0, 1, 1, 1, 1, 2]
    ])
    M = TernaryMatroid(A, 'abcdefgh')
    M.rename('R8: ' + repr(M))
    return M


def F8():
    """
    Return the matroid `F_8`, represented as circuit closures.

    The matroid `F_8` is a 8-element matroid of rank-4.
    It is a smallest non-representable matroid. See [Oxl2011]_, p. 647.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.named_matroids.F8(); M
        F8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'c', 'e', 'g'}, {'a', 'd', 'g', 'h'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'},
             {'c', 'd', 'e', 'h'}, {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: D = get_nonisomorphic_matroids([M.contract(i)
        ....:                                         for i in M.groundset()])
        sage: len(D)
        3
        sage: [N.is_isomorphic(matroids.named_matroids.Fano()) for N in D]
        [...True...]
        sage: [N.is_isomorphic(matroids.named_matroids.NonFano()) for N in D]
        [...True...]
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefgh'
    CC = {
        3: ['abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch', 'abed',
            'cfgh', 'bcef', 'adgh', 'acdf', 'aceg'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('F8: ' + repr(M))
    return M


def Q8():
    """
    Return the matroid `Q_8`, represented as circuit closures.

    The matroid `Q_8` is a 8-element matroid of rank-4.
    It is a smallest non-representable matroid. See [Oxl2011]_, p. 647.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.Q8(); M
        Q8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
             {'a', 'c', 'd', 'f'}, {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'},
             {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'c', 'd', 'e', 'h'},
             {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.flats(3))
        [{'a', 'b', 'c', 'h'}, {'a', 'b', 'd', 'e'}, {'a', 'b', 'f', 'g'},
         {'a', 'c', 'd', 'f'}, {'a', 'c', 'e'}, {'a', 'c', 'g'},
         {'a', 'd', 'g', 'h'}, {'a', 'e', 'f', 'h'}, {'a', 'e', 'g'},
         {'b', 'c', 'd', 'g'}, {'b', 'c', 'e', 'f'}, {'b', 'd', 'f'},
         {'b', 'd', 'h'}, {'b', 'e', 'g'}, {'b', 'e', 'h'}, {'b', 'f', 'h'},
         {'b', 'g', 'h'}, {'c', 'd', 'e', 'h'}, {'c', 'e', 'g'},
         {'c', 'f', 'g', 'h'}, {'d', 'e', 'f', 'g'}, {'d', 'f', 'h'},
         {'e', 'g', 'h'}]
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefgh'
    CC = {
        3: ['abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch', 'abed',
            'cfgh', 'bcef', 'adgh', 'acdf'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('Q8: ' + repr(M))
    return M


def L8():
    """
    Return the matroid `L_8`, represented as circuit closures.

    The matroid `L_8` is a 8-element matroid of rank-4.
    It is representable over all fields with at least five elements.
    It is a cube, yet it is not a tipless spike. See [Oxl2011]_, p. 648.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.L8(); M
        L8: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'h'}, {'a', 'b', 'f', 'g'}, {'a', 'c', 'e', 'g'},
             {'a', 'e', 'f', 'h'}, {'b', 'c', 'd', 'g'}, {'b', 'd', 'f', 'h'},
             {'c', 'd', 'e', 'h'}, {'d', 'e', 'f', 'g'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.equals(M.dual())
        True
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefgh'
    CC = {
        3: ['abfg', 'bcdg', 'defg', 'cdeh', 'aefh', 'abch', 'aceg', 'bdfh'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('L8: ' + repr(M))
    return M


def S8():
    """
    Return the matroid `S_8`, represented over `GF(2)`.

    The matroid `S_8` is a 8-element matroid of rank-4.
    It is representable over a field if and only if that field has
    characteristic two.
    It is the unique deletion of a non-tip element from the binary 4-spike.
    See [Oxl2011]_, p. 648.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.named_matroids.S8(); M
        S8: Binary matroid of rank 4 on 8 elements, type (2, 0)
        sage: M.contract('d').is_isomorphic(matroids.named_matroids.Fano())
        True
        sage: M.delete('d').is_isomorphic(
        ....:                           matroids.named_matroids.Fano().dual())
        False
        sage: M.is_graphic()
        False
        sage: D = get_nonisomorphic_matroids(
        ....:       list(matroids.named_matroids.Fano().linear_coextensions(
        ....:                                                 cosimple=True)))
        sage: len(D)
        2
        sage: [N.is_isomorphic(M) for N in D]
        [...True...]

    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 1, 1, 1, 1]
    ])
    M = BinaryMatroid(A, 'abcdefgh')
    M.rename('S8: ' + repr(M))
    return M


def Vamos():
    """
    Return the Vamos matroid, represented as circuit closures.

    The Vamos matroid, or Vamos cube, or `V_8` is a 8-element matroid of
    rank-4.
    It violates Ingleton's condition for representability over a division
    ring.
    It is not algebraic. See [Oxl2011]_, p. 649.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.Vamos(); M
        Vamos: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.nonbases())
        [{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
         {'c', 'd', 'e', 'f'}, {'e', 'f', 'g', 'h'}]
        sage: M.is_dependent(['c', 'd', 'g', 'h'])
        False
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefgh'
    CC = {
        3: ['abcd', 'abef', 'cdef', 'abgh', 'efgh'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('Vamos: ' + repr(M))
    return M


def T8():
    """
    Return the matroid `T_8`, represented over `GF(3)`.

    The matroid `T_8` is a 8-element matroid of rank-4.
    It is representable over a field if and only if that field has
    characteristic three.
    It is an excluded minor for the dyadic matroids. See [Oxl2011]_, p. 649.

    EXAMPLES::

        sage: M = matroids.named_matroids.T8(); M
        T8: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: M.truncation().is_isomorphic(matroids.Uniform(3, 8))
        True
        sage: M.contract('e').is_isomorphic(matroids.named_matroids.P7())
        True
        sage: M.has_minor(matroids.Uniform(3, 8))
        False

    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 1, 1, 1, 0]
    ])
    M = TernaryMatroid(A, 'abcdefgh')
    M.rename('T8: ' + repr(M))
    return M


def J():
    """
    Return the matroid `J`, represented over `GF(3)`.

    The matroid `J` is a 8-element matroid of rank-4.
    It is representable over a field if and only if that field has at least
    three elements. See [Oxl2011]_, p. 650.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.J(); M
        J: Ternary matroid of rank 4 on 8 elements, type 0-
        sage: setprint(M.truncation().nonbases())
        [{'a', 'b', 'f'}, {'a', 'c', 'g'}, {'a', 'd', 'h'}]
        sage: M.is_isomorphic(M.dual())
        True
        sage: M.has_minor(matroids.CompleteGraphic(4))
        False
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 1, 0, 0],
        [0, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 1, 1, 0, 0, 1]
    ])
    M = TernaryMatroid(A, 'abcdefgh')
    M.rename('J: ' + repr(M))
    return M


def P8():
    """
    Return the matroid `P_8`, represented over `GF(3)`.

    The matroid `P_8` is a 8-element matroid of rank-4.
    It is uniquely representable over all fields of characteristic other than
    two.
    It is an excluded minor for all fields of characteristic two with four or
    more elements. See [Oxl2011]_, p. 650.

    EXAMPLES::

        sage: M = matroids.named_matroids.P8(); M
        P8: Ternary matroid of rank 4 on 8 elements, type 2+
        sage: M.is_isomorphic(M.dual())
        True
        sage: Matroid(matrix=random_matrix(GF(4, 'a'), ncols=5,
        ....:                              nrows=5)).has_minor(M)
        False
        sage: M.bicycle_dimension()
        2

    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 2, 1, 1, 0],
        [0, 1, 0, 0, 1, 1, 0, 1],
        [0, 0, 1, 0, 1, 0, 1, 1],
        [0, 0, 0, 1, 0, 1, 1, 2]
    ])
    M = TernaryMatroid(A, 'abcdefgh')
    M.rename('P8: ' + repr(M))
    return M


def P8pp():
    """
    Return the matroid `P_8^=`, represented as circuit closures.

    The matroid `P_8^=` is a 8-element matroid of rank-4.
    It can be obtained from `P_8` by relaxing the unique pair of disjoint
    circuit-hyperplanes.
    It is an excluded minor for `GF(4)`-representability.
    It is representable over all fields with at least five elements.
    See [Oxl2011]_, p. 651.

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = matroids.named_matroids.P8pp(); M
        P8'': Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'f', 'h'}, {'a', 'c', 'e', 'f'}, {'a', 'c', 'g', 'h'},
             {'a', 'd', 'e', 'g'}, {'b', 'c', 'e', 'g'}, {'b', 'd', 'e', 'h'},
             {'b', 'd', 'f', 'g'}, {'c', 'd', 'f', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: M.is_isomorphic(M.dual())
        True
        sage: len(get_nonisomorphic_matroids([M.contract(i)
        ....:                                        for i in M.groundset()]))
        1
        sage: M.is_valid() # long time
        True

    """
    E = 'abcdefgh'
    CC = {3: ['abfh', 'bceg', 'cdfh', 'adeg', 'acef', 'bdfg', 'acgh', 'bdeh'],
          4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('P8\'\': ' + repr(M))
    return M


def K33dual():
    """
    Return the matroid `M*(K_{3, 3})`, represented over the regular partial
    field.

    The matroid `M*(K_{3, 3})` is a 9-element matroid of rank-4.
    It is an excluded minor for the class of graphic matroids.
    It is the graft matroid of the 4-wheel with every vertex except the hub
    being coloured. See [Oxl2011]_, p. 652.

    EXAMPLES::

        sage: M = matroids.named_matroids.K33dual(); M
        M*(K3, 3): Regular matroid of rank 4 on 9 elements with 81 bases
        sage: any(N.is_3connected()
        ....:     for N in M.linear_extensions(simple=True))
        False
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefghi'
    G = graphs.CompleteBipartiteGraph(3, 3)
    M = Matroid(groundset=E, graph=G, regular=True)
    M = M.dual()
    M.rename('M*(K3, 3): ' + repr(M))
    return M

# AG(2, 3) - use AG()


def TernaryDowling3():
    """
    Return the matroid `Q_3(GF(3)^\times)`, represented over `GF(3)`.

    The matroid `Q_3(GF(3)^\times)` is a 9-element matroid of rank-3.
    It is the rank-3 ternary Dowling geometry.
    It is representable over a field if and only if that field does not have
    characteristic two. See [Oxl2011]_, p. 654.

    EXAMPLES::

        sage: M = matroids.named_matroids.TernaryDowling3(); M
        Q3(GF(3)x): Ternary matroid of rank 3 on 9 elements, type 0-
        sage: len(list(M.linear_subclasses()))
        72
        sage: M.fundamental_cycle('abc', 'd')
        {'a': 2, 'b': 1, 'd': 1}

    """
    A = Matrix(GF(3), [
        [1, 0, 0, 1, 1, 0, 0, 1, 1],
        [0, 1, 0, 2, 1, 1, 1, 0, 0],
        [0, 0, 1, 0, 0, 2, 1, 2, 1]
    ])
    M = TernaryMatroid(A, 'abcdefghi')
    M.rename('Q3(GF(3)x): ' + repr(M))
    return M

# yomcat's progress flag


def CompleteGraphic(n):
    """
    Return the cycle matroid of the complete graph on `n` vertices.

    INPUT:

    - ``n`` -- an integer, the number of vertices of the underlying complete
      graph.

    OUTPUT:

    The graphic matroid associated with the `n`-vertex complete graph.
    This matroid has rank `n - 1`.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.CompleteGraphic(5); M
        M(K5): Graphic matroid of rank 4 on 10 elements
        sage: M.has_minor(matroids.Uniform(2, 4))
        False
        sage: simplify(M.contract(randrange(0,
        ....:                 10))).is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: setprint(M.closure([0, 2, 4, 5]))
        {0, 1, 2, 4, 5, 7}
        sage: M.is_valid()
        True
    """
    M = Matroid(groundset=list(range((n * (n - 1)) // 2)),
                graph=graphs.CompleteGraph(n))
    M.rename('M(K' + str(n) + '): ' + repr(M))
    return M


def Wheel(n, field=None, ring=None):
    r"""
    Return the rank-`n` wheel.

    INPUT:

    - ``n`` -- a positive integer. The rank of the desired matroid.
    - ``ring`` -- any ring. If provided, output will be a linear matroid
      over the ring or field ``ring``. If the ring is `\ZZ`, then output
      will be a regular matroid.
    - ``field`` -- any field. Same as ``ring``, but only fields are allowed.

    OUTPUT:

    The rank-`n` wheel matroid, represented as a regular matroid.

    See [Oxl2011]_, p. 659.

    EXAMPLES::

        sage: M = matroids.Wheel(5); M
        Wheel(5): Regular matroid of rank 5 on 10 elements with 121 bases
        sage: M.tutte_polynomial()
        x^5 + y^5 + 5*x^4 + 5*x^3*y + 5*x^2*y^2 + 5*x*y^3 + 5*y^4 + 10*x^3 +
        15*x^2*y + 15*x*y^2 + 10*y^3 + 10*x^2 + 16*x*y + 10*y^2 + 4*x + 4*y
        sage: M.is_valid()
        True
        sage: M = matroids.Wheel(3)
        sage: M.is_isomorphic(matroids.CompleteGraphic(4))
        True
        sage: M.is_isomorphic(matroids.Wheel(3,field=GF(3)))
        True
        sage: M = matroids.Wheel(3,field=GF(3)); M
        Wheel(3): Ternary matroid of rank 3 on 6 elements, type 0+
    """
    base_ring = ZZ
    if field is not None and ring is not None:
        raise ValueError("only one of ring and field can be specified.")
    if field is not None:
        base_ring = field
        try:
            if not base_ring.is_field():
                raise TypeError("specified ``field`` is not a field.")
        except AttributeError:
            raise TypeError("specified ``field`` is not a field.")
    if ring is not None:
        base_ring = ring
    A = Matrix(base_ring, n, 2 * n, sparse=True)
    for i in range(n):
        A[i, i] = 1
        A[i, n + i] = 1
        if i != 0:
            A[i, i + n - 1] = -1
        else:
            A[i, 2 * n - 1] = -1
    if base_ring is ZZ:
        M = RegularMatroid(A)
    else:
        M = Matroid(A)
    M.rename('Wheel(' + str(n) + '): ' + repr(M))
    return M


def Whirl(n):
    """
    Return the rank-`n` whirl.

    INPUT:

    - ``n`` -- a positive integer. The rank of the desired matroid.

    OUTPUT:

    The rank-`n` whirl matroid, represented as a ternary matroid.

    The whirl is the unique relaxation of the wheel. See [Oxl2011]_, p. 659.

    EXAMPLES::

        sage: M = matroids.Whirl(5); M
        Whirl(5): Ternary matroid of rank 5 on 10 elements, type 0-
        sage: M.is_valid()
        True
        sage: M.tutte_polynomial()
        x^5 + y^5 + 5*x^4 + 5*x^3*y + 5*x^2*y^2 + 5*x*y^3 + 5*y^4 + 10*x^3 +
        15*x^2*y + 15*x*y^2 + 10*y^3 + 10*x^2 + 15*x*y + 10*y^2 + 5*x + 5*y
        sage: M.is_isomorphic(matroids.Wheel(5))
        False
        sage: M = matroids.Whirl(3)
        sage: M.is_isomorphic(matroids.CompleteGraphic(4))
        False

    .. TODO::

        Optional arguments ``ring`` and ``x``, such that the resulting matroid
        is represented over ``ring`` by a reduced matrix like
        ``[-1  0  x]``
        ``[ 1 -1  0]``
        ``[ 0  1 -1]``
    """
    A = Matrix(GF(3), n, 2 * n, sparse=True)
    for i in range(n):
        A[i, i] = 1
        A[i, n + i] = 1
        if i != 0:
            A[i, i + n - 1] = -1
        else:
            A[i, 2 * n - 1] = 1
    M = TernaryMatroid(A)
    M.rename('Whirl(' + str(n) + '): ' + repr(M))
    return M


def Uniform(r, n):
    """
    Return the uniform matroid of rank `r` on `n` elements.

    INPUT:

    - ``r`` -- a nonnegative integer. The rank of the uniform matroid.
    - ``n`` -- a nonnegative integer. The number of elements of the uniform
      matroid.

    OUTPUT:

    The uniform matroid `U_{r,n}`.

    All subsets of size `r` or less are independent; all larger subsets are
    dependent. Representable when the field is sufficiently large. The precise
    bound is the subject of the MDS conjecture from coding theory.
    See [Oxl2011]_, p. 660.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.Uniform(2, 5); M
        U(2, 5): Matroid of rank 2 on 5 elements with circuit-closures
        {2: {{0, 1, 2, 3, 4}}}
        sage: M.dual().is_isomorphic(matroids.Uniform(3, 5))
        True
        sage: setprint(M.hyperplanes())
        [{0}, {1}, {2}, {3}, {4}]
        sage: M.has_line_minor(6)
        False
        sage: M.is_valid()
        True

    Check that bug :trac:`15292` was fixed::

        sage: M = matroids.Uniform(4,4)
        sage: len(M.circuit_closures())
        0
    """
    E = list(range(n))
    if r < n:
        CC = {r: [E]}
    else:
        CC = {}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('U(' + str(r) + ', ' + str(n) + '): ' + repr(M))
    return M

# Other interesting matroids


def PG(n, q, x=None):
    """
    Return the projective geometry of dimension ``n`` over the finite field
    of order ``q``.

    INPUT:

    - ``n`` -- a positive integer. The dimension of the projective space. This
      is one less than the rank of the resulting matroid.
    - ``q`` -- a positive integer that is a prime power. The order of the
      finite field.
    - ``x`` -- (default: ``None``) a string. The name of the generator of a
      non-prime field, used for non-prime fields. If not supplied, ``'x'`` is
      used.

    OUTPUT:

    A linear matroid whose elements are the points of `PG(n, q)`.

    EXAMPLES::

        sage: M = matroids.PG(2, 2)
        sage: M.is_isomorphic(matroids.named_matroids.Fano())
        True
        sage: matroids.PG(5, 4, 'z').size() == (4^6 - 1) / (4 - 1)
        True
        sage: M = matroids.PG(4, 7); M
        PG(4, 7): Linear matroid of rank 5 on 2801 elements represented over
        the Finite Field of size 7
    """
    if x is None:
        x = 'x'
    F = GF(q, x)
    P = ProjectiveSpace(n, F)
    A = Matrix(F, [list(p) for p in P]).transpose()
    M = Matroid(A)
    M.rename('PG(' + str(n) + ', ' + str(q) + '): ' + repr(M))
    return M


def AG(n, q, x=None):
    r"""
    Return the affine geometry of dimension ``n`` over the finite field of
    order ``q``.

    INPUT:

    - ``n`` -- a positive integer. The dimension of the projective space. This
      is one less than the rank of the resulting matroid.
    - ``q`` -- a positive integer that is a prime power. The order of the
      finite field.
    - ``x`` -- (default: ``None``) a string. The name of the generator of a
      non-prime field, used for non-prime fields. If not supplied, ``'x'`` is
      used.

    OUTPUT:

    A linear matroid whose elements are the points of `AG(n, q)`.

    The affine geometry can be obtained from the projective geometry by
    removing a hyperplane.

    EXAMPLES::

        sage: M = matroids.AG(2, 3) \ 8
        sage: M.is_isomorphic(matroids.named_matroids.AG23minus())
        True
        sage: matroids.AG(5, 4, 'z').size() == ((4 ^ 6 - 1) / (4 - 1) -
        ....:                                             (4 ^ 5 - 1)/(4 - 1))
        True
        sage: M = matroids.AG(4, 2); M
        AG(4, 2): Binary matroid of rank 5 on 16 elements, type (5, 0)

    """
    if x is None:
        x = 'x'
    F = GF(q, x)
    P = ProjectiveSpace(n, F)
    A = Matrix(F, [list(p) for p in P if not list(p)[0] == 0]).transpose()
    M = Matroid(A)
    M.rename('AG(' + str(n) + ', ' + str(q) + '): ' + repr(M))
    return M

# Named matroids


def R10():
    """
    Return the matroid `R_{10}`, represented over the regular partial field.

    The matroid `R_{10}` is a 10-element regular matroid of rank-5.
    It is the unique splitter for the class of regular matroids.
    It is the graft matroid of `K_{3, 3}` in which every vertex is coloured.
    See [Oxl2011]_, p. 656.

    EXAMPLES::

        sage: M = matroids.named_matroids.R10(); M
        R10: Regular matroid of rank 5 on 10 elements with 162 bases
        sage: cct = []
        sage: for i in M.circuits():
        ....:      cct.append(len(i))
        sage: Set(cct)
        {4, 6}
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    Check the splitter property::

        sage: matroids.named_matroids.R10().linear_extensions(simple=True)
        []
    """
    A = Matrix(ZZ, [
        [1, 0, 0, 0, 0, -1, 1, 0, 0, 1],
        [0, 1, 0, 0, 0, 1, -1, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, -1, 1, 0],
        [0, 0, 0, 1, 0, 0, 0, 1, -1, 1],
        [0, 0, 0, 0, 1, 1, 0, 0, 1, -1]
    ])
    M = RegularMatroid(A, 'abcdefghij')
    M.rename('R10: ' + repr(M))
    return M


def R12():
    """
    Return the matroid `R_{12}`, represented over the regular partial field.

    The matroid `R_{12}` is a 12-element regular matroid of rank-6.
    It induces a 3-separation in its 3-connected majors within the class of
    regular matroids.
    An excluded minor for the class of graphic or cographic matroids.
    See [Oxl2011]_, p. 657.

    EXAMPLES::

        sage: M = matroids.named_matroids.R12(); M
        R12: Regular matroid of rank 6 on 12 elements with 441 bases
        sage: M.equals(M.dual())
        False
        sage: M.is_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True
    """
    A = Matrix(ZZ, [
        [1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0],
        [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, -1],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 1, -1, -1]
    ])
    M = RegularMatroid(A, 'abcdefghijkl')
    M.rename('R12: ' + repr(M))
    return M


def NonVamos():
    """
    Return the non-Vamos matroid.

    The non-Vamos matroid, or `V_8^+` is an 8-element matroid of rank 4. It is
    a tightening of the Vamos matroid. It is representable over some field.
    See [Oxl2011]_, p. 72, 84.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.NonVamos(); M
        NonVamos: Matroid of rank 4 on 8 elements with circuit-closures
        {3: {{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
             {'c', 'd', 'e', 'f'}, {'c', 'd', 'g', 'h'}, {'e', 'f', 'g', 'h'}},
         4: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'}}}
        sage: setprint(M.nonbases())
        [{'a', 'b', 'c', 'd'}, {'a', 'b', 'e', 'f'}, {'a', 'b', 'g', 'h'},
         {'c', 'd', 'e', 'f'}, {'c', 'd', 'g', 'h'}, {'e', 'f', 'g', 'h'}]
        sage: M.is_dependent(['c', 'd', 'g', 'h'])
        True
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefgh'
    CC = {
        3: ['abcd', 'abef', 'cdef', 'abgh', 'cdgh', 'efgh'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('NonVamos: ' + repr(M))
    return M


def Pappus():
    """
    Return the Pappus matroid.

    The Pappus matroid is a 9-element matroid of rank-3.
    It is representable over a field if and only if that field either has 4
    elements or more than 7 elements.
    It is an excluded minor for the class of GF(5)-representable matroids.
    See [Oxl2011]_, p. 655.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.Pappus(); M
        Pappus: Matroid of rank 3 on 9 elements with circuit-closures
        {2: {{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'},
             {'b', 'd', 'i'}, {'b', 'f', 'g'}, {'c', 'd', 'h'},
             {'c', 'e', 'g'}, {'d', 'e', 'f'}, {'g', 'h', 'i'}},
         3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}}}
        sage: setprint(M.nonspanning_circuits())
        [{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'}, {'b', 'd', 'i'},
         {'b', 'f', 'g'}, {'c', 'd', 'h'}, {'c', 'e', 'g'}, {'d', 'e', 'f'},
         {'g', 'h', 'i'}]
        sage: M.is_dependent(['d', 'e', 'f'])
        True
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefghi'
    CC = {
        2: ['abc', 'def', 'ceg', 'bfg', 'cdh', 'afh', 'bdi', 'aei', 'ghi'],
        3: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('Pappus: ' + repr(M))
    return M


def NonPappus():
    """
    Return the non-Pappus matroid.

    The non-Pappus matroid is a 9-element matroid of rank-3.
    It is not representable over any commutative field.
    It is the unique relaxation of the Pappus matroid. See [Oxl2011]_, p. 655.

    EXAMPLES::

        sage: from sage.matroids.advanced import setprint
        sage: M = matroids.named_matroids.NonPappus(); M
        NonPappus: Matroid of rank 3 on 9 elements with circuit-closures
        {2: {{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'},
             {'b', 'd', 'i'}, {'b', 'f', 'g'}, {'c', 'd', 'h'},
             {'c', 'e', 'g'}, {'g', 'h', 'i'}},
         3: {{'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'}}}
        sage: setprint(M.nonspanning_circuits())
        [{'a', 'b', 'c'}, {'a', 'e', 'i'}, {'a', 'f', 'h'}, {'b', 'd', 'i'},
         {'b', 'f', 'g'}, {'c', 'd', 'h'}, {'c', 'e', 'g'}, {'g', 'h', 'i'}]
        sage: M.is_dependent(['d', 'e', 'f'])
        False
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefghi'
    CC = {
        2: ['abc', 'ceg', 'bfg', 'cdh', 'afh', 'bdi', 'aei', 'ghi'],
        3: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('NonPappus: ' + repr(M))
    return M


def TicTacToe():
    """
    Return the TicTacToe matroid.

    The dual of the TicTacToe matroid is not algebraic; it is unknown whether
    the TicTacToe matroid itself is algebraic. See [Hoc]_.

    EXAMPLES::

        sage: M = matroids.named_matroids.TicTacToe()
        sage: M.is_valid() # long time
        True

    """
    E = 'abcdefghi'
    CC = {
        4: ['abcdg', 'adefg', 'abceh', 'abcfi', 'cdefi', 'adghi',
            'beghi', 'cfghi'],
        5: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('TicTacToe: ' + repr(M))
    return M


def Q10():
    r"""
    Return the matroid `Q_{10}`, represented over `\GF{4}`.

    `Q_{10}` is a 10-element, rank-5, self-dual matroid. It is representable
    over `\GF{3}` and `\GF{4}`, and hence is a sixth-roots-of-unity matroid.
    `Q_{10}` is a splitter for the class of sixth-root-of-unity matroids.

    EXAMPLES::

        sage: M = matroids.named_matroids.Q10()
        sage: M.is_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    Check the splitter property. By Seymour's Theorem, and using self-duality,
    we only need to check that all 3-connected single-element extensions have
    an excluded minor for sixth-roots-of-unity. The only excluded minors that
    are quaternary are `U_{2, 5}, U_{3, 5}, F_7, F_7^*`. As it happens, it
    suffices to check for `U_{2, 5}`:

        sage: S = matroids.named_matroids.Q10().linear_extensions(simple=True)
        sage: [M for M in S if not M.has_line_minor(5)] # long time
        []
    """
    F = GF(4, 'x')
    x = F.gens()[0]
    A = Matrix(F, [
        [1, 0, 0, 0, 0, 1, x, 0, 0, x + 1],
        [0, 1, 0, 0, 0, x + 1, 1, x, 0, 0],
        [0, 0, 1, 0, 0, 0, x + 1, 1, x, 0],
        [0, 0, 0, 1, 0, 0, 0, x + 1, 1, x],
        [0, 0, 0, 0, 1, x, 0, 0, x + 1, 1]
    ])
    M = QuaternaryMatroid(A, 'abcdefghij')
    M.rename('Q10: ' + repr(M))
    return M


def N1():
    r"""
    Return the matroid `N_1`, represented over `\GF{3}`.

    `N_1` is an excluded minor for the dyadic matroids. See [Oxl2011]_, p. 554.

    EXAMPLES::

        sage: M = matroids.named_matroids.N1()
        sage: M.is_field_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 2, 0, 0, 1, 1],
        [0, 1, 0, 0, 0, 1, 2, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 1, 2, 0, 1],
        [0, 0, 0, 1, 0, 0, 0, 1, 2, 2],
        [0, 0, 0, 0, 1, 1, 1, 1, 2, 0]
    ])
    M = TernaryMatroid(A, 'abcdefghij')
    M.rename('N1: ' + repr(M))
    return M


def N2():
    r"""
    Return the matroid `N_2`, represented over `\GF{3}`.

    `N_2` is an excluded minor for the dyadic matroids. See [Oxl2011]_, p. 554.

    EXAMPLES::

        sage: M = matroids.named_matroids.N2()
        sage: M.is_field_isomorphic(M.dual())
        True
        sage: M.is_valid()
        True

    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1],
        [0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 1, 2, 2, 1, 0, 1]
    ])
    M = TernaryMatroid(A, 'abcdefghijkl')
    M.rename('N2: ' + repr(M))
    return M


def BetsyRoss():
    """
    Return the Betsy Ross matroid, represented by circuit closures.

    An extremal golden-mean matroid.
    That is, if `M` is simple, rank 3, has the Betsy Ross matroid as a
    restriction and is a Golden Mean matroid, then `M` is the Betsy Ross
    matroid.

    EXAMPLES::

        sage: M = matroids.named_matroids.BetsyRoss()
        sage: len(M.circuit_closures()[2])
        10
        sage: M.is_valid() # long time
        True
    """
    E = 'abcdefghijk'
    CC = {
        2: ['acfg', 'bdgh', 'cehi', 'befj', 'adij', 'dfk',
            'egk', 'ahk', 'bik', 'cjk'],
        3: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('BetsyRoss: ' + repr(M))
    return M


def Block_9_4():
    """
    Return the paving matroid whose non-spanning circuits form the blocks of a
    `2-(9, 4, 3)` design.

    EXAMPLES::

        sage: M = matroids.named_matroids.Block_9_4()
        sage: M.is_valid() # long time
        True
        sage: BD = BlockDesign(M.groundset(), M.nonspanning_circuits())
        sage: BD.is_t_design(return_parameters=True)
        (True, (2, 9, 4, 3))
    """
    E = 'abcdefghi'
    CC = {
        3: ['abcd', 'acef', 'bdef', 'cdeg', 'abfg', 'adeh', 'bcfh', 'acgh',
            'begh', 'dfgh', 'abei', 'cdfi', 'bcgi', 'adgi', 'efgi', 'bdhi',
            'cehi', 'afhi'],
        4: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('Block(9, 4): ' + repr(M))
    return M


def Block_10_5():
    """
    Return the paving matroid whose non-spanning circuits form the blocks of a
    `3-(10, 5, 3)` design.

    EXAMPLES::

        sage: M = matroids.named_matroids.Block_10_5()
        sage: M.is_valid() # long time
        True
        sage: BD = BlockDesign(M.groundset(), M.nonspanning_circuits())
        sage: BD.is_t_design(return_parameters=True)
        (True, (3, 10, 5, 3))
    """

    E = 'abcdefghij'
    CC = {
        4: ['abcde', 'acdfg', 'bdefg', 'bcdfh', 'abefh', 'abcgh', 'adegh',
            'cefgh', 'bcefi', 'adefi', 'bcdgi', 'acegi', 'abfgi', 'abdhi',
            'cdehi', 'acfhi', 'beghi', 'dfghi', 'abdfj', 'acefj', 'abegj',
            'cdegj', 'bcfgj', 'acdhj', 'bcehj', 'defhj', 'bdghj', 'afghj',
            'abcij', 'bdeij', 'cdfij', 'adgij', 'efgij', 'aehij', 'bfhij',
            'cghij'],
        5: [E]
    }
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('Block(10, 5): ' + repr(M))
    return M


def ExtendedBinaryGolayCode():
    """
    Return the matroid of the extended binary Golay code.

    See
    :class:`GolayCode <sage.coding.golay_code.GolayCode>`
    documentation for more on this code.

    EXAMPLES::

        sage: M = matroids.named_matroids.ExtendedBinaryGolayCode()
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(2))) # long time
        True
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(2), [
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1],
        [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    ])
    M = BinaryMatroid(A, 'abcdefghijklmnopqrstuvwx')
    M.rename('Extended Binary Golay Code: ' + repr(M))
    return M


def ExtendedTernaryGolayCode():
    """
    Return the matroid of the extended ternary Golay code.

    See
    :class:`GolayCode <sage.coding.golay_code.GolayCode>`

    EXAMPLES::

        sage: M = matroids.named_matroids.ExtendedTernaryGolayCode()
        sage: C = LinearCode(M.representation())
        sage: C.is_permutation_equivalent(codes.GolayCode(GF(3))) # long time
        True
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 0],
        [0, 1, 0, 0, 0, 0, 1, 1, 2, 1, 0, 2],
        [0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 2],
        [0, 0, 0, 1, 0, 0, 1, 2, 0, 1, 2, 1],
        [0, 0, 0, 0, 1, 0, 1, 0, 2, 2, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1]
    ])
    M = TernaryMatroid(A, 'abcdefghijkl')
    M.rename('Extended Ternary Golay Code: ' + repr(M))
    return M


def AG23minus():
    """
    Return the ternary affine plane minus a point.

    This is a sixth-roots-of-unity matroid, and an excluded minor for the
    class of near-regular matroids.
    See [Oxl2011]_, p. 653.

    EXAMPLES::

        sage: M = matroids.named_matroids.AG23minus()
        sage: M.is_valid()
        True

    """
    E = 'abcdefgh'
    CC = {2: ['abc', 'ceh', 'fgh', 'adf', 'aeg', 'cdg', 'bdh', 'bef'],
          3: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('AG23minus: ' + repr(M))
    return M


def NotP8():
    """
    Return the matroid ``NotP8``.

    This is a matroid that is not `P_8`, found on page 512 of [Oxl1992]_ (the
    first edition).

    EXAMPLES::

        sage: M = matroids.named_matroids.P8()
        sage: N = matroids.named_matroids.NotP8()
        sage: M.is_isomorphic(N)
        False
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(3), [
        [1, 0, 0, 0, 0, 1, 1, -1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, -1, 1, 1, 1]
    ])
    M = TernaryMatroid(A, 'abcdefgh')
    M.rename('NotP8: ' + repr(M))
    return M


def D16():  # A.K.A. the Carolyn Chun Matroid
    """
    Return the matroid `D_{16}`.

    Let `M` be a 4-connected binary matroid and `N` an internally 4-connected
    proper minor of `M` with at least 7 elements. Then some element of `M` can
    be deleted or contracted preserving an `N`-minor, unless `M` is `D_{16}`.
    See [CMO2012]_.

    EXAMPLES::

        sage: M = matroids.named_matroids.D16()
        sage: M
        D16: Binary matroid of rank 8 on 16 elements, type (0, 0)
        sage: M.is_valid()
        True

    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1],
        [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0]
    ])
    M = BinaryMatroid(A, 'abcdefghijklmnop')
    M.rename('D16: ' + repr(M))
    return M


def Terrahawk():  # A.K.A. the Dillon Mayhew Matroid
    """
    Return the Terrahawk matroid.

    The Terrahawk is a binary matroid that is a sporadic exception in a chain
    theorem for internally 4-connected binary matroids. See [CMO2011]_.

    EXAMPLES::

        sage: M = matroids.named_matroids.Terrahawk()
        sage: M
        Terrahawk: Binary matroid of rank 8 on 16 elements, type (0, 4)
        sage: M.is_valid()
        True

    """
    A = Matrix(GF(2), [
        [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0]
    ])
    M = BinaryMatroid(A, 'abcdefghijklmnop')
    M.rename('Terrahawk: ' + repr(M))
    return M


def R9A():
    """
    Return the matroid `R_9^A`.

    The matroid `R_9^A` is not representable over any field, yet none of the
    cross-ratios in its Tuttegroup equal 1. It is one of the 4 matroids on at
    most 9 elements with this property, the others being `{R_9^A}^*`, `R_9^B`
    and `{R_9^B}^*`.

    EXAMPLES::

        sage: M = matroids.named_matroids.R9A()
        sage: M.is_valid() # long time
        True

    """
    E = 'abcdefghi'
    CC = {3: ['abde', 'bcdf', 'aceg', 'abch', 'aefh', 'adgh', 'acdi', 'abfi',
              'defi', 'begi', 'bdhi', 'cehi', 'fghi'],
          4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('R9A: ' + repr(M))
    return M


def R9B():
    """
    Return the matroid `R_9^B`.

    The matroid `R_9^B` is not representable over any field, yet none of the
    cross-ratios in its Tuttegroup equal 1.
    It is one of the 4 matroids on at most 9 elements with this property, the
    others being `{R_9^B}^*`, `R_9^A` and `{R_9^A}^*`.

    EXAMPLES::

        sage: M = matroids.named_matroids.R9B()
        sage: M.is_valid() # long time
        True

    """
    E = 'abcdefghi'
    CC = {3: ['abde', 'bcdf', 'aceg', 'abch', 'befh', 'cdgh', 'bcei', 'adfi',
              'abgi', 'degi', 'bdhi', 'aehi', 'fghi'],
          4: [E]}
    M = CircuitClosuresMatroid(groundset=E, circuit_closures=CC)
    M.rename('R9B: ' + repr(M))
    return M


def T12():
    """
    Return the matroid `T_{12}`.

    The edges of the Petersen graph can be labeled by the 4-circuits of
    `T_{12}` so that two edges are adjacent if and only if the corresponding
    4-circuits overlap in exactly two elements.
    Relaxing a circuit-hyperplane yields an excluded minor for the class of
    matroids that are either binary or ternary. See [Oxl2011]_, p. 658.

    EXAMPLES::

        sage: M = matroids.named_matroids.T12()
        sage: M
        T12: Binary matroid of rank 6 on 12 elements, type (2, None)
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        [0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1]
    ])
    M = BinaryMatroid(A, 'abcdefghijkl')
    M.rename('T12: ' + repr(M))
    return M


def P9():
    """
    Return the matroid `P_9`.

    This is the matroid referred to as `P_9` by Oxley in his paper "The binary
    matroids with no 4-wheel minor"

    EXAMPLES::

        sage: M = matroids.named_matroids.P9()
        sage: M
        P9: Binary matroid of rank 4 on 9 elements, type (1, 1)
        sage: M.is_valid()
        True
    """
    A = Matrix(GF(2), [
        [1, 0, 0, 0, 1, 0, 0, 1, 1],
        [0, 1, 0, 0, 1, 1, 0, 0, 1],
        [0, 0, 1, 0, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 0, 0, 1, 1, 0]
    ])
    M = BinaryMatroid(A, 'abcdefghi')
    M.rename('P9: ' + repr(M))
    return M
