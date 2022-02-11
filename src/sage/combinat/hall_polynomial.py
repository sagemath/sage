r"""
Hall Polynomials
"""
#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>,
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

from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.combinat.partition import Partition
from sage.combinat.q_analogues import q_binomial


def hall_polynomial(nu, mu, la, q=None):
    r"""
    Return the (classical) Hall polynomial `P^{\nu}_{\mu,\lambda}`
    (where `\nu`, `\mu` and `\lambda` are the inputs ``nu``, ``mu``
    and ``la``).

    Let `\nu,\mu,\lambda` be partitions. The Hall polynomial
    `P^{\nu}_{\mu,\lambda}(q)` (in the indeterminate `q`) is defined
    as follows: Specialize `q` to a prime power, and consider the
    category of `\GF{q}`-vector spaces with a distinguished
    nilpotent endomorphism. The morphisms in this category shall be
    the linear maps commuting with the distinguished endomorphisms.
    The *type* of an object in the category will be the Jordan type
    of the distinguished endomorphism; this is a partition. Now, if
    `N` is any fixed object of type `\nu` in this category, then
    the polynomial `P^{\nu}_{\mu,\lambda}(q)` specialized at the
    prime power `q` counts the number of subobjects `L` of `N` having
    type `\lambda` such that the quotient object `N / L` has type
    `\mu`. This determines the values of the polynomial
    `P^{\nu}_{\mu,\lambda}` at infinitely many points (namely, at all
    prime powers), and hence `P^{\nu}_{\mu,\lambda}` is uniquely
    determined. The degree of this polynomial is at most
    `n(\nu) - n(\lambda) - n(\mu)`, where
    `n(\kappa) = \sum_i (i-1) \kappa_i` for every partition `\kappa`.
    (If this is negative, then the polynomial is zero.)

    These are the structure coefficients of the
    :class:`(classical) Hall algebra <HallAlgebra>`.

    If `\lvert \nu \rvert \neq \lvert \mu \rvert + \lvert \lambda \rvert`,
    then we have `P^{\nu}_{\mu,\lambda} = 0`. More generally, if the
    Littlewood-Richardson coefficient `c^{\nu}_{\mu,\lambda}`
    vanishes, then `P^{\nu}_{\mu,\lambda} = 0`.

    The Hall polynomials satisfy the symmetry property
    `P^{\nu}_{\mu,\lambda} = P^{\nu}_{\lambda,\mu}`.

    ALGORITHM:

    If `\lambda = (1^r)` and
    `\lvert \nu \rvert = \lvert \mu \rvert + \lvert \lambda \rvert`,
    then we compute `P^{\nu}_{\mu,\lambda}` as follows (cf. Example 2.4
    in [Sch2006]_):

    First, write `\nu = (1^{l_1}, 2^{l_2}, \ldots, n^{l_n})`, and
    define a sequence `r = r_0 \geq r_1 \geq \cdots \geq r_n` such that

    .. MATH::

        \mu = \left( 1^{l_1 - r_0 + 2r_1 - r_2}, 2^{l_2 - r_1 + 2r_2 - r_3},
        \ldots , (n-1)^{l_{n-1} - r_{n-2} + 2r_{n-1} - r_n},
        n^{l_n - r_{n-1} + r_n} \right).

    Thus if `\mu = (1^{m_1}, \ldots, n^{m_n})`, we have the following system
    of equations:

    .. MATH::

        \begin{aligned}
        m_1 & = l_1 - r + 2r_1 - r_2,
        \\ m_2 & = l_2 - r_1 + 2r_2 - r_3,
        \\ & \vdots ,
        \\ m_{n-1} & = l_{n-1} - r_{n-2} + 2r_{n-1} - r_n,
        \\ m_n & = l_n - r_{n-1} + r_n
        \end{aligned}

    and solving for `r_i` and back substituting we obtain the equations:

    .. MATH::

        \begin{aligned}
        r_n & = r_{n-1} + m_n - l_n,
        \\ r_{n-1} & = r_{n-2} + m_{n-1} - l_{n-1} + m_n - l_n,
        \\ & \vdots ,
        \\ r_1 & = r + \sum_{k=1}^n (m_k - l_k),
        \end{aligned}

    or in general we have the recursive equation:

    .. MATH::

        r_i = r_{i-1} + \sum_{k=i}^n (m_k - l_k).

    This, combined with the condition that `r_0 = r`, determines the
    `r_i` uniquely (recursively). Next we define

    .. MATH::

        t = (r_{n-2} - r_{n-1})(l_n - r_{n-1})
        + (r_{n-3} - r_{n-2})(l_{n-1} + l_n - r_{n-2}) + \cdots
        + (r_0 - r_1)(l_2 + \cdots + l_n - r_1),

    and with these notations we have

    .. MATH::

        P^{\nu}_{\mu,(1^r)} = q^t \binom{l_n}{r_{n-1}}_q
        \binom{l_{n-1}}{r_{n-2} - r_{n-1}}_q \cdots \binom{l_1}{r_0 - r_1}_q.

    To compute `P^{\nu}_{\mu,\lambda}` in general, we compute the product
    `I_{\mu} I_{\lambda}` in the Hall algebra and return the coefficient of
    `I_{\nu}`.

    EXAMPLES::

        sage: from sage.combinat.hall_polynomial import hall_polynomial
        sage: hall_polynomial([1,1],[1],[1])
        q + 1
        sage: hall_polynomial([2],[1],[1])
        1
        sage: hall_polynomial([2,1],[2],[1])
        q
        sage: hall_polynomial([2,2,1],[2,1],[1,1])
        q^2 + q
        sage: hall_polynomial([2,2,2,1],[2,2,1],[1,1])
        q^4 + q^3 + q^2
        sage: hall_polynomial([3,2,2,1], [3,2], [2,1])
        q^6 + q^5
        sage: hall_polynomial([4,2,1,1], [3,1,1], [2,1])
        2*q^3 + q^2 - q - 1
        sage: hall_polynomial([4,2], [2,1], [2,1], 0)
        1

    TESTS::

        sage: hall_polynomial([3], [1], [1], 0)
        0
    """
    if q is None:
        q = ZZ['q'].gen()
    R = q.parent()

    # Make sure they are partitions
    nu = Partition(nu)
    mu = Partition(mu)
    la = Partition(la)

    if sum(nu) != sum(mu) + sum(la):
        return R.zero()

    if all(x == 1 for x in la):
        r = [len(la)]   # r will be [r_0, r_1, ..., r_n].
        exp_nu = nu.to_exp()  # exp_nu == [l_1, l_2, ..., l_n].
        exp_mu = mu.to_exp()  # exp_mu == [m_1, m_2, ..., m_n].
        n = max(len(exp_nu), len(exp_mu))
        for k in range(n):
            r.append(r[-1] + sum(exp_mu[k:]) - sum(exp_nu[k:]))
        # Now, r is [r_0, r_1, ..., r_n].
        exp_nu += [0]*(n - len(exp_nu)) # Pad with 0's until it has length n
        # Note that all -1 for exp_nu is due to indexing
        t = sum((r[k-2] - r[k-1])*(sum(exp_nu[k-1:]) - r[k-1]) for k in range(2,n+1))
        if t < 0:
            # This case needs short-circuiting, since otherwise q**-t
            # might throw an exception if q is non-invertible.
            return R.zero()
        return q**t * q_binomial(exp_nu[n-1], r[n-1], q) \
               * prod([q_binomial(exp_nu[k-1], r[k-1] - r[k], q)
                       for k in range(1, n)], R.one())

    from sage.algebras.hall_algebra import HallAlgebra
    H = HallAlgebra(R, q)
    return (H[mu]*H[la]).coefficient(nu)

