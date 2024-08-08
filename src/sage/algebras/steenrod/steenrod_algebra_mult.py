# sage.doctest: needs sage.rings.finite_rings
r"""
Multiplication for elements of the Steenrod algebra

AUTHORS:

- John H. Palmieri (2008-07-30: version 0.9) initial version: Milnor
  multiplication.
- John H. Palmieri (2010-06-30: version 1.0) multiplication of
  Serre-Cartan basis elements using the Adem relations.
- Simon King (2011-10-25): Fix the use of cached functions.

.. rubric:: Milnor multiplication, `p=2`

See Milnor's paper [Mil1958]_ for proofs, etc.

To multiply Milnor basis elements `\text{Sq}(r_1, r_2, ...)` and
`\text{Sq}(s_1, s_2,...)` at the prime 2, form all possible matrices
`M` with rows and columns indexed starting at 0, with position (0,0)
deleted (or ignored), with `s_i` equal to the sum of column `i` for
each `i`, and with `r_j` equal to the 'weighted' sum of row `j`.  The
weights are as follows: elements from column `i` are multiplied by
`2^i`.  For example, to multiply `\text{Sq}(2)` and `\text{Sq}(1,1)`,
form the matrices

.. MATH::

  \begin{Vmatrix}
  * & 1 & 1 \\
  2 & 0 & 0
  \end{Vmatrix}
  \quad \text{and} \quad
  \begin{Vmatrix}
  * & 0 & 1 \\
  0 & 1 & 0
  \end{Vmatrix}

(The `*` is the ignored (0,0)-entry of the matrix.)  For each such
matrix `M`, compute a multinomial coefficient, mod 2: for each
diagonal `\{m_{ij}: i+j=n\}`, compute `(\sum m_{i,j}!) / (m_{0,n}!
m_{1,n-1}!  ... m_{n,0}!)`.  Multiply these together for all `n`.  (To
compute this mod 2, view the entries of the matrix as their base 2
expansions; then this coefficient is zero if and only if there is some
diagonal containing two numbers which have a summand in common in
their base 2 expansion.  For example, if 3 and 10 are in the same
diagonal, the coefficient is zero, because `3=1+2` and `10=2+8`: they
both have a summand of 2.)

Now, for each matrix with multinomial coefficient 1, let `t_n` be
the sum of the `n`-th diagonal in the matrix; then

.. MATH::

  \text{Sq}(r_1, r_2, ...) \text{Sq}(s_1, s_2, ...) = \sum \text{Sq}(t_1, t_2, ...)

The function :func:`milnor_multiplication` takes as input two tuples
of nonnegative integers, `r` and `s`, which represent
`\text{Sq}(r)=\text{Sq}(r_1, r_2, ...)` and
`\text{Sq}(s)=\text{Sq}(s_1, s_2, ...)`; it returns as output a
dictionary whose keys are tuples `t=(t_1, t_2, ...)` of nonnegative
integers, and for each tuple the associated value is the coefficient
of `\text{Sq}(t)` in the product formula.  (Since we are working mod 2,
this coefficient is 1 -- if it is zero, the element is omitted from
the dictionary altogether).

.. rubric:: Milnor multiplication, odd primes

As for the `p=2` case, see Milnor's paper [Mil1958]_ for proofs.

Fix an odd prime `p`.  There are three steps to multiply Milnor basis
elements `Q_{f_1} Q_{f_2} ... \mathcal{P}(q_1, q_2, ...)` and
`Q_{g_1} Q_{g_2} ... \mathcal{P}(s_1, s_2,...)`: first, use the formula

.. MATH::

    \mathcal{P}(q_1, q_2, ...) Q_k = Q_k \mathcal{P}(q_1, q_2, ...)
    + Q_{k+1} \mathcal{P}(q_1 - p^k, q_2, ...)
    + Q_{k+2} \mathcal{P}(q_1, q_2 - p^k, ...)
    + ...

Second, use the fact that the `Q_k`'s form an exterior algebra: `Q_k^2 =
0` for all `k`, and if `i \neq j`, then `Q_i` and `Q_j` anticommute:
`Q_i Q_j = -Q_j Q_i`.  After these two steps, the product is a linear
combination of terms of the form

.. MATH::

    Q_{e_1} Q_{e_2} ... \mathcal{P}(r_1, r_2, ...) \mathcal{P}(s_1, s_2, ...).

Finally, use Milnor matrices to multiply the pairs of
`\mathcal{P}(...)` terms, as at the prime 2: form all possible
matrices `M` with rows and columns indexed starting at 0, with
position (0,0) deleted (or ignored), with `s_i` equal to the sum of
column `i` for each `i`, and with `r_j` equal to the weighted sum of
row `j`: elements from column `i` are multiplied by `p^i`.  For
example when `p=5`, to multiply `\mathcal{P}(5)` and
`\mathcal{P}(1,1)`, form the matrices

.. MATH::

    \begin{Vmatrix}
    * & 1 & 1 \\
    5 & 0 & 0
    \end{Vmatrix}
    \quad \text{and} \quad
    \begin{Vmatrix}
    * & 0 & 1 \\
    0 & 1 & 0
    \end{Vmatrix}

For each such matrix `M`, compute a multinomial coefficient, mod `p`:
for each diagonal `\{m_{ij}: i+j=n\}`, compute `(\sum m_{i,j}!) /
(m_{0,n}!  m_{1,n-1}!  ... m_{n,0}!)`.  Multiply these together for
all `n`.

Now, for each matrix with nonzero multinomial coefficient `b_M`, let
`t_n` be the sum of the `n`-th diagonal in the matrix; then

.. MATH::

    \mathcal{P}(r_1, r_2, ...) \mathcal{P}(s_1, s_2, ...)
    = \sum b_M \mathcal{P}(t_1, t_2, ...)

For example when `p=5`, we have

.. MATH::

    \mathcal{P}(5) \mathcal{P}(1,1) = \mathcal{P}(6,1) + 2 \mathcal{P}(0,2).

The function :func:`milnor_multiplication` takes as input two pairs of
tuples of nonnegative integers, `(g,q)` and `(f,s)`, which represent
`Q_{g_1} Q_{g_2} ... \mathcal{P}(q_1, q_2, ...)` and
`Q_{f_1} Q_{f_2} ... \mathcal{P}(s_1, s_2, ...)`.  It returns as output a
dictionary whose keys are pairs of tuples `(e,t)` of nonnegative
integers, and for each tuple the associated value is the coefficient
in the product formula.

.. rubric:: The Adem relations and admissible sequences

If `p=2`, then the mod 2 Steenrod algebra is generated by Steenrod
squares `\text{Sq}^a` for `a \geq 0` (equal to the Milnor basis element
`\text{Sq}(a)`).  The *Adem relations* are as follows: if `a < 2b`,

.. MATH::

    \text{Sq}^a \text{Sq}^b = \sum_{j=0}^{a/2} \binom{b-j-1}{a-2j} \text{Sq}^{a+b-j} \text{Sq}^j

A monomial `\text{Sq}^{i_1} \text{Sq}^{i_2} ... \text{Sq}^{i_n}` is called *admissible* if
`i_k \geq 2 i_{k+1}` for all `k`.  One can use the Adem relations to
show that the admissible monomials span the Steenrod algebra, as a
vector space; with more work, one can show that the admissible
monomials are also linearly independent.  They form the *Serre-Cartan*
basis for the Steenrod algebra.  To multiply a collection of
admissible monomials, concatenate them and see if the result is
admissible.  If it is, you're done.  If not, find the first pair `\text{Sq}^a
\text{Sq}^b` where it fails to be admissible and apply the Adem relations
there.  Repeat with the resulting terms.  One can prove that this
process terminates in a finite number of steps, and therefore gives a
procedure for multiplying elements of the Serre-Cartan basis.

At an odd prime `p`, the Steenrod algebra is generated by the `p`-th
power operations `\mathcal{P}^a` (the same as `\mathcal{P}(a)` in the
Milnor basis) and the Bockstein operation `\beta` (= `Q_0` in the
Milnor basis).  The odd primary *Adem relations* are as follows: if `a
< pb`,

.. MATH::

    \mathcal{P}^a \mathcal{P}^b = \sum_{j=0}^{a/p} (-1)^{a+j}
    \binom{(b-j)(p-1)-1}{a-pj} \mathcal{P}^{a+b-j} \mathcal{P}^j

Also, if `a \leq pb`,

.. MATH::

    \mathcal{P}^a \beta \mathcal{P}^b = \sum_{j=0}^{a/p} (-1)^{a+j}
    \binom{(b-j)(p-1)}{a-pj} \beta \mathcal{P}^{a+b-j} \mathcal{P}^j +
    \sum_{j=0}^{a/p} (-1)^{a+j-1} \binom{(b-j)(p-1)-1}{a-pj-1}
    \mathcal{P}^{a+b-j} \beta \mathcal{P}^j

The *admissible* monomials at an odd prime are products of the form

.. MATH::

    \beta^{\epsilon_0} \mathcal{P}^{s_1} \beta^{\epsilon_1}
    \mathcal{P}^{s_2} ...  \mathcal{P}^{s_n} \beta^{\epsilon_n}

where `s_k \geq \epsilon_{k+1} + p s_{k+1}` for all `k`.  As at the
prime 2, these form a basis for the Steenrod algebra.

The main function for this is :func:`make_mono_admissible`,
which converts a product of Steenrod
squares or `p`-th power operations and Bocksteins into a dictionary
representing a sum of admissible monomials.
"""

#*****************************************************************************
#  Copyright (C) 2008-2010 John H. Palmieri <palmieri@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#*****************************************************************************

from sage.misc.cachefunc import cached_function

# Milnor, p=2


def milnor_multiplication(r,s):
    r"""
    Product of Milnor basis elements r and s at the prime 2.

    INPUT:

    - ``r`` -- tuple of nonnegative integers
    - ``s`` -- tuple of nonnegative integers

    OUTPUT:

    Dictionary of terms of the form (tuple: coeff), where
    'tuple' is a tuple of nonnegative integers and 'coeff' is 1.

    This computes Milnor matrices for the product of `\text{Sq}(r)`
    and `\text{Sq}(s)`, computes their multinomial coefficients, and
    for each matrix whose coefficient is 1, add `\text{Sq}(t)` to the
    output, where `t` is the tuple formed by the diagonals sums from
    the matrix.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import milnor_multiplication
        sage: milnor_multiplication((2,), (1,)) == {(0, 1): 1, (3,): 1}
        True
        sage: sorted(milnor_multiplication((4,), (2,1)).items())
        [((0, 3), 1), ((2, 0, 1), 1), ((6, 1), 1)]
        sage: sorted(milnor_multiplication((2,4), (0,1)).items())
        [((2, 0, 0, 1), 1), ((2, 5), 1)]

    These examples correspond to the following product computations:

    .. MATH::

        \text{Sq}(2) \text{Sq}(1) = \text{Sq}(0,1) + \text{Sq}(3)

        \text{Sq}(4) \text{Sq}(2,1) = \text{Sq}(6,1) + \text{Sq}(0,3) + \text{Sq}(2,0,1)

        \text{Sq}(2,4) \text{Sq}(0,1) = \text{Sq}(2, 5) + \text{Sq}(2, 0, 0, 1)

    This uses the same algorithm Monks does in his Maple package: see
    http://mathweb.scranton.edu/monks/software/Steenrod/steen.html.
    """
    result = {}
    rows = len(r) + 1
    cols = len(s) + 1
    diags = len(r) + len(s)
    # initialize matrix
    M = list(range(rows))
    for i in range(rows):
        M[i] = [0]*cols
    for j in range(1,cols):
        M[0][j] = s[j-1]
    for i in range(1,rows):
        M[i][0] = r[i-1]
        for j in range(1,cols):
            M[i][j] = 0
    found = True
    while found:
        # check diagonals
        n = 1
        okay = 1
        diagonal = [0]*diags
        while n <= diags and okay is not None:
            nth_diagonal = [M[i][n-i] for i in range(max(0,n-cols+1), min(1+n,rows))]
            okay = multinomial(nth_diagonal)
            diagonal[n-1] = okay
            n = n + 1
        if okay is not None:
            i = diags - 1
            while i >= 0 and diagonal[i] == 0:
                i = i - 1
            t = tuple(diagonal[:i+1])
            # reduce mod two:
            if t in result:
                del result[t]
            else:
                result[t] = 1
        # now look for new matrices:
        found = False
        i = 1
        while not found and i < rows:
            sum = M[i][0]
            j = 1
            while not found and j < cols:
                # check to see if column index j is small enough
                if sum >= 2**j:
                    # now check to see if there's anything above this entry
                    # to add to it
                    temp_col_sum = 0
                    for k in range(i):
                        temp_col_sum += M[k][j]
                    if temp_col_sum != 0:
                        found = True
                        for row in range(1,i):
                            M[row][0] = r[row-1]
                            for col in range(1,cols):
                                M[0][col] = M[0][col] + M[row][col]
                                M[row][col] = 0
                        for col in range(1,j):
                            M[0][col] = M[0][col] + M[i][col]
                            M[i][col] = 0
                        M[0][j] = M[0][j] - 1
                        M[i][j] = M[i][j] + 1
                        M[i][0] = sum - 2**j
                    else:
                        sum = sum + M[i][j] * 2**j
                else:
                        sum = sum + M[i][j] * 2**j
                j = j + 1
            i = i + 1
    return result


def multinomial(list):
    r"""
    Multinomial coefficient of list, mod 2.

    INPUT:

    - ``list`` -- list of integers

    OUTPUT: none if the multinomial coefficient is 0, or sum of list if it is 1

    Given the input `[n_1, n_2, n_3, ...]`, this computes the
    multinomial coefficient `(n_1 + n_2 + n_3 + ...)! / (n_1! n_2!
    n_3! ...)`, mod 2.  The method is roughly this: expand each
    `n_i` in binary.  If there is a 1 in the same digit for any `n_i`
    and `n_j` with `i\neq j`, then the coefficient is 0; otherwise, it
    is 1.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import multinomial
        sage: multinomial([1,2,4])
        7
        sage: multinomial([1,2,5])
        sage: multinomial([1,2,12,192,256])
        463

    This function does not compute any factorials, so the following are
    actually reasonable to do::

        sage: multinomial([1,65536])
        65537
        sage: multinomial([4,65535])
        sage: multinomial([32768,65536])
        98304
    """
    old_sum = list[0]
    okay = True
    i = 1
    while okay and i < len(list):
        j = 1
        while okay and j <= min(old_sum, list[i]):
            if j & old_sum == j:
                okay = (j & list[i] == 0)
            j = j << 1
        old_sum = old_sum + list[i]
        i = i + 1
    if okay:
        return old_sum
    else:
        return None

# Milnor, p odd


def milnor_multiplication_odd(m1,m2,p):
    r"""
    Product of Milnor basis elements defined by m1 and m2 at the odd prime p.

    INPUT:

    - ``m1`` -- pair of tuples (e,r), where e is an increasing tuple of
      nonnegative integers and r is a tuple of nonnegative integers
    - ``m2`` -- pair of tuples (f,s), same format as m1
    - ``p`` -- odd prime number

    OUTPUT:

    Dictionary of terms of the form (tuple: coeff), where 'tuple' is
    a pair of tuples, as for r and s, and 'coeff' is an integer mod p.

    This computes the product of the Milnor basis elements
    `Q_{e_1} Q_{e_2} ... P(r_1, r_2, ...)` and
    `Q_{f_1} Q_{f_2} ... P(s_1, s_2, ...)`.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import milnor_multiplication_odd
        sage: sorted(milnor_multiplication_odd(((0,2),(5,)), ((1,),(1,)), 5).items())
        [(((0, 1, 2), (0, 1)), 4), (((0, 1, 2), (6,)), 4)]
        sage: milnor_multiplication_odd(((0,2,4),()), ((1,3),()), 7)
        {((0, 1, 2, 3, 4), ()): 6}
        sage: milnor_multiplication_odd(((0,2,4),()), ((1,5),()), 7)
        {((0, 1, 2, 4, 5), ()): 1}
        sage: sorted(milnor_multiplication_odd(((),(6,)), ((),(2,)), 3).items())
        [(((), (0, 2)), 1), (((), (4, 1)), 1), (((), (8,)), 1)]

    These examples correspond to the following product computations:

    .. MATH::

        p=5: \quad Q_0 Q_2 \mathcal{P}(5) Q_1 \mathcal{P}(1) = 4 Q_0 Q_1 Q_2 \mathcal{P}(0,1) + 4 Q_0 Q_1 Q_2 \mathcal{P}(6)

        p=7: \quad (Q_0 Q_2 Q_4) (Q_1 Q_3) = 6 Q_0 Q_1 Q_2 Q_3 Q_4

        p=7: \quad (Q_0 Q_2 Q_4) (Q_1 Q_5) = Q_0 Q_1 Q_2 Q_3 Q_5

        p=3: \quad \mathcal{P}(6) \mathcal{P}(2) = \mathcal{P}(0,2) + \mathcal{P}(4,1) + \mathcal{P}(8)

    The following used to fail until the trailing zeroes were
    eliminated in p_mono::

        sage: A = SteenrodAlgebra(3)
        sage: a = A.P(0,3); b = A.P(12); c = A.Q(1,2)
        sage: (a+b)*c == a*c + b*c
        True

    Test that the bug reported in :issue:`7212` has been fixed::

        sage: A.P(36,6)*A.P(27,9,81)
        2 P(13,21,83) + P(14,24,82) + P(17,20,83) + P(25,18,83) + P(26,21,82) + P(36,15,80,1) + P(49,12,83) + 2 P(50,15,82) + 2 P(53,11,83) + 2 P(63,15,81)

    Associativity once failed because of a sign error::

        sage: a,b,c = A.Q_exp(0,1), A.P(3), A.Q_exp(1,1)
        sage: (a*b)*c == a*(b*c)
        True

    This uses the same algorithm Monks does in his Maple package to
    iterate through the possible matrices: see
    http://mathweb.scranton.edu/monks/software/Steenrod/steen.html.
    """
    from sage.rings.finite_rings.finite_field_constructor import GF
    F = GF(p)
    (f,s) = m2
    # First compute Q_e0 Q_e1 ... P(r1, r2, ...) Q_f0 Q_f1 ...
    # Store results (as dictionary of pairs of tuples) in 'answer'.
    answer = {m1: F(1)}
    for k in f:
        old_answer = answer
        answer = {}
        for mono in old_answer:
            if k not in mono[0]:
                q_mono = set(mono[0])
                if q_mono:
                    ind = len(q_mono.intersection(range(k,1+max(q_mono))))
                else:
                    ind = 0
                coeff = (-1)**ind * old_answer[mono]
                lst = list(mono[0])
                if ind == 0:
                    lst.append(k)
                else:
                    lst.insert(-ind,k)
                q_mono = tuple(lst)
                p_mono = mono[1]
                answer[(q_mono, p_mono)] = F(coeff)
            for i in range(1,1+len(mono[1])):
                if (k+i not in mono[0]) and (p**k <= mono[1][i-1]):
                    q_mono = set(mono[0])
                    if q_mono:
                        ind = len(q_mono.intersection(range(k+i,1+max(q_mono))))
                    else:
                        ind = 0
                    coeff = (-1)**ind * old_answer[mono]
                    lst = list(mono[0])
                    if ind == 0:
                        lst.append(k+i)
                    else:
                        lst.insert(-ind,k+i)
                    q_mono = tuple(lst)
                    p_mono = list(mono[1])
                    p_mono[i-1] = p_mono[i-1] - p**k

                    # The next two lines were added so that p_mono will not
                    # have trailing zeros. This makes p_mono uniquely
                    # determined by P(*p_mono).

                    while p_mono and p_mono[-1] == 0:
                        p_mono.pop()

                    answer[(q_mono, tuple(p_mono))] = F(coeff)
    # Now for the Milnor matrices.  For each entry '(e,r): coeff' in answer,
    # multiply r with s.  Record coefficient for matrix and multiply by coeff.
    # Store in 'result'.
    if not s:
        result = answer
    else:
        result = {}
        for (e, r) in answer:
            old_coeff = answer[(e,r)]
            # Milnor multiplication for r and s
            rows = len(r) + 1
            cols = len(s) + 1
            diags = len(r) + len(s)
            # initialize matrix
            M = list(range(rows))
            for i in range(rows):
                M[i] = [0]*cols
            for j in range(1,cols):
                M[0][j] = s[j-1]
            for i in range(1,rows):
                M[i][0] = r[i-1]
                for j in range(1,cols):
                    M[i][j] = 0
            found = True
            while found:
                # check diagonals
                n = 1
                coeff = old_coeff
                diagonal = [0]*diags
                while n <= diags and coeff != 0:
                    nth_diagonal = [M[i][n-i] for i in range(max(0,n-cols+1), min(1+n,rows))]
                    coeff = coeff * multinomial_odd(nth_diagonal,p)
                    diagonal[n-1] = sum(nth_diagonal)
                    n = n + 1
                if F(coeff) != 0:
                    i = diags - 1
                    while i >= 0 and diagonal[i] == 0:
                        i = i - 1
                    t = tuple(diagonal[:i+1])
                    if (e,t) in result:
                        result[(e,t)] = F(coeff + result[(e,t)])
                    else:
                        result[(e,t)] = F(coeff)
                    # now look for new matrices:
                found = False
                i = 1
                while not found and i < rows:
                    temp_sum = M[i][0]
                    j = 1
                    while not found and j < cols:
                        # check to see if column index j is small enough
                        if temp_sum >= p**j:
                            # now check to see if there's anything above this entry
                            # to add to it
                            temp_col_sum = 0
                            for k in range(i):
                                temp_col_sum += M[k][j]
                            if temp_col_sum != 0:
                                found = True
                                for row in range(1,i):
                                    M[row][0] = r[row-1]
                                    for col in range(1,cols):
                                        M[0][col] = M[0][col] + M[row][col]
                                        M[row][col] = 0
                                for col in range(1,j):
                                    M[0][col] = M[0][col] + M[i][col]
                                    M[i][col] = 0
                                M[0][j] = M[0][j] - 1
                                M[i][j] = M[i][j] + 1
                                M[i][0] = temp_sum - p**j
                            else:
                                temp_sum += M[i][j] * p**j
                        else:
                            temp_sum += M[i][j] * p**j
                        j = j + 1
                    i = i + 1
    return result


def multinomial_odd(list,p):
    r"""
    Multinomial coefficient of list, mod p.

    INPUT:

    - ``list`` -- list of integers
    - ``p`` -- a prime number

    OUTPUT: associated multinomial coefficient, mod p

    Given the input `[n_1, n_2, n_3, ...]`, this computes the
    multinomial coefficient `(n_1 + n_2 + n_3 + ...)! / (n_1! n_2!
    n_3! ...)`, mod `p`.  The method is this: expand each `n_i` in
    base `p`: `n_i = \sum_j p^j n_{ij}`.  Do the same for the sum of
    the `n_i`'s, which we call `m`: `m = \sum_j p^j m_j`.  Then the
    multinomial coefficient is congruent, mod `p`, to the product of
    the multinomial coefficients `m_j! / (n_{1j}! n_{2j}! ...)`.

    Furthermore, any multinomial coefficient `m! / (n_1! n_2! ...)`
    can be computed as a product of binomial coefficients: it equals

    .. MATH::

       \binom{n_1}{n_1} \binom{n_1 + n_2}{n_2} \binom{n_1 + n_2 + n_3}{n_3} ...

    This is convenient because Sage's binomial function returns
    integers, not rational numbers (as would be produced just by
    dividing factorials).

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import multinomial_odd
        sage: multinomial_odd([1,2,4], 2)
        1
        sage: multinomial_odd([1,2,4], 7)
        0
        sage: multinomial_odd([1,2,4], 11)
        6
        sage: multinomial_odd([1,2,4], 101)
        4
        sage: multinomial_odd([1,2,4], 107)
        105
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    from sage.rings.integer import Integer
    from sage.arith.misc import binomial
    n = sum(list)
    answer = 1
    F = GF(p)
    n_expansion = Integer(n).digits(p)
    list_expansion = [Integer(k).digits(p) for k in list]
    index = 0
    while answer != 0 and index < len(n_expansion):
        multi = F(1)
        partial_sum = 0
        for exp in list_expansion:
            if index < len(exp):
                partial_sum = partial_sum + exp[index]
                multi = F(multi * binomial(partial_sum, exp[index]))
        answer = F(answer * multi)
        index += 1
    return answer

# Adem relations, Serre-Cartan basis, admissible sequences


def binomial_mod2(n,k):
    r"""
    The binomial coefficient `\binom{n}{k}`, computed mod 2.

    INPUT:

    - `n`, `k` -- integers

    OUTPUT: `n` choose `k`, mod 2

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import binomial_mod2
        sage: binomial_mod2(4,2)
        0
        sage: binomial_mod2(5,4)
        1
        sage: binomial_mod2(3 * 32768, 32768)
        1
        sage: binomial_mod2(4 * 32768, 32768)
        0
    """
    if n < k:
        return 0
    elif ((n-k) & k) == 0:
        return 1
    else:
        return 0


def binomial_modp(n,k,p):
    r"""
    The binomial coefficient `\binom{n}{k}`, computed mod `p`.

    INPUT:

    - `n`, `k` -- integers
    - ``p`` -- prime number

    OUTPUT: `n` choose `k`, mod `p`

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import binomial_modp
        sage: binomial_modp(5,2,3)
        1
        sage: binomial_modp(6,2,11)  # 6 choose 2 = 15
        4
    """
    if n < k:
        return 0
    return multinomial_odd([n-k, k], p)


@cached_function
def adem(a, b, c=0, p=2, generic=None):
    r"""
    The mod `p` Adem relations.

    INPUT:

    - `a`, `b`, `c` -- nonnegative integers (optional); corresponding
      to either `P^a P^b` or (if `c` present) to `P^a \beta^b P^c`
    - ``p`` -- positive prime number (default: 2)
    - ``generic`` -- whether to use the generic Steenrod algebra, (default: depends on prime)

    OUTPUT:

    a dictionary representing the mod `p` Adem relations
    applied to `P^a P^b` or (if `c` present) to `P^a \beta^b P^c`.

    The mod `p` Adem relations for the mod `p` Steenrod algebra are as
    follows: if `p=2`, then if `a < 2b`,

    .. MATH::

       \text{Sq}^a \text{Sq}^b = \sum_{j=0}^{a/2} \binom{b-j-1}{a-2j} \text{Sq}^{a+b-j} \text{Sq}^j

    If `p` is odd, then if `a < pb`,

    .. MATH::

       P^a P^b = \sum_{j=0}^{a/p} (-1)^{a+j} \binom{(b-j)(p-1)-1}{a-pj} P^{a+b-j} P^j

    Also for `p` odd, if `a \leq pb`,

    .. MATH::

       P^a \beta P^b = \sum_{j=0}^{a/p} (-1)^{a+j} \binom{(b-j)(p-1)}{a-pj} \beta P^{a+b-j} P^j
           + \sum_{j=0}^{a/p} (-1)^{a+j-1} \binom{(b-j)(p-1)-1}{a-pj-1} P^{a+b-j} \beta P^j

    EXAMPLES:

    If two arguments (`a` and `b`) are given, then computations are
    done mod 2.  If `a \geq 2b`, then the dictionary {(a,b): 1} is
    returned.  Otherwise, the right side of the mod 2 Adem relation
    for `\text{Sq}^a \text{Sq}^b` is returned.  For example, since
    `\text{Sq}^2 \text{Sq}^2 = \text{Sq}^3 \text{Sq}^1`, we have::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import adem
        sage: adem(2,2) # indirect doctest
        {(3, 1): 1}
        sage: adem(4,2)
        {(4, 2): 1}
        sage: adem(4,4) == {(6, 2): 1, (7, 1): 1}
        True

    If `p` is given and is odd, then with two inputs `a` and `b`, the
    Adem relation for `P^a P^b` is computed.  With three inputs `a`,
    `b`, `c`, the Adem relation for `P^a \beta^b P^c` is computed.
    In either case, the keys in the output are all tuples of odd length,
    with ``(i_1, i_2, ..., i_m)`` representing

    .. MATH::

        \beta^{i_1} P^{i_2} \beta^{i_3} P^{i_4} ... \beta^{i_m}

    For instance::

        sage: adem(3,1, p=3)
        {(0, 3, 0, 1, 0): 1}
        sage: adem(3,0,1, p=3)
        {(0, 3, 0, 1, 0): 1}
        sage: adem(1,0,1, p=7)
        {(0, 2, 0): 2}
        sage: adem(1,1,1, p=5) == {(0, 2, 1): 1, (1, 2, 0): 1}
        True
        sage: adem(1,1,2, p=5) == {(0, 3, 1): 1, (1, 3, 0): 2}
        True
    """
    if generic is None:
        generic = (p != 2)
    if not generic:
        if b == 0:
            return {(a,): 1}
        elif a == 0:
            return {(b,): 1}
        elif a >= 2*b:
            return {(a,b): 1}
        result = {}
        for c in range(1 + a//2):
            if binomial_mod2(b-c-1, a-2*c) == 1:
                if c == 0:
                    result[(a+b,)] = 1
                else:
                    result[(a+b-c,c)] = 1
        return result
    # p odd
    if a == 0 and b == 0:
            return {(c,): 1}
    if c == 0:
        bockstein = 0
        A = a
        B = b
    else:
        A = a
        B = c
        bockstein = b # should be 0 or 1
    if A == 0:
        return {(bockstein, B, 0): 1}
    if B == 0:
        return {(0, A, bockstein): 1}
    if bockstein == 0:
        if A >= p*B: # admissible
            return {(0,A,0,B,0): 1}
        result = {}
        for j in range(1 + a//p):
            coeff = (-1)**(A+j) * binomial_modp((B-j) * (p-1) - 1, A - p*j, p)
            if coeff % p != 0:
                if j == 0:
                    result[(0,A+B,0)] = coeff
                else:
                    result[(0,A+B-j,0,j,0)] = coeff
    else:
        if A >= p*B + 1: # admissible
            return {(0,A,1,B,0): 1}
        result = {}
        for j in range(1 + a//p):
            coeff = (-1)**(A+j) * binomial_modp((B-j) * (p-1), A - p*j, p)
            if coeff % p != 0:
                if j == 0:
                    result[(1,A+B,0)] = coeff
                else:
                    result[(1,A+B-j,0,j,0)] = coeff
        for j in range(1 + (a-1)//p):
            coeff = (-1)**(A+j-1) * binomial_modp((B-j) * (p-1) - 1, A - p*j - 1, p)
            if coeff % p != 0:
                if j == 0:
                    result[(0,A+B,1)] = coeff
                else:
                    result[(0,A+B-j,1,j,0)] = coeff
    return result


@cached_function
def make_mono_admissible(mono, p=2, generic=None):
    r"""
    Given a tuple ``mono``, view it as a product of Steenrod
    operations, and return a dictionary giving data equivalent to
    writing that product as a linear combination of admissible
    monomials.

    When `p=2`, the sequence (and hence the corresponding monomial)
    `(i_1, i_2, ...)` is admissible if `i_j \geq 2 i_{j+1}` for all
    `j`.

    When `p` is odd, the sequence `(e_1, i_1, e_2, i_2, ...)` is
    admissible if `i_j \geq e_{j+1} + p i_{j+1}` for all `j`.

    INPUT:

    - ``mono`` -- tuple of nonnegative integers
    - ``p`` -- prime number (default: 2)
    - ``generic`` -- whether to use the generic Steenrod algebra (default: depends on prime)

    OUTPUT:

    Dictionary of terms of the form (tuple: coeff), where
    'tuple' is an admissible tuple of nonnegative integers and
    'coeff' is its coefficient.  This corresponds to a linear
    combination of admissible monomials.  When `p` is odd, each tuple
    must have an odd length: it should be of the form `(e_1, i_1, e_2,
    i_2, ..., e_k)` where each `e_j` is either 0 or 1 and each `i_j`
    is a positive integer: this corresponds to the admissible monomial

    .. MATH::

       \beta^{e_1} \mathcal{P}^{i_2} \beta^{e_2} \mathcal{P}^{i_2} ...
       \mathcal{P}^{i_k} \beta^{e_k}

    ALGORITHM:

    Given `(i_1, i_2, i_3, ...)`, apply the Adem relations to the first
    pair (or triple when `p` is odd) where the sequence is inadmissible,
    and then apply this function recursively to each of the resulting
    tuples `(i_1, ..., i_{j-1}, NEW, i_{j+2}, ...)`, keeping track of
    the coefficients.

    EXAMPLES::

        sage: from sage.algebras.steenrod.steenrod_algebra_mult import make_mono_admissible
        sage: make_mono_admissible((12,)) # already admissible, indirect doctest
        {(12,): 1}
        sage: make_mono_admissible((2,1)) # already admissible
        {(2, 1): 1}
        sage: make_mono_admissible((2,2))
        {(3, 1): 1}
        sage: make_mono_admissible((2, 2, 2))
        {(5, 1): 1}
        sage: make_mono_admissible((0, 2, 0, 1, 0), p=7)
        {(0, 3, 0): 3}

    Test the fix from :issue:`13796`::

        sage: SteenrodAlgebra(p=2, basis='adem').Q(2) * (Sq(6) * Sq(2)) # indirect doctest
        Sq^10 Sq^4 Sq^1 + Sq^10 Sq^5 + Sq^12 Sq^3 + Sq^13 Sq^2
    """
    from sage.rings.finite_rings.finite_field_constructor import GF
    if generic is None:
        generic = p != 2
    F = GF(p)
    if len(mono) == 1:
        return {mono: 1}
    if not generic and len(mono) == 2:
        return adem(*mono, p=p, generic=generic)
    if not generic:
        # check to see if admissible:
        admissible = True
        for j in range(len(mono)-1):
            if mono[j] < 2*mono[j+1]:
                admissible = False
                break
        if admissible:
            return {mono: 1}
        # else j is the first index where admissibility fails
        ans = {}
        y = adem(mono[j], mono[j+1])
        for x in y:
            new = mono[:j] + x + mono[j+2:]
            new = make_mono_admissible(new)
            for m in new:
                if m in ans:
                    ans[m] = ans[m] + y[x] * new[m]
                    if F(ans[m]) == 0:
                        del ans[m]
                else:
                    ans[m] = y[x] * new[m]
        return ans
    # p odd
    # check to see if admissible:
    admissible = True
    for j in range(1, len(mono)-2, 2):
        if mono[j] < mono[j+1] + p*mono[j+2]:
            admissible = False
            break
    if admissible:
        return {mono: 1}
    # else j is the first index where admissibility fails
    ans = {}
    y = adem(*mono[j:j+3], p=p, generic=True)
    for x in y:
        new_x = list(x)
        new_x[0] = mono[j-1] + x[0]
        if len(mono) >= j+3:
            new_x[-1] = mono[j+3] + x[-1]
        if new_x[0] <= 1 and new_x[-1] <= 1:
            new = mono[:j-1] + tuple(new_x) + mono[j+4:]
            new = make_mono_admissible(new, p, generic=True)
            for m in new:
                if m in ans:
                    ans[m] = ans[m] + y[x] * new[m]
                    if F(ans[m]) == 0:
                        del ans[m]
                else:
                    ans[m] = y[x] * new[m]
    return ans
