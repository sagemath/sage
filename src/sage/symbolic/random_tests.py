"""
Randomized tests of GiNaC / PyNaC
"""

###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################


from sage.misc.prandom import randint, random
import operator
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from sage.symbolic.expression import symbol_table, mixed_order
from sage.symbolic.constants import (pi, e, golden_ratio, log2, euler_gamma,
                                     catalan, khinchin, twinprime, mertens)
from sage.functions.hypergeometric import hypergeometric
from sage.functions.other import (cases, element_of)

###################################################################
### Generate random expressions for doctests ######################
###################################################################

def _mk_full_functions():
    r"""
    A simple function that returns a list of all Pynac functions of known
    arity, sorted by name.

    EXAMPLES::

        sage: from sage.symbolic.random_tests import _mk_full_functions
        sage: [f for (one,f,arity) in _mk_full_functions()] # random
        [Ei, abs, arccos, arccosh, arccot, arccoth, arccsc, arccsch,
        arcsec, arcsech, arcsin, arcsinh, arctan, arctan2, arctanh,
        arg, beta, binomial, ceil, conjugate, cos, cosh, cot, coth,
        csc, csch, dickman_rho, dilog, dirac_delta, elliptic_e,
        elliptic_ec, elliptic_eu, elliptic_f, elliptic_kc,
        elliptic_pi, erf, exp, factorial, floor, heaviside, imag_part,
        integrate, kronecker_delta, log, polylog, real_part, sec,
        sech, sgn, sin, sinh, tan, tanh, unit_step, zeta, zetaderiv]

    Note that this doctest will produce different output whenever a
    symbolic function is added or removed.
    """
    excluded = [hypergeometric, cases, element_of]
    items = sorted(symbol_table['functions'].items())
    return [(1.0, f, f.number_of_arguments())
            for (name, f) in items
            if hasattr(f, 'number_of_arguments') and
               f.number_of_arguments() > 0 and
               f not in excluded]

# For creating simple expressions

fast_binary = [(0.4, operator.add), (0.1, operator.sub), (0.5, operator.mul)]
fast_unary = [(0.8, operator.neg), (0.2, operator.abs)]
fast_nodes = [(0.9, fast_binary, 2), (0.1, fast_unary, 1)]

# For creating expressions with the full power of Pynac's simple expression
# subset (with no quantifiers/operators; that is, no derivatives, integrals,
# etc.)
full_binary = [(0.3, operator.add), (0.1, operator.sub), (0.3, operator.mul), (0.2, operator.truediv), (0.1, operator.pow)]
full_unary = [(0.8, operator.neg), (0.2, operator.inv)]
full_functions = _mk_full_functions()
full_nullary = [(1.0, c) for c in [pi, e]] + [(0.05, c) for c in
        [golden_ratio, log2, euler_gamma, catalan, khinchin, twinprime,
            mertens]]
full_internal = [(0.6, full_binary, 2), (0.2, full_unary, 1),
        (0.2, full_functions)]

def normalize_prob_list(pl, extra=()):
    r"""
    INPUT:

    - ``pl`` - A list of tuples, where the first element of each tuple is
      a floating-point number (representing a relative probability).  The
      second element of each tuple may be a list or any other kind of object.

    - ``extra`` - A tuple which is to be appended to every tuple in ``pl``.

    This function takes such a list of tuples (a "probability list") and
    normalizes the probabilities so that they sum to one.  If any of the
    values are lists, then those lists are first normalized; then
    the probabilities in the list are multiplied by the main probability
    and the sublist is merged with the main list.

    For example, suppose we want to select between group A and group B with
    50% probability each.  Then within group A, we select A1 or A2 with 50%
    probability each (so the overall probability of selecting A1 is 25%);
    and within group B, we select B1, B2, or B3 with probabilities in
    a 1:2:2 ratio.

    EXAMPLES::

        sage: from sage.symbolic.random_tests import *
        sage: A = [(0.5, 'A1'), (0.5, 'A2')]
        sage: B = [(1, 'B1'), (2, 'B2'), (2, 'B3')]
        sage: top = [(50, A, 'Group A'), (50, B, 'Group B')]
        sage: normalize_prob_list(top)
        [(0.250000000000000, 'A1', 'Group A'), (0.250000000000000, 'A2', 'Group A'), (0.1, 'B1', 'Group B'), (0.2, 'B2', 'Group B'), (0.2, 'B3', 'Group B')]
    """
    if len(pl) == 0:
        return pl
    result = []
    total = sum([float(p[0]) for p in pl])
    for p in pl:
        prob = p[0]
        val = p[1]
        if len(p) > 2:
            p_extra = p[2:]
        else:
            p_extra = extra
        if isinstance(val, list):
            norm_val = normalize_prob_list(val, extra=p_extra)
            for np in norm_val:
                result.append(((prob/total)*np[0], np[1]) + np[2:])
        else:
            result.append(((prob/total), val) + p_extra)
    return result

def choose_from_prob_list(lst):
    r"""
    INPUT:

    - ``lst`` - A list of tuples, where the first element of each tuple
      is a nonnegative float (a probability), and the probabilities sum
      to one.

    OUTPUT:

    A tuple randomly selected from the list according to the given
    probabilities.

    EXAMPLES::

        sage: from sage.symbolic.random_tests import *
        sage: v = [(0.1, False), (0.9, True)]
        sage: choose_from_prob_list(v)  # random
        (0.900000000000000, True)
        sage: true_count = 0
        sage: total_count = 0
        sage: def more_samples():
        ....:     global true_count, total_count
        ....:     for _ in range(10000):
        ....:         total_count += 1.0
        ....:         if choose_from_prob_list(v)[1]:
        ....:             true_count += 1.0
        sage: more_samples()
        sage: while abs(true_count/total_count - 0.9) > 0.01:
        ....:     more_samples()
    """
    r = random()
    for i in range(len(lst)-1):
        if r < lst[i][0]:
            return lst[i]
        r -= lst[i][0]
    return lst[-1]

def random_integer_vector(n, length):
    r"""
    Give a random list of length *length*, consisting of nonnegative
    integers that sum to *n*.

    This is an approximation to IntegerVectors(n, length).random_element().
    That gives values uniformly at random, but might be slow; this
    routine is not uniform, but should always be fast.

    (This routine is uniform if ``length`` is 1 or 2; for longer vectors,
    we prefer approximately balanced vectors, where all the values
    are around `n/{length}`.)

    EXAMPLES::

        sage: from sage.symbolic.random_tests import *
        sage: a = random_integer_vector(100, 2); a  # random
        [11, 89]
        sage: len(a)
        2
        sage: sum(a)
        100

        sage: b = random_integer_vector(10000, 20)
        sage: len(b)
        20
        sage: sum(b)
        10000

    The routine is uniform if ``length`` is 2::

        sage: true_count = 0
        sage: total_count = 0
        sage: def more_samples():
        ....:     global true_count, total_count
        ....:     for _ in range(1000):
        ....:         total_count += 1.0
        ....:         if a == random_integer_vector(100, 2):
        ....:             true_count += 1.0
        sage: more_samples()
        sage: while abs(true_count/total_count - 0.01) > 0.01:
        ....:     more_samples()
    """
    if length == 0:
        return []
    elif length == 1:
        return [n]
    elif length == 2:
        v = randint(0, n)
        return [v, n-v]
    else:
        v = randint(0, 2*n//length)
        return [v] + random_integer_vector(n-v, length-1)

def random_expr_helper(n_nodes, internal, leaves, verbose):
    r"""
    Produce a random symbolic expression of size *n_nodes* (or slightly
    larger).  Internal nodes are selected from the *internal* probability
    list; leaves are selected from *leaves*.  If *verbose* is True,
    then a message is printed before creating an internal node.

    EXAMPLES::

        sage: from sage.symbolic.random_tests import *
        sage: a = random_expr_helper(9, [(0.5, operator.add, 2),
        ....:     (0.5, operator.neg, 1)], [(0.5, 1), (0.5, x)], True)
        About to apply <built-in function ...

    In small cases we will see all cases quickly::

        sage: def next_expr():
        ....:     return random_expr_helper(
        ....:         6, [(0.5, operator.add, 2), (0.5, operator.neg, 1)],
        ....:         [(0.5, 1), (0.5, x)], False)
        sage: all_exprs = set()
        sage: for a in range(-4, 5):
        ....:     for b in range(-4+abs(a), 5-abs(a)):
        ....:         if a % 2 and abs(a) + abs(b) == 4 and sign(a) != sign(b):
        ....:             continue
        ....:         all_exprs.add(a*x + b)
        sage: our_exprs = set()
        sage: while our_exprs != all_exprs:
        ....:    our_exprs.add(next_expr())

    """
    if n_nodes == 1:
        return choose_from_prob_list(leaves)[1]
    else:
        r = choose_from_prob_list(internal)
        n_nodes -= 1
        n_children = r[2]
        n_spare_nodes = n_nodes - n_children
        if n_spare_nodes <= 0:
            n_spare_nodes = 0
        nodes_per_child = random_integer_vector(n_spare_nodes, n_children)
        children = [random_expr_helper(n+1, internal, leaves, verbose) for n in nodes_per_child]
        if verbose:
            print("About to apply %r to %r" % (r[1], children))
        return r[1](*children)

def random_expr(size, nvars=1, ncoeffs=None, var_frac=0.5,
                internal=full_internal,
                nullary=full_nullary, nullary_frac=0.2,
                coeff_generator=QQ.random_element, verbose=False):
    r"""
    Produce a random symbolic expression of the given size.  By
    default, the expression involves (at most) one variable, an arbitrary
    number of coefficients, and all of the symbolic functions and constants
    (from the probability lists full_internal and full_nullary).  It is
    possible to adjust the ratio of leaves between symbolic constants,
    variables, and coefficients (var_frac gives the fraction of variables,
    and nullary_frac the fraction of symbolic constants; the remaining
    leaves are coefficients).

    The actual mix of symbolic constants and internal nodes can be modified
    by specifying different probability lists.

    To use a different type for coefficients, you can specify
    coeff_generator, which should be a function that will return
    a random coefficient every time it is called.

    This function will often raise an error because it tries to create
    an erroneous expression (such as a division by zero).

    EXAMPLES::

        sage: from sage.symbolic.random_tests import *
        sage: some_functions = [arcsinh, arctan, arctan2, arctanh,
        ....: arg, beta, binomial, ceil, conjugate, cos, cosh, cot, coth,
        ....: elliptic_pi, erf, exp, factorial, floor, heaviside, imag_part,
        ....: sech, sgn, sin, sinh, tan, tanh, unit_step, zeta, zetaderiv]
        sage: my_internal = [(0.6, full_binary, 2), (0.2, full_unary, 1),
        ....: (0.2, [(1.0,f,f.number_of_arguments()) for f in some_functions])]
        sage: set_random_seed(1)
        sage: random_expr(50, nvars=3, internal=my_internal,
        ....:   coeff_generator=CDF.random_element)  # not tested  # known bug
        (v1^(0.9713408427702117 + 0.195868299334218*I)/cot(-pi + v1^2 + v3) + tan(arctan(v2 + arctan2(-0.35859061674557324 + 0.9407509502498164*I, v3) - 0.8419115504372718 + 0.30375717982404615*I) + arctan2((0.2275357305882964 - 0.8258002386106038*I)/factorial(v2), -v3 - 0.7604559947718565 - 0.5543672548552057*I) + ceil(1/arctan2(v1, v1))))/v2
        sage: random_expr(5, verbose=True)  # not tested  # known bug
        About to apply <built-in function inv> to [31]
        About to apply sgn to [v1]
        About to apply <built-in function add> to [1/31, sgn(v1)]
        sgn(v1) + 1/31

    """
    vars = [(1.0, SR.var('v%d' % (n+1))) for n in range(nvars)]
    if ncoeffs is None:
        ncoeffs = size
    coeffs = [(1.0, coeff_generator()) for _ in range(ncoeffs)]
    leaves = [(var_frac, vars), (1.0 - var_frac - nullary_frac, coeffs), (nullary_frac, nullary)]
    leaves = normalize_prob_list(leaves)

    internal = normalize_prob_list(internal)

    return random_expr_helper(size, internal, leaves, verbose)


###################################################################
### Test the ordering of operands #################################
###################################################################

def assert_strict_weak_order(a, b, c, cmp_func):
    r"""
    Check that ``cmp_func`` is a strict weak order on the elements a,b,c.

    A strict weak order is a binary relation ``<`` such that

    * For all `x`, it is not the case that `x < x` (irreflexivity).

    * For all `x\not=y`, if `x < y` then it is not the case that `y <
      x` (asymmetry).

    * For all `x`, `y`, and `z`, if `x < y` and `y < z` then `x < z`
      (transitivity).

    * For all `x`, `y`, and `z`, if x is incomparable with `y`, and
      `y` is incomparable with `z`, then `x` is incomparable with `z`
      (transitivity of incomparability).

    INPUT:

    - ``a``, ``b``, ``c`` -- anything that can be compared by ``cmp_func``.

    - ``cmp_func`` -- function of two arguments that returns their
      comparison (i.e. either ``True`` or ``False``).

    OUTPUT:

    Does not return anything. Raises a ``ValueError`` if ``cmp_func``
    is not a strict weak order on the three given elements.

    REFERENCES:

    :wikipedia:`Strict_weak_ordering`

    EXAMPLES:

    The usual ordering of integers is a strict weak order::

        sage: from sage.symbolic.random_tests import assert_strict_weak_order
        sage: a, b, c = [randint(-10, 10) for i in range(3)]
        sage: assert_strict_weak_order(a, b, c, lambda x, y: x < y)

        sage: x = [-SR(oo), SR(0), SR(oo)]
        sage: cmp_M = matrix(3, 3, 0)
        sage: for i in range(3):
        ....:     for j in range(3):
        ....:         if x[i] < x[j]:
        ....:             cmp_M[i, j] = -1
        ....:         elif x[i] > x[j]:
        ....:             cmp_M[i, j] = 1
        sage: cmp_M
        [ 0 -1 -1]
        [ 1  0 -1]
        [ 1  1  0]
    """
    from sage.matrix.constructor import matrix
    from sage.combinat.permutation import Permutations
    x = (a, b, c)

    cmp_M = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            cmp_M[i, j] = (cmp_func(x[i], x[j]) == 1)   # or -1, doesn't matter

    msg = 'the binary relation failed to be a strict weak order on the elements \n'
    msg += ' a = {}\n b = {}\n c = {}\n'.format(a, b, c)
    msg += str(cmp_M)

    for i in range(3):
        # irreflexivity
        if cmp_M[i, i]:
            raise ValueError(msg)

        # asymmetric
        for j in range(i):
            if cmp_M[i, j] and cmp_M[j, i]:
                raise ValueError(msg)

    def incomparable(i, j):
        return not (cmp_M[i, j] or cmp_M[j, i])

    for i, j, k in Permutations([0, 1, 2]):
        # transitivity
        if cmp_M[i, j] and cmp_M[j, k] and not cmp_M[i, k]:
            raise ValueError(msg)

        # transitivity of incomparability
        if (incomparable(i, j) and incomparable(j, k) and
                not incomparable(i, k)):
            raise ValueError(msg)


def test_symbolic_expression_order(repetitions=100):
    r"""
    Tests whether the comparison of random symbolic expressions
    satisfies the strict weak order axioms.

    This is important because the C++ extension class uses
    ``std::sort()`` which requires a strict weak order. See also
    :trac:`9880`.

    EXAMPLES::

        sage: from sage.symbolic.random_tests import test_symbolic_expression_order
        sage: test_symbolic_expression_order(200)
        sage: test_symbolic_expression_order(10000)  # long time
    """
    rnd_length = 50
    nvars = 10
    ncoeffs = 10
    var_frac = 0.5
    nullary_frac = 0.05

    def coeff_generator():
        return randint(-100,100)/randint(1,100)

    def make_random_expr():
        while True:
            try:
                return random_expr(
                    rnd_length, nvars=nvars, ncoeffs=ncoeffs, var_frac=var_frac,
                    nullary_frac=nullary_frac, coeff_generator=coeff_generator,
                    internal=fast_nodes)
            except (ZeroDivisionError, ValueError):
                pass

    for rep in range(repetitions):
        a = make_random_expr()
        b = make_random_expr()
        c = make_random_expr()
        assert_strict_weak_order(a, b, c, mixed_order)
        assert_strict_weak_order(a, b, c, lambda x,y: x._cmp_add(y))
        assert_strict_weak_order(a, b, c, lambda x,y: x._cmp_mul(y))
