"""
Enumeration of totally real fields: PHC interface

AUTHORS:

- John Voight (2007-09-19): initial version
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein and John Voight
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import sage.misc.misc


def coefficients_to_power_sums(n, m, a):
    r"""
    Take the list ``a``, representing a list of initial coefficients of
    a (monic) polynomial of degree `n`, and return the power sums
    of the roots of `f` up to `(m-1)`-th powers.

    INPUT:

    - ``n`` -- integer; the degree
    - ``a`` -- list of integers; the coefficients

    OUTPUT: list of integers

    .. NOTE::

        This uses Newton's relations, which are classical.

    EXAMPLES::

        sage: from sage.rings.number_field.totallyreal_phc import coefficients_to_power_sums
        sage: coefficients_to_power_sums(3,2,[1,5,7])
        [3, -7, 39]
        sage: coefficients_to_power_sums(5,4,[1,5,7,9,8])
        [5, -8, 46, -317, 2158]
    """
    S = [n] + [0]*m
    for k in range(1,m+1):
        S[k] = -sum([a[n-i]*S[k-i] for i in range(1,k)])-k*a[n-k]
    return S


def __lagrange_bounds_phc(n, m, a, tmpfile=None):
    r"""
    This function determines the bounds on the roots in
    the enumeration of totally real fields via Lagrange multipliers.

    It is used internally by the main function
    enumerate_totallyreal_fields_prim(), which should be consulted for
    further information.

    INPUT:

    - ``k`` -- integer; the index of the next coefficient
    - ``a`` -- list of integers; the coefficients

    OUTPUT: the lower and upper bounds as real numbers

    .. NOTE::

        See Cohen [Coh2000]_ for the general idea and unpublished work of the
        author for more detail.

    EXAMPLES::

        sage: # optional - phc
        sage: from sage.rings.number_field.totallyreal_phc import __lagrange_bounds_phc
        sage: __lagrange_bounds_phc(3,5,[8,1,2,0,1])
        []
        sage: x, y = __lagrange_bounds_phc(3,2,[8,1,2,0,1])
        sage: x
        -1.3333333333333299
        sage: y < 0.00000001
        True
        sage: __lagrange_bounds_phc(3,1,[8,1,2,0,1])
        []
    """

    # Compute power sums.
    S = coefficients_to_power_sums(n,m,a)

    # Look for phc.
    fi, fo = os.popen2('which phc')
    find_phc = fo.readlines()
    fi.close()
    fo.close()
    if find_phc == []:
        raise RuntimeError("PHCpack not installed.")

    # Initialization.
    if tmpfile is None:
        tmpfile = sage.misc.misc.tmp_filename()
    f = open(tmpfile + '.phc', 'w')
    f.close()

    output_data = []

    # By the method of Lagrange multipliers, if we maximize x_n subject to
    #     S_j(x) = S[j] (j = 1, ..., m),
    # then there are at most m-1 distinct values amongst the x_i.
    # Therefore we must solve the implied equations for each partition of n-1
    # into m-1 parts.
    for P in sage.combinat.partition.Partitions(n-1,length=m-1):
        f = open(tmpfile, 'w')
        # First line: number of variables/equations
        f.write('%d' % m + '\n')
        # In the next m-1 lines, write the equation S_j(x) = S[j]
        for j in range(1,m+1):
            for i in range(m-1):
                f.write('%d' % P[i] + '*x%d' % i + '**%d' % j + ' + ')
            f.write('xn**%d' % j + ' - (%d' % S[j] + ');\n')
        f.close()

        os.remove(tmpfile + '.phc')
        os.popen('phc -b ' + tmpfile + ' ' + tmpfile + '.phc')
        f = open(tmpfile + '.phc')
        f_str = f.read()
        pos = f_str.find('= real ')
        crits = []
        while pos != -1:
            posl = f_str.rfind('xn', 0, pos)
            f_str_split = f_str[posl:pos].split()
            crits += [float(f_str_split[2])]
            pos = f_str.find('= real ', pos+1)

        if len(crits) > 0:
            output_data += [[P, min(crits), max(crits)]]

    if len(output_data) > 0:
        return [min([v[1] for v in output_data]), max([v[2] for v in output_data])]
    else:
        return []
