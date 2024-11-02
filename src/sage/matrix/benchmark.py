# sage.doctest: needs sage.flint
"""
Benchmarks for matrices

This file has many functions for computing timing benchmarks
of various methods for random matrices with given bounds for
the entries.  The systems supported are Sage and Magma.

The basic command syntax is as follows::

    sage: import sage.matrix.benchmark as b
    sage: print("starting"); import sys; sys.stdout.flush(); b.report([b.det_ZZ], 'Test', systems=['sage'])
    starting...
    ======================================================================
              Test
    ======================================================================
    ...
    ======================================================================
"""

from .constructor import random_matrix, Matrix
from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.misc.timing import cputime
from cysignals.alarm import AlarmInterrupt, alarm, cancel_alarm

lazy_import('sage.interfaces.magma', 'magma')

verbose = False

timeout = 60


def report(F, title, systems=['sage', 'magma'], **kwds):
    """
    Run benchmarks with default arguments for each function in the list F.

    INPUT:

    - ``F`` -- list of callables used for benchmarking
    - ``title`` -- string describing this report
    - ``systems`` -- list of systems (supported entries are 'sage' and 'magma')
    - ``**kwds`` -- keyword arguments passed to all functions in ``F``

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: print("starting"); import sys; sys.stdout.flush(); b.report([b.det_ZZ], 'Test', systems=['sage'])
        starting...
        ======================================================================
                  Test
        ======================================================================
        ...
        ======================================================================
    """
    import os
    if len(systems) > 2:
        raise NotImplementedError("at most two systems ('sage' or 'magma')")
    print('=' * 70)
    print(' ' * 10 + title)
    print('=' * 70)
    os.system('uname -a')
    print('\n')
    for f in F:
        print("-"*70)
        print(f.__doc__.strip())
        print(('%15s' * len(systems)) % tuple(systems))
        w = []
        for s in systems:
            alarm(timeout)
            try:
                t = f(system=s, **kwds)
            except AlarmInterrupt:
                t = -timeout
            cancel_alarm()
            w.append(float(t))
        if len(w) > 1:
            if w[1] == 0:
                w.append(0.0)
            else:
                w.append(w[0]/w[1])

        w = tuple(w)
        print(('%15.3f'*len(w)) % w)
    print('=' * 70)


#######################################################################
# Dense Benchmarks over ZZ
#######################################################################

def report_ZZ(**kwds):
    """
    Reports all the benchmarks for integer matrices and few
    rational matrices.

    INPUT:

    - ``**kwds`` -- passed through to :func:`report`

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: print("starting"); import sys; sys.stdout.flush(); b.report_ZZ(systems=['sage'])  # long time (15s on sage.math, 2012)
        starting...
        ======================================================================
        Dense benchmarks over ZZ
        ======================================================================
        ...
        ======================================================================
    """
    F = [vecmat_ZZ, rank_ZZ, rank2_ZZ, charpoly_ZZ, smithform_ZZ,
         det_ZZ, det_QQ, matrix_multiply_ZZ, matrix_add_ZZ,
         matrix_add_ZZ_2,
         nullspace_ZZ]

    title = 'Dense benchmarks over ZZ'
    report(F, title, **kwds)

# Integer Nullspace

def nullspace_ZZ(n=200, min=0, max=2**32, system='sage'):
    """
    Nullspace over ZZ:
    Given a n+1 x n matrix over ZZ with random entries
    between min and max, compute the nullspace.

    INPUT:

    - ``n`` -- matrix dimension (default: ``200``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``2**32``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_ZZ(200)
        sage: tm = b.nullspace_ZZ(200, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n+1, n, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


def charpoly_ZZ(n=100, min=0, max=9, system='sage'):
    """
    Characteristic polynomial over ZZ:
    Given a n x n matrix over ZZ with random entries between min and
    max, compute the charpoly.

    INPUT:

    - ``n`` -- matrix dimension (default: ``100``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.charpoly_ZZ(100)
        sage: tm = b.charpoly_ZZ(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        v = A.charpoly()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := CharacteristicPolynomial(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


def rank_ZZ(n=700, min=0, max=9, system='sage'):
    """
    Rank over ZZ:
    Given a n x (n+10) matrix over ZZ with random entries
    between min and max, compute the rank.

    INPUT:

    - ``n`` -- matrix dimension (default: ``700``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank_ZZ(300)
        sage: tm = b.rank_ZZ(300, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n+10, x=min, y=max+1)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(IntegerRing(), n, n+10)![Random(%s,%s) : i in [1..n*(n+10)]];
t := Cputime();
K := Rank(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)

def rank2_ZZ(n=400, min=0, max=2**64, system='sage'):
    """
    Rank 2 over ZZ:
    Given a (n + 10) x n matrix over ZZ with random entries
    between min and max, compute the rank.

    INPUT:

    - ``n`` -- matrix dimension (default: ``400``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``2**64``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank2_ZZ(300)
        sage: tm = b.rank2_ZZ(300, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n+10, n, x=min, y=max+1)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(IntegerRing(), n+10, n)![Random(%s,%s) : i in [1..n*(n+10)]];
t := Cputime();
K := Rank(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)

# Smith Form

def smithform_ZZ(n=128, min=0, max=9, system='sage'):
    """
    Smith Form over ZZ:
    Given a n x n matrix over ZZ with random entries
    between min and max, compute the Smith normal form.

    INPUT:

    - ``n`` -- matrix dimension (default: ``128``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.smithform_ZZ(100)
        sage: tm = b.smithform_ZZ(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        v = A.elementary_divisors()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := ElementaryDivisors(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


def matrix_multiply_ZZ(n=300, min=-9, max=9, system='sage', times=1):
    """
    Matrix multiplication over ZZ
    Given an n x n matrix A over ZZ with random entries
    between min and max, inclusive, compute A * (A+1).

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``min`` -- minimal value for entries of matrix (default: ``-9``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``1``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_multiply_ZZ(200)
        sage: tm = b.matrix_multiply_ZZ(200, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        B = A + 1
        t = cputime()
        for z in range(times):
            v = A * B
        return cputime(t)/times
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
""" % (n,min,max,times)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError('unknown system "%s"' % system)

def matrix_add_ZZ(n=200, min=-9, max=9, system='sage', times=50):
    """
    Matrix addition over ZZ
    Given an n x n matrix A and B over ZZ with random entries between
    ``min`` and ``max``, inclusive, compute A + B ``times`` times.

    INPUT:

    - ``n`` -- matrix dimension (default: ``200``)
    - ``min`` -- minimal value for entries of matrix (default: ``-9``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``50``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_add_ZZ(200)
        sage: tm = b.matrix_add_ZZ(200, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        B = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        for z in range(times):
            v = A + B
        return cputime(t)/times
    elif system == 'magma':
        code = """
n := %s;
min := %s;
max := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(min,max) : i in [1..n^2]];
B := MatrixAlgebra(IntegerRing(), n)![Random(min,max) : i in [1..n^2]];
t := Cputime();
for z in [1..%s] do
    K := A + B;
end for;
s := Cputime(t);
""" % (n,min,max,times)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError('unknown system "%s"' % system)

def matrix_add_ZZ_2(n=200, bits=16, system='sage', times=50):
    """
    Matrix addition over ZZ.
    Given an n x n matrix A and B over ZZ with random ``bits``-bit
    entries, compute A + B.

    INPUT:

    - ``n`` -- matrix dimension (default: ``200``)
    - ``bits`` -- bitsize of entries
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``50``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_add_ZZ_2(200)
        sage: tm = b.matrix_add_ZZ_2(200, system='magma')  # optional - magma
    """
    b = 2**bits
    return matrix_add_ZZ(n=n, min=-b, max=b,system=system, times=times)

def det_ZZ(n=200, min=1, max=100, system='sage'):
    """
    Dense integer determinant over ZZ.
    Given an n x n matrix A over ZZ with random entries
    between min and max, inclusive, compute det(A).

    INPUT:

    - ``n`` -- matrix dimension (default: ``200``)
    - ``min`` -- minimal value for entries of matrix (default: ``1``)
    - ``max`` -- maximal value for entries of matrix (default: ``100``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_ZZ(200)
        sage: tm = b.det_ZZ(200, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
d := Determinant(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


def det_QQ(n=300, num_bound=10, den_bound=10, system='sage'):
    """
    Dense rational determinant over QQ.
    Given an n x n matrix A over QQ with random entries
    with numerator bound and denominator bound, compute det(A).

    INPUT:

    - ``n`` -- matrix dimension (default: ``200``)
    - ``num_bound`` -- numerator bound, inclusive (default: ``10``)
    - ``den_bound`` -- denominator bound, inclusive (default: ``10``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_QQ(200)
        sage: ts = b.det_QQ(10, num_bound=100000, den_bound=10000)
        sage: tm = b.det_QQ(200, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(QQ, n, n, num_bound=num_bound, den_bound=den_bound)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(RationalField(), n)![Random(%s,%s)/Random(1,%s) : i in [1..n^2]];
t := Cputime();
d := Determinant(A);
s := Cputime(t);
""" % (n,-num_bound, num_bound, den_bound)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


def vecmat_ZZ(n=300, min=-9, max=9, system='sage', times=200):
    """
    Vector matrix multiplication over ZZ.

    Given an n x n  matrix A over ZZ with random entries
    between min and max, inclusive, and v the first row of A,
    compute the product v * A.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``min`` -- minimal value for entries of matrix (default: ``-9``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)
    - ``times`` -- number of runs (default: ``200``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.vecmat_ZZ(300)  # long time
        sage: tm = b.vecmat_ZZ(300, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        v = A.row(0)
        t = cputime()
        for z in range(times):
            w = v * A
        return cputime(t)/times
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
v := A[1];
t := Cputime();
for z in [1..%s] do
    K := v * A;
end for;
s := Cputime(t);
""" % (n,min,max,times)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError('unknown system "%s"' % system)


#######################################################################
# Dense Benchmarks over GF(p), for small p.
#######################################################################

def report_GF(p=16411, **kwds):
    """
    Run all the reports for finite field matrix operations, for
    prime p=16411.

    INPUT:

    - ``p`` -- ignored
    - ``**kwds`` -- passed through to :func:`report`

    .. NOTE::

        right now, even though p is an input, it is being ignored!  If
        you need to check the performance for other primes, you can
        call individual benchmark functions.

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: print("starting"); import sys; sys.stdout.flush(); b.report_GF(systems=['sage'])
        starting...
        ======================================================================
        Dense benchmarks over GF with prime 16411
        ======================================================================
        ...
        ======================================================================
    """
    F = [rank_GF, rank2_GF, nullspace_GF, charpoly_GF,
         matrix_multiply_GF, det_GF]
    title = 'Dense benchmarks over GF with prime %i' % p
    report(F, title, **kwds)

# Nullspace over GF

def nullspace_GF(n=300, p=16411, system='sage'):
    """
    Given a n+1 x n  matrix over GF(p) with random
    entries, compute the nullspace.

    INPUT:

    - ``n`` -- matrix dimension (default: 300)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_GF(300)
        sage: tm = b.nullspace_GF(300, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n+1)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(RMatrixSpace(GF(%s), n, n+1));
t := Cputime();
K := Kernel(A);
s := Cputime(t);
""" % (n,p)
        if verbose:
            print(code)
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError('unknown system "%s"' % system)


# Characteristic Polynomial over GF

def charpoly_GF(n=100, p=16411, system='sage'):
    """
    Given a n x n matrix over GF with random entries, compute the
    charpoly.

    INPUT:

    - ``n`` -- matrix dimension (default: 100)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.charpoly_GF(100)
        sage: tm = b.charpoly_GF(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n)
        t = cputime()
        v = A.charpoly()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := CharacteristicPolynomial(A);
s := Cputime(t);
""" % (n,p)
        if verbose:
            print(code)
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError('unknown system "%s"' % system)

def matrix_add_GF(n=1000, p=16411, system='sage', times=100):
    """
    Given two n x n matrix over GF(p) with random entries, add them.

    INPUT:

    - ``n`` -- matrix dimension (default: 300)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``100``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_add_GF(500, p=19)
        sage: tm = b.matrix_add_GF(500, p=19, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n)
        B = random_matrix(GF(p), n, n)
        t = cputime()
        for n in range(times):
            v = A + B
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
B := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
for z in [1..%s] do
    K := A + B;
end for;
s := Cputime(t);
""" % (n,p,p,times)
        if verbose:
            print(code)
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError('unknown system "%s"' % system)


# Matrix multiplication over GF(p)

def matrix_multiply_GF(n=100, p=16411, system='sage', times=3):
    """
    Given an n x n matrix A over GF(p) with random entries, compute
    A * (A+1).

    INPUT:

    - ``n`` -- matrix dimension (default: 100)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``3``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_multiply_GF(100, p=19)
        sage: tm = b.matrix_multiply_GF(100, p=19, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n)
        B = A + 1
        t = cputime()
        for n in range(times):
            v = A * B
        return cputime(t) / times
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
""" % (n,p,times)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError('unknown system "%s"' % system)


def rank_GF(n=500, p=16411, system='sage'):
    """
    Rank over GF(p):
    Given a n x (n+10) matrix over GF(p) with random entries, compute the rank.

    INPUT:

    - ``n`` -- matrix dimension (default: 300)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank_GF(1000)
        sage: tm = b.rank_GF(1000, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n+10)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := Rank(A);
s := Cputime(t);
""" % (n,p)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)

def rank2_GF(n=500, p=16411, system='sage'):
    """
    Rank over GF(p): Given a (n + 10) x n matrix over GF(p) with
    random entries, compute the rank.

    INPUT:

    - ``n`` -- matrix dimension (default: 300)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.rank2_GF(500)
        sage: tm = b.rank2_GF(500, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n+10, n)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := Rank(A);
s := Cputime(t);
""" % (n,p)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)

def det_GF(n=400, p=16411 , system='sage'):
    """
    Dense determinant over GF(p).
    Given an n x n matrix A over GF with random entries compute
    det(A).

    INPUT:

    - ``n`` -- matrix dimension (default: 300)
    - ``p`` -- prime number (default: ``16411``)
    - ``system`` -- either 'magma' or 'sage' (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_GF(1000)
        sage: tm = b.det_GF(1000, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
d := Determinant(A);
s := Cputime(t);
""" % (n,p)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


#######################################################################
# Dense Benchmarks over QQ
#######################################################################

def hilbert_matrix(n):
    """
    Return the Hilbert matrix of size n over rationals.

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: b.hilbert_matrix(3)
        [  1 1/2 1/3]
        [1/2 1/3 1/4]
        [1/3 1/4 1/5]
    """
    A = Matrix(QQ,n,n)
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            A[i,j] = QQ(1)/((i+1)+(j+1)-1)
    return A

# Reduced row echelon form over QQ

def echelon_QQ(n=100, min=0, max=9, system='sage'):
    """
    Given a n x (2*n) matrix over QQ with random integer entries
    between min and max, compute the reduced row echelon form.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``min`` -- minimal value for entries of matrix (default: ``-9``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.echelon_QQ(100)
        sage: tm = b.echelon_QQ(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, 2*n, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = A.echelon_form()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n, 2*n)![Random(%s,%s) : i in [1..n*2*n]];
t := Cputime();
K := EchelonForm(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)

# Invert a matrix over QQ.

def inverse_QQ(n=100, min=0, max=9, system='sage'):
    """
    Given a n x n matrix over QQ with random integer entries
    between min and max, compute the reduced row echelon form.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``min`` -- minimal value for entries of matrix (default: ``-9``)
    - ``max`` -- maximal value for entries of matrix (default: ``9``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.inverse_QQ(100)
        sage: tm = b.inverse_QQ(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = ~A
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(RationalField(), n)![Random(%s,%s) : i in [1..n*n]];
t := Cputime();
K := A^(-1);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


# Matrix multiplication over QQ
def matrix_multiply_QQ(n=100, bnd=2, system='sage', times=1):
    """
    Given an n x n matrix A over QQ with random entries
    whose numerators and denominators are bounded by bnd,
    compute A * (A+1).

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``bnd`` -- numerator and denominator bound (default: ``bnd``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``1``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.matrix_multiply_QQ(100)
        sage: tm = b.matrix_multiply_QQ(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = random_matrix(QQ, n, n, num_bound=bnd, den_bound=bnd)
        B = A + 1
        t = cputime()
        for z in range(times):
            v = A * B
        return cputime(t)/times
    elif system == 'magma':
        A = magma(random_matrix(QQ, n, n, num_bound=bnd, den_bound=bnd))
        code = """
n := %s;
A := %s;
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
""" % (n, A.name(), times)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError('unknown system "%s"' % system)


# Determinant of Hilbert matrix
def det_hilbert_QQ(n=80, system='sage'):
    """
    Run the benchmark for calculating the determinant of the hilbert
    matrix over rationals of dimension n.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.det_hilbert_QQ(50)
        sage: tm = b.det_hilbert_QQ(50, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = hilbert_matrix(n)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
h := HilbertMatrix(%s);
tinit := Cputime();
d := Determinant(h);
s := Cputime(tinit);
delete h;
""" % n
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))

# inverse of Hilbert matrix
def invert_hilbert_QQ(n=40, system='sage'):
    """
    Run the benchmark for calculating the inverse of the hilbert
    matrix over rationals of dimension n.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.invert_hilbert_QQ(30)
        sage: tm = b.invert_hilbert_QQ(30, system='magma')  # optional - magma
    """
    if system == 'sage':
        A = hilbert_matrix(n)
        t = cputime()
        d = A**(-1)
        return cputime(t)
    elif system == 'magma':
        code = """
h := HilbertMatrix(%s);
tinit := Cputime();
d := h^(-1);
s := Cputime(tinit);
delete h;
""" % n
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))

def MatrixVector_QQ(n=1000, h=100, system='sage', times=1):
    """
    Compute product of square ``n`` matrix by random vector with num and
    denom bounded by ``h`` the given number of ``times``.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``h`` -- numerator and denominator bound (default: ``bnd``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)
    - ``times`` -- number of experiments (default: ``1``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.MatrixVector_QQ(500)
        sage: tm = b.MatrixVector_QQ(500, system='magma')  # optional - magma
    """
    if system == 'sage':
        V = QQ**n
        v = V.random_element(h)
        M = random_matrix(QQ,n)
        t = cputime()
        for i in range(times):
            w = M*v
        return cputime(t)
    elif system == 'magma':
        code = """
            n:=%s;
            h:=%s;
            times:=%s;
            v:=VectorSpace(RationalField(),n)![Random(h)/(Random(h)+1) : i in [1..n]];
            M:=MatrixAlgebra(RationalField(),n)![Random(h)/(Random(h)+1) : i in [1..n^2]];
            t := Cputime();
            for z in [1..times] do
                W:=v*M;
            end for;
            s := Cputime(t);
        """ % (n,h,times)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


#######################################################################
# Dense Benchmarks over machine reals
# Note that the precision in reals for MAGMA is base 10, while in
# sage it is in base 2
#######################################################################

# Real Nullspace

def nullspace_RR(n=300, min=0, max=10, system='sage'):
    """
    Nullspace over RR:
    Given a n+1 x n matrix over RR with random entries
    between min and max, compute the nullspace.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``10``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_RR(100)
        sage: tm = b.nullspace_RR(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        from sage.rings.real_mpfr import RR
        A = random_matrix(ZZ, n+1, n, x=min, y=max+1).change_ring(RR)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RealField(16), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
""" % (n,min,max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)


def nullspace_RDF(n=300, min=0, max=10, system='sage'):
    """
    Nullspace over RDF:
    Given a n+1 x n  matrix over RDF with random entries
    between min and max, compute the nullspace.

    INPUT:

    - ``n`` -- matrix dimension (default: ``300``)
    - ``min`` -- minimal value for entries of matrix (default: ``0``)
    - ``max`` -- maximal value for entries of matrix (default: ``10``)
    - ``system`` -- either ``'sage'`` or ``'magma'`` (default: ``'sage'``)

    EXAMPLES::

        sage: import sage.matrix.benchmark as b
        sage: ts = b.nullspace_RDF(100)  # long time
        sage: tm = b.nullspace_RDF(100, system='magma')  # optional - magma
    """
    if system == 'sage':
        from sage.rings.real_double import RDF
        A = random_matrix(ZZ, n+1, n, x=min, y=max+1).change_ring(RDF)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RealField(16), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
""" % (n, min, max)
        if verbose:
            print(code)
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError('unknown system "%s"' % system)
