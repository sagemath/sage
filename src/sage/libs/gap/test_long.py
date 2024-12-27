"""
Long tests for GAP

These stress test the garbage collection inside GAP
"""
from sage.libs.gap.libgap import libgap


def check_loop_1():
    """
    EXAMPLES::

        sage: from sage.libs.gap.test_long import check_loop_1
        sage: check_loop_1()  # long time (up to 25s on sage.math, 2013)
    """
    libgap.collect()
    for i in range(10000):
        _ = libgap.CyclicGroup(2)


def check_loop_2():
    """
    EXAMPLES::

        sage: from sage.libs.gap.test_long import check_loop_2
        sage: check_loop_2()  # long time (10s on sage.math, 2013)
    """
    G = libgap.FreeGroup(2)
    a, b = G.GeneratorsOfGroup()
    for i in range(100):
        rel = libgap([a**2, b**2, a * b * a * b])
        H = G / rel
        H1 = H.GeneratorsOfGroup()[0]
        n = H1.Order()
        assert n.sage() == 2

    for i in range(300000):
        n = libgap.Order(H1)


def check_loop_3():
    """
    EXAMPLES::

        sage: from sage.libs.gap.test_long import check_loop_3
        sage: check_loop_3()  # long time (31s on sage.math, 2013)
    """
    G = libgap.FreeGroup(2)
    a, b = G.GeneratorsOfGroup()
    for i in range(300000):
        lis = libgap([])
        lis.Add(a ** 2)
        lis.Add(b ** 2)
        lis.Add(b * a)
