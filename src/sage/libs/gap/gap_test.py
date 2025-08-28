import pytest
from sage.libs.gap.libgap import libgap


def test_libgap_can_read_and_write_files(tmpfile):
    """
    Test that libgap can write to a temporary file and
    subsequently read it.

    See :issue:`16502`, :issue:`15833`.
    """
    message = "Ceci n'est pas une groupe"
    libgap.PrintTo(tmpfile.name, message)
    with open(tmpfile.name) as f:
        contents = f.read()
    assert contents == message
    SystemFile = libgap.function_factory('StringFile')
    assert SystemFile(tmpfile.name) == libgap(message)


def test_gc_loop_1():
    r"""
    Stress test for garbage collection in libgap.

    Manually runs the GAP garbage collector, and then creates 10,000
    instances of the cyclic group of order two in a loop. In each
    iteration, the python variable is overwritten, meaning that python
    is free to garbage collect the object.
    """
    libgap.collect()
    for _ in range(10000):
        G = libgap.CyclicGroup(2)
    assert True


def test_gc_loop_2():
    r"""
    Stress test for garbage collection in libgap.

    Create the free group on two elements (``a`` and ``b``) and then
    construct a quotient group of order two in a loop by specifying
    some relations. The python variables are overwritten in each
    iteration, meaning that python is free to garbage-collect them.
    (We also ensure that the quotient group has the expected order.)

    After that loop, we take one of the generators of the quotient
    group (from the final iteration), and compute its order in a
    loop. The python reference from the final iteration lives on,
    so this generator should not be collected.
    """
    G = libgap.FreeGroup(2)
    a, b = G.GeneratorsOfGroup()
    two = libgap(2)

    for _ in range(100):
        rel = libgap([a**2, b**2, a*b*a*b])
        H = G / rel
        H1 = H.GeneratorsOfGroup()[0]
        n = H1.Order()
        assert n == two

    result = True
    for i in range(300000):
        n = libgap.Order(H1)
        result &= (n == two)
    assert result


def test_gc_loop_3():
    r"""
    Stress test for garbage collection in libgap.

    Create the free group on two elements (``a`` and ``b``) and then
    add some of its elements to a list in a loop. A new, empty list is
    created at each iteration, so this mainly serves to guarantee that
    the generators of the group are not garbage-collected.
    """
    G = libgap.FreeGroup(2)
    a, b = G.GeneratorsOfGroup()
    for _ in range(300000):
        lis = libgap([])
        lis.Add(a ** 2)
        lis.Add(b ** 2)
        lis.Add(b * a)
    assert True
