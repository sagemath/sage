r"""
Abstract base classes for classes in ``sage.geometry``
"""


class LatticePolytope:
    r"""
    Abstract base class for :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.geometry.abc
        sage: P = LatticePolytope([(1,2,3), (4,5,6)])
        sage: isinstance(P, sage.geometry.abc.LatticePolytope)
        True

    By design, there is a unique direct subclass::

        sage: sage.geometry.abc.LatticePolytope.__subclasses__()
        [<class 'sage.geometry.lattice_polytope.LatticePolytopeClass'>]

        sage: len(sage.geometry.abc.Polyhedron.__subclasses__()) <= 1
        True
    """

    pass


class ConvexRationalPolyhedralCone:
    r"""
    Abstract base class for :class:`~sage.geometry.cone.ConvexRationalPolyhedralCone`

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.geometry.abc
        sage: C = cones.nonnegative_orthant(2)
        sage: isinstance(C, sage.geometry.abc.ConvexRationalPolyhedralCone)
        True

    By design, there is a unique direct subclass::

        sage: sage.geometry.abc.ConvexRationalPolyhedralCone.__subclasses__()
        [<class 'sage.geometry.cone.ConvexRationalPolyhedralCone'>]

        sage: len(sage.geometry.abc.Polyhedron.__subclasses__()) <= 1
        True
    """

    pass


class Polyhedron:
    r"""
    Abstract base class for :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES::

        sage: import sage.geometry.abc
        sage: P = polytopes.cube()
        sage: isinstance(P, sage.geometry.abc.Polyhedron)
        True

    By design, there is a unique direct subclass::

        sage: sage.geometry.abc.Polyhedron.__subclasses__()
        [<class 'sage.geometry.polyhedron.base0.Polyhedron_base0'>]

        sage: len(sage.geometry.abc.Polyhedron.__subclasses__()) <= 1
        True
    """

    pass
