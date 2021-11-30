r"""
Base class for mutable polyhedra.

Just like vectors and matrices they can be set immutable.
The constructor does this by default.
"""

from .base import Polyhedron_base


class Polyhedron_mutable(Polyhedron_base):
    """
    Base class for polyhedra that allow mutability.

    This should not be used directly.
    """

    def __hash__(self):
        r"""
        TESTS::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: set([p])
            Traceback (most recent call last):
            ...
            TypeError: mutable polyhedra are unhashable
            sage: p.set_immutable()
            sage: set([p])
            {A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex}
        """
        if self._is_mutable:
            raise TypeError("mutable polyhedra are unhashable")
        return Polyhedron_base.__hash__(self)

    def _clear_cache(self):
        r"""
        Clear the Vrepresentation and Hrepresentation data of ``self``.

        TESTS::

            sage: p = polytopes.permutahedron(4)
            sage: P = p.parent()
            sage: q = P._element_constructor_(p, mutable=True)
            sage: TestSuite(q).run()
            sage: q._clear_cache()
            sage: TestSuite(q).run()

        ::

            sage: q.set_immutable()
            sage: q._clear_cache()
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra
        """
        if not self._is_mutable:
            raise TypeError("cannot clear cache of immutable polyhedra")

        # Invalidate object pointing towards this polyhedron (faces etc.).
        for ob in self._dependent_objects:
            ob._polyhedron = None
        backend_object = self.__dict__["_" + self._backend_object_name]
        del self.__dict__
        self.__dict__["_" + self._backend_object_name] = backend_object
        self._is_mutable = True
        self._dependent_objects = []

    def _add_dependent_object(self, ob):
        r"""
        Add an object that has ``self`` has attribute ``_polyhedron``.

        When ``self`` is modified, we delete this attribute to invalidate those objects.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: class foo:
            ....:     def __init__(self, p):
            ....:         self._polyhedron = p
            ....:
            sage: a = foo(p)
            sage: a.__dict__
            {'_polyhedron': A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex}
            sage: p._add_dependent_object(a)
            sage: p._clear_cache()
            sage: a.__dict__
            {'_polyhedron': None}

        TESTS::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: n = NewtonPolygon(p)
            sage: n
            Finite Newton polygon with 1 vertex: (1, 1)
            sage: n = NewtonPolygon(p)
            sage: p._clear_cache()
            sage: n
            <repr(<sage.geometry.newton_polygon.ParentNewtonPolygon_with_category.element_class at ...>) failed: AttributeError: 'NoneType' object has no attribute 'vertices'>

        ::

            sage: f = p.faces(0)[0]; f
            A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: p._clear_cache()
            sage: f
            <repr(<sage.geometry.polyhedron.face.PolyhedronFace at ...>) failed: AttributeError: 'NoneType' object has no attribute 'parent'>

        ::

            sage: v = p.vertices()[0]
            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: v = p.Vrepresentation(0); v
            A vertex at (1, 1)
            sage: h = p.Hrepresentation(0); h
            An equation (0, 1) x - 1 == 0
            sage: p._clear_cache()
            sage: v.polyhedron() is None
            True
            sage: h.polyhedron() is None
            True

        ::

            sage: p = Polyhedron([[1, 0], [0, 1]], mutable=True)
            sage: r = p.relative_interior()
            sage: p._clear_cache()
            sage: r
            Relative interior of None
        """
        if ob._polyhedron is not self:
            raise ValueError
        self._dependent_objects.append(ob)

    def is_mutable(self):
        r"""
        Return True if the polyhedron is mutable, i.e. it can be modified in place.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_mutable()
            True
            sage: p = Polyhedron([[1, 1]], mutable=False)
            sage: p.is_mutable()
            False
        """
        return self._is_mutable

    def is_immutable(self):
        r"""
        Return True if the polyhedron is immutable, i.e. it cannot be modified in place.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_immutable()
            False
            sage: p = Polyhedron([[1, 1]], mutable=False)
            sage: p.is_immutable()
            True
        """
        return not self._is_mutable

    def set_immutable(self):
        r"""
        Make this polyhedron immutable. This operation cannot be undone.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.set_immutable(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        raise NotImplementedError("a derived class must implement this")

    def Vrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Vrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.Vrepresentation(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        # A derived class must implemented it to recalculate, if necessary.
        raise NotImplementedError("a derived class must implement this")

    def Hrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Hrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.Hrepresentation(p)
            Traceback (most recent call last):
            ...
            NotImplementedError: a derived class must implement this
        """
        # A derived class must implemented it to recalculate, if necessary.
        raise NotImplementedError("a derived class must implement this")
