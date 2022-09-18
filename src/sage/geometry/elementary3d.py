import math

from sage.misc.lazy_attribute import lazy_attribute

from sage.structure.sage_object import SageObject

from sage.plot.graphics import Graphics

from sage.modules.free_module_element import vector

from .elementary2d import GeometricFigure


class Space(GeometricFigure):
    def __init__(self):
        self._figures = []

    def _repr_(self):
        num = len(self._figures)
        return f"Euclidean Space containing {num} figures"

    def add(self, figure):
        self._figures.append(figure)

    def clear(self):
        self._figures = []

    def show(self, *args, **kwargs):
        g = Graphics()
        for figure in self._figures:
            g += figure.show()
        if 'aspect_ratio' not in kwargs:
            kwargs['aspect_ratio'] = True
        g.show(*args, **kwargs)


default_space = Space()


class Figure3D(GeometricFigure):
    def __init__(self, space=None):
        if space is None:
            space = default_space
        self._space = space
        space.add(self)

    def space(self):
        return self._space

    def point(self, t=0, s=0):
        return Point(self._pvec(t, s))

    def near_point(self, t=0, s=0, angle=None, distance=1):
        if angle is None:
            angle = math.pi / 2
        return Point(self._pvec(t) + distance * self._turn(self._tvec(t), angle))


class Point(Figure3D):
    def __init__(self, pos):
        self._pos = pos

    def __iter__(self):
        return iter(self._pos)

    def _repr_(self):
        return f"Point at {self._pos}"

    def show(self, start=0, end=1, *args, **kwargs):
        from sage.plot3d.point import point3d
        return point3d(self.pos, *args, **kwargs)


class Line(Figure3D):
    def __init__(self, start, end, space=None):
        super().__init__(space=space)
        self.pos_start = start
        self.pos_end = end

    def __init__(self, start, end, plane=None):
        super().__init__(plane=plane)
        self._start = tuple(start)
        self._end = tuple(end)

    def _repr_(self):
        return f"Line from {self._start} to {self._end}"

    @lazy_attribute
    def _svec(self):
        return vector(self._start)

    @lazy_attribute
    def _evec(self):
        return vector(self._end)

    def _pvec(self, t=0):
        v = self._evec - self._svec
        return self._svec + t*v

    def _tvec(self, t=0):
        v = self._evec - self._svec
        v = ~(v.norm()) * v
        if t >= 0:
            return v
        return -v

    def start_point(self):
        return Point(self._start)

    def end_point(self):
        return Point(self._end)

    def show(self, start=0, end=1, *args, **kwargs):
        from sage.plot3d.line import line3d
        return line3d([self._start, self._end], *args, **kwargs)



class Plane(Figure3D):
    """
    EXAMPLES::

        sage: from sage.geometry.elementary3d import Plane
        sage: pl = Plane((2,0,0),(0,2,0),(0,0,2))
    """
    def __init__(self, pt1, pt2, pt3):
        """
        """
        self._pos1 = tuple(pt1)
        self._pos2 = tuple(pt2)
        self._pos3 = tuple(pt3)

    @lazy_attribute
    def _ovec(self):
        return vector(self._pos1)

    @lazy_attribute
    def _xvec(self):
        v = vector(self._pos2) - self._ovec
        v = ~(v.norm()) * v
        return v

    @lazy_attribute
    def _yvec(self):
        v1 = self._xvec
        v2 = vector(self._pos3) - self._ovec
        yvec = -(v1.inner_product(v2))*v1 + v2
        if not yvec:
            raise(f'points {self._pos1}, {self._pos2}, and {self._pos3} are colinear')
        return yvec

    def add(self, fig):
        from sage.geometry.elementary2d import Point, Line
        if isinstance(fig, Point):
            return PlanarPoint(fig, self)
        if isinstance(fig, Line):
            return PlanarLine(fig, self)


class PlanarFigure(Figure3D):
    def __init__(self, plane, space):
        super().__init__(space=space)
        self._plane = plane

    def _pvec(self, vec2d):
        v = vec2d[0] * self._plane._xvec + vec2d[1] * self._plane._yvec
        return self._plane._ovec + v


class PlanarPoint(PlanarFigure):
    """
    EXAMPLES::

        sage: from sage.geometry import elementary as el
        sage: pl = el.plane3d((1,0,0),(0,1,0),(0,0,1))
        sage: pt1 = el.point2d((0,0))
        sage: pt2 = el.point2d((1,0))
        sage: pt3 = el.point2d((0,1))
        sage: pt4 = el.point2d((1,1))
        sage: p1 = pl.add(pt1)
        sage: p2 = pl.add(pt2)
        sage: p3 = pl.add(pt3)
        sage: p4 = pl.add(pt4)
        sage: p1.show() + p2.show() + p3.show() + p4.show()
    """
    def __init__(self, point2d, plane, space=None):
        super().__init__(plane, space)
        self._point2d = point2d

    def show(self, start=0, end=1, *args, **kwargs):
        from sage.plot.plot3d.shapes2 import point3d
        v = self._pvec(list(self._point2d))
        return point3d([v], *args, **kwargs)


class PlanarLine(PlanarFigure):
    """
    EXAMPLES::

        sage: from sage.geometry import elementary as el
        sage: pl = el.plane3d((1,0,0),(0,1,0),(0,0,1))
        sage: pt1 = el.point2d((0,0))
        sage: pt2 = el.point2d((1,0))
        sage: pt3 = el.point2d((0,1))
        sage: pt4 = el.point2d((1,1))
        sage: l1 = el.line2d(pt1, pt4)
        sage: l2 = el.line2d(pt3, pt2)
        sage: line1 = pl.add(l1)
        sage: line2 = pl.add(l2)
        sage: line1.show() + line2.show()
    """
    def __init__(self, line2d, plane, space=None):
        super().__init__(plane, space)
        self._line2d = line2d

    def show(self, start=0, end=1, *args, **kwargs):
        from sage.plot.plot3d.shapes2 import line3d
        p1 = self._pvec(list(self._line2d.point(0)))
        p2 = self._pvec(list(self._line2d.point(1)))
        return line3d([list(p1), list(p2)], *args, **kwargs)


class PlanarCircle(PlanarFigure):
    def __init__(self, circle2d, plane, space=None):
        super().__init__(plane, space)
        self._circle2d = circle2d


class Sphere(Figure3D):
    def __init__(self, pt1, pt2, pt3):
        """
        """
        self._pos1 = tuple(pt1)
        self._pos2 = tuple(pt2)
        self._pos3 = tuple(pt3)

    @lazy_attribute
    def _ovec(self):
        return vector(self._pos1)

    @lazy_attribute
    def _xvec(self):
        v = vector(self._pos2) - self._ovec
        v = ~(v.norm()) * v
        return v

    @lazy_attribute
    def _yvec(self):
        v1 = self._xvec
        v2 = vector(self._pos3) - self._ovec
        yvec = -(v1.inner_product(v2))*v1 + v2
        if not yvec:
            raise(f'points {self._pos1}, {self._pos2}, and {self._pos3} are colinear')
        return yvec

    @lazy_attribute
    def _zvec(self):
        return self._xvec.cross_product(self._yvec)

    def radius(self):
        return (vector(self._pos2) - self._ovec).norm()

    def _pvec(self, t=0, s=0):
        z = math.sin(s * math.pi)
        c = math.cos(s * math.pi)
        x = c * math.cos(t * math.pi)
        y = c * math.sin(t * math.pi)
        return self._ovec + x * self._xvec + y * self._yvec + z * self._zvec

    def _tyvec(self, t=0, s=0):
        dx = -math.sin(s * math.pi) * math.cos(t * math.pi)
        dy = -math.sin(s * math.pi) * math.sin(t * math.pi)
        dx = math.cos(s * math.pi)
        retun vector([dx, dy, dz])

    def _txvec(self, t=0, s=0):
        return self._txvec(t, s).cross_product(self._pvec(t, s) - self._ovec)

    def tangent_plane(self, t=0, s=0)
        pt1 = list(self._pvec)
        pt2 = list(self._pvec(t, s) + self._txvec(t, s))
        pt3 = list(self._pvec(t, s) + self._tyvec(t, s))
        return Plane(pt1, pt2, pt3)


