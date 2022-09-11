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


default_space = Space3D()


class Figure3D(GeometricFigure):
    def __init__(self, space=None):
        if space is None:
            space = default_space
        self._space = space
        space.add(self)

    def space(self):
        return self._space


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
    def __init__(self, pt1, pt2, pt3):
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


class PlanarPoint(Figure3D):
    def __init__(self, point2d, plane, space=None):
        super().__init__(space=space)
        self._point2d = point2d
        self._plane = plane

class PlanarLine(Figure3D):
    def __init__(self, line2d, plane, space=None):
        super().__init__(space=space)
        self._line2d = line2d
        self._plane = plane

class PlanarCircle(Figure3D):
    def __init__(self, circle2d, plane, space=None):
        super().__init__(space=space)
        self._circle2d = circle2d
        self._plane = plane
