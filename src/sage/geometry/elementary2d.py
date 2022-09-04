import math

from sage.misc.lazy_attribute import lazy_attribute

from sage.structure.sage_object import SageObject

from sage.plot.graphics import Graphics

from sage.modules.free_module_element import vector


class GeometricFigure(SageObject):
    pass


class Plane2D(GeometricFigure):
    def __init__(self):
        self._figures = []

    def _repr_(self):
        num = len(self._figures)
        return f"Euclidean Plane containing {num} figures"

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


default_plane = Plane2D()


class Figure2D(GeometricFigure):
    def __init__(self, plane=None):
        if plane is None:
            plane = default_plane
        self._plane = plane
        plane.add(self)

    def plane(self):
        return self._plane

    @staticmethod
    def _turn(vec, angle):
        theta = math.atan(vec[1] / vec[0]) + angle
        return vector([math.cos(theta), math.sin(theta)])

    def point(self, t=0):
        return Point(self._pvec(t))

    def near_point(self, t=0, angle=None, distance=1):
        if angle is None:
            angle = math.pi / 2
        return Point(self._pvec(t) + distance * self._turn(self._tvec(t), angle))


class Point2D(Figure2D):
    def __init__(self, pos, plane=None):
        super().__init__(plane=plane)
        self.pos = tuple(pos)

    def __iter__(self):
        return iter(self.pos)

    def _repr_(self):
        return f"Point at {self.pos}"

    def show(self, start=0, end=1, *args, **kwargs):
        from sage.plot.point import point2d
        return point2d(self.pos, *args, **kwargs)


class Line2D(Figure2D):
    """
    EXAMPLES::

        sage: from sage.geometry.elementary2d import Line
        sage: l = Line((1,2), (4,7))
        sage: p1 = l.near_point(t=0, distance=2)
        sage: p2 = l.near_point(t=1, distance=2)
        sage: l2 = Line(p1, p2)
        sage: l3 = Line(l.start_point(), p1)
        sage: l4 = Line(l.end_point(), p2)
        sage: l.plane()
        Euclidean Plane containing 8 figures
        sage: l.plane().show()
    """
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
        from sage.plot.line import line2d
        return line2d([self._start, self._end], *args, **kwargs)


class Circle2D(Figure2D):
    def __init__(self, center, radius, plane=None):
        super().__init__(plane=plane)
        self.center = center
        self.length = length


class Ellipse2D(Figure2D):
    def __init__(self, focus1, focus2, semimajor, plane=None):
        super().__init__(plane=plane)
        self.focus1 = focus1
        self.focus2 = focus2
        self.semimajor = semimajor


class Parabola2D(Figure2D):
    def __init__(self, focus, directrix: Line2D, plane=None):
        super().__init__(plane=plane)
        self.focus: focus
        slef.directrix: directrix


def parabola_through(pt1, pt2, pt3):
    pass

def conic_by_equation(equ):
    poly = equ.left_hand_side() - equ.right_hand_side()
    a = poly.coefficient(x,2).coefficient(y,0)
    b = poly.coefficient(x,1).coefficient(y,1)
    c = poly.coefficient(x,0).coefficient(y,2)
    d = poly.coefficient(x,1).coefficient(y,0)
    e = poly.coefficient(x,0).coefficient(y,1)
    f = poly.coefficient(x,0).coefficient(y,0)

    A = matrix(2, [a, b / 2, b / 2, c])
