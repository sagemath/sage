class GeometricFigure:
    pass


class Figure2D(GeometricFigure):
    pass


class Point(Figure2D):
    def __init__(self, pos):
        self.pos = pos


class Line(Figure2D):
    def __init__(self, start, end):
        self.start_point = start
        self.end_point = end

    def intersect(fig):
        pass

    def render(self, *args, **kwargs):
        pass


class Circle(Figure2D):
    def __init__(self, center, radius):
        self.center = center
        self.length = length


class Ellipse(Figure2D):
    def __init__(self, focus1, focus2, semimajor):
        self.focus1 = focus1
        self.focus2 = focus2
        self.semimajor = semimajor


class Parabola(Figure2D):
    def __init__(self, focus, directrix: Line):
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
