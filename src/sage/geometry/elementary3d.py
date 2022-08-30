class GeometricFigure:
    pass


class Figure3D(GeometricFigure):
    pass


class Point(Figure3D):
    def __init__(self, pos):
        self.pos = pos


class Line(Figure3D):
    def __init__(self, start, end):
        self.pos_start = start
        self.pos_end = end



