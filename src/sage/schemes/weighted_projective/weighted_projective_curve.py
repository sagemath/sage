from sage.schemes.curves.curve import Curve_generic


# TODO: Implement (some) embedding into straight projective space
class WeightedProjectiveCurve(Curve_generic):
    def __init__(self, A, X, *kwargs):
        # TODO ensure that A is the right type?
        # Something like a `is_WeightProjectiveSpace` which means making a
        # WeightProjectiveSpace class?
        super().__init__(A, X, *kwargs)

    def curve(self):
        return
