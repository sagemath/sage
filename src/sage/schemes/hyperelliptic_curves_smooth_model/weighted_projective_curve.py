from sage.schemes.curves.curve import Curve_generic
from sage.schemes.toric.toric_subscheme import AlgebraicScheme_subscheme_toric


class WeightedProjectiveCurve(Curve_generic, AlgebraicScheme_subscheme_toric):
    def __init__(self, A, X, *kwargs):
        # TODO ensure that A is the right type?
        # Something like a `is_WeightProjectiveSpace` which means making a
        # WeightProjectiveSpace class?
        super().__init__(A, X, *kwargs)
