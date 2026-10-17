from typing import assert_type

from sage.misc.latex import LatexExpr


assert_type(LatexExpr("a") + LatexExpr("b"), LatexExpr)
assert_type(LatexExpr("a") + "b", LatexExpr)
assert_type("a = " + LatexExpr("b"), LatexExpr)
