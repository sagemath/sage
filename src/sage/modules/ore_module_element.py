r"""
Elements in Ore modules

AUTHOR:

- Xavier Caruso (2024-10)
"""

# ***************************************************************************
#    Copyright (C) 2024 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.misc.latex import latex
from sage.modules.free_module_element import FreeModuleElement_generic_dense

class OreModuleElement(FreeModuleElement_generic_dense):
    r"""
    A generic element of a Ore module.
    """
    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((t+1)*X^2 + 1)

            sage: v - w  # indirect doctest
            v - w
            sage: w - v  # indirect doctest
            -v + w
            sage: X^5*v  # indirect doctest
            ((-4*t + 2)/(t^4 + 4*t^3 + 6*t^2 + 4*t + 1))*v + ((t - 5)/(t^3 + 3*t^2 + 3*t + 1))*w
        """
        parent = self.parent()
        if parent._names is None:
            return self.parent()._repr_element(self)
        else:
            rank = parent.rank()
            names = parent._names
            s = ""
            for i in range(rank):
                c = self[i]
                sc = str(c)
                if sc == "0":
                    continue
                if sc == "1":
                    s += " + %s" % names[i]
                elif sc == "-1":
                    s += " - %s" % names[i]
                elif c._is_atomic():
                    if sc[0] == "-":
                        s += " - %s*%s" % (-c, names[i])
                    else:
                        s += " + %s*%s" % (sc, names[i])
                else:
                    s += " + (%s)*%s" % (c, names[i])
            if s == "":
                return "0"
            elif s[1] == '-':
                return '-' + s[3:]
            else:
                return s[3:]

    def _latex_(self):
        r"""
        Return a LaTeX representation of this element.

        EXAMPLES::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module((t+1)*X^2 + 1)

            sage: latex(v - w)
            v - w
            sage: latex(w - v)
            - v + w
            sage: latex(X^5*v)
            \left(\frac{-4 t + 2}{t^{4} + 4 t^{3} + 6 t^{2} + 4 t + 1}\right) v + \left(\frac{t - 5}{t^{3} + 3 t^{2} + 3 t + 1}\right) w
        """
        parent = self.parent()
        if parent._names is None:
            return self.parent()._latex_element(self)
        else:
            rank = parent.rank()
            names = parent._latex_names
            s = ""
            for i in range(rank):
                c = self[i]
                sc = str(c)
                if sc == "0":
                    continue
                if sc == "1":
                    s += " + %s" % names[i]
                elif sc == "-1":
                    s += " - %s" % names[i]
                elif c._is_atomic():
                    if sc[0] == "-":
                        s += " - %s %s" % (latex(-c), names[i])
                    else:
                        s += " + %s %s" % (latex(c), names[i])
                else:
                    s += " + \\left(%s\\right) %s" % (latex(c), names[i])
            if s == "":
                return "0"
            elif s[1] == '-':
                return s[1:]
            else:
                return s[3:]

    def is_mutable(self):
        r"""
        Always return ``False`` since elements in Ore modules
        are all immutable.

        EXAMPLES::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M = S.quotient_module(X^2 + t)

            sage: v, w = M.basis()
            sage: v
            (1, 0)
            sage: v.is_mutable()
            False
            sage: v[1] = 1
            Traceback (most recent call last):
            ...
            ValueError: vectors in Ore modules are immutable
        """
        return False

    def __setitem__(self, i, v):
        r"""
        Always raise an error since elements in Ore modules are
        all immutable.

        TESTS::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module(X^2 + t)
            sage: w[1] = 0
            Traceback (most recent call last):
            ...
            ValueError: vectors in Ore modules are immutable
        """
        raise ValueError("vectors in Ore modules are immutable")

    def __hash__(self):
        r"""
        Return a hash of this element.

        TESTS::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module(X^2 + t)
            sage: hash(v)  # random
            -5164621852614943976
            sage: hash(w)  # random
            -1950498447580522560
        """
        return hash(tuple(self))

    def vector(self):
        r"""
        Return the coordinates vector of this element.

        EXAMPLES::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module(X^2 + t)
            sage: v.vector()
            (1, 0)

        We underline that this vector is not an element of the
        Ore module; it lives in `K^2`. Compare::

            sage: v.parent()
            Ore module <v, w> over Fraction Field of Univariate Polynomial Ring in t over Rational Field twisted by d/dt
            sage: v.vector().parent()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in t over Rational Field
        """
        V = self.parent().module()
        return V(self.list())

    def image(self):
        r"""
        Return the image of this element by the pseudomorphism
        defining the action of the Ore variable on this Ore module.

        EXAMPLES::

            sage: K.<t> = Frac(QQ['t'])
            sage: S.<X> = OrePolynomialRing(K, K.derivation())
            sage: M.<v,w> = S.quotient_module(X^2 + t)
            sage: v.image()
            w
            sage: w.image()
            -t*v

        TESTS:

        We check that this corresponds to the action of `X`::

            sage: x = M.random_element()
            sage: x.image() == X*x
            True
        """
        return self.parent()._pseudohom(self)
