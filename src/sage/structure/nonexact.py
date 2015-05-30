"Precision management for non-exact objects"

import sage.rings.integer

class Nonexact:
    def __init__(self, prec=20):
        self.__default_prec = sage.rings.integer.Integer(prec)

    def default_prec(self):
        r"""
        Return the default precision for self.  Use
        ``set_default_prec`` to set the default precision.
        """
        try:
            return self.__default_prec
        except AttributeError:
            self.default_prec = 20
            return self.__default_prec
        return self.__default_prec

    def set_default_prec(self, prec):
        self.__default_prec = sage.rings.integer.Integer(prec)


