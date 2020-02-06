"""
Highlevel functions for managing options of the libSingular interface

AUTHOR:

- Martin Albrecht
"""

class LibSingularGBDefaultContext:
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.standard_options import LibSingularGBDefaultContext
            sage: from sage.libs.singular.option import opt
            sage: P.<a,b,c> = PolynomialRing(QQ, 3, order='lex')
            sage: I = sage.rings.ideal.Katsura(P, 3)
            sage: I.gens()
            [a + 2*b + 2*c - 1, a^2 - a + 2*b^2 + 2*c^2, 2*a*b + 2*b*c - b]
            sage: opt['red_tail'] = False
            sage: opt['red_sb'] = False
            sage: import sage.libs.singular.function_factory
            sage: groebner = sage.libs.singular.function_factory.ff.groebner
            sage: rgb = groebner(I)
            sage: rgb
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, a + 2*b + 2*c - 1]
            sage: with LibSingularGBDefaultContext(): rgb = groebner(I)
            sage: rgb
            [84*c^4 - 40*c^3 + c^2 + c, 7*b + 210*c^3 - 79*c^2 + 3*c, 7*a - 420*c^3 + 158*c^2 + 8*c - 7]

        """
        from sage.libs.singular.option import opt_ctx
        self.libsingular_option_context = opt_ctx

    def __enter__(self):
        """
        EXAMPLES::

            sage: from sage.libs.singular.standard_options import LibSingularGBDefaultContext
            sage: from sage.libs.singular.option import opt
            sage: P.<a,b,c,d> = PolynomialRing(QQ, 4, order='lex')
            sage: I = sage.rings.ideal.Katsura(P, 4)
            sage: I.gens()
            [a + 2*b + 2*c + 2*d - 1, a^2 - a + 2*b^2 + 2*c^2 + 2*d^2, 2*a*b + 2*b*c - b + 2*c*d, 2*a*c + b^2 + 2*b*d - c]

            sage: opt['red_tail'] = False
            sage: opt['red_sb'] = False
            sage: import sage.libs.singular.function_factory
            sage: groebner = sage.libs.singular.function_factory.ff.groebner
            sage: rgb = groebner(I)
            sage: rgb
            [128304*d^8 - 93312*d^7 + 15552*d^6 + 3144*d^5 - 1120*d^4 + 36*d^3 + 15*d^2 - d,
             5913075*c + 371438283744*d^7 - 237550027104*d^6 + 22645939824*d^5 + 11520686172*d^4 - 2024910556*d^3 - 132524276*d^2 + 30947828*d,
             1044*b + 10976*c^3 + 67914*c^2*d - 9709*c^2 + 100296*c*d^2 - 32862*c*d + 2412*c + 64302*d^3 - 29483*d^2 + 2683*d,
             a + 2*b + 2*c + 2*d - 1]

            sage: with LibSingularGBDefaultContext(): rgb = groebner(I)
            sage: rgb
            [128304*d^8 - 93312*d^7 + 15552*d^6 + 3144*d^5 - 1120*d^4 + 36*d^3 + 15*d^2 - d,
             5913075*c + 371438283744*d^7 - 237550027104*d^6 + 22645939824*d^5 + 11520686172*d^4 - 2024910556*d^3 - 132524276*d^2 + 30947828*d,
             1971025*b - 97197721632*d^7 + 73975630752*d^6 - 12121915032*d^5 - 2760941496*d^4 + 814792828*d^3 - 1678512*d^2 - 9158924*d,
             5913075*a - 159690237696*d^7 + 31246269696*d^6 + 27439610544*d^5 - 6475723368*d^4 - 838935856*d^3 + 275119624*d^2 + 4884038*d - 5913075]

        """
        self.libsingular_option_context.__enter__()
        self.libsingular_option_context.opt.reset_default()
        self.libsingular_option_context.opt['red_sb'] = True
        self.libsingular_option_context.opt['red_tail'] = True
        self.libsingular_option_context.opt['deg_bound'] = 0
        self.libsingular_option_context.opt['mult_bound'] = 0

    def __exit__(self, typ, value, tb):
        """
        EXAMPLES::

            sage: from sage.libs.singular.standard_options import LibSingularGBDefaultContext
            sage: from sage.libs.singular.option import opt
            sage: P.<a,b,c,d> = PolynomialRing(GF(7), 4, order='lex')
            sage: I = sage.rings.ideal.Katsura(P, 4)
            sage: I.gens()
            [a + 2*b + 2*c + 2*d - 1, a^2 - a + 2*b^2 + 2*c^2 + 2*d^2, 2*a*b + 2*b*c - b + 2*c*d, 2*a*c + b^2 + 2*b*d - c]

            sage: opt['red_tail'] = False
            sage: opt['red_sb'] = False
            sage: import sage.libs.singular.function_factory
            sage: groebner = sage.libs.singular.function_factory.ff.groebner
            sage: rgb = groebner(I)
            sage: rgb
            [d^7 + 3*d^6 - d^5 + 3*d^4 + d^3 - d^2 + 3*d,
             c - 3*d^6 + 3*d^5 - d^4 + d^3 - 2*d,
             b + 3*c*d - 3*c + d^2 + 2*d,
             a + 2*b + 2*c + 2*d - 1]
            sage: with LibSingularGBDefaultContext(): rgb = groebner(I)
            sage: rgb
             [d^7 + 3*d^6 - d^5 + 3*d^4 + d^3 - d^2 + 3*d,
              c - 3*d^6 + 3*d^5 - d^4 + d^3 - 2*d,
              b - 3*d^6 + 2*d^4 + d^3 + 2*d^2 - 3*d,
              a - 2*d^6 + d^5 - 2*d^4 + 3*d^3 + 3*d^2 - 2*d - 1]
        """
        self.libsingular_option_context.__exit__(typ,value,tb)

def libsingular_gb_standard_options(func):
    """
    Decorator to force a reduced Singular groebner basis.

    TESTS::

        sage: P.<a,b,c,d,e> = PolynomialRing(GF(127))
        sage: J = sage.rings.ideal.Cyclic(P).homogenize()
        sage: from sage.misc.sageinspect import sage_getsource
        sage: "basis.reduced" in sage_getsource(J.interreduced_basis)
        True

    The following tests against a bug that was fixed in :trac:`11298`::

        sage: from sage.misc.sageinspect import sage_getsourcelines, sage_getargspec
        sage: P.<x,y> = QQ[]
        sage: I = P*[x,y]
        sage: sage_getargspec(I.interreduced_basis)
        ArgSpec(args=['self'], varargs=None, keywords=None, defaults=None)
        sage: sage_getsourcelines(I.interreduced_basis)
        (['    @handle_AA_and_QQbar\n',
          '    @singular_gb_standard_options\n',
          '    @libsingular_gb_standard_options\n',
          '    def interreduced_basis(self):\n',
          ...
          '        return self.basis.reduced()\n'], ...)

    .. note::

       This decorator is used automatically internally so the user
       does not need to use it manually.
    """
    from sage.misc.decorators import sage_wraps
    @sage_wraps(func)
    def wrapper(*args, **kwds):
        """
        Execute function in ``LibSingularGBDefaultContext``.
        """
        with LibSingularGBDefaultContext():
            return func(*args, **kwds)
    return wrapper
