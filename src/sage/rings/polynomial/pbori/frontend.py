# sage_setup: distribution = sagemath-brial
# sage.doctest: needs sage.rings.polynomial.pbori
# Import basic functionality
r"""
This module defines an initial ring, and patches the declare_ring to use
a given context.

EXAMPLES::

    sage: from sage.rings.polynomial.pbori.frontend import x
    sage: x(0)
    x(0)
    sage: x(0)*x(0)
    x(0)
    sage: x(0) + x(0)
    0
    sage: x(9999)
    x(9999)
    sage: x(9999)*x(9999)
    x(9999)
    sage: x(9999) + x(9999)
    0

    sage: from sage.rings.polynomial.pbori.frontend import x, polybori_start
    sage: context = dict(globals())
    sage: polybori_start(context)
    ipbori...
    sage: r = context['declare_ring']('abc')
    sage: context['a']
    a
    sage: r.variable(0)
    a
"""


from .PyPolyBoRi import Ring
from .pbori import VariableFactory
from .blocks import declare_ring as orig_declare_ring


def block_scheme_names(blocks):
    r"""
    Helper for Singular interface.
    """
    context = {}
    from .blocks import declare_block_scheme
    declare_block_scheme(blocks, context)

    return list(context.keys())


ipbname = 'ipbori'


def polybori_start(global_context):
    def declare_ring(blocks, context=None):
        if context is None:
            context = global_context

        return orig_declare_ring(blocks, context)
    declare_ring.__doc__ = orig_declare_ring.__doc__
    global_context["declare_ring"] = declare_ring

    print(ipbname + """ -- The interactive command line tool of PolyBoRi/BRiAL %s
""" % global_context.get("polybori_version", ''))


# Here come the defaults
r = Ring(10000)
x = VariableFactory(r)
