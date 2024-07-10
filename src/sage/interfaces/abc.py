r"""
Abstract base classes for interface elements
"""


class AxiomElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.axiom.AxiomElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.AxiomElement.__subclasses__()) <= 1
        True
    """
    pass


class ExpectElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.expect.ExpectElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.ExpectElement.__subclasses__()) <= 1
        True
    """
    pass


class FriCASElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.fricas.FriCASElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.FriCASElement.__subclasses__()) <= 1
        True
    """
    pass


class GapElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.gap.GapElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.GapElement.__subclasses__()) <= 1
        True
    """
    pass


class GpElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.gp.GpElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.GpElement.__subclasses__()) <= 1
        True
    """
    pass


class Macaulay2Element:
    r"""
    Abstract base class for :class:`~sage.interfaces.macaulay2.Macaulay2Element`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.Macaulay2Element.__subclasses__()) <= 1
        True
    """
    pass


class MagmaElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.magma.MagmaElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.MagmaElement.__subclasses__()) <= 1
        True
    """
    pass


class SingularElement:
    r"""
    Abstract base class for :class:`~sage.interfaces.singular.SingularElement`.

    This class is defined for the purpose of ``isinstance`` tests.  It should not be
    instantiated.

    EXAMPLES:

    By design, there is a unique direct subclass::

        sage: len(sage.interfaces.abc.SingularElement.__subclasses__()) <= 1
        True
    """
    pass
