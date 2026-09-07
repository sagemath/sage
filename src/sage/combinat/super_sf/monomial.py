r"""
Monomial Supersymmetric Functions

AUTHORS:

- Shriya M
"""
from . import super_sfa
from sage.data_structures.blas_dict import convert_remove_zeroes
from sage.libs.lrcalc import lrcalc
from sage.matrix.constructor import matrix
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.partition import Partitions
from sage.categories.tensor import tensor

class SupersymFunctionAlgebra_monomial(super_sfa.SuperSymAlgebra_generic):
    r"""
    Monomial supersymmetric functions.

    The *monomial supersymmetric function* defined on variables `\mathbf{x}` and
    `\mathbf{y}`, `m_\lambda(\mathbf{x} \mid \mathbf{y})`, is given in terms of
    monomial symmetric functions and forgotten symmetric functions:

    .. MATH::

       m_\lambda(\mathbf{x} \mid \mathbf{y}) = \sum_{\mu \cup \nu = \lambda} m_\mu(\mathbf{x}) f_\nu(\mathbf{y})

    These form a non-graded multiplicative basis for the ring of supersymmetric
    functions.

    REFERENCES:

    - [BHS25]_

    INPUT:

    - ``Supersym`` -- ring of supersymmetric functions

    EXAMPLES::

        sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
        sage: s = SuperSymmetricFunctions(QQ)
        sage: m = s.m()
        sage: m
        Supersymmetric functions over Rational Field in the Schur basis
    """
    def __init__(self, SuperSym):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.super_sf.super_sf import SuperSymmetricFunctions
            sage: s = SuperSymmetricFunctions(QQ)
            sage: m = s.m()
            sage: TestSuite(m).run()
        """
        super().__init__(SuperSym=SuperSym, graded=False, prefix='m', basis_name='monomial')