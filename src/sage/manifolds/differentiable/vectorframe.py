r"""
Vector Frames

The class :class:`VectorFrame` implements vector frames on differentiable
manifolds.
By *vector frame*, it is meant a field `e` on some
differentiable manifold `U` endowed with a differentiable map
`\Phi: U \rightarrow M` to a differentiable manifold `M` such that for each
`p\in U`, `e(p)` is a vector basis of the tangent space `T_{\Phi(p)}M`.

The standard case of a vector frame *on* `U` corresponds to `U = M`
and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
immersion and `\Phi` being a curve in `M` (`U` is then an open
interval of `\RR`).

A derived class of :class:`VectorFrame` is :class:`CoordFrame`;
it regards the vector frames associated with a chart, i.e. the
so-called *coordinate bases*.

The vector frame duals, i.e. the coframes, are implemented via the class
:class:`CoFrame`. The derived class :class:`CoordCoFrame` is devoted to
coframes deriving from a chart.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Travis Scrimshaw (2016): review tweaks
- Eric Gourgoulhon (2018): some refactoring and more functionalities in the
  choice of symbols for vector frame elements (:issue:`24792`)

REFERENCES:

- [Lee2013]_

EXAMPLES:

Introducing a chart on a manifold automatically endows it with a vector frame:
the coordinate frame associated to the chart::

    sage: M = Manifold(3, 'M')
    sage: X.<x,y,z> = M.chart()
    sage: M.frames()
    [Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))]
    sage: M.frames()[0] is X.frame()
    True

A new vector frame can be defined from a family of 3 linearly independent
vector fields::

    sage: e1 = M.vector_field(1, x, y)
    sage: e2 = M.vector_field(z, -2, x*y)
    sage: e3 = M.vector_field(1, 1, 0)
    sage: e = M.vector_frame('e', (e1, e2, e3)); e
    Vector frame (M, (e_0,e_1,e_2))
    sage: latex(e)
    \left(M, \left(e_{0},e_{1},e_{2}\right)\right)

The first frame defined on a manifold is its *default frame*; in the present
case it is the coordinate frame associated to the chart ``X``::

    sage: M.default_frame()
    Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))

The default frame can be changed via the method
:meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_default_frame`::

    sage: M.set_default_frame(e)
    sage: M.default_frame()
    Vector frame (M, (e_0,e_1,e_2))

The elements of a vector frame are vector fields on the manifold::

    sage: for vec in e:
    ....:     print(vec)
    ....:
    Vector field e_0 on the 3-dimensional differentiable manifold M
    Vector field e_1 on the 3-dimensional differentiable manifold M
    Vector field e_2 on the 3-dimensional differentiable manifold M

Each element of a vector frame can be accessed by its index::

    sage: e[0]
    Vector field e_0 on the 3-dimensional differentiable manifold M
    sage: e[0].display(X.frame())
    e_0 = ∂/∂x + x ∂/∂y + y ∂/∂z
    sage: X.frame()[1]
    Vector field ∂/∂y on the 3-dimensional differentiable manifold M
    sage: X.frame()[1].display(e)
    ∂/∂y = x/(x^2 - x + z + 2) e_0 - 1/(x^2 - x + z + 2) e_1
     - (x - z)/(x^2 - x + z + 2) e_2

The slice operator ``:`` can be used to access to more than one element::

    sage: e[0:2]
    (Vector field e_0 on the 3-dimensional differentiable manifold M,
     Vector field e_1 on the 3-dimensional differentiable manifold M)
    sage: e[:]
    (Vector field e_0 on the 3-dimensional differentiable manifold M,
     Vector field e_1 on the 3-dimensional differentiable manifold M,
     Vector field e_2 on the 3-dimensional differentiable manifold M)

Vector frames can be constructed from scratch, without any connection to
previously defined frames or vector fields (the connection can be performed
later via the method
:meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_change_of_frame`)::

    sage: f = M.vector_frame('f'); f
    Vector frame (M, (f_0,f_1,f_2))
    sage: M.frames()
    [Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z)),
     Vector frame (M, (e_0,e_1,e_2)),
     Vector frame (M, (f_0,f_1,f_2))]

The index range depends on the starting index defined on the manifold::

    sage: M = Manifold(3, 'M', start_index=1)
    sage: X.<x,y,z> = M.chart()
    sage: e = M.vector_frame('e')
    sage: [e[i] for i in M.irange()]
    [Vector field e_1 on the 3-dimensional differentiable manifold M,
     Vector field e_2 on the 3-dimensional differentiable manifold M,
     Vector field e_3 on the 3-dimensional differentiable manifold M]
    sage: e[1], e[2], e[3]
    (Vector field e_1 on the 3-dimensional differentiable manifold M,
     Vector field e_2 on the 3-dimensional differentiable manifold M,
     Vector field e_3 on the 3-dimensional differentiable manifold M)

Let us check that the vector fields ``e[i]`` are the frame vectors from
their components with respect to the frame `e`::

    sage: e[1].comp(e)[:]
    [1, 0, 0]
    sage: e[2].comp(e)[:]
    [0, 1, 0]
    sage: e[3].comp(e)[:]
    [0, 0, 1]

Defining a vector frame on a manifold automatically creates the dual
coframe, which, by default, bares the same name (here `e`)::

    sage: M.coframes()
    [Coordinate coframe (M, (dx,dy,dz)), Coframe (M, (e^1,e^2,e^3))]
    sage: f = M.coframes()[1] ; f
    Coframe (M, (e^1,e^2,e^3))
    sage: f is e.coframe()
    True

Each element of the coframe is a 1-form::

    sage: f[1], f[2], f[3]
    (1-form e^1 on the 3-dimensional differentiable manifold M,
    1-form e^2 on the 3-dimensional differentiable manifold M,
    1-form e^3 on the 3-dimensional differentiable manifold M)
    sage: latex(f[1]), latex(f[2]), latex(f[3])
    (e^{1}, e^{2}, e^{3})

Let us check that the coframe `(e^i)` is indeed the dual of the vector
frame `(e_i)`::

    sage: f[1](e[1]) # the 1-form e^1 applied to the vector field e_1
    Scalar field e^1(e_1) on the 3-dimensional differentiable manifold M
    sage: f[1](e[1]).expr() # the explicit expression of e^1(e_1)
    1
    sage: f[1](e[1]).expr(), f[1](e[2]).expr(), f[1](e[3]).expr()
    (1, 0, 0)
    sage: f[2](e[1]).expr(), f[2](e[2]).expr(), f[2](e[3]).expr()
    (0, 1, 0)
    sage: f[3](e[1]).expr(), f[3](e[2]).expr(), f[3](e[3]).expr()
    (0, 0, 1)

The coordinate frame associated to spherical coordinates of the
sphere `S^2`::

    sage: M = Manifold(2, 'S^2', start_index=1) # Part of S^2 covered by spherical coord.
    sage: c_spher.<th,ph> = M.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
    sage: b = M.default_frame() ; b
    Coordinate frame (S^2, (∂/∂th,∂/∂ph))
    sage: b[1]
    Vector field ∂/∂th on the 2-dimensional differentiable manifold S^2
    sage: b[2]
    Vector field ∂/∂ph on the 2-dimensional differentiable manifold S^2

The orthonormal frame constructed from the coordinate frame::

    sage: e = M.vector_frame('e', (b[1], b[2]/sin(th))); e
    Vector frame (S^2, (e_1,e_2))
    sage: e[1].display()
    e_1 = ∂/∂th
    sage: e[2].display()
    e_2 = 1/sin(th) ∂/∂ph

The change-of-frame automorphisms and their matrices::

    sage: M.change_of_frame(c_spher.frame(), e)
    Field of tangent-space automorphisms on the 2-dimensional
     differentiable manifold S^2
    sage: M.change_of_frame(c_spher.frame(), e)[:]
    [        1         0]
    [        0 1/sin(th)]
    sage: M.change_of_frame(e, c_spher.frame())
    Field of tangent-space automorphisms on the 2-dimensional
     differentiable manifold S^2
    sage: M.change_of_frame(e, c_spher.frame())[:]
    [      1       0]
    [      0 sin(th)]
"""

# *****************************************************************************
#       Copyright (C) 2013-2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.tensor.modules.free_module_basis import (FreeModuleBasis,
                                                   FreeModuleCoBasis)
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.misc.cachefunc import cached_method


class CoFrame(FreeModuleCoBasis):
    r"""
    Coframe on a differentiable manifold.

    By *coframe*, it is meant a field `f` on some differentiable manifold `U`
    endowed with a differentiable map `\Phi: U \rightarrow M` to a
    differentiable manifold `M` such that for each `p\in U`, `f(p)` is a basis
    of the vector space `T^*_{\Phi(p)}M` (the dual to the tangent space
    `T_{\Phi(p)}M`).

    The standard case of a coframe *on* `U` corresponds to `U = M` and
    `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    INPUT:

    - ``frame`` -- the vector frame dual to the coframe
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the 1-forms constituting the coframe, or a tuple of strings,
      representing the individual symbols of the 1-forms
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the 1-forms constituting the
      coframe, or a tuple of strings, representing the individual LaTeX symbols
      of the 1-forms; if ``None``, ``symbol`` is used in place of
      ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the 1-forms
      of the coframe; if ``None``, the indices will be generated as integers
      within the range declared on the coframe's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the 1-forms of the coframe; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

    Coframe on a 3-dimensional manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x,y,z> = M.chart()
        sage: v = M.vector_frame('v')
        sage: from sage.manifolds.differentiable.vectorframe import CoFrame
        sage: e = CoFrame(v, 'e') ; e
        Coframe (M, (e^1,e^2,e^3))

    Instead of importing CoFrame in the global namespace, the coframe can be
    obtained by means of the method
    :meth:`~sage.tensor.modules.free_module_basis.FreeModuleBasis.dual_basis`;
    the symbol is then the same as that of the frame::

        sage: a = v.dual_basis() ; a
        Coframe (M, (v^1,v^2,v^3))
        sage: a[1] == e[1]
        True
        sage: a[1] is e[1]
        False
        sage: e[1].display(v)
        e^1 = v^1

    The 1-forms composing the coframe are obtained via the operator ``[]``::

        sage: e[1], e[2], e[3]
        (1-form e^1 on the 3-dimensional differentiable manifold M,
         1-form e^2 on the 3-dimensional differentiable manifold M,
         1-form e^3 on the 3-dimensional differentiable manifold M)

    Checking that `e` is the dual of `v`::

        sage: e[1](v[1]).expr(), e[1](v[2]).expr(), e[1](v[3]).expr()
        (1, 0, 0)
        sage: e[2](v[1]).expr(), e[2](v[2]).expr(), e[2](v[3]).expr()
        (0, 1, 0)
        sage: e[3](v[1]).expr(), e[3](v[2]).expr(), e[3](v[3]).expr()
        (0, 0, 1)
    """
    def __init__(self, frame, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        Construct a coframe, dual to a given vector frame.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: from sage.manifolds.differentiable.vectorframe import CoFrame
            sage: f = CoFrame(e, 'f'); f
            Coframe (M, (f^0,f^1))
            sage: TestSuite(f).run()
        """
        self._domain = frame._domain
        self._manifold = self._domain.manifold()
        FreeModuleCoBasis.__init__(self, frame, symbol,
                                   latex_symbol=latex_symbol, indices=indices,
                                   latex_indices=latex_indices)
        # The coframe is added to the domain's set of coframes, as well as to
        # all the superdomains' sets of coframes
        for sd in self._domain.open_supersets():
            sd._coframes.append(self)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: f = e.coframe()
            sage: f._repr_()
            'Coframe (M, (e^0,e^1))'
            sage: repr(f)  # indirect doctest
            'Coframe (M, (e^0,e^1))'
            sage: f  # indirect doctest
            Coframe (M, (e^0,e^1))

        Test with a nontrivial destination map::

            sage: N = Manifold(3, 'N', start_index=1)
            sage: phi = M.diff_map(N)
            sage: h = M.vector_frame('h', dest_map=phi)
            sage: h.coframe()._repr_()
            'Coframe (M, (h^1,h^2,h^3)) with values on the 3-dimensional differentiable manifold N'
        """
        description = "Coframe " + self._name
        dest_map = self._basis.destination_map()
        if dest_map is not self._domain.identity_map():
            description += " with values on the {}".format(dest_map.codomain())
        return description

    def at(self, point):
        r"""
        Return the value of ``self`` at a given point on the manifold, this
        value being a basis of the dual of the tangent space at the point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` in the domain `U` of the coframe (denoted `f` hereafter)

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_basis.FreeModuleCoBasis`
          representing the basis `f(p)` of the vector space `T^*_{\Phi(p)} M`,
          dual to the tangent space `T_{\Phi(p)} M`, where
          `\Phi: U \to M` is the differentiable map associated with
          `f` (possibly `\Phi = \mathrm{Id}_U`)

        EXAMPLES:

        Cobasis of a tangent space on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((-1,2), name='p')
            sage: f = X.coframe() ; f
            Coordinate coframe (M, (dx,dy))
            sage: fp = f.at(p) ; fp
            Dual basis (dx,dy) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: type(fp)
            <class 'sage.tensor.modules.free_module_basis.FreeModuleCoBasis_with_category'>
            sage: fp[0]
            Linear form dx on the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: fp[1]
            Linear form dy on the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: fp is X.frame().at(p).dual_basis()
            True
        """
        return self._basis.at(point).dual_basis()

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='up',
                 include_domain=True):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the 1-forms constituting the coframe, or a list/tuple of
          strings, representing the individual symbols of the 1-forms
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the 1-forms constituting
          the coframe, or a list/tuple of strings, representing the individual
          LaTeX symbols of the 1-forms; if ``None``, ``symbol`` is used in
          place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the 1-forms of the coframe; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the 1-forms;
          if ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'up'``) determines the position
          of the indices labelling the 1-forms of the coframe; can be
          either ``'down'`` or ``'up'``
        - ``include_domain`` -- boolean (default: ``True``); determining whether
          the name of the domain is included in the beginning of the coframe
          name

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e').coframe(); e
            Coframe (M, (e^0,e^1))
            sage: e.set_name('f'); e
            Coframe (M, (f^0,f^1))
            sage: e.set_name('e', latex_symbol=r'\epsilon')
            sage: latex(e)
            \left(M, \left(\epsilon^{0},\epsilon^{1}\right)\right)
            sage: e.set_name('e', include_domain=False); e
            Coframe (e^0,e^1)
            sage: e.set_name(['a', 'b'], latex_symbol=[r'\alpha', r'\beta']); e
            Coframe (M, (a,b))
            sage: latex(e)
            \left(M, \left(\alpha,\beta\right)\right)
            sage: e.set_name('e', indices=['x','y'],
            ....:            latex_indices=[r'\xi', r'\zeta']); e
            Coframe (M, (e^x,e^y))
            sage: latex(e)
            \left(M, \left(e^{\xi},e^{\zeta}\right)\right)
        """
        super().set_name(symbol, latex_symbol=latex_symbol,
                         indices=indices,
                         latex_indices=latex_indices,
                         index_position=index_position)
        if include_domain:
            # Redefinition of the name and the LaTeX name to include the domain
            self._name = "({}, {})".format(self._domain._name, self._name)
            self._latex_name = r"\left({}, {}\right)".format(
                                    self._domain._latex_name, self._latex_name)

# ******************************************************************************


class VectorFrame(FreeModuleBasis):
    r"""
    Vector frame on a differentiable manifold.

    By *vector frame*, it is meant a field `e` on some
    differentiable manifold `U` endowed with a differentiable map
    `\Phi: U\rightarrow M` to a differentiable manifold `M` such that for
    each `p\in U`, `e(p)` is a vector basis of the tangent space
    `T_{\Phi(p)}M`.

    The standard case of a vector frame *on* `U` corresponds to `U=M`
    and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi` being an
    immersion and `\Phi` being a curve in `M` (`U` is then an open interval
    of `\RR`).

    For each instantiation of a vector frame, a coframe is automatically
    created, as an instance of the class :class:`CoFrame`. It is returned by
    the method :meth:`coframe`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathfrak{X}(U, \Phi)`
      of vector fields along `U` with values on `M \supset \Phi(U)`
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the vector fields constituting the vector frame, or a tuple
      of strings, representing the individual symbols of the vector fields
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the vector fields constituting
      the vector frame, or a tuple of strings, representing the individual
      LaTeX symbols of the vector fields; if ``None``, ``symbol`` is used in
      place of ``latex_symbol``
    - ``from_frame`` -- (default: ``None``) vector frame `\tilde e` on the
      codomain `M` of the destination map `\Phi`; the constructed frame `e`
      is then such that `\forall p \in U, e(p) = \tilde{e}(\Phi(p))`
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the vector
      fields of the frame; if ``None``, the indices will be generated as
      integers within the range declared on the vector frame's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the vector fields; if
      ``None``, ``indices`` is used instead
    - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
      dual coframe; if ``None``, ``symbol`` must be a string and is used
      for the common base of the symbols of the elements of the dual coframe
    - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
      but for the dual coframe

    EXAMPLES:

    Defining a vector frame on a 3-dimensional manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e') ; e
        Vector frame (M, (e_1,e_2,e_3))
        sage: latex(e)
        \left(M, \left(e_{1},e_{2},e_{3}\right)\right)

    The individual elements of the vector frame are accessed via square
    brackets, with the possibility to invoke the slice operator '``:``' to
    get more than a single element::

        sage: e[2]
        Vector field e_2 on the 3-dimensional differentiable manifold M
        sage: e[1:3]
        (Vector field e_1 on the 3-dimensional differentiable manifold M,
         Vector field e_2 on the 3-dimensional differentiable manifold M)
        sage: e[:]
        (Vector field e_1 on the 3-dimensional differentiable manifold M,
         Vector field e_2 on the 3-dimensional differentiable manifold M,
         Vector field e_3 on the 3-dimensional differentiable manifold M)

    The LaTeX symbol can be specified::

        sage: E = M.vector_frame('E', latex_symbol=r"\epsilon")
        sage: latex(E)
        \left(M, \left(\epsilon_{1},\epsilon_{2},\epsilon_{3}\right)\right)

    By default, the elements of the vector frame are labelled by integers
    within the range specified at the manifold declaration. It is however
    possible to fully customize the labels, via the argument ``indices``::

        sage: u = M.vector_frame('u', indices=('x', 'y', 'z')) ; u
        Vector frame (M, (u_x,u_y,u_z))
        sage: u[1]
        Vector field u_x on the 3-dimensional differentiable manifold M
        sage: u.coframe()
        Coframe (M, (u^x,u^y,u^z))

    The LaTeX format of the indices can be adjusted::

        sage: v = M.vector_frame('v', indices=('a', 'b', 'c'),
        ....:                    latex_indices=(r'\alpha', r'\beta', r'\gamma'))
        sage: v
        Vector frame (M, (v_a,v_b,v_c))
        sage: latex(v)
        \left(M, \left(v_{\alpha},v_{\beta},v_{\gamma}\right)\right)
        sage: latex(v.coframe())
        \left(M, \left(v^{\alpha},v^{\beta},v^{\gamma}\right)\right)

    The symbol of each element of the vector frame can also be freely chosen,
    by providing a tuple of symbols as the first argument of ``vector_frame``;
    it is then mandatory to specify as well some symbols for the dual coframe::

        sage: h = M.vector_frame(('a', 'b', 'c'), symbol_dual=('A', 'B', 'C'))
        sage: h
        Vector frame (M, (a,b,c))
        sage: h[1]
        Vector field a on the 3-dimensional differentiable manifold M
        sage: h.coframe()
        Coframe (M, (A,B,C))
        sage: h.coframe()[1]
        1-form A on the 3-dimensional differentiable manifold M

    Example with a non-trivial map `\Phi` (see above); a vector frame along a
    curve::

        sage: U = Manifold(1, 'U')  # open interval (-1,1) as a 1-dimensional manifold
        sage: T.<t> = U.chart('t:(-1,1)')  # canonical chart on U
        sage: Phi = U.diff_map(M, [cos(t), sin(t), t], name='Phi',
        ....:                  latex_name=r'\Phi')
        sage: Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold U
         to the 3-dimensional differentiable manifold M
        sage: f = U.vector_frame('f', dest_map=Phi) ; f
        Vector frame (U, (f_1,f_2,f_3)) with values on the 3-dimensional
         differentiable manifold M
        sage: f.domain()
        1-dimensional differentiable manifold U
        sage: f.ambient_domain()
        3-dimensional differentiable manifold M

    The value of the vector frame at a given point is a basis of the
    corresponding tangent space::

        sage: p = U((0,), name='p') ; p
        Point p on the 1-dimensional differentiable manifold U
        sage: f.at(p)
        Basis (f_1,f_2,f_3) on the Tangent space at Point Phi(p) on the
         3-dimensional differentiable manifold M

    Vector frames are bases of free modules formed by vector fields::

        sage: e.module()
        Free module X(M) of vector fields on the 3-dimensional differentiable
         manifold M
        sage: e.module().base_ring()
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M
        sage: e.module() is M.vector_field_module()
        True
        sage: e in M.vector_field_module().bases()
        True

    ::

        sage: f.module()
        Free module X(U,Phi) of vector fields along the 1-dimensional
         differentiable manifold U mapped into the 3-dimensional differentiable
         manifold M
        sage: f.module().base_ring()
        Algebra of differentiable scalar fields on the 1-dimensional
         differentiable manifold U
        sage: f.module() is U.vector_field_module(dest_map=Phi)
        True
        sage: f in U.vector_field_module(dest_map=Phi).bases()
        True
    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = CoFrame

    @staticmethod
    def __classcall_private__(cls, vector_field_module, symbol,
                              latex_symbol=None, from_frame=None, indices=None,
                              latex_indices=None, symbol_dual=None,
                              latex_symbol_dual=None):
        """
        Transform input lists into tuples for the unique representation of
        VectorFrame.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module(force_free=True)
            sage: from sage.manifolds.differentiable.vectorframe import VectorFrame
            sage: e = VectorFrame(XM, ['a', 'b'], symbol_dual=['A', 'B']); e
            Vector frame (M, (a,b))
            sage: e.dual_basis()
            Coframe (M, (A,B))
            sage: e is VectorFrame(XM, ('a', 'b'), symbol_dual=('A', 'B'))
            True
        """
        if isinstance(symbol, list):
            symbol = tuple(symbol)
        if isinstance(latex_symbol, list):
            latex_symbol = tuple(latex_symbol)
        if isinstance(indices, list):
            indices = tuple(indices)
        if isinstance(latex_indices, list):
            latex_indices = tuple(latex_indices)
        if isinstance(symbol_dual, list):
            symbol_dual = tuple(symbol_dual)
        if isinstance(latex_symbol_dual, list):
            latex_symbol_dual = tuple(latex_symbol_dual)
        return super().__classcall__(cls, vector_field_module,
                                     symbol, latex_symbol=latex_symbol,
                                     from_frame=from_frame, indices=indices,
                                     latex_indices=latex_indices,
                                     symbol_dual=symbol_dual,
                                     latex_symbol_dual=latex_symbol_dual)

    def __init__(self, vector_field_module, symbol, latex_symbol=None,
                 from_frame=None, indices=None, latex_indices=None,
                 symbol_dual=None, latex_symbol_dual=None):
        r"""
        Construct a vector frame on a parallelizable manifold.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module(force_free=True)  # makes M parallelizable
            sage: from sage.manifolds.differentiable.vectorframe import VectorFrame
            sage: e = VectorFrame(XM, 'e', latex_symbol=r'\epsilon'); e
            Vector frame (M, (e_0,e_1))
            sage: TestSuite(e).run()
        """
        from sage.manifolds.differentiable.manifold import DifferentiableManifold
        # Some sanity check:
        if not isinstance(vector_field_module, FiniteRankFreeModule):
            raise ValueError("the {} has already been constructed as a "
                             "non-free module and therefore cannot have "
                             "a basis".format(vector_field_module))
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._dest_map = vector_field_module._dest_map
        self._from_frame = from_frame
        self._manifold = self._domain.manifold()
        if from_frame is not None:
            if not from_frame._domain.is_subset(self._dest_map._codomain):
                raise ValueError("the domain of the frame 'from_frame' is " +
                                 "not included in the codomain of the " +
                                 "destination map")
        if symbol is None:
            if from_frame is None:
                raise TypeError("some frame symbol must be provided")
            symbol = from_frame._symbol
            latex_symbol = from_frame._latex_symbol
            indices = from_frame._indices
            latex_indices = from_frame._latex_indices
            symbol_dual = from_frame._symbol_dual
            latex_symbol_dual = from_frame._latex_symbol_dual
        FreeModuleBasis.__init__(self, vector_field_module,
                                 symbol, latex_symbol=latex_symbol,
                                 indices=indices, latex_indices=latex_indices,
                                 symbol_dual=symbol_dual,
                                 latex_symbol_dual=latex_symbol_dual)
        # The frame is added to the domain's set of frames, as well as to all
        # the superdomains' sets of frames; moreover the first defined frame
        # is considered as the default one
        dest_map = self._dest_map
        for sd in self._domain.open_supersets():
            sd._frames.append(self)
            sd._top_frames.append(self)
            if sd._def_frame is None:
                sd._def_frame = self
            if isinstance(sd, DifferentiableManifold):
                # Initialization of the zero elements of tensor field modules:
                if dest_map in sd._vector_field_modules:
                    xsd = sd._vector_field_modules[dest_map]
                    if not isinstance(xsd, FiniteRankFreeModule):
                        for t in xsd._tensor_modules.values():
                            t.zero()._add_comp_unsafe(self)
                            # (since new components are initialized to zero)
        if dest_map is self._domain.identity_map():
            # The frame is added to the list of the domain's covering frames:
            self._domain._set_covering_frame(self)

        # Dual coframe
        self._coframe = self.dual_basis()  # self._coframe = a shortcut for
                                           # self._dual_basis

        # Derived quantities:
        # Initialization of the set of frames that are restrictions of the
        # current frame to subdomains of the frame domain:
        self._subframes = set([self])
        # Initialization of the set of frames which the current frame is a
        # restriction of:
        self._superframes = set([self])

        self._restrictions = {} # dict. of the restrictions of self to
                               # subdomains of self._domain, with the
                               # subdomains as keys
        # NB: set(self._restrictions.values()) is identical to
        #     self._subframes

    ###### Methods that must be redefined by derived classes of ######
    ###### FreeModuleBasis                                      ######

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e._repr_()
            'Vector frame (M, (e_0,e_1))'
            sage: repr(e)  # indirect doctest
            'Vector frame (M, (e_0,e_1))'
            sage: e  # indirect doctest
            Vector frame (M, (e_0,e_1))

        Test with a nontrivial destination map::

            sage: N = Manifold(3, 'N', start_index=1)
            sage: phi = M.diff_map(N)
            sage: h = M.vector_frame('h', dest_map=phi)
            sage: h._repr_()
            'Vector frame (M, (h_1,h_2,h_3)) with values on the 3-dimensional differentiable manifold N'
        """
        description = "Vector frame " + self._name
        if self._dest_map is not self._domain.identity_map():
            description += " with values on the {}".format(self._dest_map.codomain())
        return description

    def _new_instance(self, symbol, latex_symbol=None, indices=None,
                      latex_indices=None, symbol_dual=None,
                      latex_symbol_dual=None):
        r"""
        Construct a new vector frame on the same vector field module
        as ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the vector fields constituting the vector frame, or a
          tuple of strings, representing the individual symbols of the vector
          fields
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the vector fields
          constituting the vector frame, or a tuple of strings, representing
          the individual LaTeX symbols of the vector fields; if ``None``,
          ``symbol`` is used in place of ``latex_symbol``
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde e` on the
          codomain `M` of the destination map `\Phi`; the constructed frame `e`
          is then such that `\forall p \in U, e(p) = \tilde{e}(\Phi(p))`
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the vector fields of the frame; if ``None``, the indices will be
          generated as integers within the range declared on the vector frame's
          domain
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the vector fields;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT: instance of :class:`VectorFrame`

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e._new_instance('f')
            Vector frame (M, (f_0,f_1))
        """
        return VectorFrame(self._fmodule, symbol, latex_symbol=latex_symbol,
                           indices=indices, latex_indices=latex_indices,
                           symbol_dual=symbol_dual,
                           latex_symbol_dual=latex_symbol_dual)

    ###### End of methods to be redefined by derived classes ######

    def domain(self):
        r"""
        Return the domain on which ``self`` is defined.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
          representing the domain of the vector frame

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e.domain()
            2-dimensional differentiable manifold M
            sage: U = M.open_subset('U')
            sage: f = e.restrict(U)
            sage: f.domain()
            Open subset U of the 2-dimensional differentiable manifold M
        """
        return self._domain

    def ambient_domain(self):
        r"""
        Return the differentiable manifold in which ``self`` takes its values.

        The ambient domain is the codomain `M` of the differentiable map
        `\Phi: U \rightarrow M` associated with the frame.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
          representing `M`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e.ambient_domain()
            2-dimensional differentiable manifold M

        In the present case, since `\Phi` is the identity map::

            sage: e.ambient_domain() == e.domain()
            True

        An example with a non trivial map `\Phi`::

            sage: U = Manifold(1, 'U')
            sage: T.<t> = U.chart()
            sage: X.<x,y> = M.chart()
            sage: Phi = U.diff_map(M, {(T,X): [cos(t), t]}, name='Phi',
            ....:                  latex_name=r'\Phi') ; Phi
            Differentiable map Phi from the 1-dimensional differentiable
             manifold U to the 2-dimensional differentiable manifold M
            sage: f = U.vector_frame('f', dest_map=Phi); f
            Vector frame (U, (f_0,f_1)) with values on the 2-dimensional
             differentiable manifold M
            sage: f.ambient_domain()
            2-dimensional differentiable manifold M
            sage: f.domain()
            1-dimensional differentiable manifold U
        """
        return self._ambient_domain

    def destination_map(self):
        r"""
        Return the differential map associated to this vector frame.

        Let `e` denote the vector frame; the differential map associated to
        it is the map `\Phi: U\rightarrow M` such that for each `p \in U`,
        `e(p)` is a vector basis of the tangent space `T_{\Phi(p)}M`.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the differential map `\Phi`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e.destination_map()
            Identity map Id_M of the 2-dimensional differentiable manifold M

        An example with a non trivial map `\Phi`::

            sage: U = Manifold(1, 'U')
            sage: T.<t> = U.chart()
            sage: X.<x,y> = M.chart()
            sage: Phi = U.diff_map(M, {(T,X): [cos(t), t]}, name='Phi',
            ....:                  latex_name=r'\Phi') ; Phi
            Differentiable map Phi from the 1-dimensional differentiable
             manifold U to the 2-dimensional differentiable manifold M
            sage: f = U.vector_frame('f', dest_map=Phi); f
            Vector frame (U, (f_0,f_1)) with values on the 2-dimensional
             differentiable manifold M
            sage: f.destination_map()
            Differentiable map Phi from the 1-dimensional differentiable
             manifold U to the 2-dimensional differentiable manifold M
        """
        return self._dest_map

    def coframe(self):
        r"""
        Return the coframe of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e')
            sage: e.coframe()
            Coframe (M, (e^0,e^1))
            sage: X.<x,y> = M.chart()
            sage: X.frame().coframe()
            Coordinate coframe (M, (dx,dy))
        """
        return self._coframe

    def new_frame(self, change_of_frame, symbol, latex_symbol=None,
                  indices=None, latex_indices=None, symbol_dual=None,
                  latex_symbol_dual=None):
        r"""
        Define a new vector frame from ``self``.

        The new vector frame is defined from a field of tangent-space
        automorphisms; its domain is the same as that of the current frame.

        INPUT:

        - ``change_of_frame`` --
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`;
          the field of tangent space automorphisms `P` that relates
          the current frame `(e_i)` to the new frame `(n_i)` according
          to `n_i = P(e_i)`
        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the vector fields constituting the vector frame, or a
          list/tuple of strings, representing the individual symbols of the
          vector fields
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the vector fields
          constituting the vector frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the vector fields;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the vector fields of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the vector fields;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT:

        - the new frame `(n_i)`, as an instance of :class:`VectorFrame`

        EXAMPLES:

        Frame resulting from a `\pi/3`-rotation in the Euclidean plane::

            sage: M = Manifold(2, 'R^2')
            sage: X.<x,y> = M.chart()
            sage: e = M.vector_frame('e') ; M.set_default_frame(e)
            sage: M._frame_changes
            {}
            sage: rot = M.automorphism_field()
            sage: rot[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: n = e.new_frame(rot, 'n')
            sage: n[0][:]
            [1/2*sqrt(3), 1/2]
            sage: n[1][:]
            [-1/2, 1/2*sqrt(3)]
            sage: a =  M.change_of_frame(e,n)
            sage: a[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a == rot
            True
            sage: a is rot
            False
            sage: a._components # random (dictionary output)
            {Vector frame (R^2, (e_0,e_1)): 2-indices components w.r.t.
             Vector frame (R^2, (e_0,e_1)),
             Vector frame (R^2, (n_0,n_1)): 2-indices components w.r.t.
             Vector frame (R^2, (n_0,n_1))}
            sage: a.comp(n)[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a1 = M.change_of_frame(n,e)
            sage: a1[:]
            [1/2*sqrt(3)         1/2]
            [       -1/2 1/2*sqrt(3)]
            sage: a1 == rot.inverse()
            True
            sage: a1 is rot.inverse()
            False
            sage: e[0].comp(n)[:]
            [1/2*sqrt(3), -1/2]
            sage: e[1].comp(n)[:]
            [1/2, 1/2*sqrt(3)]
        """
        the_new_frame = self.new_basis(change_of_frame, symbol,
                                       latex_symbol=latex_symbol,
                                       indices=indices,
                                       latex_indices=latex_indices,
                                       symbol_dual=symbol_dual,
                                       latex_symbol_dual=latex_symbol_dual)
        for sdom in self._domain.open_supersets():
            sdom._frame_changes[(self, the_new_frame)] = \
                            self._fmodule._basis_changes[(self, the_new_frame)]
            sdom._frame_changes[(the_new_frame, self)] = \
                            self._fmodule._basis_changes[(the_new_frame, self)]
        return the_new_frame

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some open subset of its domain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `V` of the current frame domain `U`

        OUTPUT: the restriction of the current frame to `V` as a :class:`VectorFrame`

        EXAMPLES:

        Restriction of a frame defined on `\RR^2` to the unit disk::

            sage: M = Manifold(2, 'R^2', start_index=1)
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: a = M.automorphism_field()
            sage: a[:] = [[1-y^2,0], [1+x^2, 2]]
            sage: e = c_cart.frame().new_frame(a, 'e') ; e
            Vector frame (R^2, (e_1,e_2))
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1})
            sage: e_U = e.restrict(U) ; e_U
            Vector frame (U, (e_1,e_2))

        The vectors of the restriction have the same symbols as those of the
        original frame::

            sage: e_U[1].display()
            e_1 = (-y^2 + 1) ∂/∂x + (x^2 + 1) ∂/∂y
            sage: e_U[2].display()
            e_2 = 2 ∂/∂y

        They are actually the restrictions of the original frame vectors::

            sage: e_U[1] is e[1].restrict(U)
            True
            sage: e_U[2] is e[2].restrict(U)
            True
        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the provided domain is not a subdomain of " +
                                 "the current frame's domain")
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom) and subdomain in rst._restrictions:
                    res = rst._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    res._superframes.update(self._superframes)
                    for sframe2 in self._superframes:
                        sframe2._subframes.add(res)
                    return self._restrictions[subdomain]
            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for sframe in self._superframes:
                if subdomain in sframe._restrictions:
                    res = sframe._restrictions[subdomain]
                    self._restrictions[subdomain] = res
                    res._superframes.update(self._superframes)
                    for sframe2 in self._superframes:
                        sframe2._subframes.add(res)
                    return self._restrictions[subdomain]
            # If this point is reached, the restriction has to be created
            # from scratch
            sdest_map = self._dest_map.restrict(subdomain)
            resmodule = subdomain.vector_field_module(sdest_map,
                                                      force_free=True)
            if subdomain in self._restrictions:
                # the restriction has been generated during the creation of
                # resmodule (which may happen if sdest_map is not trivial)
                return self._restrictions[subdomain]
            res = VectorFrame(resmodule,
                              self._symbol, latex_symbol=self._latex_symbol,
                              indices=self._indices,
                              latex_indices=self._latex_indices,
                              symbol_dual=self._symbol_dual,
                              latex_symbol_dual=self._latex_symbol_dual)

            res._from_frame = self._from_frame

            for dom in subdomain.open_supersets():
                if dom is not subdomain:
                    dom._top_frames.remove(res)  # since it was added by
                                                 # VectorFrame constructor
            new_vectors = list()
            for i in self._fmodule.irange():
                vrest = self[i].restrict(subdomain)
                for j in self._fmodule.irange():
                    vrest.add_comp(res)[j] = 0
                vrest.add_comp(res)[i] = 1
                new_vectors.append(vrest)
            res._vec = tuple(new_vectors)
            # Update of superframes and subframes:
            for sframe in self._subframes:
                if subdomain.is_subset(sframe.domain()):
                    res._superframes.update(sframe._superframes)
            for sframe in res._superframes:
                sframe._subframes.add(res)
                sframe._restrictions[subdomain] = res # includes sframe = self
            for dom, rst in self._restrictions.items():
                if dom.is_subset(subdomain):
                    res._restrictions.update(rst._restrictions)
                    res._subframes.update(rst._subframes)
                    rst._superframes.update(res._superframes)

        return self._restrictions[subdomain]

    @cached_method
    def structure_coeff(self):
        r"""
        Evaluate the structure coefficients associated to ``self``.

        `n` being the manifold's dimension, the structure coefficients of the
        vector frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}`
        defined by

        .. MATH::

            [e_i, e_j] = C^k_{\ \, ij} e_k

        OUTPUT:

        - the structure coefficients `C^k_{\ \, ij}`, as an instance of
          :class:`~sage.tensor.modules.comp.CompWithSym`
          with 3 indices ordered as `(k,i,j)`.

        EXAMPLES:

        Structure coefficients of the orthonormal frame associated to
        spherical coordinates in the Euclidean space `\RR^3`::

            sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)  # Part of R^3 covered by spherical coordinates
            sage: c_spher.<r,th,ph> = M.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: ch_frame = M.automorphism_field()
            sage: ch_frame[1,1], ch_frame[2,2], ch_frame[3,3] = 1, 1/r, 1/(r*sin(th))
            sage: M.frames()
            [Coordinate frame (R^3, (∂/∂r,∂/∂th,∂/∂ph))]
            sage: e = c_spher.frame().new_frame(ch_frame, 'e')
            sage: e[1][:]  # components of e_1 in the manifold's default frame (∂/∂r, ∂/∂th, ∂/∂th)
            [1, 0, 0]
            sage: e[2][:]
            [0, 1/r, 0]
            sage: e[3][:]
            [0, 0, 1/(r*sin(th))]
            sage: c = e.structure_coeff() ; c
            3-indices components w.r.t. Vector frame (R^3, (e_1,e_2,e_3)), with
             antisymmetry on the index positions (1, 2)
            sage: c[:]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, -1/r, 0], [1/r, 0, 0], [0, 0, 0]],
             [[0, 0, -1/r], [0, 0, -cos(th)/(r*sin(th))], [1/r, cos(th)/(r*sin(th)), 0]]]
            sage: c[2,1,2]  # C^2_{12}
            -1/r
            sage: c[3,1,3]  # C^3_{13}
            -1/r
            sage: c[3,2,3]  # C^3_{23}
            -cos(th)/(r*sin(th))
        """
        from sage.tensor.modules.comp import CompWithSym

        fmodule = self._fmodule
        structure_coeff = CompWithSym(self._fmodule._ring, self, 3,
                                      start_index=fmodule._sindex,
                                      output_formatter=fmodule._output_formatter,
                                      antisym=(1,2))
        si = fmodule._sindex
        nsi = si + fmodule.rank()
        for k in range(si, nsi):
            ce_k = self._coframe._vec[k-si]
            for i in range(si, nsi):
                e_i = self._vec[i-si]
                for j in range(i+1, nsi):
                    e_j = self._vec[j-si]
                    structure_coeff[[k,i,j]] = ce_k(e_j.lie_der(e_i))
        return structure_coeff

    def along(self, mapping):
        r"""
        Return the vector frame deduced from the current frame via a
        differentiable map, the codomain of which is included in the domain of
        of the current frame.

        If `e` is the current vector frame, `V` its domain and if
        `\Phi: U \rightarrow V` is a differentiable map from some
        differentiable manifold `U` to `V`, the returned object is
        a vector frame `\tilde e` along `U` with values on `V` such that

        .. MATH::

           \forall p \in U,\  \tilde e(p) = e(\Phi(p)).

        INPUT:

        - ``mapping`` -- differentiable map `\Phi: U \rightarrow V`

        OUTPUT: vector frame `\tilde e` along `U` defined above

        EXAMPLES:

        Vector frame along a curve::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R = Manifold(1, 'R')  # R as a 1-dimensional manifold
            sage: T.<t> = R.chart()  # canonical chart on R
            sage: Phi = R.diff_map(M, {(T,X): [cos(t), t]}, name='Phi',
            ....:                  latex_name=r'\Phi') ; Phi
            Differentiable map Phi from the 1-dimensional differentiable
             manifold R to the 2-dimensional differentiable manifold M
            sage: e = X.frame() ; e
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: te = e.along(Phi) ; te
            Vector frame (R, (∂/∂x,∂/∂y)) with values on the 2-dimensional
             differentiable manifold M

        Check of the formula `\tilde e(p) = e(\Phi(p))`::

            sage: p = R((pi,)) ; p
            Point on the 1-dimensional differentiable manifold R
            sage: te[0].at(p) == e[0].at(Phi(p))
            True
            sage: te[1].at(p) == e[1].at(Phi(p))
            True

        The result is cached::

            sage: te is e.along(Phi)
            True
        """
        dom = self._domain
        if mapping.codomain().is_subset(dom):
            rmapping = mapping
        else:
            rmapping = None
            for doms, rest in mapping._restrictions.items():
                if doms[1].is_subset(dom):
                    rmapping = rest
                    break
            else:
                raise ValueError("the codomain of {} is not ".format(mapping) +
                                 " included in the domain of {}".format(self))
        vmodule = rmapping.domain().vector_field_module(dest_map=rmapping)
        return vmodule.basis(from_frame=self)

    def at(self, point):
        r"""
        Return the value of ``self`` at a given point, this value being
        a basis of the tangent vector space at the point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`; point
          `p` in the domain `U` of the vector frame (denoted `e` hereafter)

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the basis `e(p)` of the tangent vector space
          `T_{\Phi(p)} M`, where `\Phi: U \to M` is the differentiable map
          associated with `e` (possibly `\Phi = \mathrm{Id}_U`)

        EXAMPLES:

        Basis of a tangent space to a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((-1,2), name='p')
            sage: e = X.frame() ; e
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: ep = e.at(p) ; ep
            Basis (∂/∂x,∂/∂y) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: type(ep)
            <class 'sage.tensor.modules.free_module_basis.FreeModuleBasis_with_category'>
            sage: ep[0]
            Tangent vector ∂/∂x at Point p on the 2-dimensional differentiable
             manifold M
            sage: ep[1]
            Tangent vector ∂/∂y at Point p on the 2-dimensional differentiable
             manifold M

        Note that the symbols used to denote the vectors are same as those
        for the vector fields of the frame. At this stage, ``ep`` is the unique
        basis on the tangent space at ``p``::

            sage: Tp = M.tangent_space(p)
            sage: Tp.bases()
            [Basis (∂/∂x,∂/∂y) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M]

        Let us consider a vector frame that is a not a coordinate one::

            sage: aut = M.automorphism_field()
            sage: aut[:] = [[1+y^2, 0], [0, 2]]
            sage: f = e.new_frame(aut, 'f') ; f
            Vector frame (M, (f_0,f_1))
            sage: fp = f.at(p) ; fp
            Basis (f_0,f_1) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M

        There are now two bases on the tangent space::

            sage: Tp.bases()
            [Basis (∂/∂x,∂/∂y) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M,
             Basis (f_0,f_1) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M]

        Moreover, the changes of bases in the tangent space have been
        computed from the known relation between the frames ``e`` and
        ``f`` (field of automorphisms ``aut`` defined above)::

            sage: Tp.change_of_basis(ep, fp)
            Automorphism of the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: Tp.change_of_basis(ep, fp).display()
            5 ∂/∂x⊗dx + 2 ∂/∂y⊗dy
            sage: Tp.change_of_basis(fp, ep)
            Automorphism of the Tangent space at Point p on the 2-dimensional
             differentiable manifold M
            sage: Tp.change_of_basis(fp, ep).display()
            1/5 ∂/∂x⊗dx + 1/2 ∂/∂y⊗dy

        The dual bases::

            sage: e.coframe()
            Coordinate coframe (M, (dx,dy))
            sage: ep.dual_basis()
            Dual basis (dx,dy) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: ep.dual_basis() is e.coframe().at(p)
            True
            sage: f.coframe()
            Coframe (M, (f^0,f^1))
            sage: fp.dual_basis()
            Dual basis (f^0,f^1) on the Tangent space at Point p on the
             2-dimensional differentiable manifold M
            sage: fp.dual_basis() is f.coframe().at(p)
            True
        """
        # Case of a non-trivial destination map
        if self._from_frame is not None:
            if self._dest_map.is_identity():  #!# probably not necessary
                raise ValueError("the destination map should not be the identity")
            ambient_point = self._dest_map(point)
            return self._from_frame.at(ambient_point)

        # Determination of the tangent space:
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))

        if self._dest_map.is_identity():
            ambient_point = point
        else:
            ambient_point = self._dest_map(point)
        ts = ambient_point._manifold.tangent_space(ambient_point)

        # If the basis has already been constructed, it is simply returned:
        ts_frame_bases = ts._frame_bases
        if self in ts_frame_bases:
            return ts_frame_bases[self]
        for frame in ts_frame_bases:
            if self in frame._subframes or self in frame._superframes:
                return ts_frame_bases[frame]

        # If this point is reached, the basis has to be constructed from
        # scratch.
        # The names of the basis vectors set to those of the frame vector
        # fields:
        basis = ts.basis(self._symbol, latex_symbol=self._latex_symbol,
                         indices=self._indices,
                         latex_indices=self._latex_indices,
                         symbol_dual=self._symbol_dual,
                         latex_symbol_dual=self._latex_symbol_dual)
        ts_frame_bases[self] = basis
        # Update of the change of bases in the tangent space:
        for frame_pair, automorph in self._domain._frame_changes.items():
            frame1 = frame_pair[0]
            frame2 = frame_pair[1]
            if frame1 is self:
                fr2 = None
                for frame in ts_frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr2 is not None:
                    basis1 = basis
                    basis2 = ts_frame_bases[fr2]
                    auto = ts.automorphism()
                    for frame, comp in automorph._components.items():
                        bas = None
                        if frame is frame1:
                            bas = basis1
                        if frame is frame2:
                            bas = basis2
                        if bas is not None:
                            cauto = auto.add_comp(bas)
                            for ind, val in comp._comp.items():
                                cauto._comp[ind] = val(point)
                    ts._basis_changes[(basis1, basis2)] = auto
            if frame2 is self:
                fr1 = None
                for frame in ts_frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                if fr1 is not None:
                    basis1 = ts_frame_bases[fr1]
                    basis2 = basis
                    auto = ts.automorphism()
                    for frame, comp in automorph._components.items():
                        bas = None
                        if frame is frame1:
                            bas = basis1
                        if frame is frame2:
                            bas = basis2
                        if bas is not None:
                            cauto = auto.add_comp(bas)
                            for ind, val in comp._comp.items():
                                cauto._comp[ind] = val(point)
                    ts._basis_changes[(basis1, basis2)] = auto
        return basis

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='down',
                 include_domain=True):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the vector fields constituting the vector frame, or a
          list/tuple of strings, representing the individual symbols of the
          vector fields
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the vector fields
          constituting the vector frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the vector fields;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the vector fields of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the vector fields;
          if ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'down'``) determines the position
          of the indices labelling the vector fields of the frame; can be
          either ``'down'`` or ``'up'``
        - ``include_domain`` -- boolean (default: ``True``); determining whether
          the name of the domain is included in the beginning of the vector
          frame name

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: e = M.vector_frame('e'); e
            Vector frame (M, (e_0,e_1))
            sage: e.set_name('f'); e
            Vector frame (M, (f_0,f_1))
            sage: e.set_name('e', include_domain=False); e
            Vector frame (e_0,e_1)
            sage: e.set_name(['a', 'b']); e
            Vector frame (M, (a,b))
            sage: e.set_name('e', indices=['x', 'y']); e
            Vector frame (M, (e_x,e_y))
            sage: e.set_name('e', latex_symbol=r'\epsilon')
            sage: latex(e)
            \left(M, \left(\epsilon_{0},\epsilon_{1}\right)\right)
            sage: e.set_name('e', latex_symbol=[r'\alpha', r'\beta'])
            sage: latex(e)
            \left(M, \left(\alpha,\beta\right)\right)
            sage: e.set_name('e', latex_symbol='E',
            ....:            latex_indices=[r'\alpha', r'\beta'])
            sage: latex(e)
            \left(M, \left(E_{\alpha},E_{\beta}\right)\right)
        """
        super().set_name(symbol, latex_symbol=latex_symbol,
                         indices=indices,
                         latex_indices=latex_indices,
                         index_position=index_position)
        if include_domain:
            # Redefinition of the name and the LaTeX name to include the domain
            self._name = "({}, {})".format(self._domain._name, self._name)
            self._latex_name = r"\left({}, {}\right)".format(
                                    self._domain._latex_name, self._latex_name)

#******************************************************************************

class CoordCoFrame(CoFrame):
    r"""
    Coordinate coframe on a differentiable manifold.

    By *coordinate coframe*, it is meant the `n`-tuple of the
    differentials of the coordinates of some chart on the manifold,
    with `n` being the manifold's dimension.

    INPUT:

    - ``coord_frame`` -- coordinate frame dual to the coordinate coframe
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the 1-forms constituting the coframe, or a tuple of strings,
      representing the individual symbols of the 1-forms
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the 1-forms constituting the
      coframe, or a tuple of strings, representing the individual LaTeX symbols
      of the 1-forms; if ``None``, ``symbol`` is used in place of
      ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the 1-forms
      of the coframe; if ``None``, the indices will be generated as
      integers within the range declared on the vector frame's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the 1-forms of the coframe; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

    Coordinate coframe on a 3-dimensional manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x,y,z> = M.chart()
        sage: M.frames()
        [Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))]
        sage: M.coframes()
        [Coordinate coframe (M, (dx,dy,dz))]
        sage: dX = M.coframes()[0] ; dX
        Coordinate coframe (M, (dx,dy,dz))

    The 1-forms composing the coframe are obtained via the operator ``[]``::

        sage: dX[1]
        1-form dx on the 3-dimensional differentiable manifold M
        sage: dX[2]
        1-form dy on the 3-dimensional differentiable manifold M
        sage: dX[3]
        1-form dz on the 3-dimensional differentiable manifold M
        sage: dX[1][:]
        [1, 0, 0]
        sage: dX[2][:]
        [0, 1, 0]
        sage: dX[3][:]
        [0, 0, 1]

    The coframe is the dual of the coordinate frame::

        sage: e = X.frame() ; e
        Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))
        sage: dX[1](e[1]).expr(), dX[1](e[2]).expr(), dX[1](e[3]).expr()
        (1, 0, 0)
        sage: dX[2](e[1]).expr(), dX[2](e[2]).expr(), dX[2](e[3]).expr()
        (0, 1, 0)
        sage: dX[3](e[1]).expr(), dX[3](e[2]).expr(), dX[3](e[3]).expr()
        (0, 0, 1)

    Each 1-form of a coordinate coframe is closed::

        sage: dX[1].exterior_derivative()
        2-form ddx on the 3-dimensional differentiable manifold M
        sage: dX[1].exterior_derivative() == 0
        True
    """
    def __init__(self, coord_frame, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        Construct a coordinate coframe.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.differentiable.vectorframe import CoordCoFrame
            sage: f = CoordCoFrame(X.frame(), 'omega'); f
            Coordinate coframe (M, (omega^0,omega^1))
            sage: TestSuite(f).run()
        """
        if not isinstance(coord_frame, CoordFrame):
            raise TypeError("the first argument must be a coordinate frame")
        CoFrame.__init__(self, coord_frame, symbol, latex_symbol=latex_symbol,
                         indices=indices, latex_indices=latex_indices)
        self._chart = coord_frame._chart

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = X.frame().coframe()
            sage: f._repr_()
            'Coordinate coframe (M, (dx,dy))'
            sage: repr(f)  # indirect doctest
            'Coordinate coframe (M, (dx,dy))'
            sage: f  # indirect doctest
            Coordinate coframe (M, (dx,dy))
        """
        return "Coordinate coframe " + self._name

#******************************************************************************

class CoordFrame(VectorFrame):
    r"""
    Coordinate frame on a differentiable manifold.

    By *coordinate frame*, it is meant a vector frame on a differentiable
    manifold `M` that is associated to a coordinate chart on `M`.

    INPUT:

    - ``chart`` -- the chart defining the coordinates

    EXAMPLES:

    The coordinate frame associated to spherical coordinates of the
    sphere `S^2`::

        sage: M = Manifold(2, 'S^2', start_index=1)  # Part of S^2 covered by spherical coord.
        sage: M.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
        Chart (S^2, (th, ph))
        sage: b = M.default_frame()
        sage: b
        Coordinate frame (S^2, (∂/∂th,∂/∂ph))
        sage: b[1]
        Vector field ∂/∂th on the 2-dimensional differentiable manifold S^2
        sage: b[2]
        Vector field ∂/∂ph on the 2-dimensional differentiable manifold S^2
        sage: latex(b)
        \left(S^2, \left(\frac{\partial}{\partial {\theta} },\frac{\partial}{\partial {\phi} }\right)\right)
    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = CoordCoFrame

    def __init__(self, chart):
        r"""
        Construct a coordinate frame.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.differentiable.vectorframe import CoordFrame
            sage: e = CoordFrame(X); e
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: TestSuite(e).run()
        """
        from sage.misc.latex import latex
        from sage.typeset.unicode_characters import unicode_partial
        from sage.manifolds.differentiable.chart import DiffChart
        if not isinstance(chart, DiffChart):
            raise TypeError("the first argument must be a chart")
        dom = chart.domain()
        # Some sanity check:
        vmodule = dom._vector_field_modules.get(dom.identity_map())
        if vmodule and not isinstance(vmodule, FiniteRankFreeModule):
            raise ValueError("the {} has already been constructed as a "
                             "non-free module, which implies that the {} is "
                             "not parallelizable and hence cannot be the "
                             "domain of a coordinate chart".format(vmodule,
                             dom))
        self._chart = chart
        coords = chart[:] # list of all coordinates
        symbol = tuple(f"{unicode_partial}/{unicode_partial}{x!s}"
                       for x in coords)
        latex_symbol = tuple(r"\frac{\partial}{\partial" + latex(x) + "}"
                             for x in coords)
        symbol_dual = tuple("d" + str(x) for x in coords)
        latex_symbol_dual = tuple(r"\mathrm{d}" + latex(x) for x in coords)
        VectorFrame.__init__(self,
                             dom.vector_field_module(force_free=True),
                             symbol=symbol, latex_symbol=latex_symbol,
                             symbol_dual=symbol_dual,
                             latex_symbol_dual=latex_symbol_dual)
        # In the above:
        # - force_free=True ensures that a free module is constructed in case
        #   it is the first call to the vector field module on chart.domain()

    ###### Methods that must be redefined by derived classes of ######
    ###### FreeModuleBasis                                      ######

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e = X.frame()
            sage: e._repr_()
            'Coordinate frame (M, (∂/∂x,∂/∂y))'
            sage: repr(e)  # indirect doctest
            'Coordinate frame (M, (∂/∂x,∂/∂y))'
            sage: e  # indirect doctest
            Coordinate frame (M, (∂/∂x,∂/∂y))
        """
        return "Coordinate frame " + self._name

    ###### End of methods redefined by derived classes ######

    def chart(self):
        r"""
        Return the chart defining this coordinate frame.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e = X.frame()
            sage: e.chart()
            Chart (M, (x, y))
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: e.restrict(U).chart()
            Chart (U, (x, y))
        """
        return self._chart

    @cached_method
    def structure_coeff(self):
        r"""
        Return the structure coefficients associated to ``self``.

        `n` being the manifold's dimension, the structure coefficients
        of the frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}`
        defined by

        .. MATH::

            [e_i, e_j] = C^k_{\ \, ij} e_k.

        In the present case, since `(e_i)` is a coordinate frame,
        `C^k_{\ \, ij}=0`.

        OUTPUT:

        - the structure coefficients `C^k_{\ \, ij}`, as a vanishing instance
          of :class:`~sage.tensor.modules.comp.CompWithSym` with 3 indices
          ordered as `(k,i,j)`

        EXAMPLES:

        Structure coefficients of the coordinate frame associated to
        spherical coordinates in the Euclidean space `\RR^3`::

            sage: M = Manifold(3, 'R^3', r'\RR^3', start_index=1)  # Part of R^3 covered by spherical coord.
            sage: c_spher = M.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: b = M.default_frame() ; b
            Coordinate frame (R^3, (∂/∂r,∂/∂th,∂/∂ph))
            sage: c = b.structure_coeff() ; c
            3-indices components w.r.t. Coordinate frame
             (R^3, (∂/∂r,∂/∂th,∂/∂ph)), with antisymmetry on the index
             positions (1, 2)
            sage: c == 0
            True
        """
        from sage.tensor.modules.comp import CompWithSym
        # A zero CompWithSym
        return CompWithSym(self._fmodule._ring, self, 3,
                           start_index=self._fmodule._sindex,
                           output_formatter=self._fmodule._output_formatter,
                           antisym=(1,2))
