r"""
Local Frames

The class :class:`LocalFrame` implements local frames on vector bundles
(see :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle` or
:class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`).

For `k=0,1,\dots`, a *local frame* on a vector bundle `E \to M` of class `C^k`
and rank `n` is a local section `(e_1,\dots,e_n):U \to E^n` of class `C^k`
defined on some subset `U` of the base space `M`, such that `e(p)` is a basis of
the fiber `E_p` for any `p \in U`.

AUTHORS:

- Michael Jung (2019): initial version

EXAMPLES:

Defining a global frame on a topological vector bundle of rank 3::

    sage: M = Manifold(3, 'M', structure='top')
    sage: E = M.vector_bundle(3, 'E')
    sage: e = E.local_frame('e'); e
    Local frame (E|_M, (e_0,e_1,e_2))

This frame is now the default frame of the corresponding section module and
saved in the vector bundle::

    sage: e in E.frames()
    True
    sage: sec_module = E.section_module(); sec_module
    Free module C^0(M;E) of sections on the 3-dimensional topological manifold M
     with values in the real vector bundle E of rank 3
    sage: sec_module.default_basis()
    Local frame (E|_M, (e_0,e_1,e_2))

However, the default frame can be changed::

    sage: sec_module.set_default_basis(e)
    sage: sec_module.default_basis()
    Local frame (E|_M, (e_0,e_1,e_2))

The elements of a local frame are local sections in the vector bundle::

    sage: for vec in e:
    ....:     print(vec)
    Section e_0 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3
     Section e_1 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3
     Section e_2 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3

Each element of a vector frame can be accessed by its index::

    sage: e[0]
    Section e_0 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3

The slice operator ``:`` can be used to access to more than one element::

    sage: e[0:2]
    (Section e_0 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3,
     Section e_1 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3)
    sage: e[:]
    (Section e_0 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3,
     Section e_1 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3,
     Section e_2 on the 3-dimensional topological manifold M with values in the
     real vector bundle E of rank 3)

The index range depends on the starting index defined on the manifold::

    sage: M = Manifold(3, 'M', structure='top', start_index=1)
    sage: c_xyz.<x,y,z> = M.chart()
    sage: U = M.open_subset('U')
    sage: c_xyz_U = c_xyz.restrict(U)
    sage: E = M.vector_bundle(3, 'E')
    sage: e = E.local_frame('e', domain=U); e
    Local frame (E|_U, (e_1,e_2,e_3))
    sage: [e[i] for i in M.irange()]
    [Section e_1 on the Open subset U of the 3-dimensional topological manifold
     M with values in the real vector bundle E of rank 3,
     Section e_2 on the Open subset U of the 3-dimensional topological manifold
     M with values in the real vector bundle E of rank 3,
     Section e_3 on the Open subset U of the 3-dimensional topological manifold
     M with values in the real vector bundle E of rank 3]
    sage: e[1], e[2], e[3]
    (Section e_1 on the Open subset U of the 3-dimensional topological manifold
     M with values in the real vector bundle E of rank 3,
     Section e_2 on the Open subset U of the 3-dimensional topological manifold
     M with values in the real vector bundle E of rank 3,
     Section e_3 on the Open subset U of the 3-dimensional topological manifold
     M with values in the real vector bundle E of rank 3)

Let us check that the local sections ``e[i]`` are indeed the frame vectors
from their components with respect to the frame `e`::

    sage: e[1].comp(e)[:]
    [1, 0, 0]
    sage: e[2].comp(e)[:]
    [0, 1, 0]
    sage: e[3].comp(e)[:]
    [0, 0, 1]

Defining a local frame on a vector bundle, the dual coframe is automatically
created, which, by default, bares the same name (here `e`)::

    sage: E.coframes()
    [Local coframe (E|_U, (e^1,e^2,e^3))]
    sage: e_dual = E.coframes()[0] ; e_dual
    Local coframe (E|_U, (e^1,e^2,e^3))
    sage: e_dual is e.coframe()
    True

Let us check that the coframe `(e^i)` is indeed the dual of the vector
frame `(e_i)`::

    sage: e_dual[1](e[1]) # linear form e^1 applied to local section e_1
    Scalar field e^1(e_1) on the Open subset U of the 3-dimensional topological
     manifold M
    sage: e_dual[1](e[1]).expr() # the explicit expression of e^1(e_1)
    1
    sage: e_dual[1](e[1]).expr(), e_dual[1](e[2]).expr(), e_dual[1](e[3]).expr()
    (1, 0, 0)
    sage: e_dual[2](e[1]).expr(), e_dual[2](e[2]).expr(), e_dual[2](e[3]).expr()
    (0, 1, 0)
    sage: e_dual[3](e[1]).expr(), e_dual[3](e[2]).expr(), e_dual[3](e[3]).expr()
    (0, 0, 1)

Via bundle automorphisms, a new frame can be created from an existing one::

    sage: sec_module_U = E.section_module(domain=U)
    sage: change_frame = sec_module_U.automorphism()
    sage: change_frame[:] = [[0,1,0],[0,0,1],[1,0,0]]
    sage: f = e.new_frame(change_frame, 'f'); f
    Local frame (E|_U, (f_1,f_2,f_3))

A copy of this automorphism and its inverse is now part of the vector bundle's
frame changes::

    sage: E.change_of_frame(e, f)
    Automorphism of the Free module C^0(U;E) of sections on the Open subset U of
     the 3-dimensional topological manifold M with values in the real vector
     bundle E of rank 3
    sage: E.change_of_frame(e, f) == change_frame
    True
    sage: E.change_of_frame(f, e) == change_frame.inverse()
    True

Let us check the components of `f` with respect to the frame `e`::

    sage: f[1].comp(e)[:]
    [0, 0, 1]
    sage: f[2].comp(e)[:]
    [1, 0, 0]
    sage: f[3].comp(e)[:]
    [0, 1, 0]

"""

#******************************************************************************
#       Copyright (C) 2013-2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_basis import (FreeModuleBasis,
                                                   FreeModuleCoBasis)
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule

class LocalCoFrame(FreeModuleCoBasis):
    r"""
    Local coframe on a vector bundle.

    A *local coframe* on a vector bundle `E \to M` of class `C^k` is a
    local section `e^*: U \to E^n` of class `C^k` on some subset `U` of the base
    space `M`, such that `e^*(p)` is a basis of the fiber `E^*_p` of the dual
    bundle for any `p \in U`.

    INPUT:

    - ``frame`` -- the local frame dual to the coframe
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the linear forms constituting the coframe, or a tuple of
      strings, representing the individual symbols of the linear forms
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the linear forms constituting
      the coframe, or a tuple of strings, representing the individual LaTeX
      symbols of the linear forms; if ``None``, ``symbol`` is used in place of
      ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the linear
      forms  of the coframe; if ``None``, the indices will be generated as
      integers within the range declared on the coframe's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the linear forms of the coframe; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

    Local coframe on a topological vector bundle of rank 3::

        sage: M = Manifold(3, 'M', structure='top', start_index=1)
        sage: X.<x,y,z> = M.chart()
        sage: E = M.vector_bundle(3, 'E')
        sage: e = E.local_frame('e')
        sage: from sage.manifolds.local_frame import LocalCoFrame
        sage: f = LocalCoFrame(e, 'f'); f
        Local coframe (E|_M, (f^1,f^2,f^3))

    The local coframe can also be obtained by using the method
    :meth:`~sage.tensor.modules.free_module_basis.FreeModuleBasis.dual_basis` or
    :meth:`~sage.manifolds.local_frame.LocalFrame.coframe`::

        sage: e_dual = e.dual_basis(); e_dual
        Local coframe (E|_M, (e^1,e^2,e^3))
        sage: e_dual is e.coframe()
        True
        sage: e_dual is f
        False
        sage: e_dual[:] == f[:]
        True
        sage: f[1].display(e)
        f^1 = e^1

    The consisted linear forms can be obtained via the operator ``[]``::

        sage: f[1], f[2], f[3]
        (Linear form f^1 on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3,
         Linear form f^2 on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3,
         Linear form f^3 on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3)

    Checking that `f` is the dual of `e`::

        sage: f[1](e[1]).expr(), f[1](e[2]).expr(), f[1](e[3]).expr()
        (1, 0, 0)
        sage: f[2](e[1]).expr(), f[2](e[2]).expr(), f[2](e[3]).expr()
        (0, 1, 0)
        sage: f[3](e[1]).expr(), f[3](e[2]).expr(), f[3](e[3]).expr()
        (0, 0, 1)

    """
    def __init__(self, frame, symbol, latex_symbol=None, indices=None,
                 latex_indices=None):
        r"""
        Construct a local coframe, dual to a given local frame.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: from sage.manifolds.local_frame import LocalCoFrame
            sage: f = LocalCoFrame(e, 'f'); f
            Local coframe (E|_M, (f^0,f^1))
            sage: TestSuite(f).run()

        """
        self._domain = frame.domain()
        self._base_space = frame.base_space()
        self._vbundle = frame.vector_bundle()
        FreeModuleCoBasis.__init__(self, frame, symbol,
                                   latex_symbol=latex_symbol, indices=indices,
                                   latex_indices=latex_indices)
        # The coframe is added to the vector bundle's set of coframes
        self._vbundle._coframes.append(self)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: f = e.coframe()
            sage: f._repr_()
            'Local coframe (E|_M, (e^0,e^1))'
            sage: repr(f)  # indirect doctest
            'Local coframe (E|_M, (e^0,e^1))'
            sage: f  # indirect doctest
            Local coframe (E|_M, (e^0,e^1))

        """
        desc = "Local coframe " + self._name
        return desc

    def at(self, point):
        r"""
        Return the value of ``self`` at a given point on the base space, this
        value being a basis of the dual vector bundle at this point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` in the domain `U` of the coframe (denoted `f` hereafter)

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_basis.FreeModuleCoBasis`
          representing the basis `f(p)` of the vector space `E^*_p`,
          dual to the vector bundle fiber `E_p`

        EXAMPLES:

        Cobasis of a vector bundle fiber::

            sage: M = Manifold(2, 'M', structure='top', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: e_dual = e.coframe(); e_dual
            Local coframe (E|_M, (e^1,e^2))
            sage: p = M.point((-1,2), name='p')
            sage: e_dual_p = e_dual.at(p) ; e_dual_p
            Dual basis (e^1,e^2) on the Fiber of E at Point p on the
            2-dimensional topological manifold M
            sage: type(e_dual_p)
            <class 'sage.tensor.modules.free_module_basis.FreeModuleCoBasis'>
            sage: e_dual_p[1]
            Linear form e^1 on the Fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: e_dual_p[2]
            Linear form e^2 on the Fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: e_dual_p is e.at(p).dual_basis()
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
          symbols of the linear forms constituting the coframe, or a list/tuple
          of strings, representing the individual symbols of the linear forms
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the linear forms
          constituting the coframe, or a list/tuple of strings, representing the
          individual LaTeX symbols of the linear forms; if ``None``, ``symbol``
          is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the linear forms of the coframe; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the linear forms;
          if ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'up'``) determines the position
          of the indices labelling the linear forms of the coframe; can be
          either ``'down'`` or ``'up'``
        - ``include_domain`` -- (default: ``True``) boolean determining whether
          the name of the domain is included in the beginning of the coframe
          name

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e').coframe(); e
            Local coframe (E|_M, (e^0,e^1))
            sage: e.set_name('f'); e
            Local coframe (E|_M, (f^0,f^1))
            sage: e.set_name('e', latex_symbol=r'\epsilon')
            sage: latex(e)
            \left(E|_{M}, \left(\epsilon^{0},\epsilon^{1}\right)\right)
            sage: e.set_name('e', include_domain=False); e
            Local coframe (e^0,e^1)
            sage: e.set_name(['a', 'b'], latex_symbol=[r'\alpha', r'\beta']); e
            Local coframe (E|_M, (a,b))
            sage: latex(e)
            \left(E|_{M}, \left(\alpha,\beta\right)\right)
            sage: e.set_name('e', indices=['x','y'],
            ....:            latex_indices=[r'\xi', r'\zeta']); e
            Local coframe (E|_M, (e^x,e^y))
            sage: latex(e)
            \left(E|_{M}, \left(e^{\xi},e^{\zeta}\right)\right)

        """
        super(LocalCoFrame, self).set_name(symbol, latex_symbol=latex_symbol,
                                      indices=indices,
                                      latex_indices=latex_indices,
                                      index_position=index_position)
        if include_domain:
            # Redefinition of the name and the LaTeX name to include the domain
            self._name = "({}|_{}, {})".format(self._vbundle._name,
                                                self._domain._name, self._name)
            self._latex_name = r"\left({}|_{{{}}}, {}\right)".format(
                                                    self._vbundle._latex_name,
                                                    self._domain._latex_name,
                                                    self._latex_name)

#******************************************************************************

class LocalFrame(FreeModuleBasis):
    r"""
    Local frame on a vector bundle.

    A *local frame* on a vector bundle `E \to M` of class `C^k` is a local
    section `(e_1,\dots,e_n):U \to E^n` of class `C^k` defined on some subset `U`
    of the base space `M`, such that `e(p)` is a basis of the fiber `E_p` for
    any `p \in U`.

    For each instantiation of a local frame, a local coframe is automatically
    created, as an instance of the class :class:`LocalCoFrame`. It is returned
    by the method :meth:`coframe`.

    INPUT:

    - ``section_module`` -- free module of local sections over `U` in the given
      vector bundle `E \to M`
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the local sections constituting the local frame, or a tuple
      of strings, representing the individual symbols of the local sections
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the local sections constituting
      the local frame, or a tuple of strings, representing the individual
      LaTeX symbols of the local sections; if ``None``, ``symbol`` is used in
      place of ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the local
      sections of the frame; if ``None``, the indices will be generated as
      integers within the range declared on the local frame's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the  local sections; if
      ``None``, ``indices`` is used instead
    - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
      dual coframe; if ``None``, ``symbol`` must be a string and is used
      for the common base of the symbols of the elements of the dual coframe
    - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
      but for the dual coframe

    EXAMPLES:

    Defining a local frame on a 3-dimensional vector bundle over a 3-dimensional
    manifold::

        sage: M = Manifold(3, 'M', start_index=1, structure='top')
        sage: E = M.vector_bundle(3, 'E')
        sage: e = E.local_frame('e'); e
        Local frame (E|_M, (e_1,e_2,e_3))
        sage: latex(e)
        \left(E|_{M}, \left(e_{1},e_{2},e_{3}\right)\right)

    The individual elements of the vector frame are accessed via square
    brackets, with the possibility to invoke the slice operator '``:``' to
    get more than a single element::

        sage: e[2]
        Section e_2 on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3
        sage: e[1:3]
        (Section e_1 on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3,
         Section e_2 on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3)
        sage: e[:]
        (Section e_1 on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3,
         Section e_2 on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3,
         Section e_3 on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3)

    The LaTeX symbol can be specified::

        sage: eps = E.local_frame('eps', latex_symbol=r'\epsilon')
        sage: latex(eps)
        \left(E|_{M}, \left(\epsilon_{1},\epsilon_{2},\epsilon_{3}\right)\right)

    By default, the elements of the local frame are labelled by integers
    within the range specified at the manifold declaration. It is however
    possible to fully customize the labels, via the argument ``indices``::

        sage: u = E.local_frame('u', indices=('x', 'y', 'z')) ; u
        Local frame (E|_M, (u_x,u_y,u_z))
        sage: u[1]
        Section u_x on the 3-dimensional topological manifold M with values in
         the real vector bundle E of rank 3
        sage: u.coframe()
        Local coframe (E|_M, (u^x,u^y,u^z))

    The LaTeX format of the indices can be adjusted::

        sage: v = E.local_frame('v', indices=('a', 'b', 'c'),
        ....:                    latex_indices=(r'\alpha', r'\beta', r'\gamma'))
        sage: v
        Local frame (E|_M, (v_a,v_b,v_c))
        sage: latex(v)
        \left(E|_{M}, \left(v_{\alpha},v_{\beta},v_{\gamma}\right)\right)
        sage: latex(v.coframe())
        \left(E|_{M}, \left(v^{\alpha},v^{\beta},v^{\gamma}\right)\right)

    The symbol of each element of the local frame can also be freely chosen,
    by providing a tuple of symbols as the first argument of ``local_frame``;
    it is then mandatory to specify as well some symbols for the dual coframe::

        sage: h = E.local_frame(('a', 'b', 'c'), symbol_dual=('A', 'B', 'C')); h
        Local frame (E|_M, (a,b,c))
        sage: h[1]
        Section a on the 3-dimensional topological manifold M with values in the
         real vector bundle E of rank 3
        sage: h.coframe()
        Local coframe (E|_M, (A,B,C))
        sage: h.coframe()[1]
        Linear form A on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3

    Local frames are bases of free modules formed by local sections::

        sage: N = Manifold(2, 'N', structure='top', start_index=1)
        sage: X.<x,y> = N.chart()
        sage: U = N.open_subset('U')
        sage: F = N.vector_bundle(2, 'F')
        sage: f = F.local_frame('f', domain=U)
        sage: f.module()
        Free module C^0(U;F) of sections on the Open subset U of the
         2-dimensional topological manifold N with values in the real vector
         bundle F of rank 2
        sage: f.module().base_ring()
        Algebra of scalar fields on the Open subset U of the 2-dimensional
         topological manifold N
        sage: f.module() is F.section_module(domain=f.domain())
        True
        sage: f in F.section_module(domain=U).bases()
        True

    The value of the local frame at a given point is a basis of the
    corresponding fiber::

        sage: X_U = X.restrict(U) # We need coordinates on the subset
        sage: p = N((0,1), name='p') ; p
        Point p on the 2-dimensional topological manifold N
        sage: f.at(p)
        Basis (f_1,f_2) on the Fiber of F at Point p on the 2-dimensional
         topological manifold N

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = LocalCoFrame

    @staticmethod
    def __classcall_private__(cls, section_module, symbol,
                              latex_symbol=None, indices=None,
                              latex_indices=None, symbol_dual=None,
                              latex_symbol_dual=None):
        """
        Transform input lists into tuples for the unique representation of
        LocalFrame.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: from sage.manifolds.local_frame import LocalFrame
            sage: e = LocalFrame(C0, ['a', 'b'], symbol_dual=['A', 'B']); e
            Local frame (E|_M, (a,b))
            sage: e.dual_basis()
            Local coframe (E|_M, (A,B))
            sage: e is LocalFrame(C0, ('a', 'b'), symbol_dual=('A', 'B'))
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
        return super(LocalFrame, cls).__classcall__(cls, section_module,
                                        symbol, latex_symbol=latex_symbol,
                                        indices=indices,
                                        latex_indices=latex_indices,
                                        symbol_dual=symbol_dual,
                                        latex_symbol_dual=latex_symbol_dual)

    def __init__(self, section_module, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, symbol_dual=None, latex_symbol_dual=None):
        r"""
        Construct a local frame on a vector bundle.

        TESTS:

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: C0 = E.section_module(force_free=True)
            sage: from sage.manifolds.local_frame import LocalFrame
            sage: e = LocalFrame(C0, 'e', latex_symbol=r'\epsilon'); e
            Local frame (E|_M, (e_0,e_1))
            sage: TestSuite(e).run()

        """
        ###
        # Some sanity check:
        if not isinstance(section_module, FiniteRankFreeModule):
            raise ValueError("the {} has already been constructed as a "
                             "non-free module and therefore cannot have "
                             "a basis".format(section_module))
        self._domain = section_module.domain()
        self._base_space = section_module.base_space()
        self._vbundle = section_module.vector_bundle()
        FreeModuleBasis.__init__(self, section_module,
                                 symbol, latex_symbol=latex_symbol,
                                 indices=indices, latex_indices=latex_indices,
                                 symbol_dual=symbol_dual,
                                 latex_symbol_dual=latex_symbol_dual)
        if self._vbundle._def_frame is None:
            self._vbundle._def_frame = self
        # The frame is added to the domain's modules of frames, as well as to
        # all the superdomain's modules of frames; moreover the first defined
        # frame is considered as the default one
        for sd in self._domain._supersets:
            if sd in self._vbundle._section_modules:
                smodule = self._vbundle._section_modules[sd]
                if smodule.default_frame() is None:
                    smodule.set_default_frame(self)
                # Initialization of the zero element of the section module:
                if not isinstance(smodule, FiniteRankFreeModule):
                    smodule(0)._add_comp_unsafe(self)
                    # (since new components are initialized to zero)
        ###
        # Add this frame to the list of frames of the overlying vector bundle:
        self._vbundle._add_local_frame(self)

        self._coframe = self.dual_basis()  # Shortcut for self._dual_basis
        ###
        # Frame restrictions:
        self._subframes = set([self]) # Set of frames which are just a
                        # restriction of self
        self._superframes = set([self]) # Set of frames for which self is a
                        # restriction of
        self._restrictions = {} # Key: subdomain of self._domain; value:
                        # restriction of self on this subdomain

    ###### Methods that must be redefined by derived classes of ######
    ###### FreeModuleBasis                                      ######

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: repr(e)  # indirect doctest
            'Local frame (E|_M, (e_0,e_1))'
            sage: e  # indirect doctest
            Local frame (E|_M, (e_0,e_1))

        """
        desc = "Local frame " + self._name
        return desc

    def _new_instance(self, symbol, latex_symbol=None, indices=None,
                      latex_indices=None, symbol_dual=None,
                      latex_symbol_dual=None):
        r"""
        Construct a new local frame on the same section module as ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the sections constituting the local frame, or a
          tuple of strings, representing the individual symbols of the local
          sections
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the local sections
          constituting the local frame, or a tuple of strings, representing
          the individual LaTeX symbols of the local sections; if ``None``,
          ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the local sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on the local frame's
          domain
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the local sections;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT:

        - instance of :class:`LocalFrame`

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: e._new_instance('f')
            Local frame (E|_M, (f_0,f_1))

        """
        return LocalFrame(self._fmodule, symbol, latex_symbol=latex_symbol,
                           indices=indices, latex_indices=latex_indices,
                           symbol_dual=symbol_dual,
                           latex_symbol_dual=latex_symbol_dual)

    ###### End of methods to be redefined by derived classes ######

    def domain(self):
        r"""
        Return the domain on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U); e
            Local frame (E|_U, (e_0,e_1))
            sage: e.domain()
            Open subset U of the 3-dimensional topological manifold M

        """
        return self._domain

    def base_space(self):
        r"""
        Return the base space on which the overlying vector bundle is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U)
            sage: e.base_space()
            3-dimensional topological manifold M

        """
        return self._base_space

    def vector_bundle(self):
        r"""
        Return the vector bundle on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e', domain=U)
            sage: e.vector_bundle()
            Topological real vector bundle E -> M of rank 2 over the base space
            3-dimensional topological manifold M
            sage: e.vector_bundle() is E
            True

        """
        return self._vbundle

    def coframe(self):
        r"""
        Return the coframe of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: e.coframe()
            Local coframe (E|_M, (e^0,e^1))

        """
        return self._coframe

    def new_frame(self, change_of_frame, symbol, latex_symbol=None,
                  indices=None, latex_indices=None, symbol_dual=None,
                  latex_symbol_dual=None):
        r"""
        Define a new local frame from ``self``.

        The new local frame is defined from vector bundle automorphisms; its
        module is the same as that of the current frame.

        INPUT:

        - ``change_of_frame`` --
          :class:`~sage.tensor.modules.free_module_automorphism.FreeModuleAutomorphism`;
          vector bundle automorphisms `P` that relates
          the current frame `(e_i)` to the new frame `(f_i)` according
          to `f_i = P(e_i)`
        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the sections constituting the local frame, or a
          list/tuple of strings, representing the individual symbols of the
          sections
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the sections
          constituting the local frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the sections;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the sections;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT:

        - the new frame `(f_i)`, as an instance of :class:`LocalFrame`

        EXAMPLES:

        Orthogonal transformation of a frame on the 2-dimensional trivial vector
        bundle over the Euclidean plane::

            sage: M = Manifold(2, 'R^2', structure='top', start_index=1)
            sage: c_cart.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_R^2, (e_1,e_2))
            sage: orth = E.section_module().automorphism()
            sage: orth[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: f = e.new_frame(orth, 'f')
            sage: f[1][:]
            [1/2*sqrt(3), 1/2]
            sage: f[2][:]
            [-1/2, 1/2*sqrt(3)]
            sage: a =  E.change_of_frame(e,f)
            sage: a[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a == orth
            True
            sage: a is orth
            False
            sage: a._components # random (dictionary output)
            {Local frame (E|_D_0, (e_1,e_2)): 2-indices components w.r.t.
             Local frame (E|_D_0, (e_1,e_2)),
             Local frame (E|_D_0, (f_1,f_2)): 2-indices components w.r.t.
             Local frame (E|_D_0, (f_1,f_2))}
            sage: a.comp(f)[:]
            [1/2*sqrt(3)        -1/2]
             [        1/2 1/2*sqrt(3)]
            sage: a1 = E.change_of_frame(f,e)
            sage: a1[:]
            [1/2*sqrt(3)         1/2]
            [       -1/2 1/2*sqrt(3)]
            sage: a1 == orth.inverse()
            True
            sage: a1 is orth.inverse()
            False
            sage: e[1].comp(f)[:]
            [1/2*sqrt(3), -1/2]
            sage: e[2].comp(f)[:]
            [1/2, 1/2*sqrt(3)]

        """
        the_new_frame = self.new_basis(change_of_frame, symbol,
                                       latex_symbol=latex_symbol,
                                       indices=indices,
                                       latex_indices=latex_indices,
                                       symbol_dual=symbol_dual,
                                       latex_symbol_dual=latex_symbol_dual)
        self._vbundle._frame_changes[(self, the_new_frame)] = \
                            self._fmodule._basis_changes[(self, the_new_frame)]
        self._vbundle._frame_changes[(the_new_frame, self)] = \
                            self._fmodule._basis_changes[(the_new_frame, self)]
        return the_new_frame

    def restrict(self, subdomain):
        r"""
        Return the restriction of ``self`` to some open subset of its domain.

        If the restriction has not been defined yet, it is constructed here.

        INPUT:

        - ``subdomain`` -- open subset `V` of the current frame domain `U`

        OUTPUT:

        - the restriction of the current frame to `V` as a :class:`LocalFrame`

        EXAMPLES:

        Restriction of a frame defined on `\RR^2` to the unit disk::

            sage: M = Manifold(2, 'R^2', structure='top', start_index=1)
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_R^2, (e_1,e_2))
            sage: a = E.section_module().automorphism()
            sage: a[:] = [[1-y^2,0], [1+x^2, 2]]
            sage: f = e.new_frame(a, 'f'); f
            Local frame (E|_R^2, (f_1,f_2))
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1})
            sage: e_U = e.restrict(U); e_U
            Local frame (E|_U, (e_1,e_2))
            sage: f_U = f.restrict(U) ; f_U
            Local frame (E|_U, (f_1,f_2))

        The vectors of the restriction have the same symbols as those of the
        original frame::

            sage: f_U[1].display()
            f_1 = (-y^2 + 1) e_1 + (x^2 + 1) e_2
            sage: f_U[2].display()
            f_2 = 2 e_2

        Actually, the components are the restrictions of the original frame
        vectors::

            sage: f_U[1] is f[1].restrict(U)
            True
            sage: f_U[2] is f[2].restrict(U)
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
            resmodule = self._vbundle.section_module(domain=subdomain,
                                                      force_free=True)
            res = LocalFrame(resmodule,
                              self._symbol, latex_symbol=self._latex_symbol,
                              indices=self._indices,
                              latex_indices=self._latex_indices,
                              symbol_dual=self._symbol_dual,
                              latex_symbol_dual=self._latex_symbol_dual)

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

    def at(self, point):
        r"""
        Return the value of ``self`` at a given point, this value being
        a basis of the vector bundle fiber at the point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`; point
          `p` in the domain `U` of the local frame (denoted `e` hereafter)

        OUTPUT:

        - :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`
          representing the basis `e(p)` of the vector bundle fiber
          `E_p`

        EXAMPLES:

        Basis of a fiber of a trivial vector bundle::

            sage: M = Manifold(2, 'M', structure='top')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: p = M.point((-1,2), name='p')
            sage: ep = e.at(p) ; ep
            Basis (e_0,e_1) on the Fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: type(ep)
            <class 'sage.tensor.modules.free_module_basis.FreeModuleBasis'>
            sage: ep[0]
            Vector e_0 in the fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: ep[1]
            Vector e_1 in the fiber of E at Point p on the 2-dimensional
             topological manifold M

        Note that the symbols used to denote the vectors are same as those
        for the vector fields of the frame. At this stage, ``ep`` is the unique
        basis on fiber at ``p``::

            sage: Ep = E.fiber(p)
            sage: Ep.bases()
            [Basis (e_0,e_1) on the Fiber of E at Point p on the 2-dimensional
             topological manifold M]

        Let us consider another local frame::

            sage: aut = E.section_module().automorphism()
            sage: aut[:] = [[1+y^2, 0], [0, 2]]
            sage: f = e.new_frame(aut, 'f') ; f
            Local frame (E|_M, (f_0,f_1))
            sage: fp = f.at(p) ; fp
            Basis (f_0,f_1) on the Fiber of E at Point p on the 2-dimensional
             topological manifold M

        There are now two bases on the fiber::

            sage: Ep.bases()
            [Basis (e_0,e_1) on the Fiber of E at Point p on the 2-dimensional
             topological manifold M,
             Basis (f_0,f_1) on the Fiber of E at Point p on the 2-dimensional
             topological manifold M]

        Moreover, the changes of bases in the tangent space have been
        computed from the known relation between the frames ``e`` and
        ``f`` (via the automorphism ``aut`` defined above)::

            sage: Ep.change_of_basis(ep, fp)
            Automorphism of the Fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: Ep.change_of_basis(ep, fp).display()
            5 e_0⊗e^0 + 2 e_1⊗e^1
            sage: Ep.change_of_basis(fp, ep)
            Automorphism of the Fiber of E at Point p on the 2-dimensional
             topological manifold M
            sage: Ep.change_of_basis(fp, ep).display()
            1/5 e_0⊗e^0 + 1/2 e_1⊗e^1

        The dual bases::

            sage: e.coframe()
            Local coframe (E|_M, (e^0,e^1))
            sage: ep.dual_basis()
            Dual basis (e^0,e^1) on the Fiber of E at Point p on the
             2-dimensional topological manifold M
            sage: ep.dual_basis() is e.coframe().at(p)
            True
            sage: f.coframe()
            Local coframe (E|_M, (f^0,f^1))
            sage: fp.dual_basis()
            Dual basis (f^0,f^1) on the Fiber of E at Point p on the
             2-dimensional topological manifold M
            sage: fp.dual_basis() is f.coframe().at(p)
            True

        """
        # Determination of the vector bundle fiber:
        if point not in self._domain:
            raise ValueError("the {} is not a point in the ".format(point) +
                             "domain of {}".format(self))
        vbf = self._vbundle.fiber(point)
        # If the basis has already been constructed, it is simply returned:
        vbf_frame_bases = vbf._frame_bases
        if self in vbf_frame_bases:
            return vbf_frame_bases[self]
        for frame in vbf_frame_bases:
            if self in frame._subframes or self in frame._superframes:
                return vbf_frame_bases[frame]
        # If this point is reached, the basis has to be constructed from
        # scratch.
        # The names of the basis vectors set to those of the frame sections:
        basis = vbf.basis(self._symbol, latex_symbol=self._latex_symbol,
                         indices=self._indices,
                         latex_indices=self._latex_indices,
                         symbol_dual=self._symbol_dual,
                         latex_symbol_dual=self._latex_symbol_dual)
        vbf_frame_bases[self] = basis
        # Update of the change of bases in the fiber:
        for frame_pair, automorph in self._vbundle._frame_changes.items():
            frame1 = frame_pair[0]
            frame2 = frame_pair[1]
            if frame1 is self:
                fr2 = None
                for frame in vbf_frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr2 is not None:
                    basis1 = basis
                    basis2 = vbf_frame_bases[fr2]
                    auto = vbf.automorphism()
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
                    vbf._basis_changes[(basis1, basis2)] = auto
            if frame2 is self:
                fr1 = None
                for frame in vbf_frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                if fr1 is not None:
                    basis1 = vbf_frame_bases[fr1]
                    basis2 = basis
                    auto = vbf.automorphism()
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
                    vbf._basis_changes[(basis1, basis2)] = auto
        return basis

    def set_name(self, symbol, latex_symbol=None, indices=None,
                 latex_indices=None, index_position='down',
                 include_domain=True):
        r"""
        Set (or change) the text name and LaTeX name of ``self``.

        INPUT:

        - ``symbol`` -- either a string, to be used as a common base for the
          symbols of the local sections constituting the local frame, or a
          list/tuple of strings, representing the individual symbols of the
          local sections
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the local sections
          constituting the local frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the local sections;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the local sections of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the local sections;
          if ``None``, ``indices`` is used instead
        - ``index_position`` -- (default: ``'down'``) determines the position
          of the indices labelling the local sections of the frame; can be
          either ``'down'`` or ``'up'``
        - ``include_domain`` -- (default: ``True``) boolean determining whether
          the name of the domain is included in the beginning of the vector
          frame name

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e'); e
            Local frame (E|_M, (e_0,e_1))
            sage: e.set_name('f'); e
            Local frame (E|_M, (f_0,f_1))
            sage: e.set_name('e', include_domain=False); e
            Local frame (e_0,e_1)
            sage: e.set_name(['a', 'b']); e
            Local frame (E|_M, (a,b))
            sage: e.set_name('e', indices=['x', 'y']); e
            Local frame (E|_M, (e_x,e_y))
            sage: e.set_name('e', latex_symbol=r'\epsilon')
            sage: latex(e)
            \left(E|_{M}, \left(\epsilon_{0},\epsilon_{1}\right)\right)
            sage: e.set_name('e', latex_symbol=[r'\alpha', r'\beta'])
            sage: latex(e)
            \left(E|_{M}, \left(\alpha,\beta\right)\right)
            sage: e.set_name('e', latex_symbol='E',
            ....:            latex_indices=[r'\alpha', r'\beta'])
            sage: latex(e)
            \left(E|_{M}, \left(E_{\alpha},E_{\beta}\right)\right)

        """
        super(LocalFrame, self).set_name(symbol, latex_symbol=latex_symbol,
                                          indices=indices,
                                          latex_indices=latex_indices,
                                          index_position=index_position)
        if include_domain:
            # Redefinition of the name and the LaTeX name to include the domain
            self._name = "({}|_{}, {})".format(self._vbundle._name,
                                                self._domain._name, self._name)
            self._latex_name = r"\left({}|_{{{}}}, {}\right)".format(
                                                    self._vbundle._latex_name,
                                                    self._domain._latex_name,
                                                    self._latex_name)

#******************************************************************************

class TrivializationCoFrame(LocalCoFrame):
    r"""
    Trivialization coframe on a vector bundle.

    A *trivialization coframe* is the coframe of the trivialization frame
    induced by a trivialization (see: :class:`~sage.manifolds.local_frame.TrivializationFrame`).

    More precisely, a *trivialization frame* on a vector bundle `E \to M` of
    class `C^k` and rank `n` over the topological field `K` and over a
    topological manifold `M` is a local coframe induced by a local
    trivialization `\varphi:E|_U \to U \times K^n` of the domain `U \in M`.
    Namely, the local dual sections

    .. MATH::

        \varphi^*e^i := \varphi(\;\cdot\;, e^i)

    on `U` induce a local frame `(\varphi^*e^1, \dots, \varphi^*e^n)`, where
    `(e^1, \dots, e^n)` is the dual of the standard basis of `K^n`.

    INPUT:

    - ``triv_frame`` -- trivialization frame dual to the trivialization coframe
    - ``symbol`` -- either a string, to be used as a common base for the
      symbols of the dual sections constituting the coframe, or a tuple of
      strings, representing the individual symbols of the dual sections
    - ``latex_symbol`` -- (default: ``None``) either a string, to be used
      as a common base for the LaTeX symbols of the dual sections constituting
      the coframe, or a tuple of strings, representing the individual LaTeX
      symbols of the dual sections; if ``None``, ``symbol`` is used in place of
      ``latex_symbol``
    - ``indices`` -- (default: ``None``; used only if ``symbol`` is a single
      string) tuple of strings representing the indices labelling the dual
      sections of the coframe; if ``None``, the indices will be generated as
      integers within the range declared on the local frame's domain
    - ``latex_indices`` -- (default: ``None``) tuple of strings representing
      the indices for the LaTeX symbols of the dual sections of the coframe; if
      ``None``, ``indices`` is used instead

    EXAMPLES:

    Trivialization coframe on a trivial vector bundle of rank 3::

        sage: M = Manifold(3, 'M', start_index=1, structure='top')
        sage: X.<x,y,z> = M.chart()
        sage: E = M.vector_bundle(3, 'E')
        sage: phi = E.trivialization('phi'); phi
        Trivialization (phi, E|_M)
        sage: E.frames()
        [Trivialization frame (E|_M, ((phi^*e_1),(phi^*e_2),(phi^*e_3)))]
        sage: E.coframes()
        [Trivialization coframe (E|_M, ((phi^*e^1),(phi^*e^2),(phi^*e^3)))]
        sage: f = E.coframes()[0] ; f
        Trivialization coframe (E|_M, ((phi^*e^1),(phi^*e^2),(phi^*e^3)))

    The linear forms composing the coframe are obtained via the operator
    ``[]``::

        sage: f[1]
        Linear form (phi^*e^1) on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3
        sage: f[2]
        Linear form (phi^*e^2) on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3
        sage: f[3]
        Linear form (phi^*e^3) on the Free module C^0(M;E) of sections on the
         3-dimensional topological manifold M with values in the real vector
         bundle E of rank 3
        sage: f[1][:]
        [1, 0, 0]
        sage: f[2][:]
        [0, 1, 0]
        sage: f[3][:]
        [0, 0, 1]

    The coframe is the dual of the trivialization frame::

        sage: e = phi.frame() ; e
        Trivialization frame (E|_M, ((phi^*e_1),(phi^*e_2),(phi^*e_3)))
        sage: f[1](e[1]).expr(), f[1](e[2]).expr(), f[1](e[3]).expr()
        (1, 0, 0)
        sage: f[2](e[1]).expr(), f[2](e[2]).expr(), f[2](e[3]).expr()
        (0, 1, 0)
        sage: f[3](e[1]).expr(), f[3](e[2]).expr(), f[3](e[3]).expr()
        (0, 0, 1)

    """
    def __init__(self, triv_frame, symbol, latex_symbol=None,
                 indices=None, latex_indices=None):
        r"""
        Construct a local coframe from a local trivialization.

        TESTS::

            sage: M = Manifold(2, 'M', structure='top')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: from sage.manifolds.local_frame import TrivializationCoFrame
            sage: f = TrivializationCoFrame(phi.frame(), 'omega'); f
            Trivialization coframe (E|_M, (omega^0,omega^1))
            sage: TestSuite(f).run()

        """
        if not isinstance(triv_frame, TrivializationFrame):
            raise TypeError("the first argument must be a local trivialization "
                            "frame")
        LocalCoFrame.__init__(self, triv_frame, symbol,
                              latex_symbol=latex_symbol, indices=indices,
                              latex_indices=latex_indices)
        self._trivialization = triv_frame._trivialization

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: e = phi.frame().coframe()
            sage: e._repr_()
            'Trivialization coframe (E|_M, ((phi^*e^1),(phi^*e^2)))'
            sage: repr(e)  # indirect doctest
            'Trivialization coframe (E|_M, ((phi^*e^1),(phi^*e^2)))'
            sage: e  # indirect doctest
            Trivialization coframe (E|_M, ((phi^*e^1),(phi^*e^2)))

        """
        return "Trivialization coframe " + self._name

#******************************************************************************

class TrivializationFrame(LocalFrame):
    r"""
    Trivialization frame on a topological vector bundle.

    A *trivialization frame* on a topological vector bundle `E \to M` of rank
    `n` over the topological field `K` and over a topological manifold `M` is a
    local frame induced by a local trivialization `\varphi:E|_U \to U \times K^n`
    of the domain `U \in M`. More precisely, the local sections

    .. MATH::

        \varphi^*e_i := \varphi(\;\cdot\;, e_i)

    on `U` induce a local frame `(\varphi^*e_1, \dots, \varphi^*e_n)`, where
    `(e_1, \dots, e_n)` is the standard basis of `K^n`.

    INPUT:

    - ``trivialization`` -- the trivialization defined on the vector bundle

    EXAMPLES::

        sage: M = Manifold(3, 'M')
        sage: U = M.open_subset('U')
        sage: E = M.vector_bundle(2, 'E')
        sage: phi_U = E.trivialization('phi_U', domain=U)
        sage: phi_U.frame()
        Trivialization frame (E|_U, ((phi_U^*e_1),(phi_U^*e_2)))
        sage: latex(phi_U.frame())
        \left(E|_{U}, \left(\left(phi_U^* e_{ 1 }\right),\left(phi_U^* e_{ 2 }\right)\right)\right)

    """

    # The following class attribute must be redefined by any derived class:
    _cobasis_class = TrivializationCoFrame

    def __init__(self, trivialization):
        r"""
        Construct a trivialization frame.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: e = phi.frame()
            sage: TestSuite(e).run()

        """
        from sage.misc.latex import latex
        from .trivialization import Trivialization
        if not isinstance(trivialization, Trivialization):
            raise TypeError("the first argument must be a trivialization")
        ###
        # Some useful variables:
        triv = trivialization
        domain = triv.domain()
        vbundle = triv.vector_bundle()
        ###
        # Some sanity check:
        smodule = vbundle._section_modules.get(domain)
        if smodule and not isinstance(smodule, FiniteRankFreeModule):
            raise ValueError("the {} has already been constructed as a "
                             "non-free module and therefore cannot have "
                             "a basis".format(smodule))
        ###
        # Set trivialization:
        self._trivialization = triv
        ###
        # Define trivialization names
        rank = vbundle.rank()
        symbol = tuple("(" + triv._name + "^*" + "e_" + str(i) + ")"
            for i in range(1, rank + 1))
        symbol_dual = tuple("(" + triv._name + "^*" + "e^" + str(i) + ")"
            for i in range(1, rank + 1))
        latex_symbol = tuple(r'\left(' + triv._latex_name + r'^* e_{' +
                             latex(i) + r'}\right)'
                             for i in range(1, rank + 1))
        latex_symbol_dual = tuple(r'\left(' + triv._latex_name + r'^* e^{' +
                             latex(i) + r'}\right)'
                             for i in range(1, rank + 1))
        LocalFrame.__init__(self,
                        vbundle.section_module(domain=domain, force_free=True),
                        symbol=symbol, latex_symbol=latex_symbol,
                        symbol_dual=symbol_dual,
                        latex_symbol_dual=latex_symbol_dual)

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi = E.trivialization('phi')
            sage: e = phi.frame()
            sage: e._repr_()
            'Trivialization frame (E|_M, ((phi^*e_1),(phi^*e_2)))'
            sage: repr(e)  # indirect doctest
            'Trivialization frame (E|_M, ((phi^*e_1),(phi^*e_2)))'
            sage: e  # indirect doctest
            Trivialization frame (E|_M, ((phi^*e_1),(phi^*e_2)))

        """
        return "Trivialization frame " + self._name

    def trivialization(self):
        r"""
        Return the underlying trivialization of ``self``.

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: E = M.vector_bundle(2, 'E')
            sage: phi_U = E.trivialization('phi_U', domain=U)
            sage: e = phi_U.frame()
            sage: e.trivialization()
            Trivialization (phi_U, E|_U)
            sage: e.trivialization() is phi_U
            True

        """
        return self._trivialization
