# sage_setup: distribution = sagemath-symbolics
r"""
Bundle Connections

Let `E \to M` be a smooth vector bundle of rank `n` over a smooth manifold `M`
and over a non-discrete topological field `K` (typically `K=\RR` or `K=\CC`). A
*bundle connection* on this vector bundle is a `K`-linear map

.. MATH::

    \nabla : C^\infty(M;E) \to C^\infty(M;E \otimes T^*M)

such that the Leibniz rule applies for each scalar field `f \in C^\infty(M)` and
section `s \in C^\infty(M;E)`:

.. MATH::

    \nabla(f \, s) = f \cdot \nabla s + s \otimes \mathrm{d}f .

If `e` is a local frame on `E`, we have

.. MATH::

    \nabla e_i = \sum^n_{j=1} e_j \otimes \omega^j_i ,

and the corresponding `n \times n`-matrix `(\omega^j_i)_{i,j}` consisting of
one forms is called *connection matrix of* `\nabla` *with respect to* `e`.

AUTHORS:

- Michael Jung (2019) : initial version
"""
# ******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.mutability import Mutability
from sage.rings.integer import Integer
from sage.manifolds.differentiable.vector_bundle import \
    DifferentiableVectorBundle


class BundleConnection(SageObject, Mutability):
    r"""
    An instance of this class represents a bundle connection `\nabla` on a
    smooth vector bundle `E \to M`.

    INPUT:

    - ``vbundle`` -- the vector bundle on which the connection is defined
      (must be an instance of class
      :class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`)
    - ``name`` -- name given to the bundle connection
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the bundle
      connection; if ``None``, it is set to ``name``

    EXAMPLES:

    Define a bundle connection on a rank 2 vector bundle over some
    3-dimensional smooth manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x,y,z> = M.chart()
        sage: E = M.vector_bundle(2, 'E')
        sage: e = E.local_frame('e') # standard frame for E
        sage: nab = E.bundle_connection('nabla'); nab
        Bundle connection nabla on the Differentiable real vector bundle E -> M
         of rank 2 over the base space 3-dimensional differentiable manifold M

    First, let us initialize all connection 1-forms w.r.t. the frame ``e`` to
    zero::

        sage: nab[e, :] = [[0, 0], [0, 0]]

    This line can be shortened by the following::

        sage: nab[e, :] = 0  # initialize to zero

    The connection 1-forms are now initialized being differential 1-forms::

        sage: nab[e, 1, 1].parent()
        Free module Omega^1(M) of 1-forms on the 3-dimensional differentiable
         manifold M
        sage: nab[e, 1, 1].display()
        connection (1,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0

    Now, we want to specify some nonzero entries::

        sage: nab[e, 1, 2][:] = [x*z, y*z, z^2]
        sage: nab[e, 2, 1][:] = [x, x^2, x^3]
        sage: nab[e, 1, 1][:] = [x+z, y-z, x*y*z]
        sage: nab.display()
        connection (1,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = (x + z) dx + (y - z) dy + x*y*z dz
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x*z dx + y*z dy + z^2 dz
        connection (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x dx + x^2 dy + x^3 dz

    Notice, when we omit the frame, the default frame of the vector bundle is
    assumed (in this case ``e``)::

        sage: nab[2, 2].display()
        connection (2,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0

    The same holds for the assignment::

        sage: nab[1, 2] = 0
        sage: nab[e, 1, 2].display()
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0

    Keep noticed that item assignments for bundle connections only copy the
    right-hand-side and never create a binding to the original instance::

        sage: omega = M.one_form('omega')
        sage: omega[:] = [x*z, y*z, z^2]
        sage: nab[1, 2] = omega
        sage: nab[1, 2] == omega
        True
        sage: nab[1, 2] is omega
        False

    Hence, this is therefore equivalent to::

        sage: nab[2, 2].copy_from(omega)

    Preferably, we use :meth:`set_connection_form` to specify the connection
    1-forms::

        sage: nab[:] = 0  # re-initialize to zero
        sage: nab.set_connection_form(1, 2)[:] = [x*z, y*z, z^2]
        sage: nab.set_connection_form(2, 1)[:] = [x, x^2, x^3]
        sage: nab[1, 2].display()
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x*z dx + y*z dy + z^2 dz
        sage: nab[2, 1].display()
        connection (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x dx + x^2 dy + x^3 dz

    .. NOTE::

        Notice that item assignments and :meth:`set_connection_form` delete
        the connection 1-forms w.r.t. other frames for consistency reasons. To
        avoid this behavior, :meth:`add_connection_form` must be used instead.

    In conclusion, the connection 1-forms of a bundle connection are mutable
    until the connection itself is set immutable::

        sage: nab.set_immutable()
        sage: nab[1, 2] = omega
        Traceback (most recent call last):
        ...
        ValueError: object is immutable; please change a copy instead

    By definition, a bundle connection acts on vector fields and sections::

        sage: v = M.vector_field((x^2,y^2,z^2), name='v'); v.display()
        v = x^2 ∂/∂x + y^2 ∂/∂y + z^2 ∂/∂z
        sage: s = E.section((x-y^2, -z), name='s'); s.display()
        s = (-y^2 + x) e_1 - z e_2
        sage: nab_vs = nab(v, s); nab_vs
        Section nabla_v(s) on the 3-dimensional differentiable manifold M with
         values in the real vector bundle E of rank 2
        sage: nab_vs.display()
        nabla_v(s) = (-x^3*z^3 - 2*y^3 + x^2 - (x^2*y^2 + x^3)*z) e_1 +
         (-(y^2 - x)*z^4 - (x^3*y^2 + y^5 - x^4 - x*y^3)*z - z^2) e_2

    The bundle connection action certainly obeys the defining formula for
    the connection 1-forms::

        sage: vframe = X.frame()
        sage: all(nab(vframe[k], e[i]) == sum(nab[e, i, j](vframe[k])*e[j]
        ....:                                 for j in E.irange())
        ....:     for i in E.irange() for k in M.irange())
        True

    The connection 1-forms are computed automatically for different frames::

        sage: f = E.local_frame('f', ((1+x^2)*e[1], e[1]-e[2]))
        sage: nab.display(frame=f)
        connection (1,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (f_1,f_2)) = ((x^3 + x)*z + 2*x)/(x^2 + 1) dx + y*z dy + z^2 dz
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (f_1,f_2)) = -(x^3 + x)*z dx - (x^2 + 1)*y*z dy -
          (x^2 + 1)*z^2 dz
        connection (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (f_1,f_2)) = (x*z - x)/(x^2 + 1) dx -
          (x^2 - y*z)/(x^2 + 1) dy - (x^3 - z^2)/(x^2 + 1) dz
        connection (2,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (f_1,f_2)) = -x*z dx - y*z dy - z^2 dz

    The new connection 1-forms obey the defining formula, too::

        sage: all(nab(vframe[k], f[i]) == sum(nab[f, i, j](vframe[k])*f[j]
        ....:                                 for j in E.irange())
        ....:     for i in E.irange() for k in M.irange())
        True

    After the connection has been specified, the curvature 2-forms can be
    derived::

        sage: Omega = nab.curvature_form
        sage: for i in E.irange():
        ....:     for j in E.irange():
        ....:         print(Omega(i ,j, e).display())
        curvature (1,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = -(x^3 - x*y)*z dx∧dy + (-x^4*z + x*z^2) dx∧dz +
          (-x^3*y*z + x^2*z^2) dy∧dz
        curvature (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = -x dx∧dz - y dy∧dz
        curvature (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 2*x dx∧dy + 3*x^2 dx∧dz
        curvature (2,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = (x^3 - x*y)*z dx∧dy + (x^4*z - x*z^2) dx∧dz +
          (x^3*y*z - x^2*z^2) dy∧dz

    The derived forms certainly obey the structure equations, see
    :meth:`curvature_form` for details::

        sage: omega = nab.connection_form
        sage: check = []
        sage: for i in E.irange():  # long time
        ....:     for j in E.irange():
        ....:         check.append(Omega(i,j,e) == \
        ....:                       omega(i,j,e).exterior_derivative() + \
        ....:         sum(omega(k,j,e).wedge(omega(i,k,e))
        ....:             for k in E.irange()))
        sage: check  # long time
        [True, True, True, True]
    """

    def __init__(self, vbundle, name, latex_name=None):
        r"""
        Construct a bundle connection.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: from sage.manifolds.differentiable.bundle_connection \
            ....:                                       import BundleConnection
            sage: nab = BundleConnection(E, 'nabla', latex_name=r'\nabla')
            sage: nab
            Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 2 over the base space 3-dimensional differentiable
             manifold M
            sage: X.<x,y,z> = M.chart()
            sage: e = E.local_frame('e')
            sage: nab[:] = 0
            sage: nab.set_connection_form(1, 0)[:] = [x*z, y*z, z^2]
            sage: TestSuite(nab).run()
        """
        if not isinstance(vbundle, DifferentiableVectorBundle):
            raise TypeError("the first argument must be a differentiable " +
                            "vector bundle")
        Mutability.__init__(self)
        self._vbundle = vbundle
        self._domain = vbundle.base_space()
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._connection_forms = {}  # dict. of con. forms, with frames as keys
        self._coefficients = self._connection_forms
        # Initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: E = M.vector_bundle(3, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._repr_()
            'Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 3 over the base space 5-dimensional differentiable
             manifold M'
            sage: repr(nab)  # indirect doctest
            'Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 3 over the base space 5-dimensional differentiable
             manifold M'
        """
        description = "Bundle connection"
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._vbundle)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: E = M.vector_bundle(3, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._latex_()
            '\\nabla'
            sage: latex(nab)  # indirect doctest
            \nabla
            sage: nab = E.bundle_connection('D')
            sage: nab._latex_()
            'D'
            sage: latex(nab)  # indirect doctest
            D
        """
        return self._latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(4, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._init_derived()
        """
        self._curvature_forms = {}  # dict. of dict. of curvature forms
        # (key: local frame)
        self._hash = -1

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(4, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._del_derived()
        """
        self._curvature_forms.clear()

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a bundle connection

        OUTPUT: ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')  # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab[:] = 0
            sage: nab[0, 1][:] = [x^2, x]
            sage: nab[1, 0][:] = [y^2, y]
            sage: nab1 = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab1[:] = 0
            sage: (nab1 == nab) or (nab == nab1)
            False
            sage: nab1[0, 1][:] = [x, x^2]
            sage: nab1[1, 0][:] = [y, y^2]
            sage: (nab1 == nab) or (nab == nab1)
            False
            sage: nab1[0, 1][:] = [x^2, x]
            sage: nab1[1, 0][:] = [y^2, y]
            sage: (nab1 == nab) and (nab == nab1)
            True
        """
        if other is self:
            return True
        if not isinstance(other, BundleConnection):
            return False
        if other._domain != self._domain:
            return False
        if self._connection_forms == {}:
            return False
        for frame in self._connection_forms:
            if frame not in other._connection_forms:
                return False
            for ind in self._connection_forms[frame]:
                if (other._connection_forms[frame][ind] !=
                        self._connection_forms[frame][ind]):
                    return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- an affine connection

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')  # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab[:] = 0
            sage: nab[0, 1][:] = [x^2, x]
            sage: nab[1, 0][:] = [y^2, y]
            sage: nab1 = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab1[:] = 0
            sage: (nab1 != nab) and (nab != nab1)
            True
            sage: nab1[0, 1][:] = [x, x^2]
            sage: nab1[1, 0][:] = [y, y^2]
            sage: (nab1 != nab) and (nab != nab1)
            True
            sage: nab1[0, 1][:] = [x^2, x]
            sage: nab1[1, 0][:] = [y^2, y]
            sage: (nab1 != nab) or (nab != nab1)
            False
        """
        return not (self == other)

    def vector_bundle(self):
        r"""
        Return the vector bundle on which the bundle connection is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`
          representing the vector bundle on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', r'\nabla')
            sage: nab.vector_bundle()
            Differentiable real vector bundle E -> M of rank 2 over the base
             space 3-dimensional differentiable manifold M
        """
        return self._vbundle

    def _new_forms(self, frame):
        r"""
        Create the connection forms w.r.t. the given frame.

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: forms = nab._new_forms(e)
            sage: [forms[k] for k in sorted(forms)]
            [1-form connection (1,1) of bundle connection nabla w.r.t. Local
              frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
              manifold M,
             1-form connection (1,2) of bundle connection nabla w.r.t. Local
              frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
              manifold M,
             1-form connection (2,1) of bundle connection nabla w.r.t. Local
              frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
              manifold M,
             1-form connection (2,2) of bundle connection nabla w.r.t. Local
              frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
              manifold M]
        """
        dom = frame._domain
        forms_dict = {}
        for i in self._vbundle.irange():
            for j in self._vbundle.irange():
                # set names:
                name = "connection ({},{}) of bundle ".format(i, j)
                name += "connection " + self._name + " w.r.t. {}".format(frame)
                latex_name = r"\omega^" + str(j) + r"_{\ \, " + str(i) + "}"
                form = dom.diff_form(1, name=name, latex_name=latex_name)
                forms_dict[(i, j)] = form
        return forms_dict

    def connection_forms(self, frame=None):
        r"""
        Return the connection forms relative to the given frame.

        If `e` is a local frame on `E`, we have

        .. MATH::

            \nabla e_i = \sum^n_{j=1} e_j \otimes \omega^j_i ,

        and the corresponding `n \times n`-matrix `(\omega^j_i)_{i,j}`
        consisting of one forms is called *connection matrix of* `\nabla` *with
        respect to* `e`.

        If the connection coefficients are not known already, they are computed
        from the above formula.

        INPUT:

        - ``frame`` -- (default: ``None``) local frame relative to which the
          connection forms are required; if none is provided, the
          vector bundle's default frame is assumed

        OUTPUT:

        - connection forms relative to the frame ``frame``, as a dictionary
          with tuples `(i, j)` as key and one forms as instances of
          :class:`~sage.manifolds.differentiable.diff_form` as value
          representing the matrix entries.

        EXAMPLES:

        Connection forms of a bundle connection on a rank 2 vector bundle
        over a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla', r'\nabla')
            sage: nab[:] = 0  # initialize curvature forms
            sage: forms = nab.connection_forms()
            sage: [forms[k] for k in sorted(forms)]
            [1-form connection (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 3-dimensional differentiable
             manifold M,
            1-form connection (1,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 3-dimensional differentiable
             manifold M,
            1-form connection (2,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 3-dimensional differentiable
             manifold M,
            1-form connection (2,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 3-dimensional differentiable
             manifold M]
        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._domain)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided")
        if frame not in self._connection_forms:
            # the connection forms must be computed
            #
            # Check whether frame is a subframe of a frame in which the
            # forms are already known:
            for oframe in self._connection_forms:
                if frame in oframe._subframes:
                    self._connection_forms[frame] = self._new_forms(frame)
                    comp_store = self._connection_forms[frame]
                    ocomp_store = self._connection_forms[oframe]
                    for ind, value in ocomp_store.items():
                        comp_store[ind] = value.restrict(frame._domain)
                    break
            else:
                # If not, the forms must be computed from scratch:
                vb = self._vbundle
                dom = frame._domain
                vframe = dom.default_frame()
                # it is important to use _new_forms instead of
                # self.set_connection_form:
                omega = self._new_forms(frame)
                for d in dom.irange():
                    for i in vb.irange():
                        sec_nab = self(vframe[d], frame[i])
                        for j in vb.irange():
                            omega[(i, j)][vframe, d] = sec_nab[[frame, j]]
                self._connection_forms[frame] = omega
        return self._connection_forms[frame]

    def connection_form(self, i, j, frame=None):
        r"""
        Return the connection 1-form corresponding to the given index and
        local frame.

        .. SEEALSO::

            Consult :meth:`connection_forms` for detailed information.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^j_i`
        - ``frame`` -- (default: ``None``) local frame relative to which the
          connection 1-forms are defined; if ``None``, the default frame of the
          vector bundle's corresponding section module is assumed.

        OUTPUT:

        - the 1-form `\omega^j_i`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab.set_connection_form(0, 1)[:] = [x^2, x]
            sage: nab.set_connection_form(1, 0)[:] = [y^2, y]
            sage: nab.connection_form(0, 1).display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy
            sage: nab.connection_form(1, 0).display()
            connection (1,0) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = y^2 dx + y dy
        """
        return self.connection_forms(frame)[(i, j)]

    def __call__(self, v, s):
        r"""
        Action of the connection on a vector field and local section.

        INPUT:

        - ``v`` -- a vector field `v` on the base space
        - ``s`` -- a local section `s`

        OUTPUT: local section `\nabla_v s`

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab[:] = 0
            sage: nab[1,2][1] = x*y
            sage: v = M.vector_field('v')
            sage: v[:] = [-y, x]
            sage: s = E.section('s')
            sage: s[:] = [y, -x]
            sage: nab.__call__(v, s)
            Section nabla_v(s) on the 2-dimensional differentiable manifold M
             with values in the real vector bundle E of rank 2
        """
        from sage.manifolds.section import TrivialSection
        from sage.tensor.modules.format_utilities import format_unop_latex
        if isinstance(s, TrivialSection):
            return self._derive_trivial(v, s)
        # Resulting section
        vb = self._vbundle
        dom = s.domain()
        if s._name is None or v._name:
            name_resu = None
        else:
            name_resu = self._name + '_' + v._name + '(' + s._name + ')'
        if s._latex_name is None or v._latex_name is None:
            latex_name_resu = None
        else:
            nab_v_latex = self._latex_name + '_{' + v._latex_name + '} '
            latex_name_resu = format_unop_latex(nab_v_latex, s._latex_name)
        resu = vb.section(domain=dom, name=name_resu,
                          latex_name=latex_name_resu)
        # gluing process
        for dom, rst in s._restrictions.items():
            # the computation is performed only if dom is not a subdomain
            # of another restriction:
            for odom in s._restrictions:
                if dom in odom._subsets and dom is not odom:
                    break
            else:
                # dom is not a subdomain and the computation is performed:
                resu._restrictions[rst._domain] = self(rst)
        return resu

    def _derive_trivial(self, v, s):
        r"""
        Action of the connection on a local section whose module is free.

        INPUT:

        - ``v`` -- a vector field `v` on the base space
        - ``s`` -- a local section `s` whose module is free

        OUTPUT: local section `\nabla_v s`

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab[:] = 0
            sage: nab[1,2][1] = x*y
            sage: v = M.vector_field('v')
            sage: v[:] = [-y, x]
            sage: s = E.section('s')
            sage: s[:] = [y, -x]
            sage: nab._derive_trivial(v, s)
            Section nabla_v(s) on the 2-dimensional differentiable manifold M
             with values in the real vector bundle E of rank 2
        """
        vb = self._vbundle
        dom = s.domain()
        # pick the first frame whose forms of self on dom are known
        for frame in self._connection_forms:
            if dom.is_subset(frame.domain()):
                frame = frame.restrict(dom)
                break
        else:
            raise ValueError("no local frame found for the computation")
        # Resulting section
        from sage.tensor.modules.format_utilities import format_unop_latex
        if s._name is None or v._name is None:
            name_resu = None
        else:
            name_resu = self._name + '_' + v._name + '(' + s._name + ')'
        if s._latex_name is None or v._latex_name is None:
            latex_name_resu = None
        else:
            nab_v_latex = self._latex_name + '_{' + v._latex_name + '} '
            latex_name_resu = format_unop_latex(nab_v_latex, s._latex_name)
        res = vb.section(domain=dom, name=name_resu,
                         latex_name=latex_name_resu)
        for j in vb.irange():
            ds_comp = s[[frame, j]].differential()
            res_comp = ds_comp(v)
            res_comp += sum(s[[frame, i]] * self[frame, i, j](v)
                            for i in vb.irange())
            res[frame, j] = res_comp
        return res

    def add_connection_form(self, i, j, frame=None):
        r"""
        Return the connection form `\omega^j_i` in a given frame for
        assignment.

        See method :meth:`connection_forms` for details about the definition of
        the connection forms.

        To delete the connection forms in other frames, use the method
        :meth:`set_connection_form` instead.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^j_i`
        - ``frame`` -- (default: ``None``) local frame in which the connection
          1-form is defined; if ``None``, the default frame of the vector
          bundle is assumed.

        .. WARNING::

            If the connection has already forms in other frames, it is the
            user's responsibility to make sure that the 1-forms to be added
            are consistent with them.

        OUTPUT:

        - connection 1-form `\omega^j_i` in the given frame, as an instance of
          the class :class:`~sage.manifolds.differentiable.diff_form.DiffForm`;
          if such connection 1-form did not exist previously, it is created.
          See method :meth:`connection_forms` for the storage convention of the
          connection 1-forms.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab.add_connection_form(0, 1, frame=e)[:] = [x^2, x]
            sage: nab[e, 0, 1].display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy

        Since ``e`` is the vector bundle's default local frame, its mention may
        be omitted::

            sage: nab.add_connection_form(1, 0)[:] = [y^2, y]
            sage: nab[1, 0].display()
            connection (1,0) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = y^2 dx + y dy

        Adding connection 1-forms w.r.t. to another local frame::

            sage: f = E.local_frame('f')
            sage: nab.add_connection_form(1, 1, frame=f)[:] = [x, y]
            sage: nab[f, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (f_0,f_1)) = x dx + y dy

        The forms w.r.t. the frame ``e`` have been kept::

            sage: nab[e, 0, 1].display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy

        To delete them, use the method :meth:`set_connection_form` instead.
        """
        self._require_mutable()
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._domain)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided")
        # Are the components already known?
        if frame not in self._connection_forms:
            if frame not in self._vbundle._frames:
                raise ValueError("the {} is not".format(frame) +
                                 " a frame on the {}".format(self._domain))
            self._connection_forms[frame] = self._new_forms(frame)
        self._del_derived()  # deletes the derived quantities
        return self._connection_forms[frame][(i, j)]

    def set_connection_form(self, i, j, frame=None):
        r"""
        Return the connection form `\omega^j_i` in a given frame for
        assignment.

        See method :meth:`connection_forms` for details about the definition of
        the connection forms.

        The connection forms with respect to other frames are deleted,
        in order to avoid any inconsistency. To keep them, use the method
        :meth:`add_connection_form` instead.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^j_i`
        - ``frame`` -- (default: ``None``) local frame in which the connection
          1-form is defined; if ``None``, the default frame of the vector
          bundle is assumed.

        OUTPUT:

        - connection 1-form `\omega^j_i` in the given frame, as an instance of
          the class :class:`~sage.manifolds.differentiable.diff_form.DiffForm`;
          if such connection 1-form did not exist previously, it is created.
          See method :meth:`connection_forms` for the storage convention of the
          connection 1-forms.

        EXAMPLES:

        Setting the connection forms of a bundle connection w.r.t. some local
        frame::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab.set_connection_form(0, 1)[:] = [x^2, x]
            sage: nab[0, 1].display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy

        Since ``e`` is the vector bundle's default local frame, its mention may
        be omitted::

            sage: nab.set_connection_form(1, 0)[:] = [y^2, y]
            sage: nab[1, 0].display()
            connection (1,0) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = y^2 dx + y dy

        Setting connection 1-forms w.r.t. to another local frame::

            sage: f = E.local_frame('f')
            sage: nab.set_connection_form(1, 1, frame=f)[:] = [x, y]
            sage: nab[f, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (f_0,f_1)) = x dx + y dy

        The forms w.r.t. the frame ``e`` have been deleted::

            sage: nab[e, 0, 1].display()
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components in
             the Local frame (E|_M, (f_0,f_1))

        To keep them, use the method :meth:`add_connection_form` instead.
        """
        self._require_mutable()
        omega = self.add_connection_form(i, j, frame=frame)
        self.del_other_forms(frame)
        return omega

    def del_other_forms(self, frame=None):
        r"""
        Delete all the connection forms but those corresponding to ``frame``.

        INPUT:

        - ``frame`` -- (default: ``None``) local frame, the connection forms
          w.r.t. which are to be kept; if ``None``, the default frame of the
          vector bundle is assumed.

        EXAMPLES:

        We first create two sets of connection forms::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: e = E.local_frame('e')
            sage: nab.set_connection_form(1, 1, frame=e)[:] = [x^2, x]
            sage: f = E.local_frame('f')
            sage: nab.add_connection_form(1, 1, frame=f)[:] = [y^2, y]
            sage: nab[e, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[f, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (f_1,f_2)) = y^2 dx + y dy

        Let us delete the connection forms w.r.t. all frames except for
        frame ``e``::

            sage: nab.del_other_forms(e)
            sage: nab[e, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy

        The connection forms w.r.t. frame ``e`` have indeed been
        deleted::

            sage: nab[f, :]
            Traceback (most recent call last):
            ...
            ValueError: no basis could be found for computing the components in
             the Local frame (E|_M, (e_1,e_2))
        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._domain)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided")
        if frame not in self._connection_forms:
            raise ValueError("the coefficients w.r.t. {}".format(frame) +
                             " have not been defined")
        to_be_deleted = []
        for other_frame in self._connection_forms:
            if other_frame != frame:
                to_be_deleted.append(other_frame)
        for other_frame in to_be_deleted:
            del self._connection_forms[other_frame]

    def curvature_form(self, i, j, frame=None):
        r"""
        Return the curvature 2-form corresponding to the given index and local
        frame.

        The *curvature 2-forms* with respect to the frame `e` are the 2-forms
        `\Omega^j_i` given by the formula

        .. MATH::

            \Omega^j_i = \mathrm{d} \omega^j_i + \sum^n_{k=1} \omega^j_k
                \wedge \omega^k_i

        INPUT:

        - ``i``, ``j`` -- indices identifying the 2-form `\Omega^j_i`
        - ``frame`` -- (default: ``None``) local frame relative to which the
          curvature 2-forms are defined; if ``None``, the default frame
          of the vector bundle is assumed.

        OUTPUT:

        - the 2-form `\Omega^j_i`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(1, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: e = E.local_frame('e')
            sage: nab.set_connection_form(1, 1)[:] = [x^2, x]
            sage: curv = nab.curvature_form(1, 1); curv
            2-form curvature (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1)) on the 2-dimensional differentiable manifold M
            sage: curv.display()
            curvature (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1)) = dx∧dy
        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._domain)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided")
        if frame not in self._curvature_forms:
            self._curvature_forms[frame] = {}
        if (i, j) not in self._curvature_forms[frame]:
            name = "curvature ({},{}) of bundle connection ".format(i, j) + \
                   self._name + " w.r.t. {}".format(frame)
            latex_name = r"\Omega^" + str(i) + r"_{\ \, " + \
                         str(j) + "}"
            omega = self.connection_form
            curv_form = omega(i, j, frame).exterior_derivative()
            curv_form += sum(omega(k, j, frame).wedge(omega(i, k, frame))
                             for k in self._vbundle.irange())
            curv_form.set_name(name=name, latex_name=latex_name)
            self._curvature_forms[frame][(i, j)] = curv_form
        return self._curvature_forms[frame][(i, j)]

    def __hash__(self):
        r"""
        Hash function.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: hash(nab)
            Traceback (most recent call last):
            ...
            ValueError: object is mutable; please make it immutable first
            sage: nab.set_immutable()
            sage: hash(nab) == nab.__hash__()
            True

        Let us check that ``nab`` can be used as a dictionary key::

            sage: {nab: 1}[nab]
            1
        """
        self._require_immutable()
        if self._hash == -1:
            self._hash = hash((type(self).__name__, self._vbundle))
        return self._hash

    def __getitem__(self, args):
        r"""
        Return a component of ``self`` with respect to some frame.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are returned

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla')
            sage: e = E.local_frame('e')
            sage: nab[:] = 0
            sage: nab[1, 2][:] = [x^2, x]
            sage: nab[1, 1][:] = [y, y]
            sage: nab[1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = y dx + y dy
            sage: nab[e, 1, 2].display()
            connection (1,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[e, :]
            [[1-form connection (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (1,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M],
             [1-form connection (2,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (2,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M]]
            sage: nab[:]
            [[1-form connection (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (1,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M],
             [1-form connection (2,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (2,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M]]
        """
        # extract frame from first index:
        vb = self._vbundle
        if isinstance(args, (int, Integer, slice)):
            smodule = vb.section_module(domain=self._domain)
            frame = smodule.default_frame()
        elif not isinstance(args[0], (int, Integer, slice)):
            frame = args[0]
            args = args[1:]
        else:
            smodule = vb.section_module(domain=self._domain)
            frame = smodule.default_frame()
        # indexing:
        if isinstance(args, slice):
            indices = args
        elif isinstance(args[0], slice):
            indices = args[0]
        else:
            indices = args
        if isinstance(indices, slice):
            if indices.start is None and indices.stop is None:
                return [[self.connection_form(i, j, frame=frame)
                         for j in vb.irange()] for i in vb.irange()]
            else:
                raise NotImplementedError("[start:stop] syntax not "
                                          "implemented")
        if len(indices) != 2:
            raise ValueError("index must be a pair of integers")
        (i, j) = indices
        if i not in vb.irange() or j not in vb.irange():
            raise IndexError("index out of range")
        form = self.connection_form(indices[0], indices[1], frame=frame)
        return form

    def __setitem__(self, args, value):
        r"""
        Set the components of ``self`` corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (usually a pair of integers); if ``[:]``
          is provided, all the components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla')
            sage: e = E.local_frame('e')
            sage: a = M.one_form([x^2, x], name='a')
            sage: a.display()
            a = x^2 dx + x dy
            sage: nab[:] = 0
            sage: nab[1, 1] = a
            sage: nab[1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[e, 2, 2] = a
            sage: nab[e, 2, 2].display()
            connection (2,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[:] = [[0, 0], [a, a]]
            sage: nab[e, 2, 1].display()
            connection (2,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[e, :] = [[a, a], [0, 0]]
            sage: nab[e, 1, 2].display()
            connection (1,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
        """
        # extract frame from first index:
        vb = self._vbundle
        if isinstance(args, (int, Integer, slice)):
            smodule = vb.section_module(domain=self._domain)
            frame = smodule.default_frame()
        elif not isinstance(args[0], (int, Integer, slice)):
            frame = args[0]
            args = args[1:]
        else:
            smodule = vb.section_module(domain=self._domain)
            frame = smodule.default_frame()
        # determine indices:
        if isinstance(args, slice):
            indices = args
        elif isinstance(args[0], slice):
            indices = args[0]
        else:
            indices = args
        # set values:
        if isinstance(indices, (list, tuple)):
            if len(indices) != 2:
                raise IndexError("a tuple of integers must be provided")
            (i, j) = indices
            if i not in vb.irange() or j not in vb.irange():
                raise IndexError("index out of range")
            omega = self.set_connection_form(i, j, frame=frame)
            dom = frame.domain()
            dmodule = dom.diff_form_module(1)
            try:
                form = dmodule(value)
            except (TypeError, ValueError):
                msg = "{} must be convertible ".format(value)
                msg += "into an element of {}".format(dmodule)
                raise TypeError(msg)
            omega.copy_from(form)  # copy all components from value
        if isinstance(indices, slice):
            if indices.start is None and indices.stop is None:
                # check types:
                if isinstance(value, (int, Integer)) and value == 0:
                    for i in vb.irange():
                        for j in vb.irange():
                            self[frame, i, j] = 0
                elif not isinstance(value, (list, tuple)):
                    raise TypeError("in case of [:] syntax, zero or a "
                                    "list/tuple as value should be provided")
                elif any(not isinstance(row, (list, tuple)) for row in value):
                    raise TypeError("in case of [:] syntax, the list/tuple "
                                    "of value must contain lists/tuples")
                else:
                    # check lengths:
                    rk = vb._rank
                    if len(value) != rk:
                        raise ValueError("value must have "
                                         "length {}".format(rk))
                    if any(len(row) != rk for row in value):
                        raise ValueError("lists in value must have length "
                                         "{}".format(rk))
                    # perform designation:
                    sind = vb._base_space._sindex
                    for i in vb.irange():
                        for j in vb.irange():
                            self[frame, i, j] = value[i - sind][j - sind]
            else:
                raise NotImplementedError("[start:stop] syntax not "
                                          "implemented")

    def display(self, frame=None, vector_frame=None, chart=None,
                only_nonzero=True):
        r"""
        Display all the connection 1-forms w.r.t. to a given local frame, one
        per line.

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        INPUT:

        - ``frame`` -- (default: ``None``) local frame of the vector bundle
          relative to which the connection 1-forms are defined; if ``None``,
          the default frame of the bundle is used
        - ``vector_frame`` -- (default: ``None``) vector frame of the manifold
          relative to which the connection 1-forms should be displayed; if
          ``None``, the default frame of the local frame's domain is used
        - ``chart`` -- (default: ``None``) chart specifying the coordinate
          expression of the connection 1-forms; if ``None``,
          the default chart of the domain of ``frame`` is used
        - ``only_nonzero`` -- boolean (default: ``True``); if ``True``, only
          nonzero connection coefficients are displayed

        EXAMPLES:

        Set connection 1-forms::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla'); nab
            Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 2 over the base space 3-dimensional differentiable
             manifold M
            sage: nab[:] = 0
            sage: nab[1, 1][:] = [x, y, z]
            sage: nab[2, 2][:] = [x^2, y^2, z^2]

        By default, only the nonzero connection coefficients are displayed::

            sage: nab.display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x dx + y dy + z dz
            connection (2,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + y^2 dy + z^2 dz
            sage: latex(nab.display())
            \begin{array}{lcl} \omega^1_{\ \, 1} = x \mathrm{d} x +
             y \mathrm{d} y + z \mathrm{d} z \\ \omega^2_{\ \, 2} = x^{2}
             \mathrm{d} x + y^{2} \mathrm{d} y + z^{2} \mathrm{d} z \end{array}

        By default, the displayed connection 1-forms are those w.r.t.
        the default frame of the vector bundle. The aforementioned is
        therefore equivalent to::

            sage: nab.display(frame=E.default_frame())
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x dx + y dy + z dz
            connection (2,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + y^2 dy + z^2 dz

        Moreover, the connection 1-forms are displayed w.r.t. the default
        vector frame on the local frame's domain, i.e.::

            sage: domain = e.domain()
            sage: nab.display(vector_frame=domain.default_frame())
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x dx + y dy + z dz
            connection (2,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + y^2 dy + z^2 dz

        By default, the parameter ``only_nonzero`` is set to ``True``.
        Otherwise, the connection 1-forms being zero are shown as well::

            sage: nab.display(only_nonzero=False)
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x dx + y dy + z dz
            connection (1,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = 0
            connection (2,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = 0
            connection (2,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + y^2 dy + z^2 dz
        """
        vb = self._vbundle
        if frame is None:
            smodule = vb.section_module(domain=self._domain)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a local frame must be provided")
        dom = frame.domain()
        if frame is None:
            vmodule = dom.vector_field_module()
            vector_frame = vmodule.default_frame()
            if vector_frame is None:
                raise ValueError("a vector frame must be provided")
        if chart is None:
            chart = dom.default_chart()
            if chart is None:
                raise ValueError("a chart must be provided")
        # create output:
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion

        rlatex = r'\begin{array}{lcl}'
        rtxt = ''
        for i in vb.irange():
            for j in vb.irange():
                omega = self[frame, i, j]
                if only_nonzero and (omega == 0):
                    continue
                omega_out = omega.display(vector_frame, chart)
                rlatex += latex(omega_out) + r' \\'
                rtxt += str(omega_out) + ' \n'
        rtxt = rtxt[:-1]  # remove the last new line
        rlatex = rlatex[:-2] + r'\end{array}'
        return FormattedExpansion(rtxt, rlatex)

    def copy(self, name, latex_name=None):
        r"""
        Return an exact copy of ``self``.

        INPUT:

        - ``name`` -- name given to the copy
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          copy; if none is provided, the LaTeX symbol is set to ``name``

        .. NOTE::

            The name and the derived quantities are not copied.

        EXAMPLES::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla')
            sage: nab.set_connection_form(1, 1)[:] = [x^2, x-z, y^3]
            sage: nab.set_connection_form(1, 2)[:] = [1, x, z*y^3]
            sage: nab.set_connection_form(2, 1)[:] = [1, 2, 3]
            sage: nab.set_connection_form(2, 2)[:] = [0, 0, 0]
            sage: nab.display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + (x - z) dy + y^3 dz
            connection (1,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = dx + x dy + y^3*z dz
            connection (2,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = dx + 2 dy + 3 dz
            sage: nab_copy = nab.copy('nablo'); nab_copy
            Bundle connection nablo on the Differentiable real vector bundle
             E -> M of rank 2 over the base space 3-dimensional differentiable
             manifold M
            sage: nab is nab_copy
            False
            sage: nab == nab_copy
            True
            sage: nab_copy.display()
            connection (1,1) of bundle connection nablo w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + (x - z) dy + y^3 dz
            connection (1,2) of bundle connection nablo w.r.t. Local frame
             (E|_M, (e_1,e_2)) = dx + x dy + y^3*z dz
            connection (2,1) of bundle connection nablo w.r.t. Local frame
             (E|_M, (e_1,e_2)) = dx + 2 dy + 3 dz
        """
        copy = type(self)(self._vbundle, name, latex_name=latex_name)
        for frame, form_dict in self._connection_forms.items():
            copy._coefficients[frame] = copy._new_forms(frame=frame)
            for ind, form in form_dict.items():
                copy._coefficients[frame][ind].copy_from(form)
        return copy

    def set_immutable(self):
        r"""
        Set ``self`` and all restrictions of ``self`` immutable.

        EXAMPLES:

        An affine connection can be set immutable::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla')
            sage: nab.set_connection_form(1, 1)[:] = [x^2, x-z, y^3]
            sage: nab.set_connection_form(1, 2)[:] = [1, x, z*y^3]
            sage: nab.set_connection_form(2, 1)[:] = [1, 2, 3]
            sage: nab.set_connection_form(2, 2)[:] = [0, 0, 0]
            sage: nab.is_immutable()
            False
            sage: nab.set_immutable()
            sage: nab.is_immutable()
            True

        The coefficients of immutable elements cannot be changed::

            sage: f = E.local_frame('f')
            sage: nab.add_connection_form(1, 1, frame=f)[:] = [x, y, z]
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead
        """
        for form_dict in self._connection_forms.values():
            for form in form_dict.values():
                form.set_immutable()
        Mutability.set_immutable(self)
