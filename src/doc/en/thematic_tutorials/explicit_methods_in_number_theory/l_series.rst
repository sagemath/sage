:math:`L`-series
================

:math:`L`-series of :math:`\Delta`
----------------------------------

Thanks to wrapping work of Jennifer Balakrishnan of M.I.T., we can
compute explicitly with the :math:`L`-series of the modular form
:math:`\Delta`. Like for elliptic curves, behind the scene
this uses either Dokchitsers :math:`L`-functions calculation Pari program
or Pari own implementation.

::

    sage: L = delta_lseries(); L
    L-series of conductor 1 and weight 12
    sage: L(1)
    0.0374412812685155

:math:`L`-series of a Cusp Form
-------------------------------

In some cases we can also compute with
:math:`L`-series attached to a cusp form.

::

     sage: f = CuspForms(2,8).newforms()[0]
     sage: L = f.lseries(); L
     L-series associated to the cusp form ...
     sage: L(1)
     0.0884317737041015

:math:`L`-series of a General Newform
-------------------------------------
Computing with the :math:`L`-series of a general newform is also possible.

::

    sage: S = CuspForms(23,2); S
    Cuspidal subspace of dimension 2 of Modular Forms space of
    dimension 3 for Congruence Subgroup Gamma0(23) of weight
    2 over Rational Field
    sage: f = S.newforms('a')[0]; f
    q + a0*q^2 + (-2*a0 - 1)*q^3 + (-a0 - 1)*q^4 + 2*a0*q^5 + O(q^6)
    sage: L = f.lseries(); L
    L-series associated to the cusp form ...
    sage: L(2.1)
    0.743131662917287
