#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

This document has the following goals:
  -- to explain the file structure and class heirarchy for p-adics in Sage
  -- to serve as a location to collect priorities for action items on p-adics
  -- to serve as a collection point for tips about altering and adding to these files.

The parent heirarchy uses multiple inheritance.  The numbers before each class name give the superclasses immediately above that class.

LocalGeneric                                          [1] (local_generic.py)
1  CappedAbsoluteGeneric                              [2] (capped_absolute_generic.py)
1  CappedRelativeGeneric                              [3] (capped_relative_generic.py)
3     CappedRelativeFieldGeneric                      [4] (capped_relative_field_generic.py)
3     CappedRelativeRingGeneric                       [5] (capped_relative_ring_generic.py)
1  FixedModGeneric                                    [6] (fixed_mod_generic.py)
1  LazyGeneric                                        [7] (lazy_generic.py)
7     LazyFieldGeneric                                [8] (lazy_field_generic.py)
7     LazyRingGeneric                                 [9] (lazy_ring_generic.py)

1  pAdicGeneric                                       [A] (padic_generic.py)
A     pAdicRingGeneric                                [B] (padic_ring_generic.py)
B, 5     pAdicCappedRelativeRingGeneric               [C] (padic_capped_relative_ring_generic.py)
B, 2     pAdicCappedAbsoluteRingGeneric               [D] (padic_capped_absolute_ring_generic.py)
B, 6     pAdicFixedModRingGeneric                     [E] (padic_fixed_mod_ring_generic.py)
B, 9     pAdicLazyRingGeneric                         [F] (padic_lazy_ring_generic.py)
A     pAdicFieldGeneric                               [G] (padic_field_generic.py)
G, 4     pAdicCappedRelativeFieldGeneric              [H] (padic_capped_relative_field_generic.py)
G, 8     pAdicLazyFieldGeneric                        [I] (padic_lazy_field_generic.py)

A     pAdicBaseGeneric                                [J] (padic_base_generic.py)
B, J     pAdicRingBaseGeneric                         [K] (padic_ring_base_generic.py)
K, C        pAdicRingCappedRelative                   [L] (padic_ring_capped_relative.py)
K, D        pAdicRingCappedAbsolute                   [M] (padic_ring_capped_absolute.py)
K, E        pAdicRingFixedMod                         [N] (padic_ring_fixed_mod.py)
K, F        pAdicRingLazy                             [O] (padic_ring_lazy.py)
G, J     pAdicFieldBaseGeneric                        [P] (padic_field_base_generic.py)
P, H        pAdicFieldCappedRelative                  [Q] (padic_field_capped_relative.py)
P, I        pAdicFieldLazy                            [R] (padic_field_lazy.py)

A     pAdicExtensionGeneric                           [S] (padic_extension_generic.py)
S        EisensteinExtensionGeneric                   [T] (eisenstein_extension_generic.py)
T, C        EisensteinExtensionRingCappedRelative     [U] (padic_extension_leaves.py)
T, H        EisensteinExtensionFieldCappedRelative    [V] (padic_extension_leaves.py)
T, D        EisensteinExtensionRingCappedAbsolute     [W] (padic_extension_leaves.py)
T, E        EisensteinExtensionRingFixedMod           [X] (padic_extension_leaves.py)
T, F        EisensteinExtensionRingLazy               [Y] (padic_extension_leaves.py, does not yet exist)
T, I        EisensteinExtensionFieldLazy              [Z] (padic_extension_leaves.py, does not yet exist)
S        UnramifiedExtensionGeneric                   [a] (unramified_extension_generic.py)
a, C        UnramifiedExtensionRingCappedRelative     [b] (padic_extension_leaves.py)
a, H        UnramifiedExtensionFieldCappedRelative    [c] (padic_extension_leaves.py)
a, D        UnramifiedExtensionRingCappedAbsolute     [d] (padic_extension_leaves.py)
a, E        UnramifiedExtensionRingFixedMod           [e] (padic_extension_leaves.py)
a, F        UnramifiedExtensionRingLazy               [f] (padic_extension_leaves.py, does not yet exist)
a, I        UnramifiedExtensionFieldLazy              [g] (padic_extension_leaves.py, does not yet exist)
S        TwoStepExtensionGeneric                      [h] (two_step_extension_generic.py, does not yet exist)
h, C        TwoStepExtensionRingCappedRelative        [i] (padic_extension_leaves.py, does not yet exist)
h, H        TwoStepExtensionFieldCappedRelative       [j] (padic_extension_leaves.py, does not yet exist)
h, D        TwoStepExtensionRingCappedAbsolute        [k] (padic_extension_leaves.py, does not yet exist)
h, E        TwoStepExtensionRingFixedMod              [l] (padic_extension_leaves.py, does not yet exist)
h, F        TwoStepExtensionRingLazy                  [m] (padic_extension_leaves.py, does not yet exist)
h, I        TwoStepExtensionFieldLazy                 [n] (padic_extension_leaves.py, does not yet exist)
S        pAdicRelativeExtensionGeneric                [h] (padic_relative_extension_generic.py, does not yet exist)
h, C        pAdicRelativeExtensionRingCappedRelative  [i] (padic_extension_leaves.py, does not yet exist)
h, H        pAdicRelativeExtensionFieldCappedRelative [j] (padic_extension_leaves.py, does not yet exist)
h, D        pAdicRelativeExtensionRingCappedAbsolute  [k] (padic_extension_leaves.py, does not yet exist)
h, E        pAdicRelativeExtensionRingFixedMod        [l] (padic_extension_leaves.py, does not yet exist)
h, F        pAdicRelativeExtensionRingLazy            [m] (padic_extension_leaves.py, does not yet exist)
h, I        pAdicRelativeExtensionFieldLazy           [n] (padic_extension_leaves.py, does not yet exist)


The element heirarchy is simpler (it's a tree), because Cython does not allow multiple inheritance.
Hopefully it will eventually include power series rings.

LocalGenericElement (local_generic_element.pyx)                                   :
  pAdicGenericElement (padic_generic_element.pyx)                                 :
    pAdicBaseGenericElement (padic_base_generic_element.pyx)                      : the base elements are implemented using GMP mpz_t
      pAdicCappedRelativeElement (padic_capped_relative_element.pyx)              :
      pAdicCappedAbsoluteElement (padic_capped_absolute_element.pyx)              :
      pAdicFixedModElement (padic_fixed_mod_element.pyx)                          :
      pAdicLazyElement (padic_lazy_element.py, currently not suppored)            :
    pAdicExtElement (padic_ext_element.pyx)                                       :
      pAdicZZpXElement (padic_ZZ_pX_element.pyx)                                  : unramified and eisenstein extensions of Qp and Zp, using ntl ZZ_pX
        pAdicZZpXCRElement (padic_ZZ_pX_CR_element.pyx)                           :
        pAdicZZpXCAElement (padic_ZZ_pX_CA_element.pyx)                           :
        pAdicZZpXFMElement (padic_ZZ_pX_FM_element.pyx)                           :
        pAdicZZpXLElement (padic_ZZ_pX_L_element.pyx, does not yet exist)         :
      pAdicZZpEXElement (padic_ZZ_pEX_element.pyx, does not yet exist)            : generic absolute extensions of Qp and Zp, using ntl ZZ_pEX
        pAdicZZpEXCRElement (padic_ZZ_pEX_CR_element.pyx, does not yet exist)     :
        pAdicZZpEXCAElement (padic_ZZ_pEX_CA_element.pyx, does not yet exist)     :
        pAdicZZpEXFMElement (padic_ZZ_pEX_FM_element.pyx, does not yet exist)     :
        pAdicZZpEXLElement (padic_ZZ_pEX_L_element.pyx, does not yet exist)       :
      pAdicRelExtElement (padic_rel_ext_element.pyx, does not yet exist)          : relative extensions of base and extension rings and fields.
        pAdicRelExtCRElement (padic_rel_ext_CR_element.pyx, does not yet exist)   :
        pAdicRelExtCAElement (padic_rel_ext_CA_element.pyx, does not yet exist)   :
        pAdicRelExtFMElement (padic_rel_ext_FM_element.pyx, does not yet exist)   :
        pAdicRelExtLElement (padic_rel_ext_L_element.pyx, does not yet exist)     :

To Do:
  Residue Fields
  log in high ramification
  exp
  residue systems
  sqrt
  hensel lifting
  norms and traces
  degree
  Documentation (file, priority, assigned to):
    padic_printing.pyx (10/10)
    padic_capped_relative_element.pyx (8/10)
    padic_capped_absolute_element.pyx (8/10)
    padic_fixed_mod_element.pyx (8/10)
    factory.py (7/10)
    padic_generic_element.pyx (7/10)
    padic_ZZ_pX_CR_element.pyx (5/10)
    padic_ZZ_pX_CA_element.pyx (5/10)
    padic_ZZ_pX_FM_element.pyx (5/10)
    padic_ZZ_pX_element.pyx (5/10)
    padic_ext_element.pyx (5/10)
    padic_extension_generic.py (5/10)
    local_generic.py (5/10)
    padic_generic.py (5/10)
    local_generic_element.pyx (4/10)
    padic_extension_leaves.py (4/10)
    pow_computer.pyx (4/10)
    pow_computer_ext.pyx (4/10)
    eisenstein_extension_generic.py (3/10)
    unramified_extension_generic.py (3/10)
    padic_field_generic.py (2/10)
    padic_ring_generic.py (2/10)
    padic_base_generic_element.pyx (2/10)
    padic_base_generic.py (2/10)
    padic_field_capped_relative.py (2/10)
    padic_printing_defaults.py (2/10)
    padic_ring_base_generic.py (1/10)
    padic_ring_capped_absolute.py (1/10)
    padic_ring_capped_relative.py (1/10)
    padic_ring_fixed_mod.py (1/10)
    capped_absolute_generic.py (1/10)
    capped_relative_generic.py (1/10)
    fixed_mod_generic.py (1/10)
    padic_field_lazy.py (0/10)
    lazy_generic.py (0/10)
    padic_lazy_element.py (0/10)
    padic_lazy_field_generic.py (0/10)
    padic_lazy_generic.py (0/10)
    padic_lazy_ring_generic.py (0/10)
    padic_ring_lazy.py (0/10)
    valuation.py (0/10)
