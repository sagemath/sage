#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed.math@gmail.com>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

This document aims at explaining the (somehow complicated) file structure and
class hierarchy for `p`-adics in Sage.


Below is a representation of the parent hierarchy.
The representation is simplified is the sense that multiple inheritance does not
appear; we assume however that these relationships are obvious from the names of
the classes (e.g. pAdicRingBaseGeneric derives from pAdicBaseGeneric).
The file indicated between brackets is the file where the class is defined.


pAdicGeneric ......................................................... [padic_generic.py]

    pAdicRingGeneric ................................................. [generic_nodes.py]
    pAdicFieldGeneric ................................................ [generic_nodes.py]
    pAdicLatticeGeneric .............................................. [generic_nodes.py]
    pAdicRelaxedGeneric .............................................. [generic_nodes.py]

    pAdicBaseGeneric ................................................. [padic_base_generic.py]

        pAdicRingBaseGeneric ......................................... [generic_nodes.py]
            pAdicRingCappedRelative .................................. [padic_base_leaves.py]
            pAdicRingCappedAbsolute .................................. [padic_base_leaves.py]
            pAdicRingFloatingPoint ................................... [padic_base_leaves.py]
            pAdicRingFixedMod ........................................ [padic_base_leaves.py]
            pAdicRingLattice ......................................... [padic_base_leaves.py]
            pAdicRingRelaxed ......................................... [padic_base_leaves.py]

        pAdicFieldBaseGeneric ........................................ [generic_nodes.py]
            pAdicFieldCappedRelative ................................. [padic_base_leaves.py]
            pAdicFieldFloatingPoint .................................. [padic_base_leaves.py]
            pAdicFieldLattice ........................................ [padic_base_leaves.py]
            pAdicFieldRelaxed ........................................ [padic_base_leaves.py]

    pAdicExtensionGeneric ............................................ [padic_extension_generic.py]

        UnramifiedExtensionGeneric ................................... [unramified_extension_generic.py]
            UnramifiedExtensionRingCappedRelative .................... [padic_extension_leaves.py]
            UnramifiedExtensionRingCappedAbsolute .................... [padic_extension_leaves.py]
            UnramifiedExtensionRingFixedMod .......................... [padic_extension_leaves.py]
            UnramifiedExtensionRingFloatingPoint ..................... [padic_extension_leaves.py]
            UnramifiedExtensionFieldCappedRelative ................... [padic_extension_leaves.py]
            UnramifiedExtensionFieldFloatingPoint .................... [padic_extension_leaves.py]

        EisensteinExtensionGeneric ................................... [eisenstein_extension_generic.py]
            EisensteinExtensionRingCappedRelative .................... [padic_extension_leaves.py]
            EisensteinExtensionRingCappedAbsolute .................... [padic_extension_leaves.py]
            EisensteinExtensionRingFixedMod .......................... [padic_extension_leaves.py]
            EisensteinExtensionFieldCappedRelative ................... [padic_extension_leaves.py]
            # Eisenstein extensions on top of an unramified extension
            RelativeRamifiedExtensionRingCappedRelative .............. [relative_extension_leaves.py]
            RelativeRamifiedExtensionRingCappedAbsolute .............. [relative_extension_leaves.py]
            RelativeRamifiedExtensionRingFixedMod .................... [relative_extension_leaves.py]
            RelativeRamifiedExtensionRingFloatingPoint ............... [relative_extension_leaves.py]
            RelativeRamifiedExtensionFieldCappedRelative ............. [relative_extension_leaves.py]
            RelativeRamifiedExtensionFieldFloatingPoint .............. [relative_extension_leaves.py]


Below is the hierarchy for element classes.
Here all dependencies are represented since Cython does not allow multiple inheritance.


pAdicGenericElement .................................................. [padic_generic_element.pyx]
    pAdicTemplateElement ............................................. [padic_template_element.pxi]
        CAElement .................................................... [CA_template.pxi]
            pAdicCappedAbsoluteElement ............................... [padic_capped_absolute_element.pyx]
            qAdicCappedAbsoluteElement ............................... [qadic_flint_CA.pyx]
            RelativeRamifiedCappedAbsoluteElement .................... [relative_ramified_CA.pyx]

        CRElement .................................................... [CR_template.pxi]
            pAdicCappedRelativeElement ............................... [padic_capped_relative_element.pyx]
            qAdicCappedRelativeElement ............................... [qadic_flint_CR.pyx]
            RelativeRamifiedCappedRelativeElement .................... [relative_ramified_CR.pyx]

        FMElement .................................................... [FM_template.pxi]
            pAdicFixedModElement ..................................... [padic_fixed_mod_element.pyx]
            qAdicFixedModElement ..................................... [qadic_flint_FM.pyx]
            RelativeRamifiedFixedModElement .......................... [relative_ramified_FM.pyx]

        FPElement .................................................... [FP_template.pxi]
            pAdicFloatingPointElement ................................ [padic_floating_point_element.pyx]
            qAdicFloatingPointElement ................................ [qadic_flint_FP.pyx]
            RelativeRamifiedFloatingPointElement ..................... [relative_ramified_FP.pyx]

    pAdicExtElement .................................................. [padic_ext_element.pyx]
        pAdicZZpXElement ............................................. [padic_ZZ_pX_element.pyx]
            pAdicZZpXCAElement ....................................... [padic_ZZ_pX_CA_element.pyx]
            pAdicZZpXCRElement ....................................... [padic_ZZ_pX_CR_element.pyx]
            pAdicZZpXFMElement ....................................... [padic_ZZ_pX_FM_element.pyx]

    pAdicLatticeElement .............................................. [padic_lattice_element.py]
        pAdicLatticeCapElement ....................................... [padic_lattice_element.py]
        pAdicLatticeFloatElement ..................................... [padic_lattice_element.py]

    RelaxedElement ................................................... [relaxed_template.pxi]
        RelaxedElement_* ............................................. [relaxed_template.pxi]
            pAdicRelaxedElement_* .................................... [padic_relaxed_element.pxd]
