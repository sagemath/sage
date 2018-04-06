"""
Shimura Mass
"""
######################################################
## Routines to compute the mass of a quadratic form ##
######################################################

## Import all general mass finding routines
from sage.quadratic_forms.quadratic_form__mass__Siegel_densities import \
        mass__by_Siegel_densities, \
        Pall_mass_density_at_odd_prime, \
        Watson_mass_at_2, \
        Kitaoka_mass_at_2, \
        mass_at_two_by_counting_mod_power
from sage.quadratic_forms.quadratic_form__mass__Conway_Sloane_masses import \
        parity, \
        is_even, \
        is_odd, \
        conway_species_list_at_odd_prime, \
        conway_species_list_at_2, \
        conway_octane_of_this_unimodular_Jordan_block_at_2, \
        conway_diagonal_factor, \
        conway_cross_product_doubled_power, \
        conway_type_factor, \
        conway_p_mass, \
        conway_standard_p_mass, \
        conway_standard_mass, \
        conway_mass
#        conway_generic_mass, \
#        conway_p_mass_adjustment

###################################################


def shimura_mass__maximal(self,):
    """
    Use Shimura's exact mass formula to compute the mass of a maximal
    quadratic lattice. This works for any totally real number field,
    but has a small technical restriction when `n` is odd.

    INPUT:

        none

    OUTPUT:

        a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.shimura_mass__maximal()

    """
    pass



def GHY_mass__maximal(self):
    """
    Use the GHY formula to compute the mass of a (maximal?) quadratic
    lattice. This works for any number field.

    Reference:  See [GHY, Prop 7.4 and 7.5, p121] and [GY, Thrm 10.20, p25].

    INPUT:

        none

    OUTPUT:

        a rational number

    EXAMPLES::

        sage: Q = DiagonalQuadraticForm(ZZ, [1,1,1])
        sage: Q.GHY_mass__maximal()

    """
    pass
