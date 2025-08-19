# sage.doctest: optional - p_group_cohomology
r"""
Tests for the optional ``p_group_cohomology`` package.

AUTHOR:

- Simon King

TESTS::

    sage: from pGroupCohomology import CohomologyRing

Computation of a modular cohomology ring of a prime power group in
characteristic 2, and comparison with stored result in a database::

    sage: CohomologyRing.set_workspace(tmp_dir())
    sage: H = CohomologyRing(64,14,from_scratch=True)
    sage: H.make()
    sage: CohomologyRing.set_workspace(tmp_dir())
    sage: H0 = CohomologyRing(64,14)
    sage: H.is_isomorphic(H0)
    ('1*a_2_1', '1*c_2_2', '1*c_4_4', '1*a_1_0', '1*a_1_1', '1*a_3_3')

Computation of a modular cohomology ring of a prime power group in odd
characteristic, and some algebraic constructions in the cohomology
ring::

    sage: H = CohomologyRing(27,4)
    sage: H.make()
    sage: print(H)
    <BLANKLINE>
    Cohomology ring of Extraspecial 3-group of order 27
    and exponent 9 with coefficients in GF(3)
    <BLANKLINE>
    Computation complete
    Minimal list of generators:
    [b_2_1: 2-Cocycle in H^*(M27; GF(3)),
     c_6_2: 6-Cocycle in H^*(M27; GF(3)),
     a_1_0: 1-Cocycle in H^*(M27; GF(3)),
     a_1_1: 1-Cocycle in H^*(M27; GF(3)),
     a_3_1: 3-Cocycle in H^*(M27; GF(3)),
     a_5_1: 5-Cocycle in H^*(M27; GF(3))]
    Minimal list of algebraic relations:
    [b_2_1*a_1_0,
     a_1_0*a_3_1,
     b_2_1*a_3_1,
     a_1_0*a_5_1,
     a_3_1*a_5_1]
    sage: H.5.massey_power()
    <a_3_1; 1>: 8-Cocycle in H^*(M27; GF(3))
    sage: H.5.massey_power().as_polynomial()
    '-c_6_2*a_1_0*a_1_1'
    sage: H.essential_ideal()
    a_1_0*a_1_1,
    a_1_1*a_3_1
    sage: ascii_art(H.bar_code('LowerCentralSeries')[2])        # known bug
        *
      *-*
      *-*
    *

Computation of a modular cohomology ring of a non prime power group in
characteristic 2::

    sage: H = CohomologyRing(libgap.AlternatingGroup(6),
    ....:                    GroupName='A(6)', prime=2,
    ....:                    from_scratch=True)
    sage: H.make()
    sage: print(H)
    Cohomology ring of A(6) with coefficients in GF(2)
    <BLANKLINE>
    Computation complete
    Minimal list of generators:
    [c_2_0: 2-Cocycle in H^*(A(6); GF(2)),
     b_3_0: 3-Cocycle in H^*(A(6); GF(2)),
     b_3_1: 3-Cocycle in H^*(A(6); GF(2))]
    Minimal list of algebraic relations:
    [b_3_0*b_3_1]
"""
