# sage_setup: distribution = sagemath-schemes
"""
Periods of modular forms
"""

# The following idea just occurred to me.
# We could use that $\langle T_n(g), x\rangle = \langle g, T_n(x)\rangle$
# for any Hecke operator $T_n$, so that we only need to compute
# the period integrals $\langle g, x_i\rangle$.  Then we obtain all pairings
# $\langle T_n(g), x_i \rangle = \langle g , T_n(x_i) \rangle$.
# Since the $T_n(g)$ span the simple $\T$-module $S_k(\Gamma;\Q)[I]$,
# this must give all pairings.   However, it requires computing
# only $2d$ pairings instead of $2d^2$ pairings, which is potentially
# a huge savings when $d$ is large.
