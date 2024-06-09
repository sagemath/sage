# sage_setup: distribution = sagemath-repl
# sage.doctest: needs sage.graphs sage.modules
r"""
This file contains doctests of the article ::

    Automata and Transducers
    in the Computer Algebra System Sage

by Clemens Heuberger, Daniel Krenn, and Sara Kropf, :arxiv:`1404.7458`.

IF IT BECOMES NECESSARY TO CHANGE ANY TESTS IN THIS FILE, THERE
NEEDS TO BE A ONE-YEAR DEPRECATION PERIOD. ALSO, PLEASE INFORM:

- Clemens Heuberger <clemens.heuberger@aau.at>,
- Daniel Krenn <devel@danielkrenn.at>, AND
- Sara Kropf <sara.kropf@aau.at>

IN THIS CASE REGARDING THE CHANGES!
"""

r"""

Sage example in fsm-in-sage.tex, line 376::

    sage: NAF_Abb = Transducer([(-1, 0, 0, None), (-1, 1, 1, None), (0, 0, 0, 0),
    ....:                       (0, 1, 1, 0), (1, 0, 0, 1), (1, 2, 1, -1),
    ....:                       (2, 1, 0, 0), (2, 2, 1, 0)],
    ....:                       initial_states=[-1], final_states=[0],
    ....:                       input_alphabet=[0, 1])
    sage: NAF_Abb.state(-1).format_label=lambda: r'\mathcal{I}'
    sage: NAF_Abb.latex_options(
    ....:     coordinates = {-1: (1.5, 3),
    ....:                    0: (0, 0),
    ....:                    1: (3, 0),
    ....:                    2: (6, 0)},
    ....:     initial_where = {-1: 'above'},
    ....:     format_letter=NAF_Abb.format_letter_negative,
    ....:     loop_where = lambda x: 'below')


Sage example in fsm-in-sage.tex, line 377::

    sage: str(latex(NAF_Abb))
    '\\begin{tikzpicture}[auto, initial text=, >=latex]\n\\node[state, initial, initial where=above] (v0) at (1.500000, 3.000000) {$\\mathcal{I}$};\n\\node[state, accepting] (v1) at (0.000000, 0.000000) {$0$};\n\\node[state] (v2) at (3.000000, 0.000000) {$1$};\n\\node[state] (v3) at (6.000000, 0.000000) {$2$};\n\\path[->] (v0) edge node[rotate=63.43, anchor=south] {$0\\mid \\varepsilon$} (v1);\n\\path[->] (v0) edge node[rotate=-63.43, anchor=south] {$1\\mid \\varepsilon$} (v2);\n\\path[->] (v1) edge[loop below] node {$0\\mid 0$} ();\n\\path[->] (v1.5.00) edge node[rotate=0.00, anchor=south] {$1\\mid 0$} (v2.175.00);\n\\path[->] (v2.185.00) edge node[rotate=360.00, anchor=north] {$0\\mid 1$} (v1.355.00);\n\\path[->] (v2.5.00) edge node[rotate=0.00, anchor=south] {$1\\mid \\overline{1}$} (v3.175.00);\n\\path[->] (v3.185.00) edge node[rotate=360.00, anchor=north] {$0\\mid 0$} (v2.355.00);\n\\path[->] (v3) edge[loop below] node {$1\\mid 0$} ();\n\\end{tikzpicture}'


Sage example in fsm-in-sage.tex, line 395::

    sage: NAF1 = Transducer([('I', 0, 0, None), ('I', 1, 1, None),
    ....:                    (0, 0, 0, 0), (0, 1, 1, 0),
    ....:                    (1, 0, 0, 1), (1, 2, 1, -1),
    ....:                    (2, 1, 0, 0), (2, 2, 1, 0)],
    ....:                   initial_states=['I'], final_states=[0],
    ....:                   input_alphabet=[0, 1])


Sage example in fsm-in-sage.tex, line 422::

    sage: NAF = NAF1


Sage example in fsm-in-sage.tex, line 434::

    sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False


Sage example in fsm-in-sage.tex, line 455::

    sage: str(12.digits(base=2))
    '[0, 0, 1, 1]'


Sage example in fsm-in-sage.tex, line 460::

    sage: NAF.process(12.digits(base=2)),
    ((False, 2, [0, 0, -1]),)


Sage example in fsm-in-sage.tex, line 461::

    sage: str(NAF.process(12.digits(base=2)))
    '(False, 2, [0, 0, -1])'


Sage example in fsm-in-sage.tex, line 471::

    sage: NAF_of_12 = NAF(12.digits(base=2) + [0, 0, 0])


Sage example in fsm-in-sage.tex, line 472::

    sage: str(NAF_of_12)
    '[0, 0, -1, 0, 1, 0]'


Sage example in fsm-in-sage.tex, line 482::

    sage: NAF = NAF.with_final_word_out(0)


Sage example in fsm-in-sage.tex, line 494::

    sage: NAF_of_12 = NAF(12.digits(base=2))


Sage example in fsm-in-sage.tex, line 496::

    sage: str(NAF_of_12)
    '[0, 0, -1, 0, 1]'


Sage example in fsm-in-sage.tex, line 534::

    sage: def NAF_transition(state_from, read):
    ....:     if state_from == 'I':
    ....:         write = None
    ....:         state_to = read
    ....:         return (state_to, write)
    ....:     current = 2*read + state_from
    ....:     if current % 2 == 0:
    ....:         write = 0
    ....:     elif current % 4 == 1:
    ....:         write = 1
    ....:     else:
    ....:         write = -1
    ....:     state_to = (current - write) / 2
    ....:     return (state_to, write)


Sage example in fsm-in-sage.tex, line 544::

    sage: NAF2 = Transducer(NAF_transition,
    ....:                   initial_states=['I'],
    ....:                   final_states=[0],
    ....:                   input_alphabet=[0, 1]).with_final_word_out(0)


Sage example in fsm-in-sage.tex, line 548::

    sage: NAF == NAF2
    True


Sage example in fsm-in-sage.tex, line 549::

    sage: str(NAF==NAF2)
    'True'


Sage example in fsm-in-sage.tex, line 579::

    sage: def f(state_from, read):
    ....:     current = 3*read + state_from
    ....:     write = current % 2
    ....:     state_to = (current - write) / 2
    ....:     return (state_to, write)


Sage example in fsm-in-sage.tex, line 588::

    sage: Triple = Transducer(f, input_alphabet=[0, 1],
    ....:                     initial_states=[0],
    ....:                     final_states=[0]).with_final_word_out(0)


Sage example in fsm-in-sage.tex, line 592::

    sage: three_times_four = Triple(4.digits(base=2))


Sage example in fsm-in-sage.tex, line 593::

    sage: str(three_times_four)
    '[0, 0, 1, 1]'


Sage example in fsm-in-sage.tex, line 602::

    sage: Id = Transducer([(0, 0, 0, 0), (0, 0, 1, 1)],
    ....:                 initial_states=[0], final_states=[0],
    ....:                 input_alphabet=[0, 1])


Sage example in fsm-in-sage.tex, line 608::

    sage: prebuiltId = transducers.Identity([0, 1])


Sage example in fsm-in-sage.tex, line 620::

    sage: sage.combinat.finite_state_machine.\
    ....:     FSMOldCodeTransducerCartesianProduct = False
    sage: Combined_3n_n = Triple.cartesian_product(Id).relabeled()


Sage example in fsm-in-sage.tex, line 630::

    sage: twelve_and_four = Combined_3n_n(4.digits(base=2))


Sage example in fsm-in-sage.tex, line 631::

    sage: str(twelve_and_four)
    '[(0, 0), (0, 0), (1, 1), (1, None)]'


Sage example in fsm-in-sage.tex, line 639::

    sage: def g(read0, read1):
    ....:     return ZZ(read0) - ZZ(read1)


Sage example in fsm-in-sage.tex, line 643::

    sage: Minus = transducers.operator(g, input_alphabet=[None, -1, 0, 1])


Sage example in fsm-in-sage.tex, line 644::

    sage: latex(ZZ(None))
    0


Sage example in fsm-in-sage.tex, line 650::

    sage: prebuiltMinus = transducers.sub([-1, 0, 1])


Sage example in fsm-in-sage.tex, line 654::

    sage: latex(Combined_3n_n.state(1))
    1


Sage example in fsm-in-sage.tex, line 657::

    sage: final_word_out = Combined_3n_n.state(1).final_word_out


Sage example in fsm-in-sage.tex, line 658::

    sage: str(final_word_out)
    '[(1, None)]'


Sage example in fsm-in-sage.tex, line 663::

    sage: NAF3 = Minus(Combined_3n_n).relabeled()


Sage example in fsm-in-sage.tex, line 672::

    sage: NAF_of_12 = NAF3(12.digits(base=2))


Sage example in fsm-in-sage.tex, line 673::

    sage: str(NAF_of_12)
    '[0, 0, 0, -1, 0, 1]'


Sage example in fsm-in-sage.tex, line 736::

    sage: NAF = NAF3


Sage example in fsm-in-sage.tex, line 741::

    sage: NAF3n = NAF(Triple)


Sage example in fsm-in-sage.tex, line 749::

    sage: Combined_NAF_3n_n = NAF3n.cartesian_product(NAF).relabeled()


Sage example in fsm-in-sage.tex, line 757::

    sage: T = Minus(Combined_NAF_3n_n).relabeled()


Sage example in fsm-in-sage.tex, line 762::

    sage: str(T)
    'Transducer with 9 states'


Sage example in fsm-in-sage.tex, line 769::

    sage: expansion_of_12 = T(12.digits(base=2))


Sage example in fsm-in-sage.tex, line 772::

    sage: str(expansion_of_12)
    '[0, 0, 0, 2, 0, -1, 1]'


Sage example in fsm-in-sage.tex, line 806::

    sage: def minus(trans1, trans2):
    ....:     if trans1.word_in == trans2.word_in:
    ....:         return (trans1.word_in,
    ....:                 trans1.word_out[0] - trans2.word_out[0])
    ....:     else:
    ....:         raise LookupError


Sage example in fsm-in-sage.tex, line 815::

    sage: from itertools import zip_longest
    sage: def final_minus(state1, state2):
    ....:     return [x - y for x, y in
    ....:         zip_longest(state1.final_word_out,
    ....:                     state2.final_word_out,
    ....:                     fillvalue=0)]


Sage example in fsm-in-sage.tex, line 829::

    sage: Talternative = NAF3n.product_FiniteStateMachine(
    ....:                            NAF, minus,
    ....:                            final_function=final_minus).relabeled()


Sage example in fsm-in-sage.tex, line 845::

    sage: Talternative == T
    True


Sage example in fsm-in-sage.tex, line 846::

    sage: str(Talternative==T)
    'True'


Sage example in fsm-in-sage.tex, line 854::

    sage: for t in T.iter_states():
    ....:     other = Talternative.state(t.label())
    ....:     assert t.is_final == other.is_final
    ....:     if t.is_final:
    ....:         assert t.final_word_out == other.final_word_out


Sage example in fsm-in-sage.tex, line 872::

    sage: sage.combinat.finite_state_machine.setup_latex_preamble()


Sage example in fsm-in-sage.tex, line 888::

    sage: T.set_coordinates({
    ....:     0: (-2, 0.75),
    ....:     1: (0, -1),
    ....:     2: (-6, -1),
    ....:     3: (6, -1),
    ....:     4: (-4, 2.5),
    ....:     5: (-6, 5),
    ....:     6: (6, 5),
    ....:     7: (4, 2.5),
    ....:     8: (2, 0.75)})


Sage example in fsm-in-sage.tex, line 905::

    sage: T.latex_options(format_letter=T.format_letter_negative,
    ....:                 accepting_where={
    ....:                   0: 'right',
    ....:                   1: 'below',
    ....:                   2: 'below',
    ....:                   3: 'below',
    ....:                   4: 60,
    ....:                   5: 'above',
    ....:                   6: 'above',
    ....:                   7: 120,
    ....:                   8: 'left'},
    ....:                 accepting_show_empty=True)


Sage example in fsm-in-sage.tex, line 919::

    sage: str(latex(T))
    '\\begin{tikzpicture}[auto, initial text=, >=latex, accepting text=, accepting/.style=accepting by arrow, accepting distance=7ex]\n\\node[state, initial] (v0) at (-2.000000, 0.750000) {$0$};\n\\path[->] (v0.0.00) edge node[rotate=0.00, anchor=south] {$\\$ \\mid \\varepsilon$} ++(0.00:7ex);\n\\node[state] (v1) at (0.000000, -1.000000) {$1$};\n\\path[->] (v1.270.00) edge node[rotate=450.00, anchor=south] {$\\$ \\mid \\overline{2} 0 1$} ++(270.00:7ex);\n\\node[state] (v2) at (-6.000000, -1.000000) {$2$};\n\\path[->] (v2.270.00) edge node[rotate=450.00, anchor=south] {$\\$ \\mid 0 1$} ++(270.00:7ex);\n\\node[state] (v3) at (6.000000, -1.000000) {$3$};\n\\path[->] (v3.270.00) edge node[rotate=450.00, anchor=south] {$\\$ \\mid 0 \\overline{1} 1$} ++(270.00:7ex);\n\\node[state] (v4) at (-4.000000, 2.500000) {$4$};\n\\path[->] (v4.60.00) edge node[rotate=60.00, anchor=south] {$\\$ \\mid 1$} ++(60.00:7ex);\n\\node[state] (v5) at (-6.000000, 5.000000) {$5$};\n\\path[->] (v5.90.00) edge node[rotate=90.00, anchor=south] {$\\$ \\mid \\overline{1} 0 1$} ++(90.00:7ex);\n\\node[state] (v6) at (6.000000, 5.000000) {$6$};\n\\path[->] (v6.90.00) edge node[rotate=90.00, anchor=south] {$\\$ \\mid \\overline{1} 1$} ++(90.00:7ex);\n\\node[state] (v7) at (4.000000, 2.500000) {$7$};\n\\path[->] (v7.120.00) edge node[rotate=300.00, anchor=south] {$\\$ \\mid 1 \\overline{1} 1$} ++(120.00:7ex);\n\\node[state] (v8) at (2.000000, 0.750000) {$8$};\n\\path[->] (v8.180.00) edge node[rotate=360.00, anchor=south] {$\\$ \\mid 0 \\overline{2} 0 1$} ++(180.00:7ex);\n\\path[->] (v0) edge[loop above] node {$0\\mid 0$} ();\n\\path[->] (v0) edge node[rotate=-41.19, anchor=south] {$1\\mid 0$} (v1);\n\\path[->] (v1) edge node[rotate=360.00, anchor=south] {$0\\mid \\overline{2}$} (v2);\n\\path[->] (v1) edge node[rotate=0.00, anchor=south] {$1\\mid 2$} (v3);\n\\path[->] (v2) edge node[rotate=60.26, anchor=south] {$0\\mid 0$} (v4);\n\\path[->] (v2.95.00) edge node[rotate=90.00, anchor=south] {$1\\mid 0$} (v5.265.00);\n\\path[->] (v3.95.00) edge node[rotate=90.00, anchor=south] {$0\\mid 0$} (v6.265.00);\n\\path[->] (v3) edge node[rotate=299.74, anchor=south] {$1\\mid 0$} (v7);\n\\path[->] (v4) edge node[rotate=-41.19, anchor=south] {$0\\mid 1$} (v0);\n\\path[->] (v4) edge node[rotate=308.66, anchor=south] {$1\\mid \\overline{1}$} (v5);\n\\path[->] (v5.-85.00) edge node[rotate=90.00, anchor=north] {$0\\mid \\overline{1}$} (v2.85.00);\n\\path[->] (v5) edge node[rotate=-14.04, anchor=south] {$1\\mid 1$} (v7);\n\\path[->] (v6.-85.00) edge node[rotate=90.00, anchor=north] {$1\\mid 1$} (v3.85.00);\n\\path[->] (v6) edge node[rotate=14.04, anchor=south] {$0\\mid \\overline{1}$} (v4);\n\\path[->] (v7) edge node[rotate=51.34, anchor=south] {$0\\mid 1$} (v6);\n\\path[->] (v7) edge node[rotate=41.19, anchor=south] {$1\\mid \\overline{1}$} (v8);\n\\path[->] (v8) edge node[rotate=41.19, anchor=south] {$0\\mid 0$} (v1);\n\\path[->] (v8) edge[loop above] node {$1\\mid 0$} ();\n\\end{tikzpicture}'


Sage example in fsm-in-sage.tex, line 946::

    sage: R = T.output_projection()


Sage example in fsm-in-sage.tex, line 951::

    sage: latex(len(R.states()))
    10


Sage example in fsm-in-sage.tex, line 955::

    sage: R = R.split_transitions()


Sage example in fsm-in-sage.tex, line 956::

    sage: latex(len(R.states()))
    23


Sage example in fsm-in-sage.tex, line 959::

    sage: str(R.is_deterministic())
    'False'


Sage example in fsm-in-sage.tex, line 964::

    sage: Rdet = R.determinisation()


Sage example in fsm-in-sage.tex, line 967::

    sage: latex(len(Rdet.states()))
    22


Sage example in fsm-in-sage.tex, line 974::

    sage: Rdet12 = Rdet(expansion_of_12)


Sage example in fsm-in-sage.tex, line 977::

    sage: str(Rdet12)
    'True'


Sage example in fsm-in-sage.tex, line 984::

    sage: Rdet1 = Rdet.minimization()


Sage example in fsm-in-sage.tex, line 986::

    sage: latex(len(Rdet1.states()))
    17


Sage example in fsm-in-sage.tex, line 999::

    sage: Rdet2 = R.minimization(algorithm='Brzozowski')


Sage example in fsm-in-sage.tex, line 1001::

    sage: latex(len(Rdet2.states()))
    17


Sage example in fsm-in-sage.tex, line 1034::

    sage: def weight(state_from, read):
    ....:     write = ZZ(read != 0)
    ....:     return (0, write)
    sage: Weight = Transducer(weight, input_alphabet=srange(-2, 2+1),
    ....:                     initial_states=[0], final_states=[0])


Sage example in fsm-in-sage.tex, line 1044::

    sage: prebuiltWeight = transducers.weight(srange(-2, 2+1))


Sage example in fsm-in-sage.tex, line 1050::

    sage: W = Weight(T)


Sage example in fsm-in-sage.tex, line 1051::

    sage: latex(len(W.states()))
    9


Sage example in fsm-in-sage.tex, line 1056::

    sage: W(12.digits(base=2))
    [0, 0, 0, 1, 0, 1, 1]


Sage example in fsm-in-sage.tex, line 1057::

    sage: str(W(12.digits(base=2)))
    '[0, 0, 0, 1, 0, 1, 1]'


Sage example in fsm-in-sage.tex, line 1058::

    sage: latex(add(W(12.digits(base=2))))
    3


Sage example in fsm-in-sage.tex, line 1064::

    sage: W.prepone_output()


Sage example in fsm-in-sage.tex, line 1091::

    sage: var('y')                                                                      # needs sage.symbolic
    y
    sage: def am_entry(trans):
    ....:     return y^add(trans.word_out) / 2
    sage: A = W.adjacency_matrix(entry=am_entry)                                        # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1097::

    sage: latex.matrix_column_alignment('c')


Sage example in fsm-in-sage.tex, line 1099::

    sage: latex(A)                                                                      # needs sage.symbolic
    \left(\begin{array}{ccccccccc}
    \frac{1}{2} & \frac{1}{2} \, y^{2} & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & \frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{2} & \frac{1}{2} & 0 \\
    \frac{1}{2} & 0 & 0 & 0 & 0 & \frac{1}{2} \, y & 0 & 0 & 0 \\
    0 & 0 & \frac{1}{2} \, y & 0 & 0 & 0 & 0 & \frac{1}{2} \, y & 0 \\
    0 & 0 & 0 & \frac{1}{2} \, y & \frac{1}{2} \, y & 0 & 0 & 0 & 0 \\
    0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{2} \, y & 0 & \frac{1}{2} \\
    0 & \frac{1}{2} \, y^{2} & 0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{2}
    \end{array}\right)


Sage example in fsm-in-sage.tex, line 1109::

    sage: A1mI = (A.subs(y=1) - A.parent().identity_matrix())                           # needs sage.symbolic
    sage: (pi_not_normalized,) = A1mI.left_kernel().basis()                             # needs sage.symbolic
    sage: pi = pi_not_normalized / pi_not_normalized.norm(p=1)                          # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1110::

    sage: str(pi)                                                                       # needs sage.symbolic
    '(1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9, 1/9)'


Sage example in fsm-in-sage.tex, line 1117::

    sage: expected_output = derivative(A, y).subs(y=1) * vector(len(W.states())*[1])    # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1118::

    sage: latex(expected_output)                                                        # needs sage.symbolic
    \left(1,\,0,\,0,\,0,\,\frac{1}{2},\,1,\,1,\,\frac{1}{2},\,1\right)


Sage example in fsm-in-sage.tex, line 1126::

    sage: pi * expected_output                                                          # needs sage.symbolic
    5/9


Sage example in fsm-in-sage.tex, line 1127::

    sage: latex(pi * expected_output)                                                   # needs sage.symbolic
    \frac{5}{9}


Sage example in fsm-in-sage.tex, line 1129::

    sage: latex(pi * expected_output)                                                   # needs sage.symbolic
    \frac{5}{9}


Sage example in fsm-in-sage.tex, line 1145::

    sage: var('k')                                                                      # needs sage.symbolic
    k
    sage: moments = W.asymptotic_moments(k)                                             # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1155::

    sage: latex(moments['expectation'])                                                 # needs sage.symbolic
    \frac{5}{9} \, k + \mathcal{O}\left(1\right)


Sage example in fsm-in-sage.tex, line 1162::

    sage: latex(moments['variance'])                                                    # needs sage.symbolic
    \frac{44}{243} \, k + \mathcal{O}\left(1\right)


Sage example in fsm-in-sage.tex, line 1192::

    sage: expectation_binary = Id.asymptotic_moments(k)['expectation']                  # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1195::

    sage: latex(expectation_binary)                                                     # needs sage.symbolic
    \frac{1}{2} \, k + \mathcal{O}\left(1\right)


Sage example in fsm-in-sage.tex, line 1202::

    sage: expectation_NAF = Weight(NAF).asymptotic_moments(k)['expectation']            # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1205::

    sage: latex(expectation_NAF)                                                        # needs sage.symbolic
    \frac{1}{3} \, k + \mathcal{O}\left(1\right)


Sage example in fsm-in-sage.tex, line 1211::

    sage: Abs = transducers.abs([-1, 0, 1])                                             # needs sage.symbolic


Sage example in fsm-in-sage.tex, line 1216::

    sage: latex(moments['expectation'])                                                 # needs sage.symbolic
    \frac{5}{9} \, k + \mathcal{O}\left(1\right)
"""
