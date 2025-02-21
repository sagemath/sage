import ipywidgets as widgets
from sage.misc.latex import latex
from sage.repl.rich_output.pretty_print import pretty_print
from IPython.display import clear_output


def cluster_interact(self, fig_size=1, circular=True, kind='seed'):
    r"""
    Start an interactive window for cluster seed mutations.

    Only in *Jupyter notebook mode*.

    Not to be called directly. Use the :meth:`interact` methods
    of :class:`ClusterSeed` and :class:`ClusterQuiver` instead.

    INPUT:

    - ``fig_size`` -- (default: 1) factor by which the size of the
      plot is multiplied

    - ``circular`` -- boolean (default: ``True``); if ``True``, the circular
      plot is chosen, otherwise >>spring<< is used

    - ``kind`` -- either ``'seed'`` (default) or ``'quiver'``

    TESTS::

        sage: S = ClusterSeed(['A',4])                                                  # needs sage.graphs sage.modules
        sage: S.interact()   # indirect doctest                                         # needs sage.graphs sage.modules sage.symbolic
        ...VBox(children=...
    """
    if kind not in ['seed', 'quiver']:
        raise ValueError('kind must be "seed" or "quiver"')

    show_seq = widgets.Checkbox(value=True,
                                description="Display mutation sequence")

    show_vars = widgets.Checkbox(value=True,
                                 description="Display cluster variables")

    show_matrix = widgets.Checkbox(value=True,
                                   description="Display B-matrix")

    show_lastmutation = widgets.Checkbox(value=True,
                                         description="Show last mutation vertex")

    mut_buttons = widgets.ToggleButtons(options=list(range(self._n)),
                                       style={'button_width':'initial'},
                                       description='Mutate at: ')

    which_plot = widgets.Dropdown(options=['circular', 'spring'],
                                  value='circular' if circular else "spring",
                                  description='Display style:')

    out = widgets.Output()

    seq = []

    def print_data():
        if show_seq.value:
            pretty_print("Mutation sequence: ", seq)

        if show_vars.value and kind == 'seed':
            pretty_print("Cluster variables:")
            table = "\\begin{align*}\n"
            for i in range(self._n):
                table += "\tv_{%s} &= " % i + latex(self.cluster_variable(i)) + "\\\\ \\\\\n"
            table += "\\end{align*}"
            pretty_print(table)

        if show_matrix.value:
            pretty_print("B-Matrix: ", self._M)

    def refresh(w):
        k = mut_buttons.value
        circular = bool(which_plot.value == "circular")
        with out:
            clear_output(wait=True)
            if show_lastmutation.value:
                self.show(fig_size=fig_size, circular=circular, mark=k)
            else:
                self.show(fig_size=fig_size, circular=circular)
            print_data()

    def do_mutation(*args, **kwds):
        k = mut_buttons.value
        circular = bool(which_plot.value == "circular")
        self.mutate(k)
        if seq and seq[-1] == k:
            seq.pop()
        else:
            seq.append(k)
        with out:
            clear_output(wait=True)
            if show_lastmutation.value:
                self.show(fig_size=fig_size, circular=circular, mark=k)
            else:
                self.show(fig_size=fig_size, circular=circular)
            print_data()

    mut_buttons.on_msg(do_mutation)

    show_seq.observe(refresh, 'value')
    if kind == 'seed':
        show_vars.observe(refresh, 'value')
    show_matrix.observe(refresh, 'value')
    show_lastmutation.observe(refresh, 'value')
    which_plot.observe(refresh, 'value')

    mut_buttons.on_widget_constructed(refresh)

    if kind == 'seed':
        top = widgets.HBox([show_seq, show_vars])
    else:
        top = widgets.HBox([show_seq])

    return widgets.VBox([which_plot,
                         top,
                         widgets.HBox([show_matrix, show_lastmutation]),
                         mut_buttons, out])
