r"""
LaTeX macros

AUTHORS:

- John H. Palmieri (2009-03)

The code here sets up LaTeX macro definitions for use in the
documentation. To add a macro, modify the list ``macros``, near the
end of this file, and then run 'sage -b'. The entries in this list are
used to produce ``sage_latex_macros``, a list of strings of the form
'\\newcommand...', and ``sage_mathjax_macros``, a list of strings
suitable for parsing by MathJax.  The LaTeX macros are produced using
the ``_latex_`` method for each Sage object listed in ``macros``, and
the MathJax macros are produced from the LaTeX macros.  The list of
LaTeX macros is used in the file
``sage_docbuild.conf`` to add to the preambles of
both the LaTeX file used to build the PDF version of the documentation
and the LaTeX file used to build the HTML version.

Any macro defined here may be used in docstrings or in the tutorial
(or other pieces of documentation).  In a docstring, for example,
"\ZZ" in backquotes (demarking math mode) will appear as "ZZ" in
interactive help, but will be typeset as "\\Bold{Z}" in the
reference manual.

More details on the list ``macros``: the entries are lists or tuples
of the form ``[name]`` or ``[name, arguments]``, where ``name`` is a
string and ``arguments`` consists of valid arguments for the Sage
object named ``name``.  For example, ``["ZZ"]`` and ``["GF", 2]``
produce the LaTeX macros '\\newcommand{\\ZZ}{\\Bold{Z}}' and
'\\newcommand{\\GF}[1]{\\Bold{F}_{#1}}', respectively.  (For the
second of these, ``latex(GF(2))`` is called and the string '2' gets
replaced by '#1', so ``["GF", 17]`` would have worked just as well.
``["GF", p]`` would have raised an error, though, because ``p`` is not
defined, and ``["GF", 4]`` would have raised an error, because to
define the field with four elements in Sage, you also need to specify
the name of a generator.)

To see evidence of the results of the code here, run ``sage --docbuild
tutorial latex`` (for example), and look at the resulting LaTeX file in
``SAGE_DOC/latex/en/tutorial/``.  The preamble should
contain '\newcommand' lines for each of the entries in ``macros``.
"""


def produce_latex_macro(name, *sample_args):
    r"""
    Produce a string defining a LaTeX macro.

    INPUT:

    - ``name`` -- name of macro to be defined, also name of corresponding Sage object

    - ``sample_args`` -- (optional) sample arguments for this Sage object

    EXAMPLES::

        sage: from sage.misc.latex_macros import produce_latex_macro
        sage: produce_latex_macro('ZZ')
        '\\newcommand{\\ZZ}{\\Bold{Z}}'

    If the Sage object takes arguments, then the LaTeX macro will
    accept arguments as well. You must pass valid arguments, which
    will then be converted to #1, #2, etc. in the macro
    definition. The following allows the use of "\GF{p^n}", for
    example::

         sage: produce_latex_macro('GF', 37)
         '\\newcommand{\\GF}[1]{\\Bold{F}_{#1}}'

    If the Sage object is not in the global name space, describe it
    like so::

         sage: produce_latex_macro('sage.rings.finite_rings.finite_field_constructor.FiniteField', 3)
         '\\newcommand{\\FiniteField}[1]{\\Bold{F}_{#1}}'
    """
    from sage.misc.latex import LatexCall  # type: ignore
    # this import is used inside a string below
    names_split = name.rsplit('.', 1)
    if len(names_split) == 1:
        module = 'sage.all'
        real_name = names_split[0]
    else:
        module, real_name = names_split
    newcommand = '\\newcommand{\\' + real_name + '}'
    count = 0
    args = "("
    for x in sample_args:
        count += 1
        args += str(x) + ','
    args += ')'
    exec('from ' + module + ' import ' + real_name)
    if count:
        defn = '[' + str(count) + ']{'
        defn += eval('str(LatexCall()(' + real_name + args + '))') + '}'
    else:
        defn = '{' + eval('str(LatexCall()(' + real_name + '))') + '}'
    count = 0
    for x in sample_args:
        count += 1
        defn = defn.replace(str(x), "#" + str(count))
    return newcommand + defn


def convert_latex_macro_to_mathjax(macro):
    r"""
    This converts a LaTeX macro definition (\newcommand...) to a
    MathJax macro definition (MathJax.Macro...).

    INPUT:

    - ``macro`` -- LaTeX macro definition

    See the web page
    https://docs.mathjax.org/en/latest/input/tex/macros.html for a
    description of the format for MathJax macros.

    EXAMPLES::

        sage: from sage.misc.latex_macros import convert_latex_macro_to_mathjax
        sage: convert_latex_macro_to_mathjax('\\newcommand{\\ZZ}{\\Bold{Z}}')
        ('ZZ', '\\Bold{Z}')
        sage: convert_latex_macro_to_mathjax('\\newcommand{\\GF}[1]{\\Bold{F}_{#1}}')
        ('GF', ['\\Bold{F}_{#1}', 1])
    """
    left_bracket = macro.find('[')
    right_bracket = macro.find('[')
    if left_bracket >= 0:
        right_bracket = macro.find(']')
        num_args = int(macro[left_bracket + 1: right_bracket])
    else:
        num_args = 0
    start_name = macro.find('{') + 1  # add one to go past the backslash
    end_name = macro.find('}')
    name = macro[start_name + 1: end_name]
    start_defn = macro.find('{', end_name)
    end_defn = macro.rfind('}')
    defn = macro[start_defn + 1: end_defn]
    if num_args == 0:
        return name, defn
    else:
        return name, [defn, num_args]


# To add a new macro for use in the Sage documentation, add a list or
# tuple to the following list.  Each list (or tuple) should have the
# form [name, arguments], which will be passed to the function
# produce_latex_macro: see that for more documentation.
#
# To see the results of this, run 'sage --docbuild tutorial latex' (for
# example -- you could replace 'tutorial' with your favorite piece of
# documentation), and look at the resulting tex file in
# SAGE_DOC/latex/en/tutorial.  The preamble should contain
# \newcommand's for each of the entries here.
macros = [["ZZ"],
          ["NN"],
          ["RR"],
          ["CC"],
          ["QQ"],
          ["QQbar"],
          ["GF", 2],
          ["Zp", 2],
          ["Qp", 2],
          ["Zmod", 2],
          ["CDF"],
          ["CIF"],
          ["CLF"],
          ["RDF"],
          ["RIF"],
          ["RLF"],
          ]

# Use this list to define additional latex macros for sage documentation
latex_macros = [r"\newcommand{\SL}{\mathrm{SL}}",
                r"\newcommand{\PSL}{\mathrm{PSL}}",
                r"\newcommand{\lcm}{\mathop{\operatorname{lcm}}}",
                r"\newcommand{\dist}{\mathrm{dist}}",
                r"\newcommand{\im}{\mathop{\operatorname{im}}}",
                r"\newcommand{\rank}{\mathop{\operatorname{rank}}}",
                r"\newcommand{\PP}{\mathbf{P}}",
                r"\newcommand{\OO}{\mathcal{O}}",
                ]

# The following is to allow customization of typesetting of rings:
# mathbf vs mathbb.  See latex.py for more information.
sage_configurable_latex_macros = [r"\newcommand{\Bold}[1]{\mathbf{#1}}"]


def sage_latex_macros():
    r"""
    Return list of LaTeX macros for Sage. This just runs the function
    :func:`produce_latex_macro` on the list ``macros`` defined in this
    file, and appends ``sage_configurable_latex_macros``. To add a new
    macro for permanent use in Sage, modify ``macros``.

    EXAMPLES::

        sage: from sage.misc.latex_macros import sage_latex_macros
        sage: sage_latex_macros()
        ['\\newcommand{\\ZZ}{\\Bold{Z}}', '\\newcommand{\\NN}{\\Bold{N}}', ...
    """
    return [produce_latex_macro(*x) for x in macros] + latex_macros + sage_configurable_latex_macros


def sage_mathjax_macros():
    r"""
    Return Sage's macro definitions for usage with MathJax.

    This feeds each item output by :func:`sage_latex_macros` to
    :func:`convert_latex_macro_to_mathjax`.

    EXAMPLES::

        sage: from sage.misc.latex_macros import sage_mathjax_macros
        sage: sage_mathjax_macros()
        {'Bold': ['\\mathbf{#1}', 1], 'CC': '\\Bold{C}', ...
    """
    return dict(convert_latex_macro_to_mathjax(m) for m in sage_latex_macros())
