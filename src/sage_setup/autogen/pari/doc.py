# coding=utf-8
"""
Handle PARI documentation for Sage
"""

from __future__ import unicode_literals
import re
import subprocess
from six import unichr


leading_ws = re.compile("^ +", re.MULTILINE)
double_space = re.compile("  +")

end_space = re.compile(r"(@\[end[a-z]*\])([A-Za-z])")

begin_verb = re.compile(r"@1")
end_verb = re.compile(r"@[23] *@\[endcode\]")
verb_loop = re.compile("^(    .*)@\[[a-z]*\]", re.MULTILINE)

dollars = re.compile(r"@\[dollar\]\s*(.*?)\s*@\[dollar\]", re.DOTALL)
doubledollars = re.compile(r"@\[doubledollar\]\s*(.*?)\s*@\[doubledollar\]", re.DOTALL)

math_loop = re.compile(r"(@\[start[A-Z]*MATH\][^@]*)@\[[a-z]*\]")
math_backslash = re.compile(r"(@\[start[A-Z]*MATH\][^@]*)=BACKSLASH=")

prototype = re.compile("^[^\n]*\n\n")
library_syntax = re.compile("The library syntax is.*", re.DOTALL)

newlines = re.compile("\n\n\n\n*")

bullet_loop = re.compile("(@BULLET(  [^\n]*\n)*)([^ \n])")
indent_math = re.compile("(@\\[startDISPLAYMATH\\].*\n(.+\n)*)(\\S)")

escape_backslash = re.compile(r"^(\S.*)[\\]", re.MULTILINE)
escape_mid = re.compile(r"^(\S.*)[|]", re.MULTILINE)
escape_percent = re.compile(r"^(\S.*)[%]", re.MULTILINE)
escape_hash = re.compile(r"^(\S.*)[#]", re.MULTILINE)

label_link = re.compile(r"(Section *)?\[@\[startbold\]Label: *(se:)?([^@]*)@\[endbold\]\]")


def sub_loop(regex, repl, text):
    """
    In ``text``, substitute ``regex`` by ``repl`` recursively. As long
    as substitution is possible, ``regex`` is substituted.

    INPUT:

    - ``regex`` -- a compiled regular expression

    - ``repl`` -- replacement text

    - ``text`` -- input text

    OUTPUT: substituted text

    EXAMPLES:

    Ensure there a space between any 2 letters ``x``::

        sage: from sage_setup.autogen.pari.doc import sub_loop
        sage: import re
        sage: sub_loop(re.compile("xx"), "x x", "xxx_xx")
        u'x x x_x x'
    """
    while True:
        text, n = regex.subn(repl, text)
        if not n:
            return text


def raw_to_rest(doc):
    """
    Convert raw PARI documentation (with ``@``-codes) to reST syntax.

    INPUT:

    - ``doc`` -- the raw PARI documentation

    OUTPUT: a unicode string

    EXAMPLES::

        sage: from sage_setup.autogen.pari.doc import raw_to_rest
        sage: print(raw_to_rest("@[startbold]hello world@[endbold]"))
        :strong:`hello world`

    TESTS::

        sage: raw_to_rest("@[invalid]")
        Traceback (most recent call last):
        ...
        SyntaxError: @ found: @[invalid]
    """
    doc = doc.decode("utf-8")

    # Work around a specific problem with doc of "component"
    doc = doc.replace("[@[dollar]@[dollar]]", "[]")

    # Work around a specific problem with doc of "algdivl"
    doc = doc.replace(r"\y@", r"\backslash y@")

    # Special characters
    doc = doc.replace("@[lt]", "<")
    doc = doc.replace("@[gt]", ">")
    doc = doc.replace("@[pm]", "±")
    doc = doc.replace("@[nbrk]", unichr(0xa0))
    doc = doc.replace("@[agrave]", "à")
    doc = doc.replace("@[aacute]", "á")
    doc = doc.replace("@[eacute]", "é")
    doc = doc.replace("@[ouml]", "ö")
    doc = doc.replace("@[uuml]", "ü")
    doc = doc.replace("\\'{a}", "á")

    # Remove leading whitespace from every line
    doc = leading_ws.sub("", doc)

    # Remove multiple spaces
    doc = double_space.sub(" ", doc)

    # Sphinx dislikes inline markup immediately followed by a letter:
    # insert a non-breaking space
    doc = end_space.sub("\\1" + unichr(0xa0) + "\\2", doc)

    # Remove links
    doc = label_link.sub("``\\3`` (in the PARI manual)", doc)

    # Bullet items
    doc = doc.replace("@3@[startbold]*@[endbold] ", "@BULLET  ")
    doc = sub_loop(bullet_loop, "\\1  \\3", doc)
    doc = doc.replace("@BULLET  ", "- ")

    # Verbatim blocks
    doc = begin_verb.sub("::\n\n@0", doc)
    doc = end_verb.sub("", doc)
    doc = doc.replace("@0", "    ")
    doc = doc.replace("@3", "")

    # Remove all further markup from within verbatim blocks
    doc = sub_loop(verb_loop, "\\1", doc)

    # Pair dollars -> beginmath/endmath
    doc = doc.replace("@[dollar]@[dollar]", "@[doubledollar]")
    doc = dollars.sub(r"@[startMATH]\1@[endMATH]", doc)
    doc = doubledollars.sub(r"@[startDISPLAYMATH]\1@[endDISPLAYMATH]", doc)

    # Replace special characters (except in verbatim blocks)
    # \ -> =BACKSLASH=
    # | -> =MID=
    # % -> =PERCENT=
    # # -> =HASH=
    doc = sub_loop(escape_backslash, "\\1=BACKSLASH=", doc)
    doc = sub_loop(escape_mid, "\\1=MID=", doc)
    doc = sub_loop(escape_percent, "\\1=PERCENT=", doc)
    doc = sub_loop(escape_hash, "\\1=HASH=", doc)

    # Math markup
    doc = doc.replace("@[obr]", "{")
    doc = doc.replace("@[cbr]", "}")
    doc = doc.replace("@[startword]", "\\")
    doc = doc.replace("@[endword]", "")
    doc = doc.replace("@[startlword]", "\\")
    doc = doc.replace("@[endlword]", "")
    doc = doc.replace("@[startbi]", "\\mathbb{")
    doc = doc.replace("@[endbi]", "}")

    # PARI TeX macros
    doc = doc.replace(r"\Cl", r"\mathrm{Cl}")
    doc = doc.replace(r"\Id", r"\mathrm{Id}")
    doc = doc.replace(r"\Norm", r"\mathrm{Norm}")
    doc = doc.replace(r"\disc", r"\mathrm{disc}")
    doc = doc.replace(r"\gcd", r"\mathrm{gcd}")
    doc = doc.replace(r"\lcm", r"\mathrm{lcm}")

    # Remove extra markup inside math blocks
    doc = sub_loop(math_loop, "\\1", doc)

    # Replace special characters by escape sequences
    # Note that =BACKSLASH= becomes an unescaped backslash in math mode
    # but an escaped backslash otherwise.
    doc = sub_loop(math_backslash, r"\1\\", doc)
    doc = doc.replace("=BACKSLASH=", r"\\")
    doc = doc.replace("=MID=", r"\|")
    doc = doc.replace("=PERCENT=", r"\%")
    doc = doc.replace("=HASH=", r"\#")

    # Handle DISPLAYMATH
    doc = doc.replace("@[endDISPLAYMATH]", "\n\n")
    doc = sub_loop(indent_math, "\\1    \\3", doc)
    doc = doc.replace("@[startDISPLAYMATH]", "\n\n.. MATH::\n\n    ")

    # Inline markup. We do use the more verbose :foo:`text` style since
    # those nest more easily.
    doc = doc.replace("@[startMATH]", ":math:`")
    doc = doc.replace("@[endMATH]", "`")
    doc = doc.replace("@[startpodcode]", "``")
    doc = doc.replace("@[endpodcode]", "``")
    doc = doc.replace("@[startcode]", ":literal:`")
    doc = doc.replace("@[endcode]", "`")
    doc = doc.replace("@[startit]", ":emphasis:`")
    doc = doc.replace("@[endit]", "`")
    doc = doc.replace("@[startbold]", ":strong:`")
    doc = doc.replace("@[endbold]", "`")

    # Remove prototype
    doc = prototype.sub("", doc)

    # Remove everything starting with "The library syntax is"
    # (this is not relevant for Sage)
    doc = library_syntax.sub("", doc)

    # Allow at most 2 consecutive newlines
    doc = newlines.sub("\n\n", doc)

    # Strip result
    doc = doc.strip()

    # Ensure no more @ remains
    try:
        i = doc.index("@")
    except ValueError:
        return doc
    ilow = max(0, i-30)
    ihigh = min(len(doc), i+30)
    raise SyntaxError("@ found: " + doc[ilow:ihigh])


def get_raw_doc(function):
    r"""
    Get the raw documentation of PARI function ``function``.

    INPUT:

    - ``function`` -- name of a PARI function

    EXAMPLES::

        sage: from sage_setup.autogen.pari.doc import get_raw_doc
        sage: get_raw_doc("cos")
        '@[startbold]cos@[dollar](x)@[dollar]:@[endbold]\n\n\n\nCosine of @[dollar]x@[dollar].\n\n\nThe library syntax is @[startcode]GEN @[startbold]gcos@[endbold](GEN x, long prec)@[endcode].\n\n\n'
        sage: get_raw_doc("abcde")
        Traceback (most recent call last):
        ...
        RuntimeError: no help found for 'abcde'
    """
    doc = subprocess.check_output(["gphelp", "-raw", function])
    if doc.endswith(b"""' not found !\n"""):
        raise RuntimeError("no help found for '{}'".format(function))
    return doc


def get_rest_doc(function):
    r"""
    Get the documentation of the PARI function ``function`` in reST
    syntax.

    INPUT:

    - ``function`` -- name of a PARI function

    EXAMPLES::

        sage: from sage_setup.autogen.pari.doc import get_rest_doc
        sage: print(get_rest_doc("teichmuller"))
        Teichmüller character of the :math:`p`-adic number :math:`x`, i.e. the unique
        :math:`(p-1)`-th root of unity congruent to :math:`x / p^{v_p(x)}` modulo :math:`p`...

    ::

        sage: print(get_rest_doc("weber"))
        One of Weber's three :math:`f` functions.
        If :math:`flag = 0`, returns
        <BLANKLINE>
        .. MATH::
        <BLANKLINE>
            f(x) = \exp(-i\Pi/24).\eta((x+1)/2)/\eta(x) {such that}
            j = (f^{24}-16)^3/f^{24},
        <BLANKLINE>
        where :math:`j` is the elliptic :math:`j`-invariant (see the function :literal:`ellj`).
        If :math:`flag = 1`, returns
        <BLANKLINE>
        .. MATH::
        <BLANKLINE>
            f_1(x) = \eta(x/2)/\eta(x) {such that}
            j = (f_1^{24}+16)^3/f_1^{24}.
        <BLANKLINE>
        Finally, if :math:`flag = 2`, returns
        <BLANKLINE>
        .. MATH::
        <BLANKLINE>
            f_2(x) = \sqrt{2}\eta(2x)/\eta(x) {such that}
            j = (f_2^{24}+16)^3/f_2^{24}.
        <BLANKLINE>
        Note the identities :math:`f^8 = f_1^8+f_2^8` and :math:`ff_1f_2 = \sqrt2`.

    ::

        sage: print(get_rest_doc("ellap"))
        Let :math:`E` be an :literal:`ell` structure as output by :literal:`ellinit`, defined over
        :math:`\mathbb{Q}` or a finite field :math:`\mathbb{F}_q`. The argument :math:`p` is best left omitted if the
        curve is defined over a finite field, and must be a prime number otherwise.
        This function computes the trace of Frobenius :math:`t` for the elliptic curve :math:`E`,
        defined by the equation :math:`\#E(\mathbb{F}_q) = q+1 - t`.
        <BLANKLINE>
        When the characteristic of the finite field is large, the availability of
        the :literal:`seadata` package will speed the computation.
        <BLANKLINE>
        If the curve is defined over :math:`\mathbb{Q}`, :math:`p` must be explicitly given and the
        function computes the trace of the reduction over :math:`\mathbb{F}_p`.
        The trace of Frobenius is also the :math:`a_p` coefficient in the curve :math:`L`-series
        :math:`L(E,s) = \sum_n a_n n^{-s}`, whence the function name. The equation must be
        integral at :math:`p` but need not be minimal at :math:`p`; of course, a minimal model
        will be more efficient.
        <BLANKLINE>
        ::
        <BLANKLINE>
            ? E = ellinit([0,1]); \\ y^2 = x^3 + 0.x + 1, defined over Q
            ? ellap(E, 7) \\ 7 necessary here
            %2 = -4 \\ #E(F_7) = 7+1-(-4) = 12
            ? ellcard(E, 7)
            %3 = 12 \\ OK
        <BLANKLINE>
            ? E = ellinit([0,1], 11); \\ defined over F_11
            ? ellap(E) \\ no need to repeat 11
            %4 = 0
            ? ellap(E, 11) \\ ... but it also works
            %5 = 0
            ? ellgroup(E, 13) \\ ouch, inconsistent input!
             *** at top-level: ellap(E,13)
             *** ^-----------
             *** ellap: inconsistent moduli in Rg_to_Fp:
             11
             13
        <BLANKLINE>
            ? Fq = ffgen(ffinit(11,3), 'a); \\ defines F_q := F_{11^3}
            ? E = ellinit([a+1,a], Fq); \\ y^2 = x^3 + (a+1)x + a, defined over F_q
            ? ellap(E)
            %8 = -3
        <BLANKLINE>
        :strong:`Algorithms used.` If :math:`E/\mathbb{F}_q` has CM by a principal imaginary
        quadratic order we use a fast explicit formula (involving essentially
        Kronecker symbols and Cornacchia's algorithm), in :math:`O(\log q)^2`.
        Otherwise, we use Shanks-Mestre's baby-step/giant-step method, which runs in
        time :math:`~{O}(q^{1/4})` using :math:`~{O}(q^{1/4})` storage, hence becomes
        unreasonable when :math:`q` has about 30 digits. Above this range, the :literal:`SEA`
        algorithm becomes available, heuristically in :math:`~{O}(\log q)^4`, and
        primes of the order of 200 digits become feasible. In small
        characteristic we use Mestre's (p = 2), Kohel's (p = 3,5,7,13), Satoh-Harley
        (all in :math:`~{O}(p^{2} n^2)`) or Kedlaya's (in :math:`~{O}(p n^3)`)
        algorithms.

    ::

        sage: print(get_rest_doc("bitor"))
        bitwise (inclusive)
        :literal:`or` of two integers :math:`x` and :math:`y`, that is the integer
        <BLANKLINE>
        .. MATH::
        <BLANKLINE>
            \sum
            (x_i or y_i) 2^i
        <BLANKLINE>
        See ``bitand`` (in the PARI manual) for the behavior for negative arguments.
    """
    raw = get_raw_doc(function)
    return raw_to_rest(raw)
