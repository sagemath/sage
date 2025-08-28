"Evaluating Python code without any preparsing"


class Python:
    """
    Allow for evaluating a chunk of code without any preparsing.
    """

    def eval(self, x, globals, locals=None):
        r"""
        Evaluate x with given globals; locals is completely ignored.

        This is specifically meant for evaluating code blocks with
        ``%python`` in the notebook.

        INPUT:

        - ``x`` -- string
        - ``globals`` -- dictionary
        - ``locals`` -- completely IGNORED

        EXAMPLES::

            sage: from sage.misc.python import Python
            sage: python = Python()
            sage: python.eval('2+2', globals())
            4
            ''

        Any variables that are set during evaluation of the block
        will propagate to the globals dictionary. ::

            sage: python.eval('a=5\nb=7\na+b', globals())
            12
            ''
            sage: b
            7

        The ``locals`` variable is ignored -- it is there only for
        completeness.  It is ignored since otherwise the following
        will not work::

            sage: python.eval("def foo():\n   return 'foo'\nprint(foo())\ndef mumble():\n    print('mumble {}'.format(foo()))\nmumble()", globals())
            foo
            mumble foo
            ''
            sage: mumble
            <function mumble at ...>
        """
        x = x.strip()
        y = x.split('\n')
        if not y:
            return ''
        s = '\n'.join(y[:-1]) + '\n'
        t = y[-1]
        try:
            z = compile(t + '\n', '', 'single')
        except SyntaxError:
            s += '\n' + t
            z = None

        eval(compile(s, '', 'exec'), globals, globals)

        if z is not None:
            eval(z, globals)
        return ''


python = Python()
