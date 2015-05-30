"""
Interactive Debugger for the Sage Notebook

Start the debugger in the Sage Notebook by running the ``debug()``
function.  You can run several debuggers at once.

AUTHOR:

- William Stein (2012)
"""
# Below all tests are done using sage0, which is a pexpect interface
# to Sage itself.  This allows us to test exploring a stack traceback
# using the doctest framework.

def test_function2(a, b):
    """
    Used for doctesting the notebook debugger.

    EXAMPLES::

        >>> sage.interacts.debugger.test_function2(2, 3)  # using normal prompt would confuse tests below.
        (5, 6, True, False)
    """
    x = a + b
    y = a * b
    return x, y, x<y, x>y   # < to ensure HTML is properly escaped

def test_function(n, m,level=10):
    """
    Used for doctesting the notebook debugger.

    EXAMPLES::

        >>> sage.interacts.debugger.test_function(2, 3)
        (5, 6, True, False)
    """
    # call another function so the stack is bigger
    if level > 0:
        return test_function(n,m,level=level-1)
    else:
        return test_function2(m, n)


class Debug:
    """
    Create a debugger for the most recent stack trace.

    NOTES:

    - Input is not preparsed.
    - You can define and work with many debug interacts at the same time.

    TESTS::

    The current position in the stack frame is self._curframe_index::

        sage: a = sage0.eval("sage.interacts.debugger.test_function('n', 'm')")
        sage: d = sage0('sage.interacts.debugger.Debug()')
        sage: d._curframe_index
        8
    """
    def __init__(self):
        """
        Create the debugger object from the most recent traceback.

        TESTS::

            sage: a = sage0.eval("sage.interacts.debugger.test_function('n', 'm')")
            sage: sage0('sage.interacts.debugger.Debug()')
            <sage.interacts.debugger.Debug instance at 0x...>
        """
        import inspect, sys, traceback
        try:
            tb=sys.last_traceback
            #we strip off the 5 outermost frames, since those relate only to
            #the notebook, not user code
            for i in xrange(5):
                tb=tb.tb_next
            self._stack = inspect.getinnerframes(tb)
        except AttributeError:
            raise RuntimeError("no traceback has been produced; nothing to debug")
        self._curframe_index = len(self._stack) - 1

    def curframe(self):
        """
        Return the current frame object.  This defines the local and
        global variables at a point in the stack trace, and the code.

        OUTPUT:

        - frame object

        TESTS::

            sage: a = sage0.eval("sage.interacts.debugger.test_function('n', 'm')")
            sage: d = sage0('sage.interacts.debugger.Debug()')
            sage: d.curframe()
            <frame object at 0x...>
        """
        return self._stack[self._curframe_index][0]

    def evaluate(self, line):
        """
        Evaluate the input string ``line`` in the scope of the current
        position in the stack.

        INPUT:

        - ``line`` -- string; the code to exec

        OUTPUT:

        - string (the output)

        TESTS::

             sage: _ = sage0.eval("sage.interacts.debugger.test_function('n', 'm')")
             sage: _ = sage0.eval('d = sage.interacts.debugger.Debug()')
             sage: sage0.eval("d.evaluate('print a, b')")
             'm n'
       """
        locals = self.curframe().f_locals
        globals = self.curframe().f_globals
        try:
            code = compile(line + '\n', '<stdin>', 'single')
            exec(code, globals, locals)
        except Exception:
            import sys
            t, v = sys.exc_info()[:2]
            if isinstance(t, type('')):
                exc_type_name = t
            else:
                exc_type_name = t.__name__
            print '***', exc_type_name + ':', v

    def listing(self, n=5):
        """
        Return HTML display of the lines (with numbers and a pointer
        at the current line) in the code listing for the current
        frame, with `n` lines before and after of context.

        INPUT:

        - `n` -- integer (default: 5)

        OUTPUT:

        - list of strings

        TESTS::

             sage: _ = sage0.eval("sage.interacts.debugger.test_function('n', 'm')")
             sage: _ = sage0.eval('d = sage.interacts.debugger.Debug()')
             sage: print sage0("d.listing(1)")
                 2...      x = a + b
             --&gt; ...      y = a * b
                 ...      return x, y, x&lt;y, x&gt;y   # &lt; to ensure HTML is properly escaped
             <hr>> <a href="/src/interacts/debugger.py" target="_new">src/sage/interacts/debugger.py</a>
             sage: print sage0("d.listing()")
                 2...
                 ...
                 ...      x = a + b
             --&gt; ...      y = a * b
                 ...      return x, y, x&lt;y, x&gt;y   # &lt; to ensure HTML is properly escaped
                 ...
             sage: _ = sage0.eval('d._curframe_index -= 1')
             sage: print sage0("d.listing(1)")
                 4...       else:
             --&gt; ...      test_function2(m, n)
                 ...
                 <hr>> <a href="/src/interacts/debugger.py" target="_new">src/sage/interacts/debugger.py</a>
        """
        # TODO: Currently, just as with ipdb on the command line,
        # there is no support for displaying code in Cython files.
        # This shouldn't be too hard to add.
        curframe = self.curframe()
        filename = curframe.f_code.co_filename
        lineno = curframe.f_lineno
        import linecache
        w = []
        for i in range(lineno-n, lineno+n+1):
            z = linecache.getline(filename, i, curframe.f_globals)
            if z: w.append(('--> ' if i ==lineno else '    ') + '%-5s'%i + z)
        code = ''.join(w)
        if not code.strip():
            code = '(code not available)'

        # This is a hideous hack to get around how the notebook "works".
        # If the output of anything contains the string TRACEBACK then
        # it will get mangled.  So we replace TRACEBACK in our code block
        # by the harmless version with the colon missing.  This sucks.
        from sagenb.notebook.cell import TRACEBACK
        code = code.replace(TRACEBACK, TRACEBACK[:-1])

        # Create a hyperlink to the file, if possible.
        i = filename.rfind('site-packages/sage')
        if i != -1:
            fname = filename[i+len('site-packages/sage')+1:].rstrip('/')
            file = '<a href="/src/%s" target="_new">src/sage/%s</a>'%(fname,fname)
        else:
            file = filename

        import cgi
        t = """%s<hr>> %s"""%(cgi.escape(code), file)
        return t

    def interact(self):
        """
        Start the interact debugger.

        TESTS::

             sage: _ = sage0.eval("sage.interacts.debugger.test_function('n', 'm')")
             sage: _ = sage0.eval('d = sage.interacts.debugger.Debug()')
             sage: _ = sage0.eval('d.interact()')  # only works in the notebook
        """
        # We use a library_interact instead of a normal interact here,
        # since this is an interact in the library, and a normal
        # "@interact" is all mangled.

        from sage.interacts.library import library_interact
        from sagenb.notebook.interact import slider, input_box, selector

        # self._last holds the last state of all controls.  This allows
        # us to deduce which control changed to cause the update, or that
        # nothing changed, in which case we assume the user requested to
        # re-evaluate the input box (for some reason -- currently there is
        # no point in doing so).  It is a shortcoming of @interact that
        # we have to do this.
        self._last = None

        # two sliders and a box to put in commands with an evaluate button.
        @library_interact
        def dbg(frame = slider(vmin=0, vmax=len(self._stack)-1, step_size=1, default=len(self._stack)-1, label='stack frame'),
                lines = slider(vmin=3, vmax=99, step_size=2, default=11, label='lines of context'),
                command = input_box("", label="", type=str),
                button = selector(['Evaluate'], label='', buttons=True)
                ):

            if self._last is None:
                self._last = {'command':command, 'button':button, 'lines':lines, 'frame':frame}

            if self._last['lines'] != lines:
                # they dragged the number-of-lines slider, so done
                pass
            elif self._last['command'] != command and command.strip():
                # they changed the command, so evaluate that
                self.evaluate(command)
            elif self._last['frame'] != frame:
                # they dragged the frame slider.
                self._curframe_index = frame
            elif command:
                # must have hit the evaluate button
                self.evaluate(command)

            print '<html><hr>' + self.listing(lines//2) + '</html>'
            # save control state for next time around
            self._last = {'command':command, 'button':button, 'lines':lines, 'frame':frame}

        dbg()

def debug():
    """
    If you get a traceback in the Sage notebook, use the ``debug()``
    command to start an interactive debugger.  Using %debug on the
    command line.

    Using the debugger you can move up and down the stack frame and
    evaluate arbitrary code anywhere in the call stack.

    .. warning::

       If you type in a command, execute it, switch to a different
       frame and press enter again, nothing happens: it doesn't
       execute the command in the new frame. You can get the command
       to execute by hitting the Evaluate button though.  This is a
       fundamental limitation of the Sage notebook in Sage-5.0.

    EXAMPLES::

        sage: debug()        # only works in the notebook
        You should use %debug on the command line.
    """
    # "EMBEDDED_MODE" is True precisely when the Sage notebook is running.
    from sage.plot.plot import EMBEDDED_MODE
    if not EMBEDDED_MODE:
        # Must be the command line, so suggest using the IPython debugger.
        print "You should use %debug on the command line."
    else:
        # Create the Debug object and make it interactive.
        Debug().interact()
