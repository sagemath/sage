r"""
Logging of Sage sessions

.. TODO::

    Pressing "control-D" can mess up the I/O sequence because of
    a known bug.

You can create a log of your Sage session as a web page and/or as a
latex document. Just type ``log_html()`` to create an HTML log, or
``log_dvi()`` to create a dvi (LaTeX) log. Your complete session so
far up until when you type the above command will be logged, along
with any future input. Thus you can view the log system as a way to
print or view your entire session so far, along with a way to see
nicely typeset incremental updates as you work.

If ``L=log_dvi()`` or ``L=log_html()`` is a logger, you can type
``L.stop()`` and ``L.start()`` to stop and start logging.

The environment variables ``BROWSER`` and ``DVI_VIEWER`` determine
which web browser or dvi viewer is used to display your running log.

For both log systems you must have a TeX system installed on your
computer. For HTML logging, you must have the convert command, which
comes with the free ImageMagick tools.

.. note::

   The HTML output is done via LaTeX and PNG images right now,
   sort of like how latex2html works. Obviously it would be
   interesting to do something using MathML in the long run.

AUTHORS:

- William Stein (2006-02): initial version

- William Stein (2006-02-27): changed html generation so log directory
  is relocatable (no hardcoded paths).

- William Stein (2006-03-04): changed environment variable to BROWSER.

- Didier Deshommes (2006-05-06): added MathML support; refactored
  code.

- Dan Drake (2008-03-27): fix bit rotting so
  that optional directories work, dvi logging works, viewer() command
  works, remove no-longer-working MathML logger; fix off-by-one
  problems with IPython history; add text logger; improve
  documentation about viewers.
"""

# Note: there is a web browser module standard with Python.
# But it seems so dated as to be useless.

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

import os
import time

import sage.repl.interpreter as interpreter
from . import latex
from . import misc

from   sage.misc.viewer  import browser, dvi_viewer

offset = 0
loggers = []

def update():
    for X in loggers:
        X._update()

REFRESH = '<meta http-equiv="REFRESH" content="4">'

class Log:
    """
    This is the base logger class. The two classes that you actually
    instantiate are derived from this one.
    """
    def __init__(self, dir=None, debug=False, viewer=None):
        from sage.misc.misc import sage_makedirs
        if dir is None:
            dir = misc.DOT_SAGE + 'log'
        self._time = time.strftime('%Y-%m-%d-%H%M%S')
        dir = os.path.join(os.path.abspath(dir), 'log-' + self._time)
        sage_makedirs(dir)
        self._debug = debug
        self._n = 0
        self._dir = dir
        self._filename = os.path.join(dir, self._filename())
        self._output = __IPYTHON__.output_hist
        self._input  = __IPYTHON__.input_hist_raw
        self._text = ''
        self._after_output = False
        self._init()
        self._input_text  = ''
        self._stopped = False
        self._viewer = viewer
        loggers.append(self)
        print('Now logging to ' + self._filename)
        self._update()
        self.view()

    def __repr__(self):
        return "Logger"

    def _latex_(self):
        return "\\mathrm{%s}"%self

    def dir(self):
        """
        Return the directory that contains the log files.
        """
        return self._dir

    def stop(self):
        """
        Stop the logger. To restart use the start function.
        """
        self._stopped = True

    def start(self):
        """
        Start the logger. To stop use the stop function.
        """
        self._stopped = False

    def _update(self):
        """
        There is an off-by-one issue with IPython's input and output
        history; ``__IPYTHON__.input_hist_raw`` is a *list* containing
        the un-preparsed Sage commands. However,
        ``__IPYTHON__.output_hist`` is a dictionary whose keys are
        integers and whose values are outputs.  This is good because
        not every input has an output.

        **BUT**, the output from::

            __IPYTHON__.input_hist_raw[n]

        is stored in::

            __IPYTHON__.output_hist[n+1] !

        This is annoying and it may be a bug. Right now the loggers
        correct for this, but if modifying or extending this code,
        consider yourself warned.
        """
        if self._stopped:
            return
        # see note at end of this function for info about output and
        # input
        (O, I) = (self._output, self._input)
        K = O.keys()
        while self._n < max(len(I), max(K + [-1])):
            n = self._n
            followed_by_output = (n+1 in K)
            if n < len(I):
                L = I[n]
            else:
                L = "Sorry; raw input log out of sync"
            Lstrip = L.strip()
            if len(Lstrip) > 0 and Lstrip[0] != '%':
                self._write(self._get_input(n, followed_by_output))
                self._input_text += L
            #m = n + offset
            m = n
            if m+1 in K:
                self._write(self._get_output(m+1))
                # what does the following line do? Nothing is done with
                # this s. Commenting out for now.
                #s = '# ' + '\n# '.join(str(O[m]).split('\n')) + '\n\n'
            self._n += 1
        A = open(self._filename,'w')
        A.write(self._header() + '\n' + self._text + '\n' + self._footer())
        A.close()
        self._update_plain()
        self._build()

    def _write(self, lines):
        self._text += lines

    def _plain_text(self):
        return self._input_text

    def _input_log_name(self):
        return os.path.join(self._dir, 'input-' + self._time)

    def _update_plain(self):
        open(self._input_log_name(),'w').write(self._input_text)


class log_html(Log):
    r"""
    Create a running log of your Sage session as a web page.

    Easy usage: ``log_html()``

    TODO: Pressing "control-D" can mess up the I/O sequence because of
    a known bug.

    Use ``L=log_html([optional directory])`` to create an HTML
    log. Your complete session so far up until when you type the above
    command will be logged, along with any future input. Thus you can
    view the log system as a way to print or view your entire session
    so far, along with a way to see nicely typeset incremental updates
    as you work.

    If L is a logger, you can type ``L.stop()`` and
    ``L.start()`` to stop and start logging.

    The environment variable ``WEB_BROWSER`` determines which web
    browser or dvi viewer is used to display your running log. You can
    also specify a viewer when you start the logger with something
    like ``log_html([opt. dir], viewer='firefox')``.

    You must have a TeX system installed on your computer, and you
    must have the convert command, which comes with the free
    ImageMagick tools.
    """

    def _init(self):
        self._images = os.path.join(self._dir, 'images')
        if not os.path.exists(self._images):
            os.makedirs(self._images)


    def view(self):
        if not self._viewer is None:
            viewer = self._viewer
        else:
            viewer = browser()
        os.system('%s "%s"&'%(viewer, self._filename))

    def _build(self):
        return

    def __repr__(self):
        return "HTML Logger"

    def _get_input(self, n, followed_by_output):
        if n >= len(self._input):
            return
        return """<font color=darkblue>     %s %s:</font> %s"""%(
            n, interpreter._prompt, self._input[n])

    def _get_output(self, n):
        x = self._output[n]
        try:
            L = latex.latex(x)
        except Exception:
            L = "\\mbox{error TeXing object}"
        single_png = os.path.join(self._images, '%s.png' % n)
        try:
            x.png(single_png)
        except AttributeError:
            latex.png(x, single_png, debug=self._debug)
        oi = os.path.join(self._dir, 'images', 'o' + '%s.html' % n)
        open(oi,'w').write('<pre>OUTPUT:\n%s\n\n\nLATEX:\n%s</pre><img src="%s">'%(
            x, L, single_png))
        extra_img_opts = ''
        #if sage.plot.graphics.is_Graphics(x):
        #    extra_img_opts = 'width=300'
        return """<center> <table border=0 cellpadding=20 cellspacing=2
                bgcolor=lightgrey>
               <tr><td bgcolor=white>
               <a href="%s">
               <img src="%s" alt="%s" %s>
                </a>
             </td></tr></table> </center>\n<hr>\n"""%('images/o%s.html'%n,
                                                      'images/%s.png'%n, L,
                                                      extra_img_opts)

    def _filename(self):
        return 'index.html'

    def _header(self):
        T = self._title()
        inlog = os.path.split(self._input_log_name())[1]
        return '<html>\n%s\n<title>%s</title>\n<body><h1 align=center>%s</h1>\n<h2 align=center><a href="%s">%s</a></h2><pre>'%(REFRESH,
            T,T, inlog, inlog)

    def _footer(self):
        return "</pre>\n</body>\n</html>"

    def _title(self):
        return 'Sage Log %s'%self._time


class log_dvi(Log):
    """
    Create a running log of your Sage session as a nicely typeset dvi
    file.

    Easy usage: ``log_dvi()``

    TODO: Pressing "control-D" can mess up the I/O sequence because of
    a known bug.

    Use ``L=log_dvi([optional directory])`` to create a dvi log. Your
    complete session so far up until when you type the above command
    will be logged, along with any future input. Thus you can view the
    log system as a way to print or view your entire session so far,
    along with a way to see nicely typeset incremental updates as you
    work.

    If L is a logger, you can type ``L.stop()`` and ``L.start()`` to
    stop and start logging.

    The environment variable ``DVI_VIEWER`` determines which web
    browser or dvi viewer is used to display your running log. You can
    also specify a viewer when you start the logger with something
    like ``log_dvi([opt. dir], viewer='xdvi')``.

    You must have a LaTeX system installed on your computer and a dvi
    viewer.
    """
    def _init(self):
        self._in_verbatim = False

    def __repr__(self):
        return "DVI Logger"

    def _build(self):
        cmd = 'cd %s; latex \\\\nonstopmode \\\\input %s '%(
            self._dir, self._filename)
        cmd += " 2>/dev/null 1>/dev/null"
        dvifile = '%s.dvi'%os.path.splitext(self._filename)[0]
        if os.path.exists(dvifile):
            cmd += ' & '
        os.system(cmd)

    def view(self):
        if not self._viewer is None:
            viewer = self._viewer
        else:
            viewer = dvi_viewer()
        self._build()
        F = os.path.splitext(self._filename)[0] + '.dvi'
        cmd = 'cd %s; %s %s '%(
            self._dir, viewer, F)
        os.system(cmd + " 2>/dev/null 1>/dev/null &")

    def _get_input(self, n, followed_by_output):
        if n >= len(self._input):
            return
        s = ''
        if not self._in_verbatim:
            s += '\\begin{verbatim}'
            self._in_verbatim = True
        I = self._input[n]
        #print('input: %s' % I)
        s += "%s %s: %s"%(n,  interpreter._prompt, I)
        s += '\\end{verbatim}'
        if followed_by_output:
            self._in_verbatim = False
        else:
            s += '\\begin{verbatim}'
        self._after_output = False
        return s

    def _get_output(self, n):
        s = ''
        self._after_output = True
        if self._in_verbatim:
            s += '\\end{verbatim}\n'
            self._in_verbatim = False

        L = latex.latex(self._output[n])
        # If we explicitly ask for LaTeX output,
        # no need to format L
        if "latex" in self._input[n-1]:
            s+= "\\begin{verbatim}%s\\end{verbatim}" %L
        else:
            s += '\n\\begin{center}$\\displaystyle %s $\\end{center}\n'%L
        return s

    def _filename(self):
        return 'sagelog.tex'

    def _header(self):
        return """
\\documentclass{article}
\\title{%s}\\author{}
\\begin{document}
\\maketitle
"""%self._title()

    def _footer(self):
        if self._in_verbatim:
            return r"\end{verbatim}\end{document}"
        else:
            return r"\end{document}"


    def _title(self):
        return '\\SAGE Log %s'%self._time


class log_text(Log):
    """
    Create a running log of your Sage session as a plain text file.

    Easy usage: ``log_text()``

    TODO: Pressing "control-D" can mess up the I/O sequence because of
    a known bug.

    Use ``L=log_text([optional directory])`` to create a text
    log. Your complete session so far up until when you type the above
    command will be logged, along with any future input. Thus you can
    view the log system as a way to print or view your entire session
    so far, along with a way to see incremental updates as you work.

    Unlike the html and dvi loggers, this one does not automatically
    start a viewer unless you specify one; you can do that when you
    start the logger with something like
    ``log_text([opt. dir], viewer='xterm -e tail -f')``.

    If L is a logger, you can type ``L.stop()`` and
    ``L.start()`` to stop and start logging.
    """
    def _init(self):
        return

    def __repr__(self):
        return "Text Logger"

    def _build(self):
        return

    def view(self):
        if not self._viewer is None:
            viewer = self._viewer
        else:
            return
        cmd = 'cd %s; %s %s ' % (self._dir, viewer, self._filename)
        os.system(cmd + " 2>/dev/null 1>/dev/null &")

    def _get_input(self, n, followed_by_output):
        if n >= len(self._input):
            return
        else:
            return "%s %s: %s"%(n, interpreter._prompt, self._input[n])

    def _get_output(self, n):
        return '\n  ' + str(self._output[n]) + '\n\n'

    def _filename(self):
        return 'sagelog.txt'

    def _header(self):
        return self._title()

    def _footer(self):
        return ''

    def _title(self):
        return 'Sage Log %s' % self._time
