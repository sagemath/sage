info: stand-alone Info documentation reader
===========================================

Description
-----------

GNU Info is the stand-alone "info" reader that is part of the GNU
Texinfo suite of tools. Several packages (Maxima, Singular, ...)
install documentation in "info" format, which can be read either
with Emacs, the stand-alone "info" reader, and some other software.
In particular, the interactive help system of ``singular_console()``
uses the ``info`` program in environments in which a web browser is
not available; if ``info`` is not installed, it falls back to a
basic pager with limited capabilities.

Website: https://www.gnu.org/software/texinfo/manual/info-stnd/info-stnd.html


License
-------

GPL-3+ (``info/*.c`` comments in the source repository)
