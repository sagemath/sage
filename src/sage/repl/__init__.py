# sage_setup: distribution = sagemath-repl
# IPython calls this when "%load_ext sage.repl" is used.
# The Sage application loads it when starting up.
def load_ipython_extension(*args):
    import sage.repl.ipython_extension
    sage.repl.ipython_extension.load_ipython_extension(*args)


# The above used to be in sage.__init__, allowing users to use "%load_ext sage".
# But we are clearing out the __init__.py file as a preparation for
# making sage a native namespace package.
#
# So we make "%load_ext sage" work by monkey-patching the function
# into the sage package upon importing sage.repl.
import sage
sage.load_ipython_extension = load_ipython_extension
