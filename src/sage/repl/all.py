from sage.misc.lazy_import import lazy_import

from sage.repl.preparse import preparse, implicit_multiplication

lazy_import('sage.repl.interpreter', 'preparser')

lazy_import('sage.repl.attach', [
        'attach', 'detach', 'attached_files', 'load_attach_path',
        'reset_load_attach_path', 'load_attach_mode'])

from sage.repl.rich_output.pretty_print import pretty_print, show
