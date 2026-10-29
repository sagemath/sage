py.install_sources(
  '__init__.py',
  'all.py',
  'attach.py',
  'configuration.py',
  'image.py',
  'inputhook.py',
  'interface_magic.py',
  'interpreter.py',
  'ipython_extension.py',
  'ipython_tests.py',
  'load.py',
  'preparse.py',
  'prompts.py',
  'user_globals.py',
  subdir: 'sage/repl',
)

install_subdir('display', install_dir: sage_install_dir / 'repl')
install_subdir('rich_output', install_dir: sage_install_dir / 'repl')
subdir('ipython_kernel')
