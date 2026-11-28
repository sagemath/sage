# Build configuration of sage_docbuild.

# sage_docbuild is actually not installed by Meson and runs from the source tree,
# so we have to work in the source tree directly.

docbuild_conf_data = configuration_data()

mathjax_dirs = []
if sage_local != ''
  mathjax_dirs += sage_local / 'share' / 'mathjax' / 'mathjax'
endif
mathjax_dirs += [
  # conda-forge
  prefix / 'lib' / 'mathjax',
  prefix / 'Library' / 'lib' / 'mathjax',
  datadir / 'mathjax',
  # Arch, Gentoo, Void, ...
  '/usr/share/mathjax',
]
mathjax_dir = ''
foreach dir : mathjax_dirs
  if mathjax_dir == '' and fs.is_file(dir / 'tex-chtml.js')
    mathjax_dir = dir
  endif
endforeach

if mathjax_dir == ''
  message(
    'MathJax 3 not found. The documentation will load MathJax from a CDN.',
  )
else
  message('Using MathJax from ' + mathjax_dir)
endif
docbuild_conf_data.set_quoted('SAGE_MATHJAX_DIR', mathjax_dir)

docbuild_config = configure_file(
  input: 'config.py.in',
  output: 'config.py',
  configuration: docbuild_conf_data,
)
run_command(
  py,
  '-c',
  'import shutil, sys; shutil.copyfile(sys.argv[1], sys.argv[2])',
  docbuild_config,
  meson.current_source_dir() / 'config.py',
  check: true,
)
