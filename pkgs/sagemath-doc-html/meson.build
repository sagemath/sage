project('sagemath-doc-html',
  version: '10.2'
)

docdir = get_option('datadir') / 'share' / 'doc' / 'sage'

html = custom_target('HTML documentation',
  output: 'index.html',
  command: [join_paths(meson.current_source_dir(), 'doc', 'make.sh'),
            '-C', join_paths(meson.current_source_dir(), 'doc'), 'doc-html'],
  env: {'SAGE_DOC': meson.current_build_dir()},
  install: true,
  install_dir: docdir
)

install_subdir(meson.current_build_dir() / 'html',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
