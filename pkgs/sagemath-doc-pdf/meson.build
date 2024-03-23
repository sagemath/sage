project('sagemath-doc-pdf',
  version: '10.2'
)

docdir = get_option('datadir') / 'share' / 'doc' / 'sage'

pdf = custom_target('PDF documentation',
  build_by_default: true,
  output: 'doc-pdf',  # phony
  command: ['env', 'SAGE_DOC=' + meson.current_build_dir(),
            'make', '-j20', '-C', join_paths(meson.current_source_dir(), 'doc'), '@OUTPUT@'],
)

install_subdir(meson.current_build_dir() / 'pdf',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
