project('sagemath-doc-pdf',
        version: run_command('cat',
                             join_paths(meson.current_source_dir(), 'VERSION.txt'),
                             check: true).stdout().strip()
)

docdir = get_option('datadir') / 'share' / 'doc' / 'sage'

pdf = custom_target('PDF documentation',
  build_by_default: true,
  output: 'doc-pdf',  # phony
  command: [join_paths(meson.current_source_dir(), 'doc', 'make.sh'),
            '-C', join_paths(meson.current_source_dir(), 'doc'), '@OUTPUT@'],
  env: {'SAGE_DOC': meson.current_build_dir()}
)

install_subdir(meson.current_build_dir() / 'pdf',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
