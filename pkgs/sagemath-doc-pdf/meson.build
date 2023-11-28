project('sagemath-doc-pdf',
        version: files('VERSION.txt')
)

# This is a meson-python project, so this is subject to the wheel mappings
# defined in https://github.com/mesonbuild/meson-python/blob/main/mesonpy/__init__.py
docdir = get_option('datadir') / 'share' / 'doc' / 'sage'

pdf = custom_target('PDF documentation',
  build_by_default: true,
  output: 'doc-pdf',  # phony
  command: [meson.current_source_dir() / 'doc' / 'make.sh',
            '-C', meson.current_source_dir() / 'doc', '@OUTPUT@'],
  env: {'SAGE_DOC': meson.current_build_dir(),
        'SAGE_DOC_SRC': meson.current_source_dir() / 'doc'}
)

install_subdir(meson.current_build_dir() / 'pdf',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
