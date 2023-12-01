project('sagemath-doc-pdf',
        version: files('VERSION.txt')
)

# This is a meson-python project, so this is subject to the wheel mappings
# defined in https://github.com/mesonbuild/meson-python/blob/main/mesonpy/__init__.py
docdir = get_option('datadir') / 'share' / 'doc' / 'sage'

doc_src = custom_target('Writable documentation source tree',
  output: 'doc-src',  # phony
  command: [meson.current_source_dir() / 'copy_doc_src.py'],
  env: {'SAGE_DOC': meson.current_build_dir(),
        'SAGE_DOC_SRC': meson.current_source_dir() / 'doc'}
)

pdf = custom_target('PDF documentation',
  build_by_default: true,
  output: 'doc-pdf',  # phony
  depends: [doc_src],
  command: [meson.current_source_dir() / 'make.sh',
            '-C', meson.current_build_dir() / 'src', '@OUTPUT@'],
  env: {'SAGE_DOC': meson.current_build_dir(),
        'SAGE_DOC_SRC': meson.current_build_dir() / 'src'}
)

install_subdir(meson.current_build_dir() / 'pdf',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
