project('sagemath-doc-html',
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

html = custom_target('HTML documentation',
  output: 'index.html',
  depends: [doc_src],
  command: [meson.current_source_dir() / 'make.sh',
            '-C', meson.current_build_dir() / 'src', 'doc-html'],
  env: {'SAGE_DOC': meson.current_build_dir(),
        'SAGE_DOC_SRC': meson.current_build_dir() / 'src'},
  install: true,
  install_dir: docdir
)

install_subdir(meson.current_build_dir() / 'html',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
