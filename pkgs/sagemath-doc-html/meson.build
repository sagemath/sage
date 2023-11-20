project('sagemath-doc-html',
        version: files('VERSION.txt')
)

# This is a meson-python project, so this is subject to the wheel mappings
# defined in https://github.com/mesonbuild/meson-python/blob/main/mesonpy/__init__.py
docdir = get_option('datadir') / 'share' / 'doc' / 'sage'

html = custom_target('HTML documentation',
  output: 'index.html',
  command: [meson.current_source_dir() / 'doc' / 'make.sh',
            '-C', meson.current_source_dir() / 'doc', 'doc-html'],
  env: {'SAGE_DOC': meson.current_build_dir()},
  install: true,
  install_dir: docdir
)

install_subdir(meson.current_build_dir() / 'html',
               install_dir: docdir)

meson.add_dist_script('meson_dist_script.sh')
