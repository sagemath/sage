project('sagemath-doc-html',
  version: '10.2'
)

py = import('python').find_installation(pure: true)

html = custom_target('HTML documentation',
  output: 'html',
  command: ['env', 'SAGE_DOC=' + meson.current_build_dir(),
            'make', '-j20', '-C', join_paths(meson.current_source_dir(), 'doc'), 'doc-html'],
  install: true,
  install_dir: join_paths('share', 'doc', 'sage')
)
