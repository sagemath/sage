kernel_data = configuration_data()
kernel_data.set('SAGE_VERSION', meson.project_version())

configure_file(
  input: 'kernel.json.in',
  output: 'kernel.json',
  install_dir: kernel_dir,
  install: true,
  configuration: kernel_data,
)

configure_file(
  input: 'logo.svg',
  output: 'logo.svg',
  copy: true,
  install_dir: kernel_dir,
  install: true,
)

configure_file(
  input: 'logo-64x64.png',
  output: 'logo-64x64.png',
  copy: true,
  install_dir: kernel_dir,
  install: true,
)

install_data(
  'logo.svg',
  'logo-64x64.png',
  install_dir: sage_install_dir / 'repl/ipython_kernel/kernel_spec',
)
