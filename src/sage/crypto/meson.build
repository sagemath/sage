py.install_sources(
  '__init__.py',
  'all.py',
  'cipher.py',
  'classical.py',
  'classical_cipher.py',
  'cryptosystem.py',
  'lattice.py',
  'lfsr.py',
  'lwe.py',
  'sboxes.py',
  'stream.py',
  'stream_cipher.py',
  'util.py',
  subdir: 'sage/crypto',
)

extension_data = {
  'boolean_function' : files('boolean_function.pyx'),
  'sbox' : files('sbox.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/crypto',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

install_subdir('block_cipher', install_dir: sage_install_dir / 'crypto')
install_subdir('key_exchange', install_dir: sage_install_dir / 'crypto')
install_subdir('mq', install_dir: sage_install_dir / 'crypto')
install_subdir('public_key', install_dir: sage_install_dir / 'crypto')
