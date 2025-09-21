py.install_sources(
  '__init__.py',
  'abstract_code.py',
  'ag_code.py',
  'all.py',
  'bch_code.py',
  'binary_code.pxd',
  'bounds_catalog.py',
  'channel.py',
  'channels_catalog.py',
  'code_bounds.py',
  'code_constructions.py',
  'codes_catalog.py',
  'cyclic_code.py',
  'databases.py',
  'decoder.py',
  'decoders_catalog.py',
  'delsarte_bounds.py',
  'encoder.py',
  'encoders_catalog.py',
  'extended_code.py',
  'gabidulin_code.py',
  'golay_code.py',
  'goppa_code.py',
  'grs_code.py',
  'guava.py',
  'hamming_code.py',
  'information_set_decoder.py',
  'linear_code.py',
  'linear_code_no_metric.py',
  'linear_rank_metric.py',
  'parity_check_code.py',
  'punctured_code.py',
  'reed_muller_code.py',
  'self_dual_codes.py',
  'subfield_subcode.py',
  'two_weight_db.py',
  subdir: 'sage/coding',
)

extension_data = {
  'ag_code_decoders' : files('ag_code_decoders.pyx'),
  'binary_code' : files('binary_code.pyx'),
  'kasami_codes' : files('kasami_codes.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/coding',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

subdir('codecan')
install_subdir('guruswami_sudan', install_dir: sage_install_dir / 'coding')
install_subdir('source_coding', install_dir: sage_install_dir / 'coding')
