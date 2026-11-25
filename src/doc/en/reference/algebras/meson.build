doc_sources = [
  'conf.py',
  'cubic_hecke_algebra.rst',
  'fusion_rings.rst',
  'index.rst',
  'lie_algebras.rst',
  'lie_conformal_algebras.rst',
  'quantum_groups.rst',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('media')
