py.install_sources(
  'abstract_word.py',
  'all.py',
  'alphabet.py',
  'finite_word.py',
  'infinite_word.py',
  'lyndon_word.py',
  'morphic.py',
  'morphism.py',
  'paths.py',
  'shuffle_product.py',
  'suffix_trees.py',
  'word.py',
  'word_generators.py',
  'word_infinite_datatypes.py',
  'word_options.py',
  'words.py',
  subdir: 'sage/combinat/words',
)

extension_data = {
  'word_char' : files('word_char.pyx'),
  'word_datatypes' : files('word_datatypes.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/words',
    install: true,
    include_directories: [inc_cpython, inc_data_structures, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

