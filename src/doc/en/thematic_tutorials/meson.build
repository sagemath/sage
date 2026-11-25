doc_sources = [
  'algebraic_combinatorics.rst',
  'coding_theory.rst',
  'coercion_and_categories.rst',
  'conf.py',
  'cython_interface.rst',
  'functional_programming.rst',
  'geometry.rst',
  'group_theory.rst',
  'index.rst',
  'lie.rst',
  'linear_programming.rst',
  'numtheory_rsa.rst',
  'profiling.rst',
  'sandpile.rst',
  'steenrod_algebra_modules.rst',
  'structures_in_coding_theory.rst',
  'sws2rst.rst',
  'toctree.rst',
  'tutorial-comprehensions.rst',
  'tutorial-implementing-algebraic-structures.rst',
  'tutorial-objects-and-classes.rst',
  'tutorial-programming-python.rst',
  'vector_calculus.rst',
]

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('geometry')
subdir('vector_calculus')
subdir('explicit_methods_in_number_theory')
subdir('lie')
subdir('numerical_sage')
subdir('media')
subdir('algebraic_combinatorics')
