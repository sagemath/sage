py.install_sources(
  'all.py',
  'hexad.py',
  'quantumino.py',
  'sudoku.py',
  subdir: 'sage/games',
)

extension_data = {'sudoku_backtrack' : files('sudoku_backtrack.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/games',
    install: true,
    include_directories: [],
    dependencies: [py_dep, gmp],
  )
endforeach

