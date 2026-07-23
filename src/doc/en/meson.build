doc_sources = ['conf.py', 'index.rst']

foreach file : doc_sources
  doc_src += fs.copyfile(file)
endforeach

subdir('reference')
subdir('constructions')
subdir('website')
subdir('developer')
subdir('thematic_tutorials')
subdir('installation')
subdir('tutorial')
subdir('a_tour_of_sage')
subdir('faq')
subdir('prep')
