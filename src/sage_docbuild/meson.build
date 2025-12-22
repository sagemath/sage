# Check for system-wide MathJax installation
mathjax_path = ''
mathjax_found = false
foreach dir : ['/usr/share/mathjax', '/usr/local/share/mathjax']
  if not mathjax_found and fs.is_file(dir / 'tex-chtml.js')
    mathjax_path = dir / 'tex-chtml.js'
    mathjax_found = true
    message('Found MathJax at: ' + dir)
  endif
endforeach
if not mathjax_found
  mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js'
  message('Using CDN MathJax: ' + mathjax_path)
endif
docs_config_data = configuration_data()
docs_config_data.set('MATHJAX_PATH', mathjax_path)
configure_file(
  input: 'config.py.in',
  output: 'config.py',
  configuration: docs_config_data,
)
