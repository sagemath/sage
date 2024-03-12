# Workaround for https://github.com/mesonbuild/meson/issues/2546
fs = import('fs')
headers_copy = []
foreach header : headers
    headers_copy += fs.copyfile(header)
endforeach
