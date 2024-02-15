# Workaround for https://github.com/mesonbuild/meson/issues/2546
fs = import('fs')
foreach header : headers
    fs.copyfile(header)
endforeach
