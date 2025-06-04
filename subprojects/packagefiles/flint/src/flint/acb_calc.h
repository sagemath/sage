// Workaround for https://github.com/mesonbuild/meson/issues/10298
// to make sure that imports like "#include <flint/acb_calc.h>" work.

#pragma once
#include <acb_calc.h>
