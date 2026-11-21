# Resolve a cyclic import (libgap > [utils >] element > libgap)
import sage.libs.gap.element  # noqa: F401
