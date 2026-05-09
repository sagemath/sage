# Resolve a cyclic import (categories.map > structure.element > ... > categories.map)
import sage.structure.element  # noqa: F401
