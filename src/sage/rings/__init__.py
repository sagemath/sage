# Break cyclic import (integer_ring > integer > integer_ring)
import sage.rings.integer  # noqa: F401
