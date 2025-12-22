from typing import Any

from sage.categories.map import Map
from sage.categories.morphism import Morphism
from sage.rings.integer import Integer
from sage.structure.parent import Parent

class RingMap(Morphism):
    pass

class RingMap_lift(RingMap):
    S: Parent[Any]
    to_S: Map[Any, Any]

class RingHomomorphism(RingMap):
    _lift: Morphism
    _cached_methods: dict

class RingHomomorphism_im_gens(RingHomomorphism):
    _im_gens: object
    _base_map: object

class RingHomomorphism_from_base(RingHomomorphism):
    _underlying: object

class RingHomomorphism_from_fraction_field(RingHomomorphism):
    _morphism: object

class RingHomomorphism_cover(RingHomomorphism):
    pass

class RingHomomorphism_from_quotient(RingHomomorphism):
    phi: object

class FrobeniusEndomorphism_generic(RingHomomorphism):
    _p: Integer
    _q: Integer
    _power: int
