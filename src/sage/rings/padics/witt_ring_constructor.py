import sage.rings.ring as ring

from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()
from sage.sets.primes import Primes
_Primes = Primes()

from sage.rings.padics.witt_ring import WittRing_p_typical
from sage.rings.padics.witt_ring import WittRing_finite_field
from sage.rings.padics.witt_ring import WittRing_p_invertible
from sage.rings.padics.witt_ring import WittRing_non_p_typical

def WittRing(base_ring, prec=1, p=None, algorithm='auto'):
    
    if not ring.is_Ring(base_ring):
        raise TypeError(f'Base ring {base_ring} must be a ring.')
    
    if base_ring not in _CommutativeRings:
        raise TypeError(f'Cannot create Ring of Witt Vectors over {base_ring}: {base_ring} is not a commutative ring.')
    
    char = base_ring.characteristic()
    prime = None
    if p is None:
        if char not in _Primes:
            raise ValueError(f'Cannot create Ring of Witt Vectors over {base_ring}: {base_ring} has non-prime characteristic and no prime was supplied.')
        else:
            prime = char
    else:
        prime = p
    
    if algorithm is None:
        algorithm = 'none'
    elif algorithm not in ['none', 'auto', 'standard', 'finotti']:
        raise ValueError(f"'{algorithm}' is not a valid algorithm. It must be one of 'none', 'auto', 'standard', or 'finotti'.")
    
    if prime == char: # p-typical
        if base_ring.is_field() and base_ring.is_finite():
            # TODO: document that this ignores the choice of algorithm
            return WittRing_finite_field(base_ring.field(), prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'finotti'
            return WittRing_p_typical(base_ring, prec, prime, 
                                      algorithm=algorithm,
                                      category=_CommutativeRings)
    else: # non-p-typical
        if algorithm == 'finotti':
            raise ValueError(f"The 'finotti' algorithm only works for p-typical Witt Rings.")
        if base_ring(prime).is_unit():
            # TODO: document that this ignores the choice of algorithm
            return WittRing_p_invertible(base_ring, prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'standard'
            return WittRing_non_p_typical(base_ring, prec, prime, 
                                          algorithm=algorithm,
                                          category=_CommutativeRings)