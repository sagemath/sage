

def arithmetic(t=None):
    """
    Controls the default proof strategy for integer arithmetic algorithms
    (such as primality testing).

    INPUT:

    t -- boolean or ``None``

    OUTPUT:

    If t is ``True``, requires integer arithmetic operations to (by
    default) return results that are true unconditionally: the
    correctness will not depend on an algorithm with a nonzero
    probability of returning an incorrect answer or on the truth of
    any unproven conjectures.

    If t is ``False``, allows integer arithmetic operations to (by
    default) return results that may depend on unproven conjectures or
    on probabilistic algorithms.  Such algorithms often have a
    substantial speed improvement over those requiring proof.

    If t is ``None``, returns the integer arithmetic proof status.

    EXAMPLES::

        sage: proof.arithmetic()
        True
        sage: proof.arithmetic(False)
        sage: proof.arithmetic()
        False
        sage: proof.arithmetic(True)
        sage: proof.arithmetic()
        True
    """
    from .proof import _proof_prefs
    return _proof_prefs.arithmetic(t)


def elliptic_curve(t=None):
    """
    Controls the default proof strategy for elliptic curve algorithms.

    INPUT:

    t -- boolean or ``None``

    OUTPUT:

    If t is ``True``, requires elliptic curve algorithms to (by
    default) return results that are true unconditionally: the
    correctness will not depend on an algorithm with a nonzero
    probability of returning an incorrect answer or on the truth of
    any unproven conjectures.

    If t is ``False``, allows elliptic curve algorithms to (by
    default) return results that may depend on unproven conjectures or
    on probabilistic algorithms.  Such algorithms often have a
    substantial speed improvement over those requiring proof.

    If t is ``None``, returns the current elliptic curve proof status.

    EXAMPLES::

        sage: proof.elliptic_curve()
        True
        sage: proof.elliptic_curve(False)
        sage: proof.elliptic_curve()
        False
        sage: proof.elliptic_curve(True)
        sage: proof.elliptic_curve()
        True
    """
    from .proof import _proof_prefs
    return _proof_prefs.elliptic_curve(t)


def linear_algebra(t=None):
    """
    Controls the default proof strategy for linear algebra algorithms.

    INPUT:

    t -- boolean or ``None``

    OUTPUT:

    If t is ``True``, requires linear algebra algorithms to (by
    default) return results that are true unconditionally: the
    correctness will not depend on an algorithm with a nonzero
    probability of returning an incorrect answer or on the truth of
    any unproven conjectures.

    If t is ``False``, allows linear algebra algorithms to (by
    default) return results that may depend on unproven conjectures or
    on probabilistic algorithms.  Such algorithms often have a
    substantial speed improvement over those requiring proof.

    If t is ``None``, returns the current linear algebra proof status.

    EXAMPLES::

        sage: proof.linear_algebra()
        True
        sage: proof.linear_algebra(False)
        sage: proof.linear_algebra()
        False
        sage: proof.linear_algebra(True)
        sage: proof.linear_algebra()
        True
    """
    from .proof import _proof_prefs
    return _proof_prefs.linear_algebra(t)


def number_field(t=None):
    """
    Controls the default proof strategy for number field algorithms.

    INPUT:

    t -- boolean or ``None``

    OUTPUT:

    If t is ``True``, requires number field algorithms to (by default)
    return results that are true unconditionally: the correctness will
    not depend on an algorithm with a nonzero probability of returning
    an incorrect answer or on the truth of any unproven conjectures.

    If t is ``False``, allows number field algorithms to (by default)
    return results that may depend on unproven conjectures or on
    probabilistic algorithms.  Such algorithms often have a
    substantial speed improvement over those requiring proof.

    If t is ``None``, returns the current number field proof status.

    EXAMPLES::

        sage: proof.number_field()
        True
        sage: proof.number_field(False)
        sage: proof.number_field()
        False
        sage: proof.number_field(True)
        sage: proof.number_field()
        True
    """
    from .proof import _proof_prefs
    return _proof_prefs.number_field(t)


def polynomial(t=None):
    """
    Controls the default proof strategy for polynomial algorithms.

    INPUT:

    t -- boolean or ``None``

    OUTPUT:

    If t is ``True``, requires polynomial algorithms to (by default)
    return results that are true unconditionally: the correctness will
    not depend on an algorithm with a nonzero probability of returning
    an incorrect answer or on the truth of any unproven conjectures.

    If t is ``False``, allows polynomial algorithms to (by default)
    return results that may depend on unproven conjectures or on
    probabilistic algorithms.  Such algorithms often have a
    substantial speed improvement over those requiring proof.

    If t is ``None``, returns the current polynomial proof status.

    EXAMPLES::

        sage: proof.polynomial()
        True
        sage: proof.polynomial(False)
        sage: proof.polynomial()
        False
        sage: proof.polynomial(True)
        sage: proof.polynomial()
        True
    """
    from .proof import _proof_prefs
    return _proof_prefs.polynomial(t)


def all(t=None):
    """
    Controls the default proof strategy throughout Sage.

    INPUT:

    t -- boolean or ``None``

    OUTPUT:

    If t is ``True``, requires Sage algorithms to (by default) return
    results that are true unconditionally: the correctness will not
    depend on an algorithm with a nonzero probability of returning an
    incorrect answer or on the truth of any unproven conjectures.

    If t is ``False``, allows Sage algorithms to (by default) return
    results that may depend on unproven conjectures or on
    probabilistic algorithms.  Such algorithms often have a
    substantial speed improvement over those requiring proof.

    If t is ``None``, returns the current global Sage proof status.

    EXAMPLES::

        sage: proof.all()
        {'arithmetic': True,
         'elliptic_curve': True,
         'linear_algebra': True,
         'number_field': True,
         'other': True,
         'polynomial': True}
        sage: proof.number_field(False)
        sage: proof.number_field()
        False
        sage: proof.all()
        {'arithmetic': True,
         'elliptic_curve': True,
         'linear_algebra': True,
         'number_field': False,
         'other': True,
         'polynomial': True}
        sage: proof.number_field(True)
        sage: proof.number_field()
        True
    """
    from .proof import _proof_prefs
    if t is None:
        return _proof_prefs._require_proof.copy()
    for s in _proof_prefs._require_proof:
        _proof_prefs._require_proof[s] = bool(t)


from .proof import WithProof
