# sage_setup: distribution = sagemath-objects
"Global proof preferences"

from sage.structure.sage_object import SageObject


class _ProofPref(SageObject):
    """
    An object that holds global proof preferences.  For now these are merely
    boolean flags for various parts of Sage that use probabilistic
    algorithms.

    A ``True`` flag means that the subsystem (such as linear algebra or number
    fields) should return results that are true unconditionally: the
    correctness should not depend on an algorithm with a nonzero probability of
    returning an incorrect answer or on the truth of any unproven conjectures.

    A ``False`` flag means that the subsystem can use faster methods to return
    answers that have a very small probability of being wrong.
    """
    def __init__(self, proof=True):
        self._require_proof = {}
        self._require_proof["arithmetic"] = proof
        self._require_proof["elliptic_curve"] = proof
        self._require_proof["linear_algebra"] = proof
        self._require_proof["number_field"] = proof
        self._require_proof["polynomial"] = proof
        self._require_proof["other"] = proof

    def arithmetic(self, t=None):
        """
        Controls the default proof strategy for integer arithmetic algorithms (such as primality testing).

        INPUT:

        - ``t`` -- boolean or ``None``

        OUTPUT:

        - if ``t == True``, requires integer arithmetic operations to (by
          default) return results that are true unconditionally: the
          correctness will not depend on an algorithm with a nonzero
          probability of returning an incorrect answer or on the truth of any
          unproven conjectures.
        - if ``t == False``, allows integer arithmetic operations to (by
          default) return results that may depend on unproven conjectures or on
          probabilistic algorithms.  Such algorithms often have a substantial
          speed improvement over those requiring proof.
        - if ``t == None``, returns the integer arithmetic proof status.

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
        if t is None:
            return self._require_proof["arithmetic"]
        self._require_proof["arithmetic"] = bool(t)

    def elliptic_curve(self, t=None):
        """
        Controls the default proof strategy for elliptic curve algorithms.

        INPUT:

        - ``t`` -- boolean or ``None``

        OUTPUT:

        - if ``t == True``, requires elliptic curve algorithms to (by default)
          return results that are true unconditionally: the correctness will
          not depend on an algorithm with a nonzero probability of returning an
          incorrect answer or on the truth of any unproven conjectures.
        - if ``t == False``, allows elliptic curve algorithms to (by default)
          return results that may depend on unproven conjectures or on
          probabilistic algorithms.  Such algorithms often have a substantial
          speed improvement over those requiring proof.
        - if ``t == None``, returns the current elliptic curve proof status.

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
        if t is None:
            return self._require_proof["elliptic_curve"]
        self._require_proof["elliptic_curve"] = bool(t)

    def linear_algebra(self, t=None):
        """
        Controls the default proof strategy for linear algebra algorithms.

        INPUT:

        - ``t`` -- boolean or ``None``

        OUTPUT:

        - if ``t == True``, requires linear algebra algorithms to (by default)
          return results that are true unconditionally: the correctness will
          not depend on an algorithm with a nonzero probability of returning an
          incorrect answer or on the truth of any unproven conjectures.
        - if ``t == False``, allows linear algebra algorithms to (by default)
          return results that may depend on unproven conjectures or on
          probabilistic algorithms.  Such algorithms often have a substantial
          speed improvement over those requiring proof.
        - if ``t == None``, returns the current linear algebra proof status.

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
        if t is None:
            return self._require_proof["linear_algebra"]
        self._require_proof["linear_algebra"] = bool(t)

    def number_field(self, t=None):
        """
        Controls the default proof strategy for number field algorithms.

        INPUT:

        - ``t`` -- boolean or ``None``

        OUTPUT:

        - if ``t == True``, requires number field algorithms to (by default)
          return results that are true unconditionally: the correctness will
          not depend on an algorithm with a nonzero probability of returning an
          incorrect answer or on the truth of any unproven conjectures.
        - if ``t == False``, allows number field algorithms to (by default)
          return results that may depend on unproven conjectures or on
          probabilistic algorithms.  Such algorithms often have a substantial
          speed improvement over those requiring proof.
        - if ``t == None``, returns the current number field proof status.

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
        if t is None:
            return self._require_proof["number_field"]
        self._require_proof["number_field"] = bool(t)

    def polynomial(self, t=None):
        """
        Controls the default proof strategy for polynomial algorithms.

        INPUT:

        - ``t`` -- boolean or ``None``

        OUTPUT:

        - if ``t == True``, requires polynomial algorithms to (by default)
          return results that are true unconditionally: the correctness will
          not depend on an algorithm with a nonzero probability of returning an
          incorrect answer or on the truth of any unproven conjectures.
        - if ``t == False``, allows polynomial algorithms to (by default)
          return results that may depend on unproven conjectures or on
          probabilistic algorithms.  Such algorithms often have a substantial
          speed improvement over those requiring proof.
        - if ``t == None``, returns the current polynomial proof status.

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
        if t is None:
            return self._require_proof["polynomial"]
        self._require_proof["polynomial"] = bool(t)


_proof_prefs = _ProofPref(True)  # Creates the global object that stores proof preferences.


def get_flag(t=None, subsystem=None):
    """
    Used for easily determining the correct proof flag to use.

    EXAMPLES::

        sage: from sage.structure.proof.proof import get_flag
        sage: get_flag(False)
        False
        sage: get_flag(True)
        True
        sage: get_flag()
        True
        sage: proof.all(False)
        sage: get_flag()
        False
    """
    if t is None:
        if subsystem in ["arithmetic", "elliptic_curve", "linear_algebra", "number_field","polynomial"]:
            return _proof_prefs._require_proof[subsystem]
        else:
            return _proof_prefs._require_proof["other"]
    return t


class WithProof():
    """
    Use :class:`WithProof` to temporarily set the value of one of the proof
    systems for a block of code, with a guarantee that it will be set
    back to how it was before after the block is done, even if there is an error.

    EXAMPLES:

    This would hang "forever" if attempted with ``proof=True``::

        sage: proof.arithmetic(True)
        sage: with proof.WithProof('arithmetic', False):                                # needs sage.libs.pari
        ....:      print((10^1000 + 453).is_prime())
        ....:      print(1/0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: rational division by zero
        sage: proof.arithmetic()
        True
    """
    def __init__(self, subsystem, t):
        """
        TESTS::

            sage: proof.arithmetic(True)
            sage: P = proof.WithProof('arithmetic',False); P
            <sage.structure.proof.proof.WithProof object at ...>
            sage: P._subsystem
            'arithmetic'
            sage: P._t
            False
            sage: P._t_orig
            True
        """
        self._subsystem = str(subsystem)
        self._t = bool(t)
        self._t_orig = _proof_prefs._require_proof[subsystem]

    def __enter__(self):
        """
        TESTS::

            sage: proof.arithmetic(True)
            sage: P = proof.WithProof('arithmetic',False)
            sage: P.__enter__()
            sage: proof.arithmetic()
            False
            sage: proof.arithmetic(True)
        """
        _proof_prefs._require_proof[self._subsystem] = self._t

    def __exit__(self, *args):
        """
        TESTS::

            sage: proof.arithmetic(True)
            sage: P = proof.WithProof('arithmetic',False)
            sage: P.__enter__()
            sage: proof.arithmetic()
            False
            sage: P.__exit__()
            sage: proof.arithmetic()
            True
        """
        _proof_prefs._require_proof[self._subsystem] = self._t_orig
