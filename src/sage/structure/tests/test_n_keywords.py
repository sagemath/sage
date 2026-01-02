def test_n_accepts_extra_keywords():
    """
    Test that Element.n() accepts extra keyword arguments
    for backward compatibility.
    """
    from sage.all import SR
    x = SR(1)/3
    assert x.n(context=None) == x.n()
