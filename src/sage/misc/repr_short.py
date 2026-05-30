def repr_short(obj):
    """
    Return a short representation string for Sage objects.
    Falls back to repr(obj) if no _repr_short_ is defined.
    """
    if hasattr(obj, "_repr_short_"):
        return obj._repr_short_()
    return repr(obj)