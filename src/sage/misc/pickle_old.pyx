# sage_setup: distribution = sagemath-modules

from sage.misc.superseded import deprecation_cython as deprecation
import sage.rings.fraction_field_element


def make_element_old(_class, _dict, parent):
    """
    Used for unpickling old pickles of Element objects (and subclasses).

    This function is deprecated and is kept only to support old pickles.
    """
    deprecation(4698, "Your data is stored in an old format. Please use the save() function to store your data in a more recent format.")
    if _class is sage.rings.fraction_field_element.FractionFieldElement:
        return sage.rings.fraction_field_element.make_element_old(parent, _dict)
    new_object = _class.__new__(_class)
    new_object._set_parent(parent)
    new_object.__dict__ = _dict
    return new_object
