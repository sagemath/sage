try: 
    import sage.structure.element  # to break a cyclic import (misc.constant_function > structure.element > ... > misc.constant_function)
except ModuleNotFoundError:
    # Ignore, the import above may fail while collecting metadata info (at which point sage.structure is not yet installed)
    pass
