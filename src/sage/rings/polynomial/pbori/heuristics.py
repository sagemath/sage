from .PyPolyBoRi import Polynomial, gauss_on_polys


def dense_system(I):
    I = (Polynomial(p) for p in I)
    I = (p for p in I if not p.is_zero())
    for p in I:
        d = p.deg()
        if p.deg() == 1:
            continue
        else:
            try:
                if len(p) > 2 ** d + 5:
                    return True
            except OverflowError:
                return True
    return False


def gauss_on_linear(I):
    I = (Polynomial(p) for p in I)
    linear = []
    non_linear = []
    for p in I:
        if p.is_zero():
            continue
        if p.deg() <= 1:
            linear.append(p)
        else:
            non_linear.append(p)
    if not linear:
        return non_linear
    linear = list(gauss_on_polys(linear))
    return linear + non_linear
