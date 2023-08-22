from .multi_polynomial_element import MPolynomial_element
class MPolynomial_subring_element(MPolynomial_element):

    def __init__(self, parent, x):
        self._parent=parent
        super().__init__(parent, x)

    #########################
    ##!!!!!! DISCUSS !!!!!!##
    #########################
    def degree(self):
        return self.element().degree()
    
    def construction(self, *args, **kwargs):
        return self._parent.construct(self, *args, **kwargs)

    def __eq__(self, other):
        return self.element()==other.element()
    
    def as_quotientring_element(self, *args, **kwargs):
        return self._parent._quotientring_element_representation_(self, *args, **kwargs)