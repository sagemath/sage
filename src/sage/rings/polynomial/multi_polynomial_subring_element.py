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
