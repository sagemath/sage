from sage.structure.parent import Parent
from sage.categories.drinfeld_modules import DrinfeldModules
from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism
from sage.categories.homset import Homset, Hom

class DrinfeldModuleHomset(Homset):
    
    Element = DrinfeldModuleMorphism
    __contains__ = Parent.__contains__

    def __init__(self, X, Y, category=None, check=True):
        if category is None:
            category = X.category()
        if check:
            if X.category() != Y.category() \
                    or not isinstance(X.category(), DrinfeldModules):
                raise NotImplementedError('Drinfeld modules must be in the same category')
            if category != X.category():
                raise NotImplementedError('category should be DrinfeldModules')
        base = category.base()
        super().__init__(X, Y, category=category, base=base, check=check)

    def _latex_(self):
        return f'\\text{{Set{{ }}of{{ }}Drinfeld{{ }}module{{ }}morphisms' \
                f'{{ }}from}}{latex(self.domain())}\\text{{{{ }}to{{ }}}}' \
                f'{self.codomain()}'

    def _repr_(self):
        return f'Set of Drinfeld module morphisms:\n' \
                f'  From: {self.domain()}\n' \
                f'  To:   {self.codomain()}'

    def __contains__(self, x):
        try:
            x = self(x)
            return True
        except (AttributeError, ValueError, TypeError):
            return False

    def _element_constructor_(self, *args, **kwds):
        # NOTE: This used to be self.element_class(self, ...), but this
        # would call __init__ instead of __classcall_private__. This
        # seems to work, but I don't know what I'm doing.
        return DrinfeldModuleMorphism(self, *args, **kwds)
