from sage.rings.homset import RingHomset_generic
from sage.rings.ring_extension import RingExtension_class
from sage.rings.ring_extension_morphism import RingExtensionHomomorphism

class RingExtensionHomset(RingHomset_generic):
    def __call__(self, *args, **kwargs):
        return RingExtensionHomomorphism(self, *args, **kwargs)

    def _coerce_impl(self, x):
        if isinstance(x, RingExtensionHomomorphism):
            x = x._backend()
        domain = self.domain()
        if isinstance(domain, RingExtension_class):
            domain = domain._backend()
        codomain = self.codomain()
        if isinstance(codomain, RingExtension_class):
            codomain = codomain._backend()
        if domain is x.domain() and codomain is x.codomain():
            return RingExtensionHomomorphism(self, x)
        raise TypeError
