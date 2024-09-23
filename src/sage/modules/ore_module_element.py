from sage.modules.free_module_element import FreeModuleElement_generic_dense

class OreModule_element(FreeModuleElement_generic_dense):
    def _repr_(self):
        parent = self.parent()
        if parent._names is None:
            return self.parent()._repr_element(self)
        else:
            rank = parent.rank()
            names = parent._names
            s = ""
            for i in range(rank):
                c = self[i]
                if c == 0:
                    continue
                if c == 1:
                    s += " + %s" % names[i]
                else:
                    s += " + (%s)*%s" % (c, names[i])
            if s == "":
                return "0"
            else:
                return s[3:]

    def vector(self):
        M = self.parent()
        V = M.base_ring() ** M.rank()
        return V(self)

    def image(self):
        return self.parent()._pseudohom(self)
