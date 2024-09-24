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
                sc = str(c)
                if sc == "0":
                    continue
                if sc == "1":
                    s += " + %s" % names[i]
                elif sc == "-1":
                    s += " - %s" % names[i]
                elif c._is_atomic():
                    if sc[0] == "-":
                        s += " - %s*%s" % (-c, names[i])
                    else:
                        s += " + %s*%s" % (sc, names[i])
                else:
                    s += " + (%s)*%s" % (c, names[i])
            if s == "":
                return "0"
            elif s[1] == '-':
                return '-' + s[3:]
            else:
                return s[3:]

    def vector(self):
        V = self.parent().module()
        return V(self.list())

    def image(self):
        return self.parent()._pseudohom(self)
