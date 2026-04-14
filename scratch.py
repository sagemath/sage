import types
m = types.ModuleType("m")
class Lazy:
    def __get__(self, inst, owner):
        print("evaluated!")
        return 42

m.foo = Lazy()
print("getattr gives:", type(getattr(m, "foo")))
