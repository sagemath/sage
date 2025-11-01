"""
The set of prime numbers and its subsets defined by congruence conditions

AUTHORS:

 - William Stein (2005): original version
 - Florent Hivert (2009-11): adapted to the category framework
 - Xavier Caruso (2025-10): implement congruence conditions
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2009 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#                     2025 Xavier Caruso <xavier@caruso.ovh>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from .set import Set_generic
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.arith.misc import euler_phi
from sage.structure.unique_representation import UniqueRepresentation


class Primes(Set_generic, UniqueRepresentation):
    @staticmethod
    def __classcall_private__(cls, modulus=1, classes=None, exceptions=None):
        modulus = ZZ(modulus)
        if modulus == 0:
            raise ValueError("modulus must be nonzero")
        if modulus < 0:
            modulus = -modulus
        if classes is None:
            classes = [1]
        indic = modulus * [False]
        if exceptions is None:
            exceptions = {}
        for c in classes:
            indic[ZZ(c) % modulus] = True
        for c in range(modulus):
            g = modulus.gcd(c)
            if g > 1:
                if indic[c]:
                    if c == 0:
                        if modulus not in exceptions:
                            exceptions[modulus] = True
                    else:
                        if c not in exceptions:
                            exceptions[c] = True
                indic[c] = None
        for p, mult in modulus.factor():
            while mult > 0:
                m = modulus // p
                add_true = []
                add_false = []
                add_excluded = []
                for c in range(m):
                    cs = [indic[c + m*i] for i in range(p) if indic[c + m*i] is not None]
                    if not cs:
                        pass
                    elif all(cs):
                        if m.gcd(c) == 1:
                            add_true.append(c)
                        for i in range(p):
                            j = c + m*i
                            if indic[j] is None:
                                if j == 0:
                                    add_excluded.append(modulus)
                                else:
                                    add_excluded.append(j)
                    elif not any(cs):
                        if m.gcd(c) == 1:
                            add_false.append(c)
                    else:
                        mult = 0
                        break
                else:
                    for c in add_true:
                        indic[c] = True
                    for c in add_false:
                        indic[c] = False
                    for c in add_excluded:
                        if c not in exceptions:
                            exceptions[c] = False
                    modulus = m
                    mult -= 1
        classes = tuple([c for c in range(modulus) if indic[c] is True])
        excep = []
        for c, v in list(exceptions.items()):
            c = ZZ(c)
            if c.is_prime() and (v != (indic[c % modulus] is True)):
                excep.append((c, v))
        excep.sort()
        return cls.__classcall__(cls, modulus, classes, tuple(excep))

    def __init__(self, modulus, classes, exceptions):
        if classes:
            category = InfiniteEnumeratedSets()
        else:
            category = FiniteEnumeratedSets()
        super().__init__(facade=ZZ, category=category)
        self._modulus = modulus
        self._classes = set(classes)
        self._exceptions = {}
        self._included = []
        self._excluded = []
        for c, v in exceptions:
            self._exceptions[c] = v
            if v:
                self._included.append(c)
            else:
                self._excluded.append(c)
        self._included.sort()
        self._excluded.sort()

    def _repr_(self):
        classes = sorted(list(self._classes))
        sc = ", ".join([str(c) for c in classes])
        si = ", ".join([str(i) for i in self._included])
        se = ", ".join([str(e) for e in self._excluded])
        if sc == "":
            if si == "":
                return "Empty set of primes"
            else:
                return "Finite set of primes: %s" % si
        if self._modulus == 1:
            s = "Set of all prime numbers"
        else:
            s = "Set of prime numbers congruent to %s modulo %s" % (sc, self._modulus)
        if si != "":
            s += " with %s included" % si
        if se != "":
            if si == "":
                s += " with %s excluded" % se
            else:
                s += " and %s excluded" % se
        s += ": %s, ..." % (", ".join([str(c) for c in self.first_n(4)]))
        return s

    def __contains__(self, x):
        try:
            if x not in ZZ:
                return False
        except TypeError:
            return False
        x = ZZ(x)
        if not x.is_prime():
            return False
        e = self._exceptions.get(x, None)
        return (e is True) or (e is None and x % self._modulus in self._classes)

    def next(self, pr):
        pr = ZZ(pr)
        if not self._classes:
            if not (self._included and pr < self._included[-1]):
                raise ValueError("no element greater that %s in this set" % pr)
            min = 0
            max = len(self._included)
            while min < max:
                i = (min + max) // 2
                print(min, max, i)
                if self._included[i] <= pr:
                    min = i + 1
                if self._included[i] > pr:
                    max = i
            return self._included[min]
        while True:
            pr = pr.next_prime()
            e = self._exceptions.get(pr, None)
            if (e is True) or (e is None and pr % self._modulus in self._classes):
                return pr

    def first(self):
        return self.next(1)

    def first_n(self, n):
        pr = 1
        ans = []
        for _ in range(n):
            try:
                pr = self.next(pr)
            except ValueError:
                return ans
            ans.append(pr)
        return ans

    def _an_element_(self, n):
        if self.is_finite():
            if self._included:
                return self._included[0]
            raise ValueError("this set is empty")
        return self.next(42)

    def unrank(self, n):
        pr = 1
        for _ in range(n):
            try:
                pr = self.next(pr)
            except ValueError:
                raise ValueError("this set has less than %s elements" % n)
        return pr

    def is_empty(self):
        return not bool(self._classes) and not bool(self._exceptions)

    def is_finite(self):
        return not bool(self._classes)

    def is_cofinite(self):
        return self._modulus == 1 and bool(self._classes)

    def density(self):
        return len(self._classes) / euler_phi(self._modulus)

    def include(self, elements, check=True):
        if elements in ZZ:
            elements = [elements]
        exceptions = self._exceptions.copy()
        for x in elements:
            if check and not x.is_prime():
                raise ValueError("%s is not a prime number" % x)
            exceptions[x] = True
        return Primes(self._modulus, self._classes, exceptions)

    def exclude(self, elements):
        if elements in ZZ:
            elements = [elements]
        exceptions = self._exceptions.copy()
        for x in elements:
            exceptions[x] = False
        return Primes(self._modulus, self._classes, exceptions)

    def complement_in_primes(self):
        modulus = self._modulus
        classes = [c for c in range(modulus)
                   if c % self._modulus not in self._classes]
        exceptions = {c: not v for c, v in self._exceptions.items()}
        return Primes(modulus, classes, exceptions)

    def intersection(self, other):
        if other in ZZ:
            return self
        if not isinstance(other, Primes):
            raise NotImplementedError("boolean operations are only implemented with other sets of primes")
        modulus = self._modulus.lcm(other._modulus)
        classes = [c for c in range(modulus)
                   if (c % self._modulus in self._classes
                   and c % other._modulus in other._classes)]
        exceptions = {}
        for c, v in self._exceptions.items():
            if v and c in other:
                exceptions[c] = True
            if not v:
                exceptions[c] = False
        for c, v in other._exceptions.items():
            if v and c in self:
                exceptions[c] = True
            if not v:
                exceptions[c] = False
        return Primes(modulus, classes, exceptions)

    def union(self, other):
        if other is ZZ:
            return ZZ
        if not isinstance(other, Primes):
            raise NotImplementedError("boolean operations are only implemented with other sets of primes")
        modulus = self._modulus.lcm(other._modulus)
        classes = [c for c in range(modulus)
                   if (c % self._modulus in self._classes
                    or c % other._modulus in other._classes)]
        exceptions = {}
        for c, v in self._exceptions.items():
            if v:
                exceptions[c] = True
            if not v and c not in other:
                exceptions[c] = False
        for c, v in other._exceptions.items():
            if v:
                exceptions[c] = True
            if not v and c not in self:
                exceptions[c] = False
        return Primes(modulus, classes, exceptions)

    def is_subset(self, other):
        return self.intersection(other) == self

    def is_supset(self, other):
        return self.intersection(other) == self

    def is_disjoint(self, other):
        return self.intersection(other).is_empty()
