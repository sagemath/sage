import random
import re
from sage.all import Set

class SetQuestionFiller:
    def __init__(self, universe=None, min_size=3, max_size=6, force_valid_partition=False):
        if universe is None:
            universe = list(range(1, 11))
        self.universe = universe
        self.min_size = min_size
        self.max_size = max_size
        self.force_valid_partition = force_valid_partition

    def randomSet(self):
        size = random.randint(self.min_size, self.max_size)
        return Set(random.sample(self.universe, size))

    def randomSubset(self, A, allow_empty=False):
        """Return a random subset of A. Default: non-empty subset."""
        elements = list(A)
        if not elements:
            return Set([])
        min_size = 0 if allow_empty else 1
        size = random.randint(min_size, len(elements))
        return Set(random.sample(elements, size))

    def randomPartition(self, A, k=3):
        """
        Generate k nonempty disjoint subsets whose union is A.
        If k > |A|, k is reduced to |A|.
        """
        elements = list(A)
        random.shuffle(elements)
        k = min(k, max(1, len(A)))
        groups = [[] for _ in range(k)]
        for elem in elements:
            groups[random.randint(0, k-1)].append(elem)
        # ensure nonempty groups
        for i in range(k):
            if not groups[i]:
                donors = [g for g in groups if len(g) > 1]
                if donors:
                    donor = random.choice(donors)
                    groups[i].append(donor.pop())
                else:
                    # all groups size 0 or 1, move from a random non-empty
                    nonempty = [g for g in groups if g]
                    if nonempty:
                        donor = random.choice(nonempty)
                        groups[i].append(donor.pop())
        return [Set(g) for g in groups]

    def fillTemplate(self, template):
        """
        Replace <setX> placeholders with sets.

        Rules (priority):
          1. If 'partition' in text -> generate A and B,C,... as a partition of A
             (if force_valid_partition True) or as random subsets of A otherwise.
          2. Else if '⊆' appears in the text or 'subset' in the text -> generate A
             and then B,C,... as subsets of A (no partition constraints).
          3. Else -> generate independent random sets for each placeholder.

        Returns (filled_template, substitutions_dict)
        """
        # placeholders like <setA>, <setB>, etc.
        labels = set(re.findall(r"<set([A-Z])>", template))
        substitutions = {}

        lower = template.lower()

        # 1) partition detection
        if "partition" in lower:
            # generate A first
            A = self.randomSet()
            substitutions["A"] = str(A)

            # create parts for the other labels
            other_labels = sorted(labels - {"A"})
            if self.force_valid_partition:
                parts = self.randomPartition(A, k=len(other_labels) if other_labels else 1)
                # zip parts to labels (if fewer parts than labels, remaining labels get empty sets)
                for label, subset in zip(other_labels, parts):
                    substitutions[label] = str(subset)
                for label in other_labels[len(parts):]:
                    substitutions[label] = str(Set([]))
            else:
                # random subsets of A (may overlap / not cover A)
                for label in other_labels:
                    substitutions[label] = str(self.randomSubset(A, allow_empty=False))

            # replace placeholders
            for label, s in substitutions.items():
                template = template.replace(f"<set{label}>", str(s))
            return template, substitutions

        # 2) subset-of-A detection: look for the subset symbol or the word 'subset'
        if "⊆" in template or "subset" in lower:
            # ensure we actually have A as a placeholder; if not, still generate A if needed
            A = self.randomSet()
            substitutions["A"] = A
            for label in sorted(labels):
                if label == "A":
                    continue
                substitutions[label] = str(self.randomSubset(A, allow_empty=False))
            # replace placeholders
            for label, s in substitutions.items():
                template = template.replace(f"<set{label}>", str(s))
            return template, substitutions

        # 3) normal case: independent random sets for each placeholder
        for label in labels:
            s = self.randomSet()
            substitutions[label] = str(s)
            template = template.replace(f"<set{label}>", str(s))

        return template, substitutions
