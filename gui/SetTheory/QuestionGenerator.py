from SetTheory.SetQuestionFiller import SetQuestionFiller
from SetTheory.SetTheoryQuestions import randomQuestion
import re
from sage.all import Set


def generateQuestion():
    q_num, template = randomQuestion()
    filler = SetQuestionFiller()

    # Identify placeholders
    placeholders = re.findall(r'<set[A-D]>', template)
    unique_placeholders = sorted(set(placeholders), key=lambda x: placeholders.index(x))

    sets_dict = {}

    if 'partition' in template:
        # Generate A and a partition of A
        A = filler.randomSet()
        partition_sets = filler.randomPartition(A, k=len(unique_placeholders)-1)
        sets_dict[unique_placeholders[0]] = str(A)
        for i, part_set in enumerate(partition_sets):
            sets_dict[unique_placeholders[i+1]] = str(part_set)
    else:
        # Generate sets based on subset relations or randomly
        A = filler.randomSet()
        sets_dict[unique_placeholders[0]] = str(A)
        for ph in unique_placeholders[1:]:
            if '⊆' in template or 'subset' in template:
                subset = filler.randomSubset(A)
                sets_dict[ph] = str(subset)
            else:
                rand_set = filler.randomSet()
                sets_dict[ph] = str(rand_set)

    # Replace placeholders in the template
    sets_list = []
    for ph, s in sets_dict.items():
        template = template.replace(ph, str(s))
        sets_list.append(s)

    return q_num, template, sets_list