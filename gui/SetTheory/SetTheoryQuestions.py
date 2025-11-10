import random


def randomQuestion():
    Questions = {1:"Given the sets A = <setA> and B = <setB>, find A ∪ B.",
             2:"Given the sets A = <setA> and B = <setB>, find A ∩ B.",
             3:"Determine if given A = <setA> and B = <setB>, if A ⊆ A ∪ B.",
             4:"Given A = <setA>, B = <setB>, C = <setC>, find\na.) A Δ B\nb.) A ∪ (B ∩ C)\nc.) A ∩ (B ∪ C)",
             5:"Given the set A = <setA>, determine if the family of sets F = {<setB>, <setC>, <setD>} form a partition of A",
             6:"Given A = <setA> and B = <setB>, find\na.) A\\B\nb.)B\\A",
             7:"Given A = <setA> and B ⊆ A and C ⊆ A where B = <setB> and C = <setC>, find\na.) B'\nb.) C'"}
    q_num = random.choice(list(Questions.keys()))
    return q_num, Questions[q_num]