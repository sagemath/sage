********
Apêndice
********

.. _section-precedence:

Precedência de operações aritméticas binárias
=============================================

Quanto é ``3^2*4 + 2%5``? A resposta (38) é determinada pela "tabela
de precedência" abaixo. A tabela abaixo é baseada na tabela em § 5.14
do *Python Language Reference Manual* by G. Rossum and F. Drake. As
operações estão listadas aqui em ordem crescente de precedência.


==========================  =================
Operadores                  Descrição
==========================  =================
or                          "ou" booleano
and                         "e" booleano
not                         "não" booleano
in, not in                  pertence
is, is not                  teste de identidade
>, <=, >, >=, ==, !=, <>    comparação
+, -                        adição, subtração
\*, /, %                    multiplicação, divisão, resto
\*\*, ^                     exponenciação
==========================  =================

Portanto, para calcular ``3^2*4 + 2%5``, O Sage inclui parenteses de
precedência da seguinte forma: ``((3^2)*4) + (2%5)``. Logo, primeiro
calcula ``3^2``, que é ``9``, então calcula ``(3^2)*4`` e ``2%5``, e
finalmente soma os dois.
