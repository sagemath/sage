********
Posfacio
********

Por quê o Python?
=================

Vantagens do Python
-------------------

A primeira linguagem de implementação do Sage é o Python (veja [Py]_),
embora rotinas que precisam ser muito rápidas são implementadas em uma
linguagem compilada. O Python possui várias vantagens:

-  **Salvar objetos** é bem suportado em Python. Existe suporte
   extenso em Python para salvar (na grande maioria dos casos) objetos
   arbitrários em arquivos em disco ou em uma base de dados.

-  Suporte excelente para **documentação** de funções e pacotes no
   código fonte, incluindo extração automática de documentação e teste
   automático de todos os exemplos. Os exemplos são automaticamente
   testados regularmente para garantir que funcionam como indicado.

-  **Gerenciamento de memória:** O Python agora possui um sistema de
   gerenciamento de memória e "garbage collector" muito bem pensados e
   robustos que lidam corretamente com referências circulares, e
   permitem variáveis locais em arquivos.

-  O Python possui **diversos pacotes** disponíveis que podem ser de
   grande interesse para os usuários do Sage: análise numérica e
   álgebra linear, visualização 2D e 3D, comunicação via rede (para
   computação distribuída e servidores, por exemplo, via twisted),
   suporte a base de dados, etc.

-  **Portabilidade:** O Python é fácil de compilar a partir do código
   fonte em poucos minutos na maioria das arquiteturas.

-  **Manuseamento de exceções:** O Python possui um sofisticado e bem
   pensado sistema de manuseamento de exceções, através do qual
   programas podem facilmente se recuperar mesmo se ocorrerem erros no
   código que está sendo executado.

-  **Debugador:** O Python inclui um debugador, de modo que quando
   alguma rotina falha por algum motivo, o usuário pode acessar
   extensiva informação sobre a pilha de cálculos e inspecionar o
   estado de todas as variáveis relevantes.

-  **Profiler:** Existe um profiler para o Python, o qual executa
   programas e cria um relatório detalhando quantas vezes e por quando
   tempo cada função foi executada.

-  **Uma Linguagem:** Em vez de escrever uma **nova linguagem** para
   matemática como foi feito para o Magma, Maple, Mathematica, Matlab,
   GP/PARI, GAP, Macaulay 2, Simath, etc., nós usamos a linguagem
   Python, que é uma linguagem de programação popular que está
   sendo desenvolvida e otimizada ativamente por centenas de
   engenheiros de software qualificados. O Python é uma grande
   história de sucesso em desenvolvimento com código aberto com um
   processo de desenvolvimento maduro (veja [PyDev]_).

.. _section-mathannoy:

O Pré-Processador: Diferenças entre o Sage e o Python
-----------------------------------------------------

Alguns aspectos matemáticos do Python podem ser confusos, logo o Sage
se comporta diferentemente do Python em diversas situações.

-  **Notação para exponenciação:** ``**`` versus ``^``. Em Python,
   ``^`` significa "xor", não exponenciação, logo em Python temos

   ::

       >>> 2^8
       10
       >>> 3^2
       1
       >>> 3**2
       9

   Esse uso de ``^`` pode parecer estranho, e é ineficiente para
   pesquisa em matemática pura, pois a função "ou exclusivo" é
   raramente usada. Por conveniência, o Sage pre-processa todos as
   linhas de comandos antes de passá-las para o Python, substituindo
   ocorrências de ``^`` que não estão em strings por ``**``:

   ::

       sage: 2^8
       256
       sage: 3^2
       9
       sage: "3^2"
       '3^2'

-  **Divisão por inteiros:** A expressão em Python ``2/3`` não se
   comporta da forma que um matemático esperaria. Em Python 2, se ``m``
   e ``n`` são inteiros (int), então ``m/n`` também é um inteiro
   (int), a saber, o quociente de ``m`` dividido por ``n``. Portanto
   ``2/3=0``. Tem havido discussões na comunidade do Python para
   modificar o Python de modo que ``2/3`` retorne um número de
   precisão flutuante (float) ``0.6666...``, e ``2//3`` retorne ``0``.

   Nós lidamos com isso no interpretador Sage, encapsulando inteiros
   literais em ``Integer()`` e fazendo a divisão um construtor para
   números racionais. Por exemplo:

   ::

       sage: 2/3
       2/3
       sage: (2/3).parent()
       Rational Field
       sage: 2//3
       0

-  **Inteiros longos:** O Python possui suporte nativo para inteiros
   com precisão arbitrária, além de int's do C. Esses são
   significantemente mais lentos do que os fornecidos pela biblioteca
   GMP, e têm a propriedade que eles são impressos com o sufixo ``L``
   para distingui-los de int's (e isso não será modificado no futuro
   próximo). O Sage implementa inteiros com precisão arbitrária usando
   a biblioteca C do GMP, e esses são impressos sem o sufixo ``L``.

Em vez de modificar o interpretador Python (como algumas pessoas
fizeram para projetos internos), nós usamos a linguagem Python
exatamente com ela é, e escrevemos um pré-processador para o IPython de
modo que o comportamento da linha de comando seja o que um matemático
espera. Isso significa que qualquer programa existente em Python pode
ser usado no Sage. Todavia, deve-se obedecer as regras padrão do
Python para escrever programas que serão importados no Sage.

(Para instalar uma biblioteca do Python, por exemplo uma que você
tenha encontrado na internet, siga as instruções, mas execute ``sage
-python`` em vez de ``python``. Frequentemente isso significa digitar
``sage -python setup.py install``.)

Eu gostaria de contribuir de alguma forma. Como eu posso?
=========================================================

Se você quiser contribuir para o Sage, a sua ajuda será muito bem
vinda! Ela pode variar desde substancial quantidade de código, até
contribuições com respeito à documentação ou notificação de defeitos
(bugs).

Explore a página na web do Sage para informações para desenvolvedores;
entre outras coisas, você pode encontrar uma lista longa de projetos
relacionados ao Sage ordenados por prioridade e categoria. O `Guia
para desenvolvedores do Sage
<http://doc.sagemath.org/html/en/developer/>`_ (em inglês) também possui
informações úteis, e você pode também visitar o grupo de discussões
``sage-devel`` no Google Groups.

Como eu faço referência ao Sage?
================================

Se você escrever um artigo usando o Sage, por favor faça referência
aos cálculos feitos com o Sage incluindo

::

    [Sage] SageMath, the Sage Mathematics Software System (Version 8.7),
           The Sage Developers, 2019, https://www.sagemath.org.

na sua bibliografia (substituindo 8.7 pela versão do Sage que você
está usando). Além disso, procure observar quais componentes do Sage
você está usando em seus cálculos, por exemplo, PARI, Singular, GAP,
Maxima, e também site esses sistemas. Se você está em dúvida sobre
qual software está sendo usado em seus cálculos, fique à vontade para
perguntar no grupo ``sage-devel`` do Google Groups. Veja
:ref:`section-univariate` para mais discussões sobre esse aspecto.

------------

Se por acaso você leu este tutorial do começo ao fim em uma só vez, e
faz idéia de quanto tempo você levou, por favor nos informe no grupo
``sage-devel`` do Google Groups.

Divirta-se com o Sage!
