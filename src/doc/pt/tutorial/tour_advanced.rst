Um Pouco Mais de Matemática Avançada
====================================

Geometria Algébrica
-------------------

Você pode definir variedades algébricas arbitrárias no Sage, mas as
vezes alguma funcionalidade não-trivial é limitada a anéis sobre
:math:`\QQ` ou corpos finitos. Por exemplo, vamos calcular a união de
duas curvas planas afim, e então recuperar as curvas como as
componentes irredutíveis da união.

::

    sage: x, y = AffineSpace(2, QQ, 'xy').gens()
    sage: C2 = Curve(x^2 + y^2 - 1)
    sage: C3 = Curve(x^3 + y^3 - 1)
    sage: D = C2 + C3
    sage: D
    Affine Plane Curve over Rational Field defined by 
       x^5 + x^3*y^2 + x^2*y^3 + y^5 - x^3 - y^3 - x^2 - y^2 + 1
    sage: D.irreducible_components()
    [
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x^2 + y^2 - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x^3 + y^3 - 1
    ]

Você também pode encontrar todos os pontos de interseção das duas
curvas, intersectando-as, e então calculando as componentes
irredutíveis.

.. link

::

    sage: V = C2.intersection(C3)
    sage: V.irreducible_components()
    [
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y - 1,
      x,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      y,
      x - 1,
    Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
      x + y + 2,
      2*y^2 + 4*y + 3
    ]

Portanto, por exemplo, :math:`(1,0)` e :math:`(0,1)` estão em ambas as
curvas (o que é claramente visível), como também estão certos pontos
(quadráticos) cuja coordenada :math:`y` satisfaz :math:`2y^2 + 4y +
3=0`.

O Sage pode calcular o ideal toroidal da cúbica torcida no espaço-3
projetivo:

::

    sage: R.<a,b,c,d> = PolynomialRing(QQ, 4)
    sage: I = ideal(b^2-a*c, c^2-b*d, a*d-b*c)
    sage: F = I.groebner_fan(); F
    Groebner fan of the ideal:
    Ideal (b^2 - a*c, c^2 - b*d, -b*c + a*d) of Multivariate Polynomial Ring
    in a, b, c, d over Rational Field
    sage: F.reduced_groebner_bases ()
    [[-c^2 + b*d, -b*c + a*d, -b^2 + a*c],
     [-b*c + a*d, -c^2 + b*d, b^2 - a*c],
     [-c^3 + a*d^2, -c^2 + b*d, b*c - a*d, b^2 - a*c],
     [-c^2 + b*d, b^2 - a*c, b*c - a*d, c^3 - a*d^2],
     [-b*c + a*d, -b^2 + a*c, c^2 - b*d],
     [-b^3 + a^2*d, -b^2 + a*c, c^2 - b*d, b*c - a*d],
     [-b^2 + a*c, c^2 - b*d, b*c - a*d, b^3 - a^2*d],
     [c^2 - b*d, b*c - a*d, b^2 - a*c]]
    sage: F.polyhedralfan()
    Polyhedral fan in 4 dimensions of dimension 4

Curvas Elípticas
----------------

A funcionalidade para curvas elípticas inclui a maior parte da
funcionalidade para curvas elípticas do PARI, acesso aos dados da base
de dados Cremona (isso requer um pacote adicional), os recursos do
mwrank, isto é, "2-descends" com cálculos do grupo de Mordell-Weil
completo, o algoritmo SEA (sigla em inglês), cálculo de todas as
isogenias, bastante código novo para curvas sobre :math:`\QQ`, e parte
do software "algebraic descent" de Denis Simons.

O comando ``EllipticCurve`` para criar curvas elípticas possui várias
formas:


-  EllipticCurve([:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]):
   Fornece a curva elíptica

   .. math::  y^2+a_1xy+a_3y=x^3+a_2x^2+a_4x+a_6,


   onde os :math:`a_i`'s são coagidos para a família de :math:`a_1`.
   Se todos os :math:`a_i` possuem parente :math:`\ZZ`, então eles são
   coagidos para :math:`\QQ`.

-  EllipticCurve([:math:`a_4`, :math:`a_6`]): Conforme acima, mas
   :math:`a_1=a_2=a_3=0`.

-  EllipticCurve(label): Fornece a curva elíptica da base de dados
   Cremona com o "label" (novo) dado. O label é uma string, tal como
   ``"11a"`` ou ``"37b2"``. As letras devem ser minúsculas (para
   distinguir dos labels antigos).

-  EllipticCurve(j): Fornece uma curva elíptica com invariante
   :math:`j`.

-  EllipticCurve(R,
   [:math:`a_1`, :math:`a_2`, :math:`a_3`, :math:`a_4`, :math:`a_6`]):
   Cria uma curva elíptica sobre um anel :math:`R` com os
   :math:`a_i`'s.


Agora ilustramos cada uma dessas construções:

::

    sage: EllipticCurve([0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve([GF(5)(0),0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

    sage: EllipticCurve([1,2])
    Elliptic Curve defined by y^2  = x^3 + x + 2 over Rational Field

    sage: EllipticCurve('37a')
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field

    sage: EllipticCurve_from_j(1)
    Elliptic Curve defined by y^2 + x*y = x^3 + 36*x + 3455 over Rational Field

    sage: EllipticCurve(GF(5), [0,0,1,-1,0])
    Elliptic Curve defined by y^2 + y = x^3 + 4*x over Finite Field of size 5

O par :math:`(0,0)` é um ponto na curva elíptica :math:`E` definida
por :math:`y^2 + y = x^3 - x`. Para criar esse ponto digite
``E([0,0])``. O Sage pode somar pontos em uma curva elíptica
(lembre-se que é possível definir uma estrutura de grupo aditivo em
curvas elípticas onde o ponto no infinito é o elemento nulo, e a some
de três pontos colineares sobre a curva é zero):

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
    sage: P = E([0,0])
    sage: P + P
    (1 : 0 : 1)
    sage: 10*P
    (161/16 : -2065/64 : 1)
    sage: 20*P
    (683916417/264517696 : -18784454671297/4302115807744 : 1)
    sage: E.conductor()
    37

As curvas elípticas sobre os números complexos são parametrizadas
pelo invariante :math:`j`. O Sage calcula o invariante :math:`j` da
seguinte forma:

::

    sage: E = EllipticCurve([0,0,0,-4,2]); E
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: E.conductor()
    2368
    sage: E.j_invariant()
    110592/37

Se criarmos uma curva com o mesmo invariante :math:`j` que a curva
:math:`E`, ela não precisa ser isomórfica a :math:`E`. No seguinte
exemplo, as curvas não são isomórficas porque os seus condutores são
diferentes.

::

    sage: F = EllipticCurve_from_j(110592/37)
    sage: F.conductor()
    37

Todavia, uma torção de :math:`F` por um fator 2 resulta em uma curva
isomórfica.

.. link

::

    sage: G = F.quadratic_twist(2); G
    Elliptic Curve defined by y^2 = x^3 - 4*x + 2 over Rational Field
    sage: G.conductor()
    2368
    sage: G.j_invariant()
    110592/37

Nós podemos calcular os coeficientes :math:`a_n` de uma
série-:math:`L` ou forma modular :math:`\sum_{n=0}^\infty
a_nq^n` associada à curva elíptica. Esse cálculo usa a biblioteca C do
PARI.

::

    sage: E = EllipticCurve([0,0,1,-1,0])
    sage: E.anlist(30)
    [0, 1, -2, -3, 2, -2, 6, -1, 0, 6, 4, -5, -6, -2, 2, 6, -4, 0, -12, 0, -4,
     3, 10, 2, 0, -1, 4, -9, -2, 6, -12]
    sage: v = E.anlist(10000)

Leva apenas um segundo para calcular todos os :math:`a_n` para
:math:`n\leq 10^5`:

.. skip

::

    sage: %time v = E.anlist(100000)
    CPU times: user 0.98 s, sys: 0.06 s, total: 1.04 s
    Wall time: 1.06

Curvas elípticas podem ser construídas usando o "label" da base de
dados Cremona. Isso importa a curva elíptica com informações prévias
sobre o seu posto, números de Tomagawa, regulador, etc.

::

    sage: E = EllipticCurve("37b2")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 1873*x - 31833 over Rational
    Field
    sage: E = EllipticCurve("389a")
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x  over Rational Field
    sage: E.rank()
    2
    sage: E = EllipticCurve("5077a")
    sage: E.rank()
    3

Nós também podemos acessar a base de dados Cremona diretamente.

::

    sage: db = sage.databases.cremona.CremonaDatabase()
    sage: db.curves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1], 'b1': [[0, 1, 1, -23, -50], 0, 3]}
    sage: db.allcurves(37)
    {'a1': [[0, 0, 1, -1, 0], 1, 1],
     'b1': [[0, 1, 1, -23, -50], 0, 3],
     'b2': [[0, 1, 1, -1873, -31833], 0, 1],
     'b3': [[0, 1, 1, -3, 1], 0, 3]}

Os objetos obtidos pela base de dados não são do tipo
``EllipticCurve``. Eles são elementos de uma base de dados e possuem
alguns campos, e apenas isso. Existe uma versão básica da base de
dados Cremona, que já é distribuída na versão padrão do Sage, e contém
informações limitadas sobre curvas elípticas de condutor :math:`\leq
10000`. Existe também uma versão estendida opcional, que contém
informações extensas sobre curvas elípticas de condutor :math:`\leq
120000` (em outubro de 2005). Por fim, existe ainda uma versão (2GB)
opcional de uma base de dados para o Sage que contém centenas de
milhares de curvas elípticas na base de dados Stein-Watkins.

Caracteres de Dirichlet
-----------------------

Um *caractere de Dirichlet* é a extensão de um homomorfismo
:math:`(\ZZ/N\ZZ)* \to R^*`, para algum anel :math:`R`, para o mapa
:math:`\ZZ \to R` obtido mapeando os inteiros :math:`x` tais que
:math:`\gcd(N,x)>1` em 0.

::

    sage: G = DirichletGroup(12)
    sage: G.list()
    [Dirichlet character modulo 12 of conductor 1 mapping 7 |--> 1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1,
    Dirichlet character modulo 12 of conductor 12 mapping 7 |--> -1, 5 |--> -1]
    sage: G.gens()
    (Dirichlet character modulo 12 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
    Dirichlet character modulo 12 of conductor 3 mapping 7 |--> 1, 5 |--> -1)
    sage: len(G)
    4

Tendo criado o grupo, a seguir calculamos um elemento e fazemos
cálculos com ele.

.. link

::

    sage: G = DirichletGroup(21)
    sage: chi = G.1; chi
    Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6
    sage: chi.values()
    [0, 1, zeta6 - 1, 0, -zeta6, -zeta6 + 1, 0, 0, 1, 0, zeta6, -zeta6, 0, -1,
     0, 0, zeta6 - 1, zeta6, 0, -zeta6 + 1, -1]
    sage: chi.conductor()
    7
    sage: chi.modulus()
    21
    sage: chi.order()
    6
    sage: chi(19)
    -zeta6 + 1
    sage: chi(40)
    -zeta6 + 1

É também possível calcular a ação do grupo de Galois
:math:`\text{Gal}(\QQ(\zeta_N)/\QQ)` sobre esses caracteres, bem como
a decomposição em produto direto correspondente à fatorização do
módulo.

.. link

::

    sage: chi.galois_orbit()
    [Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> -zeta6 + 1,
     Dirichlet character modulo 21 of conductor 7 mapping 8 |--> 1, 10 |--> zeta6]

    sage: go = G.galois_orbits()
    sage: [len(orbit) for orbit in go]
    [1, 2, 2, 1, 1, 2, 2, 1]

    sage: G.decomposition()
    [
    Group of Dirichlet characters modulo 3 with values in Cyclotomic Field of order 6 and degree 2,
    Group of Dirichlet characters modulo 7 with values in Cyclotomic Field of order 6 and degree 2
    ]

A seguir, construímos o grupo de caracteres de Dirichlet mod 20, mas
com valores em :math:`\QQ(i)`:

::

    sage: K.<i> = NumberField(x^2+1)
    sage: G = DirichletGroup(20,K)
    sage: G
    Group of Dirichlet characters modulo 20 with values in Number Field in i with defining polynomial x^2 + 1

Agora calculamos diversos invariantes de ``G``:

.. link

::

    sage: G.gens()
    (Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1,
    Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> i)

    sage: G.unit_gens()
    (11, 17)
    sage: G.zeta()
    i
    sage: G.zeta_order()
    4

No próximo exemplo criamos um caractere de Dirichlet com valores em um
corpo numérico. Nós especificamos explicitamente a escolha da raiz da
unidade no terceiro argumento do comando ``DirichletGroup`` abaixo.

::

    sage: x = polygen(QQ, 'x')
    sage: K = NumberField(x^4 + 1, 'a'); a = K.0
    sage: b = K.gen(); a == b
    True
    sage: K
    Number Field in a with defining polynomial x^4 + 1
    sage: G = DirichletGroup(5, K, a); G
    Group of Dirichlet characters modulo 5 with values in the group of order 8 generated by a in Number Field in a with defining polynomial x^4 + 1
    sage: chi = G.0; chi
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> a^2
    sage: [(chi^i)(2) for i in range(4)]
    [1, a^2, -1, -a^2]

Aqui ``NumberField(x^4 + 1, 'a')`` diz para o Sage usar o símbolo "a"
quando imprimir o que é ``K`` (um corpo numérico definido pelo
polinômio :math:`x^4 + 1`). O nome "a" não está declarado até então.
Uma vez que ``a = K.0`` (ou equivalentemente ``a = K.gen()``) é
calculado, o símbolo "a" representa a raiz do polinômio gerador
:math:`x^4+1`.

Formas Modulares
----------------

O Sage pode fazer alguns cálculos relacionados a formas modulares,
incluindo dimensões, calcular espaços de símbolos modulares,
operadores de Hecke, e decomposições.

Existem várias funções disponíveis para calcular dimensões de espaços
de formas modulares. Por exemplo,

::

    sage: from sage.modular.dims import dimension_cusp_forms
    sage: dimension_cusp_forms(Gamma0(11),2)
    1
    sage: dimension_cusp_forms(Gamma0(1),12)
    1
    sage: dimension_cusp_forms(Gamma1(389),2)
    6112

A seguir ilustramos o cálculo dos operadores de Hecke em um espaço de
símbolos modulares de nível :math:`1` e peso :math:`12`.

::

    sage: M = ModularSymbols(1,12)
    sage: M.basis()
    ([X^8*Y^2,(0,0)], [X^9*Y,(0,0)], [X^10,(0,0)])
    sage: t2 = M.T(2)
    sage: t2
    Hecke operator T_2 on Modular Symbols space of dimension 3 for Gamma_0(1)
    of weight 12 with sign 0 over Rational Field
    sage: t2.matrix()
    [ -24    0    0]
    [   0  -24    0]
    [4860    0 2049]
    sage: f = t2.charpoly('x'); f
    x^3 - 2001*x^2 - 97776*x - 1180224
    sage: factor(f)
    (x - 2049) * (x + 24)^2
    sage: M.T(11).charpoly('x').factor()
    (x - 285311670612) * (x - 534612)^2

Podemos também criar espaços para :math:`\Gamma_0(N)` e
:math:`\Gamma_1(N)`.


::

    sage: ModularSymbols(11,2)
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
     0 over Rational Field
    sage: ModularSymbols(Gamma1(11),2)
    Modular Symbols space of dimension 11 for Gamma_1(11) of weight 2 with
    sign 0 over Rational Field

Vamos calcular alguns polinômios característicos e expansões
:math:`q`.

::

    sage: M = ModularSymbols(Gamma1(11),2)
    sage: M.T(2).charpoly('x')
    x^11 - 8*x^10 + 20*x^9 + 10*x^8 - 145*x^7 + 229*x^6 + 58*x^5 - 360*x^4
         + 70*x^3 - 515*x^2 + 1804*x - 1452
    sage: M.T(2).charpoly('x').factor()
    (x - 3) * (x + 2)^2 * (x^4 - 7*x^3 + 19*x^2 - 23*x + 11)
            * (x^4 - 2*x^3 + 4*x^2 + 2*x + 11)
    sage: S = M.cuspidal_submodule()
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.q_expansion_basis(10)
    [
        q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10)
    ]

Podemos até mesmo calcular espaços de símbolos modulares com carácter.

::

    sage: G = DirichletGroup(13)
    sage: e = G.0^2
    sage: M = ModularSymbols(e,2); M
    Modular Symbols space of dimension 4 and level 13, weight 2, character
    [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
    sage: M.T(2).charpoly('x').factor()
    (x - zeta6 - 2) * (x - 2*zeta6 - 1) * (x + zeta6 + 1)^2
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of
    dimension 4 and level 13, weight 2, character [zeta6], sign 0, over
    Cyclotomic Field of order 6 and degree 2
    sage: S.T(2).charpoly('x').factor()
    (x + zeta6 + 1)^2
    sage: S.q_expansion_basis(10)
    [
    q + (-zeta6 - 1)*q^2 + (2*zeta6 - 2)*q^3 + zeta6*q^4 + (-2*zeta6 + 1)*q^5
      + (-2*zeta6 + 4)*q^6 + (2*zeta6 - 1)*q^8 - zeta6*q^9 + O(q^10)
    ]

Aqui está um outro exemplo de como o Sage pode calcular a ação de
operadores de Hecke em um espaço de formas modulares.

::

    sage: T = ModularForms(Gamma0(11),2)
    sage: T
    Modular Forms space of dimension 2 for Congruence Subgroup Gamma0(11) of
    weight 2 over Rational Field
    sage: T.degree()
    2
    sage: T.level()
    11
    sage: T.group()
    Congruence Subgroup Gamma0(11)
    sage: T.dimension()
    2
    sage: T.cuspidal_subspace()
    Cuspidal subspace of dimension 1 of Modular Forms space of dimension 2 for
    Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: T.eisenstein_subspace()
    Eisenstein subspace of dimension 1 of Modular Forms space of dimension 2
    for Congruence Subgroup Gamma0(11) of weight 2 over Rational Field
    sage: M = ModularSymbols(11); M
    Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign
    0 over Rational Field
    sage: M.weight()
    2
    sage: M.basis()
    ((1,0), (1,8), (1,9))
    sage: M.sign()
    0

Denote por :math:`T_p` os operadores de Hecke usuais (:math:`p`
primo).  Como os operadores de Hecke :math:`T_2`, :math:`T_3`,
e :math:`T_5` agem sobre o espaço de símbolos modulares?


.. link

::

    sage: M.T(2).matrix()
    [ 3  0 -1]
    [ 0 -2  0]
    [ 0  0 -2]
    sage: M.T(3).matrix()
    [ 4  0 -1]
    [ 0 -1  0]
    [ 0  0 -1]
    sage: M.T(5).matrix()
    [ 6  0 -1]
    [ 0  1  0]
    [ 0  0  1]
