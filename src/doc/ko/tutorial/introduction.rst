****
소개
****

Sage의 첫인상을 보고 싶다면, 여기에 있습니다. Sage는 주로 Python을 기반으로
하지만 사용을 위해 Python의 지식은 필요하지 않습니다. 예를 들면:

::

    sage: 2 + 2
    4
    sage: factor(-2007)
    -1 * 3^2 * 223

    sage: A = matrix(4, 4, range(16)); A
    [ 0  1  2  3]
    [ 4  5  6  7]
    [ 8  9 10 11]
    [12 13 14 15]

    sage: factor(A.charpoly())
    x^2 * (x^2 - 30*x - 80)

    sage: m = matrix(ZZ, 2, range(4))
    sage: m[0,0] = m[0,0] - 3
    sage: m
    [-3  1]
    [ 2  3]

    sage: E = EllipticCurve([1, 2, 3, 4, 5]);  # Ελλειπτική καμπύλη
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5
    over Rational Field
    sage: E.anlist(10)
    [0, 1, 1, 0, -1, -3, 0, -1, -3, -3, -3]
    sage: E.rank()
    1

    sage: k = 1/(sqrt(3)*I + 3/4 + sqrt(73)*5/9); k
    36/(20*sqrt(73) + 36*I*sqrt(3) + 27)
    sage: N(k)
    0.165495678130644 - 0.0521492082074256*I
    sage: N(k, 30)  # 30 "bits"
    0.16549568 - 0.052149208*I
    sage: latex(k)
    \frac{36}{20 \, \sqrt{73} + 36 i \, \sqrt{3} + 27}

.. _installation:

설치
====

Sage `Installation Guide
<https://doc.sagemath.org/html/en/installation/index.html>`_ 를 참조하세요 (`이
링크 <../../en/installation/index.html>`_ 는 설치 가이드의 로컬 복사본으로
이동합니다).

컴퓨터에 Sage가 설치되어 있지 않은 경우, https://sagecell.sagemath.org 에서
시도해 볼 수 있습니다.

Sage는 "배터리 포함"입니다. 즉, Sage는 Python, IPython, PARI, GAP, Singular,
Maxima, NTL, GMP 등을 사용하지만 별도로 설치할 필요가 없습니다. 이들은 Sage
배포본에 포함되어 있습니다. 그러나 Sage의 일부 기능을 사용하려면 Macaulay 또는
KASH 등 관련 프로그램이 이미 컴퓨터에 설치되어 있어야 합니다.


Sage 사용 방법
==============

Sage를 여러 가지 방법으로 사용할 수 있습니다.


- **대화형 셸:** :ref:`chapter-interactive_shell` 참조,

- **노트북 그래픽 환경:** ``sage -n jupyter`` 를 실행하세요.· 확인하세요
   `the Jupyter documentation on-line <https://jupyter-notebook.readthedocs.io/en/latest/notebook.html>`_,

- **프로그램:** 해석되고 컴파일된 프로그램을 작성하거나, Sage 라이브러리를
  사용하는 Python 스크립트를 작성합니다

Sage의 장기적 목표
==================

- **유용한:** Sage의 대상은 수학 학생, 교사 및 연구 수학자입니다. 목표는
  대수학, 기하학, 수론, 미적분학, 수치 계산 등 수학적 구조를 탐색하고 실험할 수
  있는 소프트웨어를 제공하는 것입니다. Sage는 수학 객체를 대화식으로 실험하는
  것을 용이하게 합니다.

- **효율적인:** Sage는 GMP, PARI, GAP 및 NTL과 같은 최적화되고 성숙한
  소프트웨어를 사용하므로 특정 기능에서 매우 빠릅니다.

- **무료 및 오픈 소스:** 소스 코드는 자유롭게 이용 가능하고 읽을 수 있어야
  합니다. 이렇게 함으로써 사용자들은 시스템이 실제로 무엇을 하는지 이해하고
  보다 쉽게 확장할 수 있습니다. 수학자들이 정리를 주의 깊게 읽어야만 정리를
  보다 깊게 이해할 수 있는 것처럼, 계산을 하는 사람들도 소스 코드를 읽음으로써
  컴퓨팅 프로세스를 더 깊이 이해할 수 있습니다. Sage를 논문에 사용하는 경우,
  독자들이 항상 Sage 및 해당 소스 코드에 무료로 액세스 할 수 있도록 하고 사용한
  Sage 버전을 아카이브하고 재배포할 수 있습니다.

- **사용하기 쉬운:** Sage는 Linux, OS X 및 Windows용 소스 코드에서 쉽게 빌드될
  수 있어야 합니다. 이렇게 함으로써 사용자들은 시스템을 수정할 수 있는 유연성을
  가질 수 있습니다.

- **협력적인:** 대부분의 다른 컴퓨터 대수 시스템에 강력한 인터페이스를
  제공하며, PARI, GAP, Singular, Maxima, KASH, Magma, Maple 및 Mathematica를
  포함합니다. Sage는 기존의 수학 소프트웨어를 통합하고 확장하는 것을 목표로
  합니다.

- **종합적인 안내서:** 세미나, 프로그래밍 가이드, 참조 매뉴얼 및 팁으로 수학적
  배경에 대한 많은 예제와 토론을 제공합니다.

- **확장 가능한:** 새로운 데이터 유형을 정의하거나 내장된 유형에서 유추하고
  다양한 언어로 작성된 코드를 사용할 수 있어야 합니다.

- **사용하기 편리한**: 각 구성 요소에 대한 제공된 기능을 이해하기 쉽고 문서 및
  소스 코드의 표시. 사용자 지원 수준이 높습니다
