.. _chapter-interactive_shell:

*********
대화형 셸
*********

``sage`` 명령으로 Sage 인터프리터를 시작하면 많은 함수와 클래스가 도입된 수정된
버전의 ``IPython`` 셸이 시작됩니다.

.. CODE-BLOCK:: text

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version X.X, Release Date: YYYY-MM-DD                     │
    │ Using Python X.X.X. Type "help()" for help.                        │
    └────────────────────────────────────────────────────────────────────┘
    sage:

Sage를 종료하려면 :kbd:`Ctrl+d` (또는 ``quit`` 또는 ``exit``)를 입력하세요.

(시그널을 통한 종료인 ``kill -9`` 와 같은 방법을 피하십시오. 이 방법은 일부
프로세스를 종료하지 못하거나 임시 파일을 정리하지 못할 수 있습니다.)

유용한 바로 가기 키와 명령어
============================

- ``?`` -- 함수에 대한 도움말 표시

- ``??`` -- 함수의 소스 코드 표시

- ``help()`` -- 명령 또는 클래스 매뉴얼 표시 (``q`` 종료를 위해)

- ``!`` -- 어떤 UNIX 명령어를 입력하기 전에 사용하세요

- ``_`` -- 마지막 비어 있지 않은 출력

- ``__`` -- 끝에서 두 번째 비어 있지 않은 출력

- ``ih``, ``_oh`` -- 세션 명령어의 입력과 출력 목록

- :kbd:`Ctrl+r` -- 과거 명령어에서 검색하기

- ``%history`` 또는 ``%hist`` -- 현재 세션의 모든 이전 명령어 목록을 반환합니다

- ``%time`` --명령어 앞에 오고 실행 시간을 표시합니다

- ``timeit()`` -- 반복 실행 후 명령어 실행 시간 정보를 표시합니다. 명령어는
  문자열로 입력됩니다 (``""`` 내부)

- ``cputime()`` -- 처리 시간을 표시합니다; 예시 [1]_

- ``walltime()`` -- ``cputime()`` 와 유사한 명령어이지만 '실제' 시간을
  측정합니다. 벽시계가 측정하는 것처럼

- ``%macro`` -- 매크로 명령어 생성 (여러 명령어의 단축키); 예제 [2]_.

- ``logstart`` -- 세션 입력 명령어 기록 시작; ``load()`` 명령어를 사용하여
  로드하고 다시 실행합니다

- ``save_session()`` -- 세션 저장

- ``load_session()`` -- 세션 불러오기

일부 예시
=========

.. [1] ``cputime()``

.. skip

::

    sage: t = cputime()
    sage: a = int(1938) ^ int(99484)
    sage: b = 1938 ^ 99484
    sage: c = pari(1938) ^ pari(99484)
    sage: cputime(t)
    0.11

.. [2] ``%macro``

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    E = EllipticCurve([1,2,3,4,5])
    M = ModularSymbols(37)
    %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).
    === Macro contents: ===
    E = EllipticCurve([Integer(1),Integer(2),Integer(3),Integer(4),Integer(5)])
    M = ModularSymbols(Integer(37))
