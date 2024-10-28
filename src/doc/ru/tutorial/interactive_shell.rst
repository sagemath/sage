.. _chapter-interactive_shell:

**********************
Интерактивная оболочка
**********************

Почти всегда в этом руководстве мы предполагаем, что интерпретатор
Sage был запущен командой ``sage``. Она запустит специальную версию
консоли IPython и импортирует множество функций и классов, так что
они готовы для использования в командной строке. Более тонкая настройка
производится редактированием файла ``$SAGE_ROOT/ipythonrc``. При запуске
Sage вы увидите вывод, похожий на следующий:

.. skip

::

    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage:

Чтобы выйти из Sage, нажмите Ctrl-D или введите ``quit`` или ``exit``.

.. skip

::

    sage: quit
    Exiting Sage (CPU time 0m0.00s, Wall time 0m0.89s)

Wall time — это прошедшее время. Это значение верно, потому как в "CPU time"
не входит время, использованное субпроцессами вроде GAP или Singular.

(Постарайтесь не убивать процесс Sage командой ``kill -9`` из терминала,
потому что Sage может не убить дочерние процессы, такие как Maple, или может
не очистить временные файлы из директории ``$HOME/.sage/tmp``.)

Ваша сессия Sage
================

Сессия — это последовательность вводов и выводов начиная с запуска программы
и заканчивая выходом из нее. Sage заносит всю историю вводов в log-файл,
используя IPython. Если вы используете интерактивную оболочку (не веб-интерфейс
Notebook), то вы можете ввести ``%hist``, чтобы вывести список всех введенных
команд. Вы можете ввести ``?`` в командной строке Sage, чтобы получить больше
информации о IPython, например,
"IPython предоставляет пронумерованные командные строки... с кешированием ввода и вывода. Все введенные данные сохраняются и могут быть использованы как переменные (помимо обычного вызова с помощью стрелок). Следующие глобальные переменные присутствуют всегда (не перезаписывайте их!)":

::

      _:  previous input (interactive shell and notebook)
      __: next previous input (interactive shell only)
      _oh : list of all inputs (interactive shell only)

Пример:

.. skip

::

    sage: factor(100)
     _1 = 2^2 * 5^2
    sage: kronecker_symbol(3,5)
     _2 = -1
    sage: %hist   # Работает только в интерактивной оболочке, но не в Sage notebook.
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    sage: _oh
     _4 = {1: 2^2 * 5^2, 2: -1}
    sage: _i1
     _5 = 'factor(ZZ(100))\n'
    sage: eval(_i1)
     _6 = 2^2 * 5^2
    sage: %hist
    1: factor(100)
    2: kronecker_symbol(3,5)
    3: %hist
    4: _oh
    5: _i1
    6: eval(_i1)
    7: %hist

Мы не включаем номера строк в этом учебном пособии и в другой документации Sage.

Вы также можете хранить список введенных команд сессии в виде макроса для сессии.

.. skip

::

    sage: E = EllipticCurve([1,2,3,4,5])
    sage: M = ModularSymbols(37)
    sage: %hist
    1: E = EllipticCurve([1,2,3,4,5])
    2: M = ModularSymbols(37)
    3: %hist
    sage: %macro em 1-2
    Macro `em` created. To execute, type its name (without quotes).


.. skip

::

    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field
    sage: E = 5
    sage: M = None
    sage: em
    Executing Macro...
    sage: E
    Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over
    Rational Field

При использовании интерактивной оболочки Sage, любая UNIX-команда может быть
запущена с помощью префикса ``!``. Например

.. skip

::

    sage: !ls
    auto  example.sage glossary.tex  t  tmp  tut.log  tut.tex

возвращает содержание текущей директории.

Переменная ``PATH`` сожержит директорию bin (бинарные файлы) в самом начале,
так что если вы запускаете ``gp``, ``gap``, ``singular``, ``maxima``, и т.д.
вы получаете версии, включенные в Sage.

.. skip

::

    sage: !gp
    Reading GPRC: /etc/gprc ...Done.

                               GP/PARI CALCULATOR Version 2.2.11 (alpha)
                      i686 running linux (ix86/GMP-4.1.4 kernel) 32-bit version
    ...
    sage: !singular
                         SINGULAR                             /  Development
     A Computer Algebra System for Polynomial Computations   /   version 3-0-1
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   October 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \

Журналирование ввода и вывода
=============================

Журналирование сессии Sage это не то же самое, что сохрнанение
сессии (см. :ref:`section-save` для этого). Для журналирования ввода (и, опционально,
вывода), используйте команду ``logstart``. Введите ``logstart?`` для подробностей.
Вы можете использовать эту команду для журналирования всего, что вы вводите, всего
вывода, и даже можете воспроизвести введенные данные в будущей сессии (просто
загрузив log-файл).

.. skip

::

    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: logstart setup
    Activating auto-logging. Current session state plus future input saved.
    Filename       : setup
    Mode           : backup
    Output logging : False
    Timestamping   : False
    State          : active
    sage: E = EllipticCurve([1,2,3,4,5]).minimal_model()
    sage: F = QQ^3
    sage: x,y = QQ['x,y'].gens()
    sage: G = E.gens()
    sage:
    Exiting Sage (CPU time 0m0.61s, Wall time 0m50.39s).
    was@form:~$ sage
    ┌────────────────────────────────────────────────────────────────────┐
    │ SageMath version 9.7, Release Date: 2022-01-10                     │
    │ Using Python 3.10.4. Type "help()" for help.                       │
    └────────────────────────────────────────────────────────────────────┘

    sage: load("setup")
    Loading log file <setup> one line at a time...
    Finished replaying log file <setup>
    sage: E
    Elliptic Curve defined by y^2 + x*y  = x^3 - x^2 + 4*x + 3 over Rational
    Field
    sage: x*y
    x*y
    sage: G
    [(2 : 3 : 1)]

Если вы используете Sage в ``konsole`` — терминале среды KDE в GNU/Linux —
тогда вы можете сохранить сессию следующим образом: после запуска Sage в
``konsole``, выберите "settings", потом "history...", потом "set unlimited".
Когда вы готовы сохранить сессию, выберите "edit" и "save history as..." и
введите имя файла для сохранения. После этого вы можете воспользоваться
любым текстовым редактором, например xemacs, для чтения файла.

Вставка игнорирует приглашение
==============================

Допустим, вы читаете сессию Sage или вычисления Python, и хотите скопировать
их в Sage. Но есть одна проблема: знаки ``>>>`` или ``sage:``. На самом деле
вы можете копировать и вставлять примеры, которые включают эти знаки. Дргуими
словами, Sage игнорирует символы ``>>>`` или ``sage:`` перед отправкой команд
в Python. Например,

.. skip

::

    sage: 2^10
    1024
    sage: sage: sage: 2^10
    1024
    sage: >>> 2^10
    1024

Команды измерения времени
=========================

Если вы введете команду ``%time`` в начале строки ввода, то время,
затраченное на выполнение операции, будет выведено на экран. Например, вы
можете измерить время выполнения операции возведения в степень несколькими
путями. Показания ниже будут отличаться от ваших; они могут отличаться даже
в разных версиях Sage. Чистый Python:

.. skip

::

    sage: %time a = int(1938)^int(99484)
    CPU times: user 0.66 s, sys: 0.00 s, total: 0.66 s
    Wall time: 0.66

Это означает что 0.66 секунд было затрачено в сумме, а "Wall time",
(прошедшее время), тоже 0.66 секунд. Если ваш компьютер сильно загружен
другими процессами, то "Wall time" может сильно отличаться от процессорного
времени.

Далее мы посчитаем время возведения в степень с использованием встроенного в
Sage целочисленного типа данных, реализованного (в Cython) с использованием
библиотеки GMP:

.. skip

::

    sage: %time a = 1938^99484
    CPU times: user 0.04 s, sys: 0.00 s, total: 0.04 s
    Wall time: 0.04

Используя интерфейс PARI из библиотеки C:

.. skip

::

    sage: %time a = pari(1938)^pari(99484)
    CPU times: user 0.05 s, sys: 0.00 s, total: 0.05 s
    Wall time: 0.05

GMP ведет себя лучше, но только немного (как и ожидалось, ведь версия PARI,
встроенная в Sage, использует GMP для работы с целыми числами).

Вы также можете замерить время выполнения блока команд с помощью ``cputime``,
как показано ниже:

::

    sage: t = cputime()
    sage: a = int(1938)^int(99484)
    sage: b = 1938^99484
    sage: c = pari(1938)^pari(99484)
    sage: cputime(t)                       # random output
    0.64

.. skip

::

    sage: cputime?
    ...
        Return the time in CPU second since Sage started, or with optional
        argument t, return the time since time t.
        INPUT:
            t -- (optional) float, time in CPU seconds
        OUTPUT:
            float -- time in CPU seconds

Команда ``walltime`` ведет себя так же, как ``cputime``, но она измеряет
настоящее время.

Мы также можем возвести число в степень, используя системы компьютерной
алгебры, включённые в Sage. В каждом случае мы запускаем простую команду
в системе чтобы запустить сервер для этой программы. Самое точное - время это
Wall time. Однако, если существует существенная разница между этим значением
и процессорным временем (CPU time), то, возможно, есть смысл проверить систему
на наличие проблем производительности.

.. skip

::

    sage: time 1938^99484;
    CPU times: user 0.01 s, sys: 0.00 s, total: 0.01 s
    Wall time: 0.01
    sage: gp(0)
    0
    sage: time g = gp('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: maxima(0)
    0
    sage: time g = maxima('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.30
    sage: kash(0)
    0
    sage: time g = kash('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.04
    sage: mathematica(0)
            0
    sage: time g = mathematica('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.03
    sage: maple(0)
    0
    sage: time g = maple('1938^99484')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.11
    sage: libgap(0)
    0
    sage: time g = libgap.eval('1938^99484;')
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 1.02

Заметьте, что GAP и Maxima являются самыми медленными в этом тесте (тест
был проведен на машине ``sage.math.washington.edu``). Так как они работают
с другим интерфейсом, надстроенным над ними, судить об абсолютной
производительности этих систем не стоит.

Ошибки и исключения
===================

Когда что-то идет не так, обычно можно увидеть исключение Python (Python
"exception"). Python даже попытается предположить, что вызвало ошибку. Часто
вы можете видеть имя исключения, например, :class:`NameError` или
:class:`ValueError` (см. Python Library Reference [PyLR]_ для полного списка
исключений). Например,

.. skip

::

    sage: 3_2
    ------------------------------------------------------------
       File "<console>", line 1
         ZZ(3)_2
               ^
    SyntaxError: invalid ...

    sage: EllipticCurve([0,infinity])
    ------------------------------------------------------------
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce Infinity (<class 'sage...Infinity'>) to Rational

Интерактивный отладчик может быть полезным для понимая того, что пошло не так.
Отладчик можно включать или выключать командой ``%pdb`` (по умолчанию он
выключен). Приглашение командной строки ``ipdb>`` появляется на экране,
если случилось исключение и отладчик был включен. Из отладчика вы можете
вывести на экран состояние любой локальной переменной и двигаться вверх и вниз
по стеку (execution stack). Например,

.. skip

::

    sage: %pdb
    Automatic pdb calling has been turned ON
    sage: EllipticCurve([1,infinity])
    ---------------------------------------------------------------------------
    <class 'exceptions.TypeError'>             Traceback (most recent call last)
    ...

    ipdb>

Для получения списка команд отладчика введите ``?`` в командной строке ``ipdb>``:

::

    ipdb> ?

    Documented commands (type help <topic>):
    ========================================
    EOF    break  commands   debug    h       l     pdef   quit    tbreak
    a      bt     condition  disable  help    list  pdoc   r       u
    alias  c      cont       down     ignore  n     pinfo  return  unalias
    args   cl     continue   enable   j       next  pp     s       up
    b      clear  d          exit     jump    p     q      step    w
    whatis where

    Miscellaneous help topics:
    ==========================
    exec  pdb

    Undocumented commands:
    ======================
    retval  rv

Нажмите Ctrl-D или введите ``quit`` чтобы вернуться в Sage.

.. _section-tabcompletion:

Обратный поиск и автодополнение
===============================

Сначала создадим трехмерное векторное пространство :math:`V=\QQ^3` следующим
образом:

::

    sage: V = VectorSpace(QQ,3)
    sage: V
    Vector space of dimension 3 over Rational Field

Можно использовать сокращенное обозначение:

::

    sage: V = QQ^3

Введите начало команды, потом нажмите ``Ctrl-p`` (или просто нажмите
стрелку вверх на клавиатуре) чтобы вернуться к любой из строк, которые
вы вводили, начинающейся с таких же символов. Это работает даже если вы
полность вышли из Sage и перезапустили его позже. Можно использовать и
обратный поиск по истории команд с помощью ``Ctrl-r``. Все эти возможности
используют пакет ``readline`` который доступен почти на всех разновидностях
GNU/Linux.

Можно с легкостью вывести список всех функций для :math:`V`, используя
автодополнение. Просто введите ``V.``, потом нажмите ``[TAB]`` на своей
клавиатуре:

.. skip

::

    sage: V.[tab key]
    V._VectorSpace_generic__base_field
    ...
    V.ambient_space
    V.base_field
    V.base_ring
    V.basis
    V.coordinates
    ...
    V.zero_vector

Если вы введете первые несколько символов команды, а потом нажмёте ``[TAB]``,
вы получите функции, которые начинаются с этих символов.

.. skip

::

    sage: V.i[tab key]
    V.is_ambient  V.is_dense    V.is_full     V.is_sparse

Если вам интересно, что делает какая-нибудь функция, например coordinates,
введите ``V.coordinates?`` для получения справки или ``V.coordinates??`` для
получения исходного кода (объясняется в следующем разделе).


Встроенная справочная система
=============================

Sage обладает встроенной справочной системой. Введите название функции со
знаком ? для доступа к документации по этой функции.

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates?
    Type:           instancemethod
    Base Class:     <class 'instancemethod'>
    String Form:    <bound method FreeModule_ambient_field.coordinates of Vector
    space of dimension 3 over Rational Field>
    Namespace:      Interactive
    File:           /home/was/s/local/lib/python2.4/site-packages/sage/modules/f
    ree_module.py
    Definition:     V.coordinates(self, v)
    Docstring:
        Write v in terms of the basis for self.

        Returns a list c such that if B is the basis for self, then

                sum c_i B_i = v.

        If v is not in self, raises an ArithmeticError exception.

        EXAMPLES:
            sage: M = FreeModule(IntegerRing(), 2); M0,M1=M.gens()
            sage: W = M.submodule([M0 + M1, M0 - 2*M1])
            sage: W.coordinates(2*M0-M1)
            [2, -1]

Как показано выше, вывод показывает тип объекта, файл, в котором он
определен и полезное описание функции с примерами, которые можно вставить
в вашу текущую сессию. Почти все примеры подвергаются регулярной
автоматической проверке на предмет работоспособности и наличия требуемого
поведения.

Другая возможность хорошо отражает дух открытого программного обеспечения:
если ``f`` это функция Python'а, то ``f??`` выведет исходный код, который
определяет ``f``. Например,

.. skip

::

    sage: V = QQ^3
    sage: V.coordinates??
    Type:           instancemethod
    ...
    Source:
    def coordinates(self, v):
            """
            Write $v$ in terms of the basis for self.
            ...
            """
            return self.coordinate_vector(v).list()

Отсюда мы знаем, что все, что делает функция ``coordinates``, это вызов функции
``coordinate_vector`` и превращает результат в список. Что делает функция
``coordinate_vector?``

.. skip

::

    sage: V = QQ^3
    sage: V.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            ...
            return self.ambient_vector_space()(v)

Функция ``coordinate_vector`` удерживает введенные значения во внешнем
пространстве, что позволяет добиться такого же эффекта, как при вычислении
вектора коэффициентов переменной :math:`v` с точки зрения :math:`V`.
Пространство :math:`V` уже внешнее, так как оно является :math:`\QQ^3`.
Существует также функция ``coordinate_vector`` для подпространств, и она
ведет себя по-иному. Мы создим подпространство и посмотрим:

.. skip

::

    sage: V = QQ^3; W = V.span_of_basis([V.0, V.1])
    sage: W.coordinate_vector??
    ...
    def coordinate_vector(self, v):
            """
             ...
            """
            # First find the coordinates of v wrt echelon basis.
            w = self.echelon_coordinate_vector(v)
            # Next use transformation matrix from echelon basis to
            # user basis.
            T = self.echelon_to_user_matrix()
            return T.linear_combination_of_rows(w)

(Если вы считаете, что существующая реализация неэффективна, пожалуйста,
зарегистрируйтесь и помогите оптимизировать линейную алгебру.)

Вы также можете ввести ``help(имя_команды)`` или ``help(класс)`` для
получения справки о классах или функциях в стиле man-страниц.

.. skip

::

    sage: help(VectorSpace)
    Help on function VectorSpace in module sage.modules.free_module:

    VectorSpace(K, dimension_or_basis_keys=None, sparse=False, inner_product_matrix=None, *,
                with_basis='standard', dimension=None, basis_keys=None, **args)
    EXAMPLES:

    The base can be complicated, as long as it is a field.

    ::

        sage: V = VectorSpace(FractionField(PolynomialRing(ZZ,'x')),3)
        sage: V
        Vector space of dimension 3 over Fraction Field of Univariate Polynomial Ring in x
         over Integer Ring
        sage: V.basis()
        [
        (1, 0, 0),
        (0, 1, 0),
    --More--

Когда вы вводите ``q`` для выхода из справочной системы, ваша сессия
находится в том же состоянии, что и до этого. Справка не захламляет ваш
экран, в отличие от формы ``function_name?``, которая иногда может
оставлять информацию в вашей сессии. Особенно полезно использовать
``help(module_name)``. Например, векторные пространства описаны в
``sage.modules.free_module``, поэтому введите ``help(sage.modules.free_module)``
для документации обо всем модуле. Когда вы просматриваете документацию в
справочной системе, вы можете осуществлять поиск с помощью ``/`` и в обратном
порядке с помощью ``?``.

Сохранение и загрузка отдельных объектов
========================================

Допустим вы вычислили матрицу или хуже: сложное пространство модулярных
символов, и хотите сохранить его для работы в будущем. Как это сделать? Есть
несколько способов, которыми компьютерные алгебры пользуются для сохранения
объектов.

#. **Сохранить игру:** Поддерживается сохранение и загрузка только полных сессий
   (например, GAP, Magma).

#. **Унифицированный ввод/вывод:** Вывод объектов на экран в таком виде, в
   котором они могут быть считаны позже. (GP/PARI).

#. **Eval**: Легкий способ запуска любого кода в интерпретаторе (например,
   Singular, PARI).


Так как Sage построен на Python'е, он использует иной подход: каждый объект
может быть превращен в строку, из которой в последствии можно восстановить объект.
Способ схож со способом унификации ввода и вывода, как в PARI, но в случае с
Sage нет необходимости выводить объект на экран в самой неудобной форме.
Также, поддержка сохранения и загрузки (в большинстве случаев) полностью
автоматична, не требует дополнительного программирования; это просто возможность
Python'а, которая была включена в язык с самого начала.

Почти любой объект x может быть сохранен в сжатой форме на диск при помощи
команды ''save(x, filename)'' (или во многих случаях ''x.save(filename)'').
Для загрузки объекта введите ''load(filename)''.

.. skip

::

    sage: A = MatrixSpace(QQ,3)(range(9))^2
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]
    sage: save(A, 'A')

Теперь выйдите из Sage и перезапустите. Теперь вы можете получить ''A'' обратно:

.. skip

::

    sage: A = load('A')
    sage: A
    [ 15  18  21]
    [ 42  54  66]
    [ 69  90 111]

То же самое можно делать и с более сложными объектами, например эллиптическими
кривыми. Вся информация об объекте (которая находится в кеше) сохраняется вместе
с объектом. Например,

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: v = E.anlist(100000)              # требует некоторого времени...
    sage: save(E, 'E')
    sage: quit

Сохраненная версия ``E`` занимает 153 килобита, так как в нем содержатся первые
100000 :math:`a_n`.

.. skip

::

    ~/tmp$ ls -l E.sobj
    -rw-r--r--  1 was was 153500 2006-01-28 19:23 E.sobj
    ~/tmp$ sage [...]
    sage: E = load('E')
    sage: v = E.anlist(100000)              # моментально!

(В Python, сохранение и загрузка осуществляется модулем ``cPickle``. Объект
Sage ``x`` может быть сохранен с помощью ``cPickle.dumps(x, 2)``. Обратите
внимание на ``2``!)

Sage не может сохранять и загружать объекты, созданные в других системах
компьютерной алгебры, таких как GAP, Singular, Maxima и пр. Они загружаются
в состоянии, которое помечено как "invalid". Хотя, в GAP многие объекты выводятся
в форме, из которой их потом можно восстановить, но многие не выводятся в
такой форме, поэтому их восстановление из такого вида нарочно запрещено.

.. skip

::

    sage: a = libgap(2)
    sage: a.save('a')
    sage: load('a')
    Traceback (most recent call last):
    ...
    ValueError: The session in which this object was defined is no longer
    running.

Объекты GP/PARI могут быть сохранены и загружены, так как их вид при выводе
на экран достаточен для восстановления объекта.

.. skip

::

    sage: a = gp(2)
    sage: a.save('a')
    sage: load('a')
    2

Сохраненные объекты могут быть загружены позже на компьютерах с другой
архитектурой или операционной системой, например, вы можете сохранить
огромную матрицу в 32-битной OS X и загрузить ее в 64-битную GNU/Linux,
привести к ступенчатой форме и переместить обратно. Также во многих случаях
вы можете загружать объекты в версии Sage, отличные от версии, на которой
они были сохранены. Все атрибуты объекта сохраняются вместе с классом (но не
включая исходный код), который описывает объект. Если класс более не существует
в новой версии Sage, тогда объект не может быть загружен в эту новую версию.
Но если вы загрузите ее на версию ниже, получите словарь объектов (с помощью
``x.__dict__``) и сохраните словарь, то сможете загрузить его в новую версию.

Сохранение в виде текста
------------------------

Вы также можете сохранять объекты в виде набора ASCII символов в простой
текстовый файл простым открытием файла и сохранением строки, которая выражает
(описывает) объект (вы можете записывать несколько объектов). Не забудьте
закрыть файл после добавления данных.

.. skip

::

    sage: R.<x,y> = PolynomialRing(QQ,2)
    sage: f = (x+y)^7
    sage: o = open('file.txt','w')
    sage: o.write(str(f))
    sage: o.close()

.. _section-save:

Сохранение и загрузка полных сессий
===================================

Sage обладает очень гибкими возможностями сохранения и загрузки полных сессий.

Команда ``save_session(sessionname)`` сохраняет все переменные, которые
вы задали в текущей сессии в виде словаря в заданном ``sessionname``. (В редком
случае, когда объект не поддерживает сохранения, он просто не будет включен
в словарь.) В результате будет создан файл с расширением ``.sobj`` и может быть
загружен как любой другой объект. Когда вы загружаете сохраненные объекты в
сессию, вы получаете словарь, ключами которого являются имена переменных, а
значениями — объекты.

Вы можете использовать команду ``load_session(sessionname)``, чтобы загрузить
переменные, описанные в ``sessionname``, в текущую сессию. Заметьте, что это
не удаляет переменные, заданные в этой сессии. Вместо этого, две сессии
объединяются.

Для начала запустим Sage и зададим несколько переменных.

.. skip

::

    sage: E = EllipticCurve('11a')
    sage: M = ModularSymbols(37)
    sage: a = 389
    sage: t = M.T(2003).matrix(); t.charpoly().factor()
     _4 = (x - 2004) * (x - 12)^2 * (x + 54)^2

Далее, сохраним нашу сессию, что включит в себя сохранение всех заданных
выше переменных в файл. Потом мы проверим информацию о файле. Его размер —
3 килобайта.

.. skip

::

    sage: save_session('misc')
    Saving a
    Saving M
    Saving t
    Saving E
    sage: quit
    was@form:~/tmp$ ls -l misc.sobj
    -rw-r--r--  1 was was 2979 2006-01-28 19:47 misc.sobj

Наконец, мы перезапустим Sage, зададим дополнительную переменную и загрузим
сохраненную сессию.

.. skip

::

    sage: b = 19
    sage: load_session('misc')
    Loading a
    Loading M
    Loading E
    Loading t

Каждая сохраненная переменная снова является переменной. Кроме того, переменная
``b`` не была перезаписана.

.. skip

::

    sage: M
    Full Modular Symbols space for Gamma_0(37) of weight 2 with sign 0
    and dimension 5 over Rational Field
    sage: E
    Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational
    Field
    sage: b
    19
    sage: a
    389

