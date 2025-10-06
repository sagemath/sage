*********************
Использование SageTeX
*********************

Встроенный в Sage пакет SageTeX позволяет внедрять результаты вычислений в
документ типа LaTeX. Для того, чтобы использовать данный пакет, понадобится
"установить" его в локальную систему TeX (под "установкой" подразумевается
копирование одного файла). См. :ref:`installation`, а также раздел "Make
SageTeX known to TeX" `Руководства по установке Sage
<http://doc.sagemath.org/html/en/installation/index.html>`_ (`данная ссылка
<../../en/installation/index.html>`_ ведет к локальному размещению копии руководства
по установке).

В этом уроке показан небольшой пример использования SageTeX. Полная документация
находится в :file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`, где
``SAGE_ROOT`` - это директория, в которой установлен Sage. Эта папка содержит
документацию, файл с примером и полезные скрипты Python.

Для начала работы с SageTeX следуйте указаниям по установке (в :ref:`installation`)
и вставьте следующия текст в файл названный, скажем, ``st_example.tex``:

.. warning::

  Нижеследующий текст может содержать несколько сообщений об ошибках, связанных
  с неизвестными управляющими последовательностями, если Вы используете
  интерактивную помощь. Откройте статическую версию, чтобы увидеть правильный
  текст.

.. code-block:: latex

    \documentclass{article}
    \usepackage{sagetex}

    \begin{document}

    Using Sage\TeX, one can use Sage to compute things and put them into
    your \LaTeX{} document. For example, there are
    $\sage{number_of_partitions(1269)}$ integer partitions of $1269$.
    You don't need to compute the number yourself, or even cut and paste
    it from somewhere.

    Here's some Sage code:

    \begin{sageblock}
        f(x) = exp(x) * sin(2*x)
    \end{sageblock}

    The second derivative of $f$ is

    \[
      \frac{\mathrm{d}^{2}}{\mathrm{d}x^{2}} \sage{f(x)} =
      \sage{diff(f, x, 2)(x)}.
    \]

    Here's a plot of $f$ from $-1$ to $1$:

    \sageplot{plot(f, -1, 1)}

    \end{document}

Запустите LaTeX для ``st_example.tex``. Заметьте, что LaTeX будет жаловаться
на некоторые вещи, как-то:

    Package sagetex Warning: Graphics file
    sage-plots-for-st_example.tex/plot-0.eps on page 1 does not exist. Plot
    command is on input line 25.

    Package sagetex Warning: There were undefined Sage formulas and/or
    plots. Run Sage on st_example.sage, and then run LaTeX on
    st_example.tex again.

Среди файлов, сгенерированных после запуска LaTeX, есть файл ``st_example.sage``,
являющийся скриптом Sage. Сообщение, показанное выше, предлагало запустить
``st_example.sage``, поэтому стоит так и сделать. Затем будет предложено
запустить LaTeX для ``st_example.tex`` еще раз; перед этим будет создан файл
``st_example.sout``. Этот файл содержит результаты вычислений в Sage в формате,
удобном для LaTeX. Новая папка, содержащая EPS-файл с графиком, также будет
создана автоматически. Запустите LaTeX еще раз: все, что было вычислено в Sage,
теперь включено в Ваш документ.

Перечислим шаги:

    - запустите LaTeX для .tex файла;
    - запустите Sage для сгенерированного .sage файла;
    - запустите LaTeX еще раз.

Пункт с запуском Sage можно пропустить, если никакие изменения не были
применены к командам Sage в документе.

SageTeX предлагает много возможностей, и так как Sage и LaTeX являются
мощными инструментами, то стоит изучить
:file:`SAGE_ROOT/venv/share/texmf/tex/latex/sagetex`.
