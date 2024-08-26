.. _chapter-faq-contribute:

=========================
FAQ: Contributing to Sage
=========================


How can I start contributing to Sage?
"""""""""""""""""""""""""""""""""""""

The first step is to use Sage and encourage your friends to use
Sage. If you find bugs or confusing documentation along the way,
please report your problems!

Two popular ways to contribute to Sage are to write code and to
create documentation or tutorials. Some steps in each direction
are described below.


I want to contribute code to Sage. How do I get started?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Take a look at the
`Sage Developer Guide <https://doc.sagemath.org/html/en/developer>`_.
At a minimum, the first chapter in that guide is required
reading for any Sage developer. Also pay special attention to the
`GitHub guidelines <https://doc.sagemath.org/html/en/developer/github.html>`_.
You can also join the
`sage-devel <https://groups.google.com/group/sage-devel>`_
mailing list or hang around on the
`SageMath Zulip chat channel <https://sagemath.zulipchat.com/>`_.
While you are getting to know the community, grab a copy of the Sage
source and familiarize yourself with the
`git <https://git-scm.com>`_ version control system.

The best way to become familiar with the Sage development process is
to choose an issue from
`the Sage repository on GitHub <https://github.com/sagemath/sage/issues>`_
and review the proposed changes contained in that issue. If you want
to implement something, it is a good practice to discuss your ideas on
the ``sage-devel`` mailing list first, so that other developers have a
chance to comment on your ideas/proposals. They are pretty open to new
ideas, too, as all mathematicians should be.

Sage's main programming language is
`Python <https://www.python.org>`_.
Some parts of Sage may be written in other languages, especially the
components that do the heavy number crunching, but most native
functionality is done using Python, including "glue code". One of the
good aspects of Python that Sage inherits is that working code is
considered more valuable than merely fast code. Fast code is valuable,
but clean, readable code is important. In the mathematics community,
inaccurate results are unacceptable. Correctness comes before
optimization. In the following paper

* D. Knuth. Structured Programming with go to Statements.
  *ACM Journal Computing Surveys*, 6(4), 1974.

Don Knuth observes that: "We should forget about small efficiencies,
say about 97% of the time: premature optimization is the root of all
evil."

If you do not know Python, you should start learning that language. A
good place to start is the
`Python Official Tutorial <https://docs.python.org/3/tutorial>`_
and other documents in the
`Python standard documentation <https://docs.python.org>`_.
Another good place to take a look at is
`Dive Into Python <https://diveintopython3.net>`_
by Mark Pilgrim, which may be pretty helpful on some specific topics
such as test-driven development. The book
`Building Skills in Python <http://itmaybeahack.com/homepage/books/python.html>`_
by Steven F. Lott is suitable for anyone who is already comfortable
with programming.

If you want, you can
try to learn Python by using Sage. However,
it is helpful to know what is pure Python and when Sage is doing its
"magic". There are many things that work in Python but not in Sage,
and vice versa. Furthermore, even when the syntax is identical, many
programming concepts are explained more thoroughly in Python-centered
resources than in Sage-centered resources; in the latter,
mathematics is usually the priority.


I am not a programmer. Is there another way I can help out?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Yes. As with any free open source software project, there are numerous
ways in which you could help out within the Sage community, and
programming is only one of many ways to contribute.

Many people like writing technical tutorials. One of the joys of doing
so is that you also learn something new in the process. At the same
time, you communicate your knowledge to beginners, a skill which is
useful in fields other than technical writing itself. A main point
about technical writing is that you communicate a technical subject to
beginners, so keep technical jargon to a minimum. Darrell Anderson
has written some
`tips on technical writing <http://web.archive.org/web/20130128102724/http://humanreadable.nfshost.com:80/howtos/technical_writing_tips.htm>`_,
which we highly recommend.

For the graphic designers or the artistically creative, you can
help out with improving the design of the Sage website.

If you can speak, read,
and write in another (natural) language, there are many ways in which
your contribution would be very valuable to the whole Sage
community. Say you know Italian. Then you can write a Sage tutorial in
Italian, or help out with translating the official Sage tutorial to
Italian.

The above is a very short list. There are many, many more ways in
which you can help out. Feel free to send an email to the
`sage-devel <https://groups.google.com/group/sage-devel>`_ mailing list
to ask about possible ways in which you could help out, or to suggest a
project idea.


Where can I find resources on Python or Cython?
"""""""""""""""""""""""""""""""""""""""""""""""

If you are new to Python, you can start with the `Official Python Tutorial <https://docs.python.org/3/tutorial/>`_ or one of numerous
free tutorials and courses out on the web.

To learn about Cython, start with the `Basic Cython Tutorial <https://cython.readthedocs.io/en/latest/src/tutorial/cython_tutorial.html>`_.

It is a good idea to learn about good developer tools
alongside with learning about programming in these languages.
Our developer guide explains how to :ref:`write testable examples
<section-doctest-writing>` and has information on the
:ref:`development and testing tools <chapter-tools>` that
are in use in Sage development.


Are there any coding conventions I need to follow?
""""""""""""""""""""""""""""""""""""""""""""""""""

See the Sage Developer's Guide, especially the chapter
:ref:`chapter-code-basics`.


I submitted a bug fix to the GitHub Sage repo several weeks ago. Why is it being ignored?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

There are several possible reasons why a pull request might not get attention:

Most people who work on Sage do so in their free time.
There are many open pull requests that keep contributors busy.
Understanding the changes requires knowledge that only very few people have.

It is thus important to make the job for potential reviewers as easy as possible.
Here are some tips on making your PR easy to review:

* Clearly describe the problem your PR is trying to
  solve.
* Provide background information relevant to the problem
  that your PR is solving. Such information can include links to
  online resources and any relevant papers, books and reference
  materials.
* Clearly describe how your PR solves the problem under
  consideration.
* Clearly describe how to test the changes in your PR.
* List any Issues or Pull Requests that your PR depends on.
* Make sure your PR is based on a recent (preferably, the latest) Sage beta version.
* Follow the :ref:`relevant coding conventions <chapter-code-basics>`
  as documented in the Developer Guide.


When and how might I remind the Sage community of a PR I care about?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

It is appropriate to join the
`sage-devel <https://groups.google.com/group/sage-devel>`_
mailing list and to post there about the PRs that you care about.
Although it may at times feel intimidating to post there, many
developers are eager to join the discusssion and help.

You can also try to find developers interested in reviewing your PRs
in our `SageMath Zulip chat channel <https://sagemath.zulipchat.com/>`_.
