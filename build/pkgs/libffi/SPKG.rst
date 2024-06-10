libffi: A portable foreign-function interface library
=====================================================

Description
-----------

Compilers for high level languages generate code that follow certain
conventions. These conventions are necessary, in part, for separate
compilation to work. One such convention is the "calling convention".
The "calling convention" is essentially a set of assumptions made by the
compiler about where function arguments will be found on entry to a
function. A "calling convention" also specifies where the return value
for a function is found.

Some programs may not know at the time of compilation what arguments are
to be passed to a function. For instance, an interpreter may be told at
run-time about the number and types of arguments used to call a given
function. Libffi can be used in such programs to provide a bridge from
the interpreter program to compiled code.

The libffi library provides a portable, high level programming interface
to various calling conventions. This allows a programmer to call any
function specified by a call interface description at run time.

FFI stands for Foreign Function Interface. A foreign function interface
is the popular name for the interface that allows code written in one
language to call code written in another language. The libffi library
really only provides the lowest, machine dependent layer of a fully
featured foreign function interface. A layer must exist above libffi
that handles type conversions for values passed between the two
languages.

License
-------

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Upstream Contact
----------------

- https://sourceware.org/libffi/
- https://github.com/libffi/libffi
