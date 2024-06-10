
bernmm: an implementation of the algorithm described in "A multimodular
        algorithm for computing Bernoulli numbers", by David Harvey, 2008.

version 1.1

Copyright (C) 2008, 2009, David Harvey

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

===============================================================================

To build bernmm, you need to edit the makefile to point to your GMP and NTL
install directories, and then just do "make". This will create an executable
"bernmm-test". Run it with no arguments to see the available options.

In the makefile there is also a THREAD_STACK_SIZE option; the default is
4096 KB. If this option is not set, the system default thread stack size will
be used. This proved to be a problem on OS X, where GMP's mpz_invert routine
smashed the default 512 KB stack. For more info see
https://github.com/sagemath/sage/issues/6304
http://gmplib.org/list-archives/gmp-devel/2009-August/000929.html

===============================================================================
