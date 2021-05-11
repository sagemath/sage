yasm: An assembler for the x86 and AMD64 instruction sets
=========================================================

Description
-----------

Yasm is a complete rewrite of the NASM assembler under the “new” BSD
License (some portions are under other licenses, see COPYING for
details).

Yasm currently supports the x86 and AMD64 instruction sets, accepts NASM
and GAS assembler syntaxes, outputs binary, ELF32, ELF64, 32 and 64-bit
Mach-O, RDOFF2, COFF, Win32, and Win64 object formats, and generates
source debugging information in STABS, DWARF 2, and CodeView 8 formats.

Yasm can be easily integrated into Visual Studio 2005/2008 and 2010 for
assembly of NASM or GAS syntax code into Win32 or Win64 object files.

See https://yasm.tortall.net

License
-------

Yasm is licensed under the 2-clause and 3-clause “revised” BSD licenses,
with one exception: the Bit::Vector module used by the mainline version
of Yasm to implement its large integer and machine-independent floating
point support is triple-licensed under the Artistic license, GPL, and
LGPL. The “yasm-nextgen” codebase uses a different BSD-licensed
implementation and is thus entirely under BSD-equivalent licenses. The
full text of the licenses are provided in the Yasm source distribution.


Upstream Contact
----------------

-  https://yasm.tortall.net

Dependencies
------------

-  none
