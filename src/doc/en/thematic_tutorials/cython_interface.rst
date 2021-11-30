.. nodoctest

.. _cython_interface:

========================================================
How to call a C code (or a compiled library) from Sage ?
========================================================

If you have some C/C++ code that you would like to call from Sage for your own
use, this document is for you.

- Do you want to **contribute** to Sage by adding your interface to its code? The
  (more complex) instructions are `available here
  <http://doc.sagemath.org/html/en/developer/index.html#packaging-third-party-code>`_.

.. _section-cython-interface-helloworld:

Calling "hello_world()" from hello.c
------------------------------------

Let us suppose that you have a file named ``~/my_dir/hello.c`` containing:

.. CODE-BLOCK:: c

  #include <stdio.h>

  void hello_world(){
      printf("Hello World\n");
  }

In order to call this function from Sage, you must create a Cython file (i.e. a
file whose extension is .pyx). Here, ``~/my_dir/hello_sage.pyx`` contains a
header describing the signature of the function that you want to call:

.. CODE-BLOCK:: cython

  cdef extern from "hello.c":
      void hello_world()

  def my_bridge_function():
      hello_world() # This is the C function from hello.c

You can now load this file in Sage, and call the C code though the Python
function ``my_bridge_function``::

  sage: %runfile hello_sage.pyx
  Compiling ./hello_sage.pyx...
  sage: my_bridge_function()
  Hello World

Arguments and return value
--------------------------

Calling function with more complex arguments and return values works the same
way. To learn more about the Cython language, `click here
<http://docs.cython.org/src/reference/language_basics.html>`_

The following example defines a function taking and returning ``int *``
pointers, and involves some memory allocation. The C code defines a function
whose purpose is to return the sum of two vectors as a third vector.

**The C file** (``double_vector.c``)

.. CODE-BLOCK:: c

  #include <string.h>

  int * sum_of_two_vectors(int n, int * vec1, int * vec2){
    /*
       INPUT : two arrays vec1,vec2 of n integers
       OUTPUT: an array of size n equal to vec1+vec2
    */
    int * sum = (int *) malloc(n*sizeof(int));
    int i;

    for(i=0;i<n;i++)
      sum[i] = vec1[i] + vec2[i];
    return sum;
  }

**The Cython file** (``double_vector_sage.pyx``)

.. CODE-BLOCK:: cython

  cdef extern from "double_vector.c":
      int * sum_of_two_vectors(int n, int * vec1, int * vec2)

  from libc.stdlib cimport malloc, free

  def sage_sum_of_vectors(n, list1, list2):
      cdef int * vec1 = <int *> malloc(n*sizeof(int))
      cdef int * vec2 = <int *> malloc(n*sizeof(int))

      # Fill the vectors
      for i in range(n):
          vec1[i] = list1[i]
          vec2[i] = list2[i]

      # Call the C function
      cdef int * vec3 = sum_of_two_vectors(n,vec1,vec2)

      # Save the answer in a Python object
      answer = [vec3[i] for i in range(n)]

      free(vec1)
      free(vec2)
      free(vec3)

      return answer


**Call from Sage**::

  sage: %runfile double_vector_sage.pyx
  Compiling ./double_vector_sage.pyx...
  sage: sage_sum_of_vectors(3,[1,1,1],[2,3,4])
  [3, 4, 5]

Calling code from a compiled library
------------------------------------

The procedure is very similar again. For our purposes, we build a library from
the file ``~/my_dir/hello.c``:

.. CODE-BLOCK:: c

   #include <stdio.h>

   void hello_world(){
       printf("Hello World\n");
  }

We also need a ``~/my_dir/hello.h`` header file:

.. CODE-BLOCK:: c

   void hello_world();

We can now **compile it** as a library:

.. CODE-BLOCK:: shell-session

   [user@localhost ~/my_dir/] gcc -c -Wall -Werror -fpic hello.c
   [user@localhost ~/my_dir/] gcc -shared -o libhello.so hello.o

The only files that we need now are ``hello.h`` and ``libhello.so`` (you can
remove the others if you like). We must now indicate the location of the ``.so``
and ``.h`` files in the header of our ``~/my_dir/hello_sage.pyx`` file:

.. CODE-BLOCK:: cython

   # distutils: libraries = /home/username/my_dir/hello

   cdef extern from "hello.h":
       void hello_world()

   def my_bridge_function():
       hello_world() # This is the C function from hello.c

.. NOTE::

   The instruction ``# distutils: libraries = /home/username/my_dir/hello``
   indicates that the library is actually named ``/home/username/my_dir/hello``.
   Change it according to your needs.
   For more information about these instructions, see
   http://cython.readthedocs.io/en/latest/src/reference/compilation.html#configuring-the-c-build

We can now **load** this file in Sage and **call** the function::

   sage: %runfile hello_sage.pyx
   Compiling ./hello_sage.pyx...
   sage: my_bridge_function()
   Hello World
