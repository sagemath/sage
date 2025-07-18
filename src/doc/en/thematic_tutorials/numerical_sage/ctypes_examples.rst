More complicated ctypes example
===============================

Here we will look at a more complicated example. First consider the
following C code.

.. code-block:: c

    #include <stdio.h>
    #include <stdlib.h>

    struct double_row_element_t {
      double value;
      int  col_index;
      struct double_row_element_t * next_element;
    };

    typedef struct double_row_element_t double_row_element;

    typedef struct  {
      int nrows;
      int ncols;
      int nnz;
      double_row_element** rows;
    } double_sparse_matrix;


    double_sparse_matrix * initialize_matrix(int nrows, int ncols)
    {
      int i;
      double_sparse_matrix* new_matrix;
      new_matrix = (double_sparse_matrix *) malloc(sizeof(double_sparse_matrix));
      new_matrix->rows= (double_row_element **) malloc(sizeof(double_row_element *)*nrows);
      for(i=0;i<nrows;i++)
        {
          (new_matrix->rows)[i]=(double_row_element *) malloc(sizeof(double_row_element));
          (new_matrix->rows)[i]->value=0;
          (new_matrix->rows)[i]->col_index=0;
          (new_matrix->rows)[i]->next_element = 0;
        }
      new_matrix->nrows=nrows;
      new_matrix->ncols=ncols;
      new_matrix->nnz=0;
      return new_matrix;
    }


    int free_matrix(double_sparse_matrix * matrix)
    {
      int i;
      double_row_element* next_element;
      double_row_element* current_element;
      for(i=0;i<matrix->nrows;i++)
        {
          current_element = (matrix->rows)[i];
          while(current_element->next_element!=0)
        {
          next_element=current_element->next_element;
          free(current_element);
          current_element=next_element;
        }
          free(current_element);
        }
      free(matrix->rows);
      free(matrix);
      return 1;
    }

    int set_value(double_sparse_matrix * matrix,int row, int col, double value)
    {

      int i;
      i=0;
      double_row_element* current_element;
      double_row_element* new_element;

      if(row> matrix->nrows || col > matrix->ncols || row <0 || col <0)
        return 1;

      current_element = (matrix->rows)[row];
      while(1)
        {
          if(current_element->col_index==col)
        {
        current_element->value=value;
        return 0;
        }

          else
        if(current_element->next_element!=0)
        {
              if(current_element->next_element->col_index <=col)
            current_element = current_element->next_element;
          else
            if(current_element->next_element->col_index > col)
              {
            new_element = (double_row_element *) malloc(sizeof(double_row_element));
            new_element->value=value;
            new_element->col_index=col;
            new_element->next_element=current_element->next_element;
            current_element->next_element=new_element;
            return 0;
              }
        }
        else
          {
            new_element = (double_row_element *) malloc(sizeof(double_row_element));
            new_element->value=value;
            new_element->col_index=col;
            new_element->next_element=0;
            current_element->next_element=new_element;
            break;
          }

        }

      return 0;
    }


    double get_value(double_sparse_matrix* matrix,int row, int col)
    {
      int i;
      double_row_element * current_element;
      if(row> matrix->nrows || col > matrix->ncols || row <0 || col <0)
        return 0.0;

      current_element = (matrix->rows)[row];
      while(1)
        {
          if(current_element->col_index==col)
        {
          return current_element->value;
        }
          else
        {
          if(current_element->col_index<col && current_element->next_element !=0)
            current_element=current_element->next_element;
          else
            if(current_element->col_index >col || current_element ->next_element==0)
               return 0;
        }
        }

    }

Put it in a file called linked_list_sparse.c and compile it using

.. CODE-BLOCK:: shell-session

    $ gcc -c linked_list_sparse.c
    $ gcc -shared -o linked_list_sparse.so linked_list_sparse.o

Next consider the following python helper code.

.. CODE-BLOCK:: python

    from ctypes import *

    class double_row_element(Structure):
        pass

    double_row_element._fields_=[("value",c_double),("col_index",c_int),("next_element",POINTER(double_row_element) )]


    class double_sparse_matrix(Structure):
        _fields_=[("nrows",c_int),("ncols",c_int),("nnz",c_int),("rows",POINTER(POINTER(double_row_element)))]


    double_sparse_pointer=POINTER(double_sparse_matrix)
    sparse_library=CDLL("/home/jkantor/linked_list_sparse.so")
    initialize_matrix=sparse_library.initialize_matrix
    initialize_matrix.restype=double_sparse_pointer
    set_value=sparse_library.set_value
    get_value=sparse_library.get_value
    get_value.restype=c_double
    free_matrix=sparse_library.free_matrix

Let's discuss the above code. The original C code stored a sparse
matrix as a linked list. The python code uses the ctypes Structure
class to create structures mirroring the structs in the C code. To
create python object representing a C struct, simply create class that
derives from Structure. The _fields_ attribute of the class must be set
to a list of tuples of field names and values. Note that in case you
need to refer to a struct before it is completely defined (as in the
linked list) you can first declare it with "Pass", and then specify
the field contents as above. Also note the POINTER operator which
creates a pointer out of any ctypes type. We are able to directly call
our library as follows.

.. CODE-BLOCK:: python

    m=double_sparse_pointer()
    m=initialize_matrix(c_int(10),c_int(10))
    set_value(m,c_int(4),c_int(4),c_double(5.0))
    a=get_value(m,c_int(4),c_int(4))
    print("%f"%a)
    free_matrix(m)

Note that you can access the contents of a structure just by
(struct_object).field name. However for pointers, there is a contents
attribute. So, in the above, m.contents.nrows would let you access the
nrows field.  In fact you can manually walk along the linked list as
follows.

.. CODE-BLOCK:: python

    m=double_sparse_pointer()
    m=initialize_matrix(c_int(10),c_int(10))
    set_value(m,c_int(4),c_int(4),c_double(5.0))
    a=m.contents.rows[4]
    b=a.contents.next_element
    b.contents.value
    free_matrix(m)
