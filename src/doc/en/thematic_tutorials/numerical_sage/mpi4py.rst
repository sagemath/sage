mpi4py
======

MPI, which stands for Message Passing Interface, is a common library
for parallel programming. There is a package ``mpi4py`` that builds on
the top of MPI, and lets arbitrary python objects be passed between
different processes. These packages are not available from the
Sage distribution. Install ``openmpi`` using your distribution's
package manager. Then install ``mpi4py`` using

.. skip

::

    sage: !pip install mpi4py

Now, the way that MPI works is you start a group of MPI processes,
all of the processes run the same code. Each process has a rank,
that is a number that identifies it. The following pseudocode
indicates the general format of MPI programs.

.. CODE-BLOCK:: text

       ....

    if my rank is n:
       do some computation ...
       send some stuff to the process of rank j
       receive some data from the process of rank k

    else if my rank is n+1:
       ....

Each process looks for what it's supposed to do (specified by its
rank), and processes can send data and receive data. Let's give an
example. Create a script with the following code in a file ``mpi_1.py``

.. CODE-BLOCK:: python

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    print("hello world")
    print(f"my rank is: {comm.rank}")

To run it you can do (from the command line in your Sage
directory)

.. CODE-BLOCK:: shell-session

    mpirun -np 5 ./sage -python mpi_1.py

The command ``mpirun -np 5`` starts 5 copies of a program under MPI. In
this case we have 5 copies of Sage in pure Python mode running the
script ``mpi_1.py``. The result should be 5 "hello worlds" plus 5 distinct ranks.

The two most important MPI operations are sending and receiving.
Consider the following example which you should put in a script ``mpi_2.py``

.. CODE-BLOCK:: python

    from mpi4py import MPI
    import numpy
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    v = numpy.array([rank] * 5, dtype=float)
    comm.send(v, dest=(rank+1) % size)
    data = comm.recv(source=(rank-1) % size)
    print(f"my rank is: {rank}")
    print("I received this:")
    print(data)

The same command as above with ``mpi_1.py`` replaced by ``mpi_2.py`` will
produce 5 outputs. Each process will create an array and pass
it to the next process, where the last process passes to the
first. Note that ``MPI.size`` is the total number of MPI
processes. ``MPI.COMM_WORLD`` is the communication world.

There are some subtleties regarding MPI to be aware of. Small sends
are buffered. This means if a process sends a small object it will
be stored by openmpi and that process will continue its execution
and the object it sent will be received whenever the destination
executes a receive. However, if an object is large, a process will
hang until its destination executes a corresponding receive. In
fact, the above code will hang if ``[rank]*5`` is replaced by
``[rank]*500``. It would be better to do

.. CODE-BLOCK:: python

    from mpi4py import MPI
    import numpy
    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size
    v = numpy.array([rank] * 500, dtype=float)
    if comm.rank == 0:
       comm.send(v, dest=(rank+1) % size)
    if comm.rank > 0:
        data = comm.recv(source=(rank-1) % size)
        comm.send(v, dest=(rank+1) % size)
    if comm.rank == 0:
        data = comm.recv(source=size - 1)

    print(f"my rank is: {rank}")
    print("I received this:")
    print(data)

Now, process 0 sends the data to process 1, then waits to receive from
process ``MPI.size - 1``.  Simultaneously, process 1 will send the
data to process 2, then receives the data from process 0.  This will
not lock even if the array transmitted is huge.

A common idiom is to have one process, usually the one with rank 0,
act as a leader. That process sends data out to the other
processes, compute on the results, and decides how much further
computation should proceed. Consider the following code

.. CODE-BLOCK:: python

    from mpi4py import MPI
    import numpy
    sendbuf = []
    root = 0
    comm = MPI.COMM_WORLD
    if comm.rank == 0:
        m = numpy.random.randn(comm.size, comm.size)
        print(m)
        sendbuf=m

    v = comm.scatter(sendbuf, root)

    print("I got this array:")
    print(v)

The ``scatter`` command takes a list and evenly divides it amongst all
the processes. Here the root process creates a matrix (which is
viewed as a list of rows) and then scatters it to everybody (root's
``sendbuf`` is divided equally amongst the processes). Each process
prints the row it got. Note that the ``scatter`` command is executed by
everyone, but when root executes it, it acts as a ``send`` and a
``receive`` (root gets one row from itself), while for everyone else it
is just a ``receive``.

There is a complementary ``gather`` command that collects results from
all the processes into a list. The next example uses ``scatter`` and
``gather`` together. Now the root process scatters the rows of a
matrix. Each process squares the elements of the row it receives.
The root process then gathers the rows into a new matrix.

.. CODE-BLOCK:: python

    from mpi4py import MPI
    import numpy
    comm = MPI.COMM_WORLD
    sendbuf = []
    root = 0
    if comm.rank == 0:
        m = numpy.array(range(comm.size * comm.size), dtype=float)
        m.shape = (comm.size, comm.size)
        print(m)
        sendbuf = m

    v = comm.scatter(sendbuf, root)
    print("I got this array:")
    print(v)
    v = v*v
    recvbuf = comm.gather(v, root)
    if comm.rank == 0:
        print(numpy.array(recvbuf))

There is also a ``broadcast`` command that sends a single object to
every process. Consider the following small extension. This is the
same as before, but now at the end, the root process sends everyone
the string "done", which is printed out.

.. CODE-BLOCK:: python

    v = MPI.COMM_WORLD.scatter(sendbuf, root)
    print("I got this array:")
    print(v)
    v = v*v
    recvbuf = MPI.COMM_WORLD.gather(v, root)
    if MPI.COMM_WORLD.rank == 0:
        print(numpy.array(recvbuf))

    if MPI.COMM_WORLD.rank == 0:
        sendbuf = "done"
    recvbuf = MPI.COMM_WORLD.bcast(sendbuf,root)
    print(recvbuf)
