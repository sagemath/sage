mpi4py
======

MPI which stands for message passing interface is a common library
for parallel programming. There is a package mpi4py that builds on
the top of mpi, and lets arbitrary python objects be passed between
different processes. These packages are not part of the default
sage install. To install them do

.. skip

::

    sage: optional_packages()

Find the package name openmpi-\* and mpi4py-\*and do

.. skip

::

    sage: install_package('openmpi-*')
    sage: install_package('mpi4py-*')

Note that openmpi takes a while to compile (15-20 minutes or so).
Openmpi can be run on a cluster, however this requires some set up
so that processes on different machines can communicate (though if
you are on a cluster this is probably already set up). The simplest
case is if you are on a shared memory or multicore system where
openmpi will just work with no configuration from you. To be
honest, I have never tried to run mpi4py on a cluster, though there
is much information about these topics online.

Now, the way that mpi works is you start a group of mpi processes,
all of the processes run the same code. Each process has a rank,
that is a number that identifies it. The following pseudocode
indicates the general format of MPI programs.

::

       ....

    if my rank is n:
       do some somputation ...
       send some stuff to the process of rank j
       receive some data from the process of rank k

    else if my rank is n+1:
       ....

Each processes looks for what it's supposed to do (specified by its
rank) and processes can send data and receive data. Lets give an
example. Create a script with the following code in a file mpi_1.py

::

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    print("hello world")
    print("my rank is: %d"%comm.rank)

To run it you can do (from the command line in your sage
directory)

::

    ./local/bin/mpirun -np 5 ./sage -python mpi_1.py

The command mpirun -np 5 starts 5 copies of a program under mpi. In
this case we have 5 copies of sage in pure python mode run the
script mpi_1.py. The result should be 5 "hello worlds" plus 5 distinct ranks.
The two most important mpi operations are sending and receiving.
Consider the following example which you should put in a script mpi_2.py

::

    from mpi4py import MPI
    import numpy
    comm = MPI.COMM_WORLD
    rank=comm.rank
    size=comm.size
    v=numpy.array([rank]*5,dtype=float)
    comm.send(v,dest=(rank+1)%size)
    data=comm.recv(source=(rank-1)%size)
    print("my rank is %d"%rank)
    print("I received this:")
    print(data)

The same command as above with mpi_1.py replaced by mpi_2.py will
produce 5 outputs and you will see each process creates an array and
then passes it to the next guy (where the last guy passes to the
first.) Note that MPI.size is the total number of mpi
processes. MPI.COMM WORLD is the communication world.

There are some subtleties regarding MPI to be aware of. Small sends
are buffered. This means if a process sends a small object it will
be stored by openmpi and that process will continue its execution
and the object it sent will be received whenever the destination
executes a receive. However, if an object is large a process will
hang until its destination executes a corresponding receive. In
fact the above code will hang if [rank]\*5 is replaced by
[rank]\*500. It would be better to do

::

    from mpi4py import MPI
    import numpy
    comm = MPI.COMM_WORLD
    rank=comm.rank
    size=comm.size
    v=numpy.array([rank]*500,dtype=float)
    if comm.rank==0:
       comm.send(v,dest=(rank+1)%size)
    if comm.rank > 0:
        data=comm.recv(source=(rank-1)%size)
        comm.send(v,dest=(rank+1)%size)
    if comm.rank==0:
        data=comm.recv(source=size-1)

    print("my rank is %d"%rank)
    print("I received this:")
    print(data)

Now the first process initiates a send, and then process 1 will be
ready to receive and then he will send and process 2 will be
waiting to receive, etc. This will not lock regardless of how large
of an array we pass.

A common idiom is to have one process, usually the one with rank 0
act as a leader. That processes sends data out to the other
processes and processes the results and decides how further
computation should proceed. Consider the following code

::

    from mpi4py import MPI
    import numpy
    sendbuf=[]
    root=0
    comm = MPI.COMM_WORLD
    if comm.rank==0:
        m=numpy.random.randn(comm.size,comm.size)
        print(m)
        sendbuf=m

    v=comm.scatter(sendbuf,root)

    print("I got this array:")
    print(v)

The scatter command takes a list and evenly divides it amongst all
the processes. Here the root process creates a matrix (which is
viewed as a list of rows) and then scatters it to everybody (roots
sendbuf is divided equally amongst the processes). Each process
prints the row it got. Note that the scatter command is executed by
everyone, but when root executes it, it acts as a send and a
receive (root gets one row from itself), while for everyone else it
is just a receive.

There is a complementary gather command that collects results from
all the processes into a list. The next example uses scatter and
gather together. Now the root process scatters the rows of a
matrix, each process then squares the elements of the row it gets.
Then the rows are all gathered up again by the root process who
collects them into a new matrix.

::

    from mpi4py import MPI
    import numpy
    comm = MPI.COMM_WORLD
    sendbuf=[]
    root=0
    if comm.rank==0:
        m=numpy.array(range(comm.size*comm.size),dtype=float)
        m.shape=(comm.size,comm.size)
        print(m)
        sendbuf=m

    v=comm.scatter(sendbuf,root)
    print("I got this array:")
    print(v)
    v=v*v
    recvbuf=comm.gather(v,root)
    if comm.rank==0:
        print(numpy.array(recvbuf))

There is also a broadcast command that sends a single object to
every process. Consider the following small extension. This is the
same as before, but now at the end the root process sends everyone
the string "done", which is printed out.

::

    v=MPI.COMM_WORLD.scatter(sendbuf,root)
    print("I got this array:")
    print(v)
    v=v*v
    recvbuf=MPI.COMM_WORLD.gather(v,root)
    if MPI.COMM_WORLD.rank==0:
        print(numpy.array(recvbuf))

    if MPI.COMM_WORLD.rank==0:
        sendbuf="done"
    recvbuf=MPI.COMM_WORLD.bcast(sendbuf,root)
    print(recvbuf)

MPI programming is difficult. It is "schizophrenic programming" in
that you are writing a single programming with multiple threads of
execution "many voices in one head".
