
.. highlight:: bash

Installation
============


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

Installation steps 
------------------

#. Download the sources of the solver from github:: 
 
     $ git clone https://bitbucket.org/Wentzell/ctint.git ctint.src
 
#. Create an empty build directory where you will compile the code:: 
 
     $ mkdir ctint.build && cd ctint.build 
 
#. In the build directory call cmake (take care that triqsvars.sh was sourced)
 
     $ cmake ../src
 
#. Compile the code, run the tests and install the application:: 
 
     $ make 
     $ make test 
     $ make install 
 
Version compatibility 
--------------------- 
 
Be careful that the version of the TRIQS library and of the solver must be 
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
If you want to use a version of 
the solver that is not the latest one, go into the directory with the sources 
and look at all available versions:: 
 
     $ cd src && git tag 
 
Checkout the version of the code that you want:: 
 
     $ git co 1.0.0 
 
Then follow the steps 2 to 4 described above to compile the code. 
