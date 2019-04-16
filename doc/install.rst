.. highlight:: bash

.. _install:

Compiling CTINT from source
===========================


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` library, see :ref:`TRIQS installation instruction <triqslibs:installation>`.
   In the following, we assume that Triqs is installed in the ``path_to_triqs`` directory.

Installation steps 
------------------

#. Download the source code of the latest stable version by cloning the ``TRIQS/ctint`` repository from GitHub::

     $ git clone https://github.com/TRIQS/ctint ctint.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir ctint.build && cd ctint.build
 
#. Ensure that your shell contains the TRIQS environment variables by sourcing the ``triqsvars.sh`` file from your TRIQS installation::

     $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake ../ctint.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Version compatibility
---------------------

Be careful that the version of the TRIQS library and of the solver must be
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
In particular you should make sure that the Major and Minor Version number
of the application and TRIQS agree.
If you want to use a particular version of the solver, go into the directory with the sources
and look at all available versions::

     $ cd ctint.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2.1.0

Then follow the steps 2 to 4 described above to compile the code.

Custom CMake options
--------------------

Functionality of ``ctint`` can be tweaked using extra compile-time options passed to CMake::

    cmake ../ctint.src -DOPTION1=value1 -DOPTION2=value2 ... ../ctint.src

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path other than path_to_triqs           | -DCMAKE_INSTALL_PREFIX=path_to_ctint          |
+-----------------------------------------------------------------+-----------------------------------------------+
| Allow the non-interacting Green function G0(tau) to be complex  | -DGTAU_IS_COMPLEX=ON                          |
+-----------------------------------------------------------------+-----------------------------------------------+
| Allow the interaction term to be complex                        | -DINTERACTION_IS_COMPLEX=ON                   |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
