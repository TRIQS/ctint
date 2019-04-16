.. highlight:: bash

.. _install:

Installation
============


Prerequisite
-------------------

#. The :ref:`TRIQS <triqslibs:welcome>` toolbox (see :ref:`TRIQS installation instruction <triqslibs:installation>`).
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

Installation steps 
------------------

#. Download the sources of the solver from github:: 
 
     $ git clone https://github.com/Wentzell/ctint.git ctint.src
 
#. Create an empty build directory where you will compile the code:: 
 
     $ mkdir ctint.build && cd ctint.build 
 
#. Make sure that you have added the TRIQS and Cpp2Py installation to your environment variables::

     $ source path_to_triqs/share/triqsvarsh.sh

#. In the build directory call cmake::

     $ cmake ../ctint.src
 
#. Compile the code, run the tests and install the application:: 
 
     $ make 
     $ make test 
     $ make install 
 
Version compatibility
---------------------
 
Be careful that the version of the TRIQS library and of the solver must be 
compatible (more information on the :ref:`TRIQS website <triqslibs:versions>`).
As the ctint solver is currently still in alpha phase (unstable), it can
only be compiled against the unstable branch of triqs.

Custom CMake options
--------------------

Functionality of ``ctint`` can be tweaked using extra compile-time options passed to CMake::

    cmake -DOPTION1=value1 -DOPTION2=value2 ... ../ctint.src

+-----------------------------------------------------------------------+-----------------------------------------------+
| Options                                                               | Syntax                                        |
+=======================================================================+===============================================+
| Specify an installation path other than path_to_triqs                 | -DCMAKE_INSTALL_PREFIX=path_to_ctint          |
+-----------------------------------------------------------------------+-----------------------------------------------+
| Allow the non-interacting Green function G0(tau) to be complex        | -DGTAU_IS_COMPLEX=ON                          |
+-----------------------------------------------------------------------+-----------------------------------------------+
| Allow the interaction term to be complex                              | -DINTERACTION_IS_COMPLEX=ON                   |
+-----------------------------------------------------------------------+-----------------------------------------------+
| Build the documentation locally                                       | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------------+-----------------------------------------------+
| Build the ctint solver in debugging mode (*developers only*)          | -DDEBUG_CTINT=ON                              |
+-----------------------------------------------------------------------+-----------------------------------------------+
