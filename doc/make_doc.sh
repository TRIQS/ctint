TRIQS_PATH=$HOME/Work/TRIQS_1.4/install/bin
rm -rf triqs
$TRIQS_PATH/c++2doc.py -N triqs -N applications -N impurity_solvers -N ctint common.hpp
