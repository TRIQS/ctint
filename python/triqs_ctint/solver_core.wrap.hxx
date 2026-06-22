#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_solver_core_GUARDS
#define C2PY_HXX_DECLARATION_solver_core_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_ctint::constr_params_t>     = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctint::constr_params_t> = "triqs_ctint.solver_core.ConstrParamsT";
template <> constexpr bool c2py::is_wrapped<triqs_ctint::solve_params_t>      = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctint::solve_params_t>  = "triqs_ctint.solver_core.SolveParamsT";
template <> constexpr bool c2py::is_wrapped<triqs_ctint::solver_core>         = true;
template <> inline constexpr auto c2py::tp_name<triqs_ctint::solver_core>     = "triqs_ctint.solver_core.SolverCore";
#endif