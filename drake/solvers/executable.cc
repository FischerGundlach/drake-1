/* This file is my own test environment. It is added to git, but will
be removed from git when the branch is commited*/

#include <iostream>
#include "mathematical_program.h"
#include "mosek_solver.h"
//#include "drake/common/symbolic_expression.h"

using namespace drake::solvers;
using namespace drake::symbolic;

namespace drake {
namespace solvers {
namespace executable {
namespace {

/*
int sos_decomposition1() {
  std::cout << "sos_decomposition1()" << std::endl;

  MathematicalProgram prog;
  auto Q = prog.NewSymmetricContinuousVariables<3>("Q");
  auto psd_constraint = prog.AddPositiveSemidefiniteConstraint(Q);

  prog.AddLinearConstraint(Q(0,0), 0, 0);
  prog.AddLinearConstraint(Q(0,1), 0, 0);
  prog.AddLinearConstraint(Q(1,1)+2*Q(1,2), 1, 1);
  prog.AddLinearConstraint(Q(0,2), 0, 0);
  prog.AddLinearConstraint(Q(2,2), 1, 1);

  prog.Solve();
  std::cout << "The Gram Matrix is: " << prog.GetSolution(Q) << std::endl;

  return 0;
}
*/

int sos_decomposition2() {
  std::cout << "sos_decomposition1()" << std::endl;

  MathematicalProgram prog;
  auto x = prog.NewIndeterminates(1,"x");
  Expression x0{x(0)};  // TODO(FischerGundlach) add pow(const var&, const int&)
  Expression f = pow(x0,4) + pow(x0,2);

  auto Q = prog.NewSymmetricContinuousVariables<3>("Q");
  auto psd_constraint = prog.AddPositiveSemidefiniteConstraint(Q);
  Expression sos_poly = Q(0,0) + 2*Q(0,1)*x0 + (Q(1,1)+2*Q(1,2))*pow(x0,2) + Q(2,0)*pow(x0,3) + Q(2,2)*pow(x0,4);

  Variables x_vars{x(0)};
  MonomialAsExpressionToCoefficientMap DecomposePolynomialIntoExpression()


      (f-sos_poly,{x(0)});

/*  prog.AddLinearConstraint(Q(0,0), 0, 0);
  prog.AddLinearConstraint(Q(1,0), 0, 0);
  prog.AddLinearConstraint(Q(1,1)+2*Q(1,2), 1, 1);
  prog.AddLinearConstraint(Q(2,0), 0, 0);
  prog.AddLinearConstraint(Q(2,2), 1, 1);*/

  prog.Solve();
  std::cout << "The Gram Matrix is: " << prog.GetSolution(Q) << std::endl;

//  auto x = prog.NewIndeterminates(3,"x");
//  std::cout << x.transpose()*(Q_value*x) << std::endl;
//
//  auto mono = Monomial(Q(0,0));
//  auto monoSquare = mono*mono;
//  auto monoExpression = mono*Q(0,0);

  return 0;
}

}  // namespace
}  // namespace executable
}  // namespace solvers
}  // namespace drake

int main() {
  std::cout << "runs executable.cc" << std::endl;

  drake::solvers::executable::sos_decomposition2();

  return 0;
}
