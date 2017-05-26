/* This file is my own test environment. It is added to git, but will
be removed from git when the branch is commited*/

#include <iostream>
#include "mathematical_program.h"
//#include "drake/common/symbolic_expression.h"

using namespace drake::solvers;
using namespace drake::symbolic;

namespace drake {
namespace solvers {
namespace executable {
namespace {

int psd_example() {
  std::cout << "runs psd_examples" << std::endl;

  // Create a MathematicalProgram object.
  MathematicalProgram prog;
  // Add a 2 x 2 symmetric matrix S to optimization program as new decision
  // variables.
  auto S = prog.NewSymmetricContinuousVariables<2>("S");
  // Impose a positive semidefinite constraint on S.
  auto psd_constraint = prog.AddPositiveSemidefiniteConstraint(S);
  // Add more constraints to make the program more interesting,
  // but this is not needed.
  // Add the constraint that S(1, 0) = 1.
  prog.AddBoundingBoxConstraint(1, 1, S(1, 0));
  // Minimize S(0, 0) + S(1, 1).
  prog.AddLinearCost(S(0, 0) + S(1, 1));
  // Now solve the program.
  prog.Solve();
  // Retrieve the solution of matrix S.
  auto S_value = prog.GetSolution(S);
  // Compute the eigen values of the solution, to see if they are
  // all non-negative.
  //  Eigen::Vector4d S_stacked;
  //  S_stacked << S_value.col(0), S_value.col(1);
  //  Eigen::VectorXd S_eigen_values;
  //  psd_constraint->Eval(S_stacked, S_eigen_values);
  std::cout << "S solution is: " << S_value << std::endl;
  //  std::cout << "The eigen values of S are " << S_eigen_values << std::endl;

  return 0;
}

}  // namespace
}  // namespace executable
}  // namespace solvers
}  // namespace drake

int main() {
  std::cout << "runs executable.cc" << std::endl;


  //  MathematicalProgram prog;
  //  auto z = prog.NewIndeterminates(3, "z");
  //  auto x = prog.NewContinuousVariables(4, "x");
  //  drake::symbolic::Expression expression = x(0) + 2 * x(1);
  //  auto map = DecomposePolynomialIntoExpression(expression);

  drake::solvers::executable::psd_example();

  return 0;
}
