#pragma once

#include <Eigen/Dense>

namespace drake {
namespace math {

/**
 * Computes a unique solution to the continuouse-time Lyapunov equation:
 *
 * @verbatim
 * A'X + XA + Q = 0
 * @endverbatim
 *
 * where @param Q is a Hermitian matrix with size equal to @param A.
 *
 * Limitations: Given the Eigenvalues of @param A as lambda_1, ..., lambda_n,
 * then the solution is unique if and only if lambda_i +
 * complexConjugate(lambda_j) != 0 for all i, j.
 * Further, if all lambda_i have negative real parts, and if @param Q is
 * semi-positive definite, then X is also semi-positive definite [1].
 * TODO(FischerGundlach) Check if there is some implication if Q is positive
 * definite, then X is also positive definite. Otherwise, we may have trivial
 * solutions, i.e. X = 0.
 *
 * The implementation is based on SLICOT routine SB03MD [2].
 *
 * [1] Hammarling, S.J., "Numerical solution of the stable, non-negative
 * definite Lyapunov equation," IMA J. Num. Anal., Vol. 2, pp. 303–325, 1982.
 * [2] http://slicot.org/objects/software/shared/doc/SB03MD.html
 *
 * */

Eigen::MatrixXd ContinuousLyapunovEquation(
    const Eigen::Ref<const Eigen::MatrixXd>& A,
    const Eigen::Ref<const Eigen::MatrixXd>& Q);

namespace internal {

// Subroutines which help special cases.

Eigen::Vector1d Solve1By1ContinuousLyapunovEquation(const Eigen::Ref<const Eigen::Vector1d> &A,
                                                    const Eigen::Ref<const Eigen::Vector1d> &Q);

Eigen::Matrix2d Solve2By2ContinuousLyapunovEquation(const Eigen::Ref<const Eigen::Matrix2d> &A,
                                                    const Eigen::Ref<const Eigen::Matrix2d> &Q);

Eigen::MatrixXd SolveReducedContinuousLyapunovEquation(const Eigen::Ref<const Eigen::MatrixXd> &A,
                                                    const Eigen::Ref<const Eigen::MatrixXd> &Q);

}  // namespace internal
}  // namespace math
}  // namespace drake
