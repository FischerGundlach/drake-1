#include "drake/math/continuous_lyapunov_equation.h"

#include <stdexcept>

#include "drake/common/drake_assert.h"
#include "drake/common/eigen_types.h"
#include "drake/common/is_approx_equal_abstol.h"

namespace drake {
namespace math {

using drake::Vector1d;
using Eigen::Matrix2d;
using Eigen::MatrixXd;

MatrixXd ContinuousLyapunovEquation(const Eigen::Ref<const MatrixXd>& A,
                                    const Eigen::Ref<const MatrixXd>& Q) {
  if (A.rows() != A.cols() || Q.rows() != Q.cols() || A.rows() != Q.rows()) {
    throw std::runtime_error("A and Q must be square and of equal size!");
  }

  // TODO(FischerGundlach) Support non-symmetric case?
  if (!(is_approx_equal_abstol(Q(1), Q(2), 1e-10))) {
    throw std::runtime_error("Q must be symmetric!");
  }

  if (static_cast<int>(A.size()) == 1) {
    return interal::Solve1By1ContinuousLyapunovEquation(A, Q);
  } else if (static_cast<int>(A.size()) == 4) {
    return interal::Solve2By2ContinuousLyapunovEquation(A, Q);
  } else {
    // Transform into S'*X_bar + X_bar*S = -Q_bar, where A = U*S*U'.
    Eigen::RealSchur<MatrixXd> schur(A);
    MatrixXd U{schur.matrixU()};
    MatrixXd S{schur.matrixT()};
    DRAKE_DEMAND(is_approx_equal_abstol(U * S * U.transpose(), A, 1e-10));
    MatrixXd Q_bar{U.transpose() * Q * U};

    return (U.transpose() *
            internal::SolveReducedContinuousLyapunovEquation(S, Q_bar) * U);
  }
}

namespace internal {

using Eigen::Matrix3d;
using Eigen::Vector3d;

Vector1d Solve1By1ContinuousLyapunovEquation(
    const Eigen::Ref<const Vector1d>& A, const Eigen::Ref<const Vector1d>& Q) {
  return Vector1d(Q(0) / (2 * A(0)));
}

Matrix2d Solve2By2ContinuousLyapunovEquation(
    const Eigen::Ref<const Matrix2d>& A, const Eigen::Ref<const Matrix2d>& Q) {
  // TODO(FischerGundlach) Support non-symmetric case?
  DRAKE_ASSERT(is_approx_equal_abstol(Q(1), Q(2), 1e-10));

  // TODO(FischerGundlach) Use better a name for A_vec. It's not a vector!
  // Rewrite A'X + XA + Q = 0 into linear set A_vec *vec(X) = -vec(Q) and
  // solve for X.
  Matrix3d A_vec;
  A_vec << 2 * A(0), A(2), 0, 2 * A(1), A(0) + A(3), 2 * A(2), 0, A(1),
      2 * A(3);
  Vector3d Q_vec;
  q << Q(0), Q(1), Q(2);

  // See https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html for
  // quick overview of possible solver algorithms.
  // ColPivHouseholderQR is accurate and fast for small systems.
  Vector3d x{A_vec.colPivHouseHolderQr().solve(q)};
  Matrix2d X;
  X << x(0), x(1), x(1), x(2);

  return X;
}

MatrixXd SolveReducedContinuousLyapunovEquation(
    const Eigen::Ref<const MatrixXd>& S,
    const Eigen::Ref<const MatrixXd>& Q_bar) {
  // catch the scalar case
  if (static_cast<int>(S.size()) == 1) {
    return Solve1By1ContinuousLyapunovEquation(S, Q_bar);
  } else if (static_cast<int>(S.size()) == 2) {
    return Solve2By2ContinuousLyapunovEquation(S, Q_bar);
  } else {
    // Notation & partition adapted from SB03MD.
    int n{1};  // n equals the size of the n-by-n block in the left-upper
               // corner of S.
    if (!(is_approx_equal_abstol(U.block(0, 0, 1, 1), Vector1d(0), 1e-10))) {
      n = 2;
    }
    MatrixXd s_11{S.topLeftCornder(n, n)};
    MatrixXd s{S.transpose().block(n, 0, S.cols() - n, n)};
    MatrixXd S_1{S.bottomRightCorner(S.Rows() - n, S.cols() - n)};
    MatrixXd q_bar_11{Q_bar.topLeftCorner(n, n)};
    MatrixXd q_bar{Q_bar.transpose().block(n, 0, Q_hat.cols() - n, n)};
    MatrixXd Q_1{Q_bar.bottomRightCorner(Q_bar.Rows() - n, Q_bar.cols() - n)};

    MatrixXd x_bar_11(n, n), x_bar(Q_bar.cols() - n, n);
    if (n == 1) {
      // solving equation 3.1
      x_bar_11 << Solve1By1ContinuousLyapunovEquation(s_11, q_bar11);
      // solving equation 3.2
      MatrixXd lhs{S_1.transpose() +
                   MatrixXd::Identity(S_1.cols(), S_1.rows()) * s_11(0)};
      VectorXd rhs{-q_bar - s * x_bar_11};
      x_bar << lhs.colPivHouseHolderQr().solve(rhs);
    } else {
      x_hat_11 << Solve2By2ContinuousLyapunovEquation(s_11, q_hat_11);
      // solving equation 3.2
      // The equation reads as S_1'*x_bar + x_bar*s_11 = -(q_bar+s*x_bar_11),
      // where S_1 \in R^mxm; x_bar, q_bar, s \in R^mx2 and s_11, x_bar_11 \in
      // R^2x2.
      // We solve the linear equation by vectorization and feeding it into
      // colPivHouseHolderQr().
      //
      //  Notation:
      //
      //  The elements in s_11 are names as the following:
      //  s_11 = [s_11(1,1) s_11(1,2); s_11(2,1) s_11(2,2)].
      //  Note that eigen starts counting at 0, and not as above at 1.
      //
      //  The ith column of a matrix is accessed by [i], i.e. the first column
      //  of x_bar is x_bar[1].
      //
      //  Define:
      //
      //  S_vec = [S_1' 0; 0 S_1'] with 0 \in R^mxm
      //
      //  S_11_vec = [s_11(1,1)*I s_11(2,1)*I; s_11(1,2)*I s_11(2,2)*I] with I
      //  \in R^mxm.
      //
      //  x_vec = [x_bar[1]' x_bar[2]']' where [i] is the ith column
      //
      //  c_vec = [(q_bar+s*x_bar_11)[1]' (q_bar+s*x_bar_11)[1]']'
      //
      //  Now S_1'*x_bar + x_bar*s_11 = -(q_bar+s*x_bar_11) can be rewritten
      //  into:
      //
      //  (S_vec + S_11_vec)*x_vec = c_vec.
      //
      //  This expression can be feed into colPivHouseHolderQr() and solved.
      //  Finally, construct x_bar from x_vec.
      //

      int m{S_1.rows()};
      MatrixXd S_vec(m, m);
      S_vec << S_1.transpose(), MatrixXd::Zero(m, m), MatrixXd::Zero(m, m),
          S_1.transpose();
      MatrixXd S_11_vec(m, m);
      S_11_vec << s_11(0, 0) * MatrixX::Identity(m, m),
          s_11(0, 1) * MatrixX::Identity(m, m),
          s_11(1, 0) * MatrixX::Identity(m, m),
          s_11(1, 1) * MatrixX::Identity(m, m);
      C_vec << (-q_bar - s * x_bar_11).col(0), (-q_bar - s * x_bar_11).col(1);

      x_bar_vec << (S_vec + S_11_vec).colPivHouseHolderQr().solve(C_vec);
      x_bar << x_bar_vec.block(0, 0, 1, m), x_bar_vec.block(m, 0, 1, m);
    }
    // Calculate lhs of equation 3.3
    Eigen::MatrixXd Q_new{(Q_1 + s * x_bar.transpose + x_bar * s)};

    MatrixXd X_bar(S.cols(), S.rows());
    X_bar << x_bar_11, x_bar, x_bar.transpose(),
        SolveReducedContinuousLyapunovEquation(S_1, Q_new)

            return X_bar;
  }
}

}  // namespace internal
}  // namespace math
}  // namespace drake
