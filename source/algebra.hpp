#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <iostream>

# define M_PI           3.14159265358979323846  /* pi */

/** \namespace Algebra
 * \brief Is Linear Algebra library of required mathematics for,
 * Kinetics, and Kinematics mechanics.
 *
 * Build on top of Kevin M. Lynch's, and Frank C. Park 'Modern Robotics' Book http://hades.mech.northwestern.edu/index.php/Modern_Robotics, and courses from Princeton.
 */
namespace Algebra {

  /**
   * Find if the value is negligible, or close to zero
   * @param val double value to be checked.
   * @return Boolean of true-ignore or false-can't ignore.
   */
  inline bool Epsilon(const double val) {
    return (std::abs(val) < 1e-6);
  }

  /**
   * Get the skew symmetric matrix representation
   of a vector \f$v=\begin{vmatrix}v_1&v_2&v_3\end{vmatrix}\f$
   \f$skew(v)=\begin{vmatrix}0&-v_3&v_2\\v_3&0&-v_1\\-v_2&v_1&0\end{vmatrix}\f$
   * @param omg Eigen::Vector3d 3x1 angular velocity vector
   * @see so3ToVec
   * @return Eigen::MatrixXd 3x3 skew symmetric matrix in so3
   */
  inline Eigen::Matrix3d VecToso3(const Eigen::Vector3d& omg) {
    Eigen::Matrix3d m_ret;
    m_ret << 0, -omg(2), omg(1),
      omg(2), 0, -omg(0),
      -omg(1), omg(0), 0;
    return m_ret;
  }

  /**
   * Returns  vector represented by the skew symmetric matrix.
   \f$v=\begin{vmatrix}v_1&v_2&v_3\end{vmatrix}\f$
   \f$skew(v)=\begin{vmatrix}0&-v_3&v_2\\v_3&0&1v_1\\-v_2&v_1&0\end{vmatrix}\f$
   * @param so3mat Eigen::MatrixXd 3x3 skew symmetric matrix in so3
   * @see VecToso3
   * @return v_ret Eigen::Vector3d 3x1 vector
   */
  inline Eigen::Vector3d so3ToVec(const Eigen::MatrixXd& so3mat) {
    Eigen::Vector3d v_ret;
    v_ret << so3mat(2, 1), so3mat(0, 2), so3mat(1, 0);
    return v_ret;
  }

  /**
   * Combines a rotation matrix and position vector
   *           into a single Special Euclidian Group (SE3)
   *           homogeneous transformation matrix
   *
   * @param R Rotation Matrix
   * @param p Position Vector
   * @see TransToRp
   * @return Transformation Matrix \f$\begin{vmatrix}R&p\\0&1\end{vmatrix}\f$.
   */
  inline Eigen::MatrixXd
  RpToTrans(const Eigen::Matrix3d& R, const Eigen::Vector3d& p) {
    Eigen::MatrixXd m_ret(4, 4);
    m_ret << R, p,
      0, 0, 0, 1;
    return m_ret;
  }

  /**
   * Separates the rotation matrix and position vector
   *           from the transfomation matrix representation
   * @param T Homogeneous transformation matrix \f$\begin{vmatrix}R&p\\0&1\end{vmatrix}\f$.
   * @see RpToTrans
   * @return Rp_ret std::vector of rotation matrix (R), and position vector (p).
   */
  inline std::vector<Eigen::MatrixXd>
  TransToRp(const Eigen::MatrixXd& T) {
    std::vector<Eigen::MatrixXd> Rp_ret;
    Eigen::Matrix3d R_ret;
    // Get top left 3x3 corner
    R_ret = T.block<3, 3>(0, 0);

    Eigen::Vector3d p_ret(T(0, 3), T(1, 3), T(2, 3));

    Rp_ret.push_back(R_ret);
    Rp_ret.push_back(p_ret);

    return Rp_ret;
  }

  /**
   * Translates a spatial velocity vector into a transformation matrix
   * @param V Spatial velocity vector \f$\begin{vmatrix}w\\v\end{vmatrix}\f$
   * @see se3ToVec
   * @return linear representation of Twist in special Euclidian group se3 \f$\begin{vmatrix}[w]&v\\0&0\end{vmatrix}\f$
   */
  inline Eigen::MatrixXd VecTose3(const Eigen::VectorXd& V) {
    Eigen::Vector3d exp(V(0), V(1), V(2));
    Eigen::Vector3d linear(V(3), V(4), V(5));
    Eigen::MatrixXd m_ret(4, 4);
    m_ret << VecToso3(exp), linear,
      0, 0, 0, 0;
    return m_ret;
  }

  /**
   * Translates a transformation matrix into a spatial
   *           velocity vector
   * @param T Transformation matrix \f$\begin{vmatrix}[w]&v\\0&0\end{vmatrix}\f$
   * @see VecTose3
   * @return Spatial velocity vector \f$\begin{vmatrix}w\\v\end{vmatrix}\f$
   */
  inline Eigen::VectorXd se3ToVec(const Eigen::MatrixXd& T) {
    Eigen::VectorXd m_ret(6);
    m_ret << T(2, 1), T(0, 2), T(1, 0), T(0, 3), T(1, 3), T(2, 3);
    return m_ret;
  }

  /**
   * Returns a normalized version of the input vector
   * @param V Eigen::MatrixXd
   * @return Eigen::MatrixXd normalized V
   * Note: MatrixXd is used instead of VectorXd for the case of
   *      row vectors
   * 	  Requires a copy Useful because of the MatrixXd casting
   */
  inline Eigen::MatrixXd Normalize(Eigen::MatrixXd V) {
    V.operatorNorm();
    return V;
  }

  /*
   * Translates an exponential rotation into it's individual components
   * @param expc3 Exponential rotation (rotation matrix in terms of a rotation axis and the angle of rotation) \f$[w]\theta\f$
   * @see AxisAng6
   * @return The axis and angle of rotation as [x, y, z, theta]
   */
  inline Eigen::Vector4d AxisAng3(const Eigen::Vector3d& expc3) {
    Eigen::Vector4d v_ret;
    float theta = expc3.norm();
    v_ret << expc3/theta, theta;
    //v_ret << Normalize(expc3)/theta, theta;
    return v_ret;
  }

  /**
   * Translates an exponential rotation into a rotation
   *           matrix using rodrigues formula.
   * @param so3mat exponenential representation of a rotation in so3 \f$[w]\theta\f$.
   * @return Rotation matrix \f$e^{[w]\theta}\f$.
   */
  inline Eigen::Matrix3d MatrixExp3(const Eigen::Matrix3d& so3mat) {
    Eigen::Vector3d omgtheta = so3ToVec(so3mat);

    Eigen::Matrix3d m_ret = Eigen::Matrix3d::Identity();
    if (Epsilon(so3mat.norm())) {
      return m_ret;
    }
    else {
      double theta = (AxisAng3(omgtheta))(3);
      Eigen::Matrix3d omgmat = so3mat * (1 / theta);
      //Rodrigues Formula.
      //how to write power using Eigen?
      return m_ret + std::sin(theta) * omgmat +
        ((1 - std::cos(theta)) * (omgmat * omgmat));
    }
  }

  /**
   * Function: Computes the matrix logarithm of a rotation matrix
   * @param R Rotation matrix \f$[w]\theta\f$.
   * @return matrix logarithm of a rotation \f$e^{[w]\theta}\f$.
   */
  inline Eigen::Matrix3d MatrixLog3(const Eigen::Matrix3d& R) {
    double acosinput = (R.trace() - 1) / 2.0;
    double theta = std::acos(acosinput);
    Eigen::MatrixXd m_ret = Eigen::MatrixXd::Zero(3, 3);
    if (acosinput >= 1)
      return m_ret;
    else if (acosinput <= -1) {
      Eigen::Vector3d omg;
      if (!Epsilon(1 + R(2, 2)))
        omg = (1.0 / std::sqrt(2 * (1 + R(2, 2))))*
            Eigen::Vector3d(R(0, 2), R(1, 2), 1 + R(2, 2));
      else if (!Epsilon(1 + R(1, 1)))
        omg = (1.0 / std::sqrt(2 * (1 + R(1, 1))))*
            Eigen::Vector3d(R(0, 1), 1 + R(1, 1), R(2, 1));
      else
        omg = (1.0 / std::sqrt(2 * (1 + R(0, 0))))*
            Eigen::Vector3d(1 + R(0, 0), R(1, 0), R(2, 0));
      m_ret = VecToso3(M_PI * omg);
    }
    else {
        m_ret = theta / 2.0 / sin(theta) * (R.array() - R.transpose().array());;
    }
    return  m_ret;
  }

  /**
   * Provides the adjoint representation of a transformation matrix Used to change the frame
   *           of reference for spatial velocity vectors
   * @param T 4x4 Transformation matrix SE(3) \f$\begin{vmatrix}R&p\\0&1\end{vmatrix}\f$.
   * @return 6x6 Adjoint Representation of the matrix T \f$adj(T)=\begin{vmatrix}R&0 \\ [p]R&R\end{vmatrix}\f$
   */
  inline Eigen::MatrixXd Adjoint(const Eigen::MatrixXd& T) {
    std::vector<Eigen::MatrixXd> R = TransToRp(T);
    Eigen::MatrixXd ad_ret(6, 6);
    ad_ret = Eigen::MatrixXd::Zero(6, 6);
    Eigen::MatrixXd zeroes = Eigen::MatrixXd::Zero(3, 3);
    ad_ret << R[0], zeroes,
      VecToso3(R[1]) * R[0], R[0];
    return ad_ret;
  }

  /**
   * Rotation expanded for screw axis
   * @param se3mat se3 matrix representation of exponential coordinates \f$[s]\theta\f$
   * @return 6x6 Matrix exponential Transformation \f$T=e^{[s]\theta}\f$
   */
  inline Eigen::MatrixXd
  MatrixExp6(const Eigen::MatrixXd& se3mat) {
    // Extract the angular velocity vector from the transformation matrix
    Eigen::Matrix3d se3mat_cut = se3mat.block<3, 3>(0, 0);
    Eigen::Vector3d omgtheta = so3ToVec(se3mat_cut);

    Eigen::MatrixXd m_ret(4, 4);

    // If negligible rotation, m_Ret = [[Identity, angular velocty ]]
    //                                 [[0       , 1]]
    if (Epsilon(omgtheta.norm())) {
      // Reuse previous variables that have our required size
      se3mat_cut = Eigen::MatrixXd::Identity(3, 3);
      omgtheta << se3mat(0, 3), se3mat(1, 3), se3mat(2, 3);
      m_ret << se3mat_cut, omgtheta,
        0, 0, 0, 1;
      return m_ret;
    }
    // If not negligible, MR page 105
    else {
      double theta = (AxisAng3(omgtheta))(3);
      Eigen::Matrix3d omgmat = se3mat.block<3, 3>(0, 0) / theta;
      Eigen::Matrix3d expExpand = Eigen::MatrixXd::Identity(3, 3) * theta + (1 - std::cos(theta)) * omgmat + ((theta - std::sin(theta)) * (omgmat * omgmat));
      Eigen::Vector3d linear(se3mat(0, 3), se3mat(1, 3), se3mat(2, 3));
      Eigen::Vector3d GThetaV = (expExpand*linear) / theta;
      m_ret << MatrixExp3(se3mat_cut), GThetaV,
        0, 0, 0, 1;
      return m_ret;
    }

  }
  /**
   * Computes the matrix logarithm of a homogeneous transformation matrix
   *
   * @param T: A matrix in SE3.
   * @return The matrix logarithm of \f$e^{[s]\theta}\f$.
   */
    inline Eigen::MatrixXd MatrixLog6(const Eigen::MatrixXd& T) {
    Eigen::MatrixXd m_ret(4, 4);
    auto rp = TransToRp(T);
    Eigen::Matrix3d omgmat = MatrixLog3(rp.at(0));
    Eigen::Matrix3d zeros3d = Eigen::Matrix3d::Zero(3, 3);
    if (Epsilon(omgmat.norm())) {
      m_ret << zeros3d, rp.at(1),
        0, 0, 0, 0;
    }
    else {
      double theta = std::acos((rp.at(0).trace() - 1) / 2.0);
      Eigen::Matrix3d logExpand1 = Eigen::MatrixXd::Identity(3, 3) - omgmat / 2.0;
      Eigen::Matrix3d logExpand2 = (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2)*omgmat*omgmat / theta;
      Eigen::Matrix3d logExpand = logExpand1 + logExpand2;
      m_ret << omgmat, logExpand*rp.at(1),
        0, 0, 0, 0;
		}
    return m_ret;
  }

  /**
   * Inverts a homogeneous transformation matrix
    @param T: A homogeneous transformation matrix
    @return The inverse of T Uses the structure of transformation matrices to avoid taking a matrix inverse, for efficiency.
  */
  inline Eigen::MatrixXd TransInv(const Eigen::MatrixXd& transform) {

    auto rp = TransToRp(transform);
    auto Rt = rp.at(0).transpose();
    auto t = -(Rt * rp.at(1));
    Eigen::MatrixXd inv(4, 4);
    inv = Eigen::MatrixXd::Zero(4,4);
    inv.block(0, 0, 3, 3) = Rt;
    inv.block(0, 3, 3, 1) = t;
    inv(3, 3) = 1;
    return inv;
  }
  /**
   *Inverts a rotation matrix

   @param R: A rotation matrix
   @return: The inverse of R
  */
  inline Eigen::MatrixXd RotInv(const Eigen::MatrixXd& rotMatrix) {
    return rotMatrix.transpose();
  }

  /**
     Takes a parametric description of a screw axis and converts it to a
     normalized screw axis

     @param q A point lying on the screw axis
     @param s A unit vector in the direction of the screw axis
     @param h The pitch of the screw axis
     @return A normalized screw axis described by the inputs
  */
  inline Eigen::VectorXd ScrewToAxis(Eigen::Vector3d q, Eigen::Vector3d s, double h) {
    Eigen::VectorXd axis(6);
    axis.segment(0, 3) = s;
    axis.segment(3, 3) = q.cross(s) + (h * s);
    return axis;
  }

  /**
     Converts a 6-vector of exponential coordinates into screw axis-angle form
     @param expc6 A 6-vector of exponential coordinates for rigid-body motion
     S*theta
     @return S The corresponding normalized screw axis
     @return theta The distance traveled along/about S
  */
  inline Eigen::VectorXd AxisAng6(const Eigen::VectorXd& expc6) {
    Eigen::VectorXd v_ret(7);
    double theta = Eigen::Vector3d(expc6(0), expc6(1), expc6(2)).norm();
    if (Epsilon(theta))
      theta = Eigen::Vector3d(expc6(3), expc6(4), expc6(5)).norm();
    v_ret << expc6 / theta, theta;
    return v_ret;
  }

  /**
     Returns a projection of mat into SO(3)

     @param M A matrix near SO(3) to project to SO(3)
     @return The closest matrix to R that is in SO(3)
     Projects a matrix mat to the closest matrix in SO(3) using singular-value
     decomposition (see http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
     This function is only appropriate for matrices close to SO(3).
   */
  inline Eigen::MatrixXd ProjectToSO3(const Eigen::MatrixXd& M) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd R = svd.matrixU() * svd.matrixV().transpose();
    if (R.determinant() < 0)
      // In this case the result may be far from M; reverse sign of 3rd column
      R.col(2) *= -1;
    return R;
  }

  /**
     Returns a projection of mat into SE(3)

     @param M A 4x4 matrix to project to SE(3)
     @return The closest matrix to T that is in SE(3)
     Projects a matrix mat to the closest matrix in SE(3) using singular-value
     decomposition (see
     http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
     This function is only appropriate for matrices close to SE(3).
   */
  inline Eigen::MatrixXd ProjectToSE3(const Eigen::MatrixXd& M) {
    Eigen::Matrix3d R = M.block<3, 3>(0, 0);
    Eigen::Vector3d t = M.block<3, 1>(0, 3);
    Eigen::MatrixXd T = RpToTrans(ProjectToSO3(R), t);
    return T;
  }
  /**
     Returns the Frobenius norm to describe the distance of mat from the SO(3) manifold

     @param M A 3x3 matrix
     @return A quantity describing the distance of mat from the SO(3) manifold
     Computes the distance from mat to the SO(3) manifold using the following
     method:
     If det(mat) <= 0, return a large number.
     If det(mat) > 0, return norm(mat^T.mat - I).
  */
  inline double DistanceToSO3(const Eigen::Matrix3d& M) {
    if (M.determinant() > 0)
      return (M.transpose() * M - Eigen::Matrix3d::Identity()).norm();
    else
      return 1.0e9;
  }
  /**
     Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold

     @param mat A 4x4 matrix
     @return A quantity describing the distance of mat from the SE(3) manifold
     Computes the distance from mat to the SE(3) manifold using the following
     method:
     Compute the determinant of matR, the top 3x3 submatrix of mat.
     If det(matR) <= 0, return a large number.
     If det(matR) > 0, replace the top 3x3 submatrix of mat with matR^T.matR,
     and set the first three entries of the fourth column of mat to zero. Then return norm(mat - I).
  */
  inline double DistanceToSE3(const Eigen::Matrix4d& T) {
    Eigen::Matrix3d matR = T.block<3, 3>(0, 0);
    if (matR.determinant() > 0) {
      Eigen::Matrix4d m_ret;
      m_ret << matR.transpose()*matR, Eigen::Vector3d::Zero(3),
          T.row(3);
      m_ret = m_ret - Eigen::Matrix4d::Identity();
      return m_ret.norm();
    }
    else
        return 1.0e9;
  }

    inline bool TestIfSO3(const Eigen::Matrix3d& M) {
        return std::abs(DistanceToSO3(M)) < 1e-3;
    }

    inline bool TestIfSE3(const Eigen::Matrix4d& T) {
        return std::abs(DistanceToSE3(T)) < 1e-3;
    }

    /**
     * Calculate the 6x6 matrix of the given
     *           6-vector V
     * @param V Eigen::VectorXd (6x1), V=[omg, v]
     * @return Eigen::MatrixXd (6x6)
     * Note: Can be used to calculate the Lie bracket [V1, V2] =
     *       [adV1]V2
     */
    inline Eigen::MatrixXd ad(Eigen::VectorXd V) {
        Eigen::Matrix3d omgmat =
            VecToso3(Eigen::Vector3d(V(0), V(1), V(2)));
        Eigen::MatrixXd result(6, 6);
        result.topLeftCorner<3, 3>() = omgmat;
        result.topRightCorner<3, 3>() = Eigen::Matrix3d::Zero(3, 3);
        result.bottomLeftCorner<3, 3>() =
            VecToso3(Eigen::Vector3d(V(3), V(4), V(5)));
        result.bottomRightCorner<3, 3>() = omgmat;
        return result;
    }
}
