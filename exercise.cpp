#include <iostream>
#include "algebra.hpp"

using namespace std;
using namespace Algebra;

/* 
 *this is just a working exercise to for testing the api.
 *full tests are writen (see ./tests/)
 */
int main(int argc, char **args) {
  //two frames {a,b} with respect to frame s.
  Eigen::Matrix3d Ra;
  Ra << 0, -1, 0,
    0, 0, -1,
    1, 0, 0;
  Eigen::Vector3d pa(0,0,1);
  Eigen::Matrix3d Rb;
  Rb << 1,0,0,
    0,0,1,
    0,-1,0;
  //
  Eigen::Vector3d pb(0,2,0);
  //Q1
  Eigen::Matrix4d Ta=Algebra::RpToTrans(Ra,pa);
  std::cout << "Q1: " << Ta << std::endl;
  //Q2
  Eigen::Matrix4d Tb=Algebra::RpToTrans(Rb,pb);
  std::cout << "Q2: " << Algebra::TransInv(Tb) << std::endl;
  //Q3
  std::cout << "Q3: " << Algebra::TransInv(Ta)*Tb << std::endl;
  //Q5
  auto T=Tb*Tb;
  Eigen::Vector3d pt_b(1,2,3);
  Eigen::Vector4d point_b4(pt_b(0), pt_b(1), pt_b(2), 1);
  auto Q5 = Tb*point_b4;
  std::cout << "Q5: " << Q5 << std::endl;
  //Q7
  Eigen::VectorXd Vs(6);
  Vs << 3,2,1,-1,-2,-3;
  auto Va=Algebra::Adjoint(Algebra::TransInv(Ta))*Vs;
  std::cout << "Q7: " << Va << std::endl;
  //Q8
  auto Q8=Algebra::AxisAng6(Algebra::se3ToVec(Algebra::MatrixLog6(Ta)))[6];
  std::cout << "Q8: " << Q8 << std::endl;
  //Q9
  Eigen::VectorXd S_t(6);
  S_t << 0,1,2,3,0,0;
  std::cout << "Q9: " << Algebra::MatrixExp6(Algebra::VecTose3(S_t)) << std::endl;
  //Q10
  Eigen::VectorXd Fb(6);
  Fb << 1,0,0,2,1,0;
  auto Q10 = Algebra::Adjoint(Tb)*Fb;
  std::cout << "Q10: " << Q10 << std::endl;
  //Q11
  Eigen::Matrix4d trans;
  trans << 0, -1, 0, 3,
    1,0,0,0,
    0,0,1,1,
    0,0,0,1;
  auto Q11=Algebra::TransInv(trans);
  std::cout << "Q11: " << Q11 << std::endl;
  //Q12
  Eigen::VectorXd V(6);
  V << 1,0,0,0,2,3;
  auto Q12=Algebra::VecTose3(V);
  std::cout << "Q12: " << Q12 << std::endl;
  //Q13
  float h=1;
  Eigen::Vector3d s(1,0,0);
  Eigen::Vector3d p(0,0,2);
  auto Q13=Algebra::ScrewToAxis(p,s,h);
  std::cout << "Q13: " << Q13 << std::endl;
  //Q14
  Eigen::Matrix4d so3;
  so3 << 0,-1.5708,0,2.3562,
    1.5708,0,0,-2.3562,
    0,0,0,1,
    0,0,0,0;
  auto Q14=Algebra::MatrixExp6(so3);
  std::cout << "Q14: " << Q14 << std::endl;
  //Q15
  Eigen::Matrix4d Tt;
  Tt << 0,-1,0,3,
    1,0,0,0,
    0,0,1,1,
    0,0,0,1;
  auto Q15=Algebra::MatrixLog6(Tt);
  std::cout << "Q15: " << Q15 << std::endl;

  //
  // w3
  //
  Eigen::Matrix3d ra;
  ra << 0,1,0,
    0,0,1,
    1,0,0;
  std::cout << "Q8: " << Algebra::AxisAng3(Algebra::so3ToVec(Algebra::MatrixLog3(ra)))(3) << std::endl;
  Eigen::Vector3d w_t(1,2,0);
  std::cout << "Q9: " << Algebra::MatrixExp3(Algebra::VecToso3(w_t)) << std::endl;
  Eigen::Matrix3d w_t_so3;
  w_t_so3 << 0,0.5,-1,-0.5,0,2,1,-2,0;
  std::cout << "Q10: " << Algebra::MatrixExp3(w_t_so3) << std::endl;
  Eigen::Matrix3d r;
  r << 0,0,1,-1,0,0,0,-1,0;
  std::cout << "Q11: " << Algebra::MatrixLog3(r) << std::endl;
  return 0;
}
