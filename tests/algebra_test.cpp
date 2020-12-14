#include <iostream>
#include <Eigen/Dense>
#include <iostream>
#include "gmock/gmock.h"
#include "../source/algebra.hpp"

using namespace testing;

TEST(ALGEBRA, Epsilon) {
  ASSERT_THAT(Algebra::Epsilon(0.0003), false);
}

TEST(ALGEBRA, So3ToVecTest) {
  Eigen::Matrix3d so3(3, 3);
  so3 <<
    0, -3, 2,
    3, 0, -1,
    -2, 1, 0;
  Eigen::Vector3d vec(1, 2, 3);
  EXPECT_EQ(vec, Algebra::so3ToVec(so3));
}

TEST(ALGEBRA, VecToSO3Test) {
  Eigen::Vector3d vec(1, 2, 3);
  Eigen::Matrix3d so3(3, 3);
  so3 <<
    0, -3, 2,
    3, 0, -1,
    -2, 1, 0;
  EXPECT_EQ(so3, Algebra::VecToso3(vec));
}

TEST(ALGEBRA, RpToTrans) {
  Eigen::Matrix3d R;
  R << 1, 2, 3,
    4, 5, 6,
    7, 8, 9;
  Eigen::Vector3d p;
  p << 1, 2, 3;
  Eigen::Matrix4d Rp;
  Rp << 1, 2, 3, 1,
    4, 5, 6, 2,
    7, 8, 9, 3,
    0, 0, 0, 1;
  EXPECT_EQ(Rp, Algebra::RpToTrans(R, p));
}

TEST(ALGEBRA, TransToRp) {
  Eigen::Matrix3d R;
  R << 1, 2, 3,
    4, 5, 6,
    7, 8, 9;
  Eigen::Vector3d p;
  p << 1, 2, 3;
  Eigen::Matrix4d Rp;
  Rp << 1, 2, 3, 1,
    4, 5, 6, 2,
    7, 8, 9, 3,
    0, 0, 0, 1;
  EXPECT_EQ(R, Algebra::TransToRp(Rp)[0]);
  EXPECT_EQ(p, Algebra::TransToRp(Rp)[1]);
}

TEST(ALGEBRA, VecToso3) {
  //spatial twist vector(angular, and linear) velocity
  Eigen::VectorXd V(6);
  V << 1, 2, 3, 4, 5, 6;
  //spatial twist matrix
  Eigen::Matrix4d Tm;
  Tm << 0, -3, 2, 4,
    3, 0, -1, 5,
    -2, 1, 0, 6,
    0, 0, 0, 0;
  EXPECT_EQ(Tm, Algebra::VecTose3(V));
}

TEST(ALGEBRA, so3ToVec) {
  //spatial twist vector(angular, and linear) velocity
  Eigen::VectorXd V(6);
  V << 1, 2, 3, 4, 5, 6;
  //spatial twist matrix
  Eigen::Matrix4d Tm;
  Tm << 0, -3, 2, 4,
    3, 0, -1, 5,
    -2, 1, 0, 6,
    0, 0, 0, 0;
  EXPECT_EQ(V, Algebra::se3ToVec(Tm));
}

TEST(ALGEBRA, AxisAng3) {
  float theta=0.53607;
  //note this is valid iff omg is unite length.
  Eigen::Vector3d omg(0, 0.6, 0.8);
  omg.transpose()*theta;
  Eigen::Vector4d alg_v4= Algebra::AxisAng3(omg*theta);
  //
  EXPECT_EQ(alg_v4(3), theta);
  Eigen::Vector3d alg_omg(alg_v4(0), alg_v4(1), alg_v4(2));
  ASSERT_TRUE(alg_omg.isApprox(omg, 4));
}

TEST(ALGEBRA, Adjoint) {
  Eigen::Matrix3d R;
  R << 1, 2, 3,
    4, 5, 6,
    7, 8, 9;
  Eigen::Vector3d p;
  p << 1, 2, 3;
  Eigen::Matrix4d Rp = Algebra::RpToTrans(R, p);
  Eigen::MatrixXd adj(6,6);
  Eigen::Matrix3d zeros=Eigen::MatrixXd::Zero(3,3);
  adj << R, zeros,
    Algebra::VecToso3(p)*R, R;
  EXPECT_EQ(adj, Algebra::Adjoint(Rp));
}

TEST(ALGEBRA, Adjoint2) {
  Eigen::Matrix4d T;
  T << 1, 0, 0, 0,
    0, 0, -1, 0,
    0, 1, 0, 3,
    0, 0, 0, 1;
  Eigen::MatrixXd result(6, 6);
  result <<
    1, 0, 0, 0, 0, 0,
    0, 0, -1, 0, 0, 0,
    0, 1, 0, 0, 0, 0,
    0, 0, 3, 1, 0, 0,
    3, 0, 0, 0, 0, -1,
    0, 0, 0, 0, 1, 0;

  ASSERT_TRUE(Algebra::Adjoint(T).isApprox(result, 4));
}

// the implementation is the same of python package
// but i according to the book there is missing theta need to be mutiplyed by the output.
// secondly theta calculated isn't actually theta passed to R!!!


TEST(ALGEBRA, MatrixLog3) {
  //Eigen::Vector3d omg(0, 0.8, 0.6);
  Eigen::Vector3d omg(1,2,3);
  //note this tests words only for theta < M_PI there are multiples of solutios see section 3.2 of the book
  float theta=0.5367;
  auto so3=Algebra::VecToso3(omg*theta);
  Eigen::Matrix3d R = Algebra::MatrixExp3(so3);
  Eigen::Matrix3d logR = Algebra::MatrixLog3(R);
  ASSERT_TRUE(so3.isApprox(logR, 4));
}

TEST(ALGEBRA, TransInv) {
  Eigen::MatrixXd input(4, 4);
  input << 1, 0, 0, 0,
    0, 0, -1, 0,
    0, 1, 0, 3,
    0, 0, 0, 1;
  Eigen::MatrixXd result(4, 4);
  result << 1, 0, 0, 0,
    0, 0, 1, -3,
    0, -1, 0, 0,
    0, 0, 0, 1;
  
  auto inv = Algebra::TransInv(input);
  ASSERT_TRUE(inv.isApprox(result, 4));
}


TEST(ALGEBRA, RotInv) {
  Eigen::MatrixXd input(3, 3);
  input << 0, 0, 1,
    1, 0, 0,
    0, 1, 0;
  Eigen::MatrixXd result(3, 3);
  result << 0, 1, 0,
    0, 0, 1,
    1, 0, 0;
  auto inv = Algebra::RotInv(input);
  ASSERT_TRUE(inv.isApprox(result, 4));
}

TEST(ALGEBRA, ScrewToAxis) {
  Eigen::Vector3d q, s;
  q << 3, 0, 1;
  s << 0, 0, 1;
  double h = 2;
  
  Eigen::VectorXd axis = Algebra::ScrewToAxis(q, s, h);
  Eigen::VectorXd result(6);
  result << 0, 0, 1, 0, -3, 2;
  
  ASSERT_TRUE(axis.isApprox(result, 4));
}

TEST(ALGEBRA, AxisAng6) {
  Eigen::VectorXd input(6);
  Eigen::VectorXd result(7);
  input << 1.0, 0.0, 0.0, 1.0, 2.0, 3.0;
  result << 1.0, 0.0, 0.0, 1.0, 2.0, 3.0, 1.0;
  
  Eigen::VectorXd output = Algebra::AxisAng6(input);
  ASSERT_TRUE(output.isApprox(result, 4));
}

TEST(ALGEBRA, MatrixLog6) {
  Eigen::MatrixXd Tinput(4, 4);
  Eigen::MatrixXd result(4, 4);
  Tinput << 1, 0, 0, 0,
    0, 0, -1, 0,
    0, 1, 0, 3,
    0, 0, 0, 1;
  
  result << 0, 0, 0, 0,
    0, 0, -1.57079633, 2.35619449,
    0, 1.57079633, 0, 2.35619449,
    0, 0, 0, 0;
  
  Eigen::MatrixXd Toutput = Algebra::MatrixLog6(Tinput);
  ASSERT_TRUE(Toutput.isApprox(result, 4));
}


TEST(ALGEBRA, DistanceToSO3Test) {
	Eigen::Matrix3d input;
	double result = 0.088353;
	input << 1.0, 0.0, 0.0,
		0.0, 0.1, -0.95,
		0.0, 1.0, 0.1;
	EXPECT_NEAR(result, Algebra::DistanceToSO3(input), 3);
}

TEST(ALGEBRA, DistanceToSE3Test) {
	Eigen::Matrix4d input;
	double result = 0.134931;
	input << 1.0, 0.0, 0.0, 1.2,
		0.0, 0.1, -0.95, 1.5,
		0.0, 1.0, 0.1, -0.9,
		0.0, 0.0, 0.1, 0.98;
	EXPECT_NEAR(result, Algebra::DistanceToSE3(input), 3);
}

TEST(ALGEBRA, TestIfSO3Test) {
	Eigen::Matrix3d input;
	bool result = false;
	input << 1.0, 0.0, 0.0,
          0.0, 0.1, -0.95,
		0.0, 1.0, 0.1;
	ASSERT_EQ(result, Algebra::TestIfSO3(input));
}

TEST(ALGEBRA, TestIfSE3Test) {
	Eigen::Matrix4d input;
	bool result = false;
	input << 1.0, 0.0, 0.0, 1.2,
		0.0, 0.1, -0.95, 1.5,
		0.0, 1.0, 0.1, -0.9,
		0.0, 0.0, 0.1, 0.98;
	ASSERT_EQ(result, Algebra::TestIfSE3(input));
}

TEST(ALGEBRA, adTest) {
  Eigen::VectorXd V(6);
  V << 1, 2, 3, 4, 5, 6;
  
  Eigen::MatrixXd result(6, 6);
  result << 0, -3, 2, 0, 0, 0,
    3, 0, -1, 0, 0, 0,
    -2, 1, 0, 0, 0, 0,
    0, -6, 5, 0, -3, 2,
    6, 0, -4, 3, 0, -1,
    -5, 4, 0, -2, 1, 0;
  
  ASSERT_TRUE(Algebra::ad(V).isApprox(result, 4));
}

int main(int argc, char **argv) {
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
