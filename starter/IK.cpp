#include "IK.h"
#include "FK.h"
#include "minivectorTemplate.h"
#include <Eigen/Dense>
#include <adolc/adolc.h>
#include <cassert>
#if defined(_WIN32) || defined(WIN32)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

namespace
{

// Converts degrees to radians.
template<typename real>
inline real deg2rad(real deg) { return deg * M_PI / 180.0; }

template<typename real>
Mat3<real> Euler2Rotation(const real angle[3], RotateOrder order)
{
  Mat3<real> RX = Mat3<real>::getElementRotationMatrix(0, deg2rad(angle[0]));
  Mat3<real> RY = Mat3<real>::getElementRotationMatrix(1, deg2rad(angle[1]));
  Mat3<real> RZ = Mat3<real>::getElementRotationMatrix(2, deg2rad(angle[2]));

  switch(order)
  {
    case RotateOrder::XYZ:
      return RZ * RY * RX;
    case RotateOrder::YZX:
      return RX * RZ * RY;
    case RotateOrder::ZXY:
      return RY * RX * RZ;
    case RotateOrder::XZY:
      return RY * RZ * RX;
    case RotateOrder::YXZ:
      return RZ * RX * RY;
    case RotateOrder::ZYX:
      return RX * RY * RZ;
  }
  assert(0);
}

// Performs forward kinematics, using the provided "fk" class.
// This is the function whose Jacobian matrix will be computed using adolc.
// numIKJoints and IKJointIDs specify which joints serve as handles for IK:
//   IKJointIDs is an array of integers of length "numIKJoints"
// Input: numIKJoints, IKJointIDs, fk, eulerAngles (of all joints)
// Output: handlePositions (world-coordinate positions of all the IK joints; length is 3 * numIKJoints)
template<typename real>
void forwardKinematicsFunction(
    int numIKJoints, const int * IKJointIDs, const FK & fk,
    const std::vector<real> & eulerAngles, std::vector<real> & handlePositions)
{
  // Students should implement this.
  // The implementation of this function is very similar to function computeLocalAndGlobalTransforms in the FK class.
  // The recommended approach is to first implement FK::computeLocalAndGlobalTransforms.
  // Then, implement the same algorithm into this function. To do so,
  // you can use fk.getJointUpdateOrder(), fk.getJointRestTranslation(), and fk.getJointRotateOrder() functions.
  // Also useful is the multiplyAffineTransform4ds function in minivectorTemplate.h .
  // It would be in principle possible to unify this "forwardKinematicsFunction" and FK::computeLocalAndGlobalTransforms(),
  // so that code is only written once. We considered this; but it is actually not easily doable.
  // If you find a good approach, feel free to document it in the README file, for extra credit.

    //local transform
    int numFKJoints = fk.getNumJoints();
    
    std::vector< Mat3<real> > handleRotate(numFKJoints);
    std::vector< Vec3<real> > handleTranslate(numFKJoints);

    for (int i = 0; i < numFKJoints; i++)
    {
        int updateOrder = fk.getJointUpdateOrder(i);
        Vec3d restTranslation = fk.getJointRestTranslation(i);
        RotateOrder rotateOrder = fk.getJointRotateOrder(i);
       int parent = fk.getJointParent(i);

        real angles[3] = {eulerAngles[3 * i + 0],eulerAngles[3 * i + 1], eulerAngles[3 * i + 2]};
        Mat3<real> eulerToMatrix = Euler2Rotation(angles, rotateOrder);

        const Vec3d& orient = fk.getJointOrient(i);
        real orientArr[3] = { (real)orient.data()[0], (real)orient.data()[1], (real)orient.data()[2]};
        Mat3<real> jointOrientToMatrix = Euler2Rotation(orientArr, XYZ);

        Mat3<real> localRot = jointOrientToMatrix * eulerToMatrix;

        Vec3<real> localRestTranslation((real)restTranslation.data()[0], (real)restTranslation.data()[1],
            (real)restTranslation.data()[2]);

        if (parent < 0)
        {
            // Root
            handleRotate[i] = localRot;
            handleTranslate[i] = localRestTranslation;
        }
        else
        {
            multiplyAffineTransform4ds(handleRotate[parent], handleTranslate[parent],
                localRot, localRestTranslation, handleRotate[i], handleTranslate[i]);
        }
    }


    for (int handlePos = 0; handlePos < numIKJoints; handlePos++)
    {
        int jointID = IKJointIDs[handlePos];

        handlePositions[3 * handlePos + 0] = handleTranslate[jointID][0];
        handlePositions[3 * handlePos + 1] = handleTranslate[jointID][1];
        handlePositions[3 * handlePos + 2] = handleTranslate[jointID][2];
    }
}

} // end anonymous namespaces

IK::IK(int numIKJoints, const int * IKJointIDs, FK * inputFK, int adolc_tagID)
{
  this->numIKJoints = numIKJoints;
  this->IKJointIDs = IKJointIDs;
  this->fk = inputFK;
  this->adolc_tagID = adolc_tagID;

  FKInputDim = fk->getNumJoints() * 3;
  FKOutputDim = numIKJoints * 3;

  train_adolc();
}

void IK::train_adolc()
{
  // Students should implement this.
  // Here, you should setup adol_c:
  //   Define adol_c inputs and outputs. 
  //   Use the "forwardKinematicsFunction" as the function that will be computed by adol_c.
  //   This will later make it possible for you to compute the gradient of this function in IK::doIK
  //   (in other words, compute the "Jacobian matrix" J).
  // See ADOLCExample.cpp .

    using adouble = adouble;

    std::vector<adouble> x(FKInputDim);
    std::vector<adouble> y(FKOutputDim);
    std::vector<double> output(FKOutputDim);

    trace_on(adolc_tagID);

    for (int i = 0; i < FKInputDim; i++)
    {
        x[i] <<= 0.0;
    }

    forwardKinematicsFunction(numIKJoints, IKJointIDs, *fk, x, y);

    for (int i = 0; i < FKOutputDim; i++)
    {
        y[i] >>= output[i];
    }

    trace_off();

}

void IK::doIK(const Vec3d* targetHandlePositions, Vec3d* jointEulerAngles)
{
    // You may find the following helpful:
    int numJoints = fk->getNumJoints(); // Note that is NOT the same as numIKJoints!

    // Students should implement this.
    // Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
    // Specifically, use ::function, and ::jacobian .
    // See ADOLCExample.cpp .
    //
    // Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
    // Note that at entry, "jointEulerAngles" contains the input Euler angles. 
    // Upon exit, jointEulerAngles should contain the new Euler angles.


    //vector of current euler angles of all joints, length is 3 * numJoints
    std::vector<double> theta(FKInputDim);
    for (int i = 0; i < numJoints; i++)
    {
        theta[3 * i + 0] = jointEulerAngles[i][0];
        theta[3 * i + 1] = jointEulerAngles[i][1];
        theta[3 * i + 2] = jointEulerAngles[i][2];
    }

    //calculate new angle positions, run the funcition w adolc
    std::vector<double> f(FKOutputDim);
    ::function(adolc_tagID, FKOutputDim, FKInputDim, theta.data(), f.data());

    //calculate change between current global angle positions and target global angle positions
    Eigen::VectorXd delta_b(FKOutputDim);
    for (int h = 0; h < numIKJoints; h++)
    {
        delta_b[3 * h + 0] = targetHandlePositions[h][0] - f[3 * h + 0];
        delta_b[3 * h + 1] = targetHandlePositions[h][1] - f[3 * h + 1];
        delta_b[3 * h + 2] = targetHandlePositions[h][2] - f[3 * h + 2];
    }

    Eigen::MatrixXd jacobianMatrix(FKOutputDim, FKInputDim);
    std::vector<double> jacobianBuffer(FKOutputDim * FKInputDim);
    std::vector<double*> jacobianRows(FKOutputDim);
    for (int i = 0; i < FKOutputDim; ++i)
    {
        jacobianRows[i] = &jacobianBuffer[i * FKInputDim];
    }
    //we want the smallest error between change of Euler angles we want to find(J delta_theta) 
    // and change of handle global positions (theta_b) (calculate derivatives)
    ::jacobian(adolc_tagID, FKOutputDim, FKInputDim, theta.data(), jacobianRows.data());

    // Copy the jacobianBuffer into the Eigen matrix J
    for (int row = 0; row < FKOutputDim; ++row) 
    {
        for (int col = 0; col < FKInputDim; ++col) 
        {
            jacobianMatrix(row, col) = jacobianBuffer[row * FKInputDim + col];
        }   
    }

    Eigen::VectorXd delta_theta;
    if (!isPseudoInverse)
    {
        double alpha = 1e-3;

        //(J^TJ + aI)delta_theta = J^T * delta_b
        Eigen::MatrixXd A = jacobianMatrix.transpose() * jacobianMatrix + alpha * Eigen::MatrixXd::Identity(FKInputDim, FKInputDim);
        Eigen::VectorXd b = jacobianMatrix.transpose() * delta_b;

         delta_theta = A.ldlt().solve(b);
    }
    else
    {
        // Pseudoinverse
        // delta_theta = J^T (J * J^T)^(-1) delta_b
        Eigen::MatrixXd JJt = jacobianMatrix * jacobianMatrix.transpose();

        // Solve for x in  (J J^T) x = delta_b
        Eigen::VectorXd x = JJt.ldlt().solve(delta_b);

        // delta_theta = J^T x
       delta_theta = jacobianMatrix.transpose() * x;
	}

  //update angles (output)
  for (int i = 0; i < numJoints; i++)
  {
      jointEulerAngles[i][0] += delta_theta[3 * i + 0];
      jointEulerAngles[i][1] += delta_theta[3 * i + 1];
      jointEulerAngles[i][2] += delta_theta[3 * i + 2];
  }
}

