//////////////////////////////////////////////////////////////////////////
// Author       :  Ngo Tien Dat
// Email        :  dat.ngo@epfl.ch
// Organization :  EPFL
// Purpose      :  Represent camera object and its operations
// Date         :  19 March 2012
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include "eigen3/Eigen/Eigen"

class Camera
{

private: 
  Eigen::Matrix<float,3,3> A;      // Intrinsic matrix  3x3
  Eigen::Matrix<float,3,4> Rt;     // Extrinsic matrix  3x4
  Eigen::Matrix<float,3,4> ARt;    // Projection matrix 3x4

public:

  // Default constructor
  Camera() {}

  // Constructor that create a camera in camera coordinate. Rt = [I|0]
  Camera(const Eigen::Matrix<float,3,3> A)
  {
    this->A = A;
    Rt << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0;

    ARt = A * Rt;
  }

  // Constructor that create a camera given A and Rt
  Camera(Eigen::Matrix<float,3,3> A, Eigen::Matrix<float,3,4> Rt)
  {
    this->A  = A;
    this->Rt = Rt;

    ARt = A * Rt;
  }

//  // Load from intrinsic and extrinsic matrix files
//  void LoadFromFile(std::string camIntrFile, std::string camExtFile)
//  {
//    A. load(camIntrFile);
//    Rt.load(camExtFile);

//    ARt = A * Rt;
//  }

//  // Load from projection matrix file
//  void LoadFromFile(std::string camProjFile)
//  {
//    ARt.load(camProjFile);

//    // TODO: Compute A, Rt using QR decomposition
//  }

  // Get intrinsic matrix A
  const Eigen::Matrix<float,3,3>& GetA() const
  {
    return A;
  }

  // Get extrinsic matrix Rt
  const Eigen::Matrix<float,3,4>& GetRt()  const
  {
    return Rt;
  }

  // Get projection matrix ART
  const Eigen::Matrix<float,3,4>& GetARt() const {
    return ARt;
  }

  // Project a list of 3D points onto the image plane
  // Input:
  //  + 3D points: size of #points * 3
  // Output:
  //  + 2D points: size of #points * 2
  Eigen::MatrixX2d ProjectPoints(const Eigen::MatrixX3d points3D);

  // Project a 3D point onto the image plane
  // Input:
  //  + A 3D point: 3*1
  // Output:
  //  + A 2D point: vector of size 2*1
  Eigen::VectorXf ProjectAPoint(Eigen::VectorXf point3D) const;

};




