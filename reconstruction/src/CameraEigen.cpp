//////////////////////////////////////////////////////////////////////////
// Author       :  Ngo Tien Dat
// Email        :  dat.ngo@epfl.ch
// Organization :  EPFL
// Purpose      :  Represent camera object and its operations
// Date         :  19 March 2012
//////////////////////////////////////////////////////////////////////////

#include "Camera.hpp"

using namespace Eigen;

MatrixX2d Camera::ProjectPoints(const MatrixX3d points3D)
{
  int nVertices = points3D.rows();

  // Projected points
  Matrix<float,nVertices,2> points2D;

  for (int i = 0; i < nVertices; i++)
  {    
    Vector3d point3D = points3D.row(i);
    points2D.row(i) = ProjectAPoint(point3D);
  }

  return points2D;
}

Vector2d Camera::ProjectAPoint(Vector3d point3D) const
{
  // Add homogeneous component
  point3D.resize(4); //?
  point3D(3) = 1;

  Vector3d pointUVW = ARt * point3D;

  double u = pointUVW(0) / pointUVW(2);
  double v = pointUVW(1) / pointUVW(2);

  Vector2d point2D;
  point2D << u << v;

  return point2D;
}
