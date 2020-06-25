//////////////////////////////////////////////////////////////////////////
// Author    :  Ngo Tien Dat
// Email    :  dat.ngo@epfl.ch
// Organization  :  EPFL
// Purpose    :  Triangle mesh class
// Date      :  12 March 2012
//////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <assert.h>
#include <armadillo>
#include <Camera.hpp>
#include <eigen3/Eigen/Eigen>
using namespace Eigen;

class TriangleMesh
{

protected:
  // Uninitialized mat and vec have the size of zero.

  MatrixX3d   vertexCoords;  // Coordinates of vertices. Size of #vertices * 3. Number 3 means (x,y,z)
  MatrixX3d    facets;        // Facets of the mesh. Size of #facets * 3. (vertexId1,vertexId2,vertexId3) vertex id starts from 0

  // ----- Extended dependent attributes ------
  MatrixX2d edges;          // Edges. Size of 2 * #edges. Each represetned by two vertex ids.

 VectorXd  edgeLengths;    // Edge lengths, a column vector
 VectorXd XYZColumn;      // Vectorized coordinates of vertices in a column vector with the
                             // order [x1,...,xN,y1,...,yN,z1,...,zN]

  MatrixX3d  facetNormals;    // Facet normals. Size of #facets * 3
  VectorXd  facetNormalMags; // Facet normal magnitudes, a column vector //??
  MatrixX3d  facetCentroids;  // Size of #facets * 3
  MatrixX3d  vertexNormals;   // Size of #vertices * 3
  MatrixX2d facetPairs;      // pairs of adjacent facets, size of #pairs * 2 (two facet ids)

public:
  // Constructor
  TriangleMesh(void) {}

  TriangleMesh( MatrixX3d coords,  MatrixX3d faces);

  // Destructor
  virtual ~TriangleMesh(void) {}

  // Load the triangle mesh from file pair ".pts" and ".trig"
  // Input:
  //  + baseFile: root file name without .pts and .trig
  void Load(std::string baseFile);

  // Transform the mesh into camera coordinate
  // Multiply vertex with Rt matrix and call set vertex function
  void TransformToCameraCoord(const Camera& worldCamera);

  // Get facets
  const  MatrixX2d GetFacets() const
  {
    return facets;
  }

  // Get edges: size of 2 * #edges. Each represetned by two vertex ids.
  const  MatrixX2d GetEdges() const
  {
    return edges;
  }

  // Get number of edges
  int GetNEdges() const
  {
    return edges.rows() ;
  }

  // Get vertex coordinates: size of #vertices * 3
  const  MatrixX3d GetVertexCoords() const
  {
    return vertexCoords;
  }

  // Set vertex coordinates. Need to re-compute dependent attributes
  virtual void SetVertexCoords(const arma::mat& vtCoords);

  // Get number of vertices
  int GetNVertices() const
  {
    return this->vertexCoords.rows();
  }

  // Get number of vertices
  int GetNFacets() const
  {
    return this->facets.rows();
  }

  // Get all pairs of adjacent facets
  const  MatrixX2d GetFacetPairs() const
  {
    return facetPairs;
  }

  // Get number of facet pairs
  int GetNFacetPairs() const
  {
    return facetPairs.rows();
  }

  // Get the vectorized coordinates of the mesh vertices with the order
  // [x1,...,xN,y1,...,yN,z1,...,zN]
  const VectorXd GetXYZColumn() const
  {
    return XYZColumn;
  }

  // Get edge lengths
  const VectorXd GetEdgeLengths() const
  {
    // Be sure that edge lengths were already computed
    assert(edgeLengths.rows() != 0); //?
    return this->edgeLengths;
  }

  // Compute edge lengths. Update edge length and return the reference as well
  const VectorXd ComputeEdgeLengths();

  const  MatrixX3d GetFacetNormals() const
  {
    return this->facetNormals;
  }

  // Compute facet normals
  void computeFacetNormalsNCentroids();

protected:

  // Compute edges from vertices and facets
  void computeEdges();

  // Compute neighborhood relations from facets
  // Return a boolean matrix: 1 mean true, 0 mean false
  arma::SpMat<char> computeVertexNeighborhoods();

  // Compute facet pairs that share an edge
  void computeFacetPairs();
};

