#include "skinning.h"
#include "vec3d.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

Skinning::Skinning(int numMeshVertices, const double * restMeshVertexPositions,
    const std::string & meshSkinningWeightsFilename)
{
  this->numMeshVertices = numMeshVertices;
  this->restMeshVertexPositions = restMeshVertexPositions;

  cout << "Loading skinning weights..." << endl;
  ifstream fin(meshSkinningWeightsFilename.c_str());
  assert(fin);
  int numWeightMatrixRows = 0, numWeightMatrixCols = 0;
  fin >> numWeightMatrixRows >> numWeightMatrixCols;
  assert(fin.fail() == false);
  assert(numWeightMatrixRows == numMeshVertices);
  int numJoints = numWeightMatrixCols;

  vector<vector<int>> weightMatrixColumnIndices(numWeightMatrixRows);
  vector<vector<double>> weightMatrixEntries(numWeightMatrixRows);
  fin >> ws;
  while(fin.eof() == false)
  {
    int rowID = 0, colID = 0;
    double w = 0.0;
    fin >> rowID >> colID >> w;
    weightMatrixColumnIndices[rowID].push_back(colID);
    weightMatrixEntries[rowID].push_back(w);
    assert(fin.fail() == false);
    fin >> ws;
  }
  fin.close();

  // Build skinning joints and weights.
  numJointsInfluencingEachVertex = 0;
  for (int i = 0; i < numMeshVertices; i++)
    numJointsInfluencingEachVertex = std::max(numJointsInfluencingEachVertex, (int)weightMatrixEntries[i].size());
  assert(numJointsInfluencingEachVertex >= 2);

  // Copy skinning weights from SparseMatrix into meshSkinningJoints and meshSkinningWeights.
  meshSkinningJoints.assign(numJointsInfluencingEachVertex * numMeshVertices, 0);
  meshSkinningWeights.assign(numJointsInfluencingEachVertex * numMeshVertices, 0.0);
  for (int vtxID = 0; vtxID < numMeshVertices; vtxID++)
  {
    vector<pair<double, int>> sortBuffer(numJointsInfluencingEachVertex);
    for (size_t j = 0; j < weightMatrixEntries[vtxID].size(); j++)
    {
      int frameID = weightMatrixColumnIndices[vtxID][j];
      double weight = weightMatrixEntries[vtxID][j];
      sortBuffer[j] = make_pair(weight, frameID);
    }
    sortBuffer.resize(weightMatrixEntries[vtxID].size());
    assert(sortBuffer.size() > 0);
    sort(sortBuffer.rbegin(), sortBuffer.rend()); // sort in descending order using reverse_iterators
    for(size_t i = 0; i < sortBuffer.size(); i++)
    {
      meshSkinningJoints[vtxID * numJointsInfluencingEachVertex + i] = sortBuffer[i].second;
      meshSkinningWeights[vtxID * numJointsInfluencingEachVertex + i] = sortBuffer[i].first;
    }

    // Note: When the number of joints used on this vertex is smaller than numJointsInfluencingEachVertex,
    // the remaining empty entries are initialized to zero due to vector::assign(XX, 0.0) .
  }
}

void Skinning::applySkinning(const RigidTransform4d * jointSkinTransforms, double * newMeshVertexPositions) const
{
  // Students should implement this
  // The following below is just a dummy implementation.
  
     // position p for mesh vertex i
      // - pi is the homogeneous coordinates of the rest(a.k.a undeformed) position of the vertex
      // j goes over all the joints that affect the vertex
      // wj is the skinning weight of joint j to the vertex
      // Msj is the joint's skinning transform matrix
  for(int i=0; i<numMeshVertices; i++)
  {
	  Vec4d newPos(0.0, 0.0, 0.0, 0.0);

      for (int joint = 0; joint < numJointsInfluencingEachVertex; joint++)
      {
		  int jointID = meshSkinningJoints[i * numJointsInfluencingEachVertex + joint];
		  double weight = meshSkinningWeights[i * numJointsInfluencingEachVertex + joint];

          Mat4d jointSkinTransformMatrix = jointSkinTransforms[jointID];
		  Vec4d restMeshVertexPosition = Vec4d(restMeshVertexPositions[3 * i + 0], 
              restMeshVertexPositions[3 * i + 1], restMeshVertexPositions[3 * i + 2], 1.0);

          newPos += weight * (jointSkinTransformMatrix * restMeshVertexPosition);
		  
      }
        newMeshVertexPositions[3 * i + 0] = newPos[0];
        newMeshVertexPositions[3 * i + 1] = newPos[1];
        newMeshVertexPositions[3 * i + 2] = newPos[2];



    //newMeshVertexPositions[3 * i + 0] = restMeshVertexPositions[3 * i + 0];
    //newMeshVertexPositions[3 * i + 1] = restMeshVertexPositions[3 * i + 1];
    //newMeshVertexPositions[3 * i + 2] = restMeshVertexPositions[3 * i + 2];
  }
}

