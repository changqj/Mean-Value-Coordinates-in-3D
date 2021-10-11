/**********************************************************************************
/// @file       MyMesh.h
/// @author     Qingjun Chang
/// @date       2019.05.04
/// @brief      This is a c++ class about mesh
/// @details    
**********************************************************************************/
#pragma once
#include <Eigen/Dense>
#include <string.h>
#include <iostream>
#include <fstream>
#include <Eigen/Geometry>
#include <stdlib.h>
#include <vector>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include "tutorial_shared_path.h"

using namespace std;
using namespace Eigen;

class MyMesh
{
private:
	Matrix<double, 3, Dynamic> vertices;	// vertex set of the mesh
	int numberOfVertices;
	Matrix<int, 3, Dynamic> faces;			// face set of the mesh
	int numberOfFaces;
public:
	MyMesh();
	MyMesh(const MyMesh &myMesh);
	MyMesh(const Matrix<double, 3, Dynamic> &v, const Matrix<int, 3, Dynamic> &f);	// overload constructor
	MyMesh(string fileName);		// load mesh from file
	Matrix<double, 3, Dynamic> getVertices();		// get method
	Matrix<int, 3, Dynamic> getFaces();				// get method
	void setVertices(Matrix<double, 3, Dynamic>&v);		// set method
	void setFaces(Matrix<int, 3, Dynamic> &f);				// set method
	int getNumberOfVertices();
	int getNumberOfFaces();
	void copyNumberFromMesh(const MyMesh &myMesh);
	VectorXd areaFaces();	// compute area of each face
	VectorXd meanValueCoordinates(Vector3d p,bool &postiveity);		// get mean value coordinates of point p
	VectorXd meanValueCoordinates(Vector3d p);		// get mean value coordinates of point p
	bool writeMeshToFile(string outputFileName);  // write the mesh to file
	bool fixMesh();	// repair the model
	MyMesh sqrt3subdivision(Matrix<double, 3, Dynamic> &v); // new face vertex is v, to do sqrt(3)-subdivision.
	MatrixXi computeFaceRing();  // compute ring of each face
	vector<int> findAdjacentFaces(int indexOfVertex);  // find all adjacent faces of v_i
	VectorXd distancePointFace(Vector3d p);
	virtual ~MyMesh();
};

