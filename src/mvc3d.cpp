#include "../include/stdafx.h"
#include "../include/MyMesh.h"
#include <time.h>
#include <igl/cotmatrix.h>
#include<Eigen/Sparse>
#include <igl/readOFF.h>
#include <omp.h>

#define PI 3.1415926535897932384626433832795

int main(int argc, char* argv[])
{
	// dos : MVC3D \bunny.off 0.2 0.5 0.1 \your-output-path
	
	//argv[1] = (char*) ("C:/Users/user/3D Objects/off/bunny_super_low.off");
	//argv[2] = (char*)"0.5";
	//argv[3] = (char*)"0.50";
	//argv[4] = (char*)"0.1";
	//argv[5] = (char*)"/res";
	
	//std::cout << "test" << std::endl;
	//std::cout << argv[1] << std::endl;
	//std::cout << argv[2] << std::endl;
	//std::cout << argv[3] << std::endl;
	//std::cout << argv[4] << std::endl;
	string model_name = argv[1];

	double x = std::stod(argv[2]);
	double y = std::stod(argv[3]);
	double z = std::stod(argv[4]);

	string outputPath = argv[5];

	// load polyhedron model

	std::cout << model_name << std::endl;

	MyMesh myMesh(model_name);

	ofstream out(outputPath);
	int nv, nf;
	nv = myMesh.getNumberOfVertices();
	nf = myMesh.getNumberOfFaces();
	cout << "Number of vertices of the model:" << nv << endl;
	cout << "Number of faces of the model:" << nf << endl;

	int isBoundary;
	IOFormat CommaInitFmt(-1, 0, "\t", "\n", "", "", "[", "]");
	IOFormat CommaInitFmt_out(-1, 0, "\t", "\n", "", "", "", "");

	Vector3d p(x, y, z);

	out << myMesh.meanValueCoordinates(p).transpose().format(CommaInitFmt_out) << endl;

	cout << "mean value coordinates: " << endl;
	cout << myMesh.meanValueCoordinates(p).transpose().format(CommaInitFmt_out) << endl;
	cout << myMesh.getVertices() * myMesh.meanValueCoordinates(p) << endl;

	
	return 0;
}

