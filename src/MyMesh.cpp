#include "../include/stdafx.h"
#include "../include/MyMesh.h"

string& trim(string &s)
{
	if (s.empty())
	{
		return s;
	}
	s.erase(0, s.find_first_not_of(" "));
	s.erase(s.find_last_not_of(" ") + 1);
	return s;
}

MyMesh::MyMesh()
{
	this->numberOfFaces = 0;
	this->numberOfVertices = 0;
}

MyMesh::MyMesh(const MyMesh & myMesh)
{
	this->faces = myMesh.faces;
	this->numberOfFaces = myMesh.numberOfFaces;
	this->vertices = myMesh.vertices;
	this->numberOfVertices = myMesh.numberOfVertices;
}

/*
/// @brief      overload constructor
/// @details    
*/
MyMesh::MyMesh(const Matrix<double, 3, Dynamic> &vertices, const Matrix<int, 3, Dynamic> &faces)
{
	this->vertices = vertices;
	this->numberOfVertices = vertices.cols();
	this->faces = faces;
	this->numberOfFaces = faces.cols();
}

/*
/// @brief      load a mesh from file
/// @details    
/// @param[in]  fileName: 
/// @return     
/// @attention  
*/
MyMesh::MyMesh(string fileName)
{
	MatrixXd V;
	MatrixXi F;
	igl::readOFF(fileName, V, F);
	this->vertices = V.transpose();
	this->numberOfVertices = V.rows();
	this->faces = F.transpose();
	this->numberOfFaces = F.rows();
	/*double x, y, z;
	//int v0, v1, v2;
	//// Container holding last line read
	//string readLine;
	//// Containers for delimiter positions
	//int delimiterPos_1, delimiterPos_2, delimiterPos_3, delimiterPos_4;
	//// Open file for reading
	//ifstream in(fileName.c_str());
	//// Check if file is in OFF format
	//getline(in, readLine);
	//if (readLine != "OFF")
	//{
	//	cout << "The file to read is not in OFF format." << endl;
	//	return;
	//}
	//// Read values for Nv and Nf
	//while (getline(in, readLine))
	//{
	//	readLine = trim(readLine);
	//	if (!readLine.empty() && readLine[0] != '#')
	//	{
	//		break;
	//	}
	//}
	//delimiterPos_1 = readLine.find(" ", 0);
	//int nv = atoi(readLine.substr(0, delimiterPos_1 + 1).c_str());
	//delimiterPos_2 = readLine.find(" ", delimiterPos_1);
	//int nf = atoi(readLine.substr(delimiterPos_1, delimiterPos_2 + 1).c_str());

	//this->vertices.resize(3, nv);
	//this->numberOfVertices = nv;
	//this->faces.resize(3, nf);
	//this->numberOfFaces = nf;

	//// Read the vertices
	//for (int n = 0; n < nv; n++)
	//{
	//	getline(in, readLine);
	//	delimiterPos_1 = readLine.find(" ", 0);
	//	x = atof(readLine.substr(0, delimiterPos_1).c_str());
	//	delimiterPos_2 = readLine.find(" ", delimiterPos_1 + 1);
	//	y = atof(readLine.substr(delimiterPos_1, delimiterPos_2).c_str());
	//	delimiterPos_3 = readLine.find(" ", delimiterPos_2 + 1);
	//	z = atof(readLine.substr(delimiterPos_2, delimiterPos_3).c_str());
	//	this->vertices(0, n) = x;
	//	this->vertices(1, n) = y;
	//	this->vertices(2, n) = z;
	//}
	//// Read the facades
	//for (int n = 0; n < nf; n++)
	//{
	//	getline(in, readLine);
	//	delimiterPos_1 = readLine.find(" ", 0);
	//	delimiterPos_2 = readLine.find(" ", delimiterPos_1 + 1);
	//	v0 = atoi(readLine.substr(delimiterPos_1, delimiterPos_2).c_str());
	//	delimiterPos_3 = readLine.find(" ", delimiterPos_2 + 1);
	//	v1 = atoi(readLine.substr(delimiterPos_2, delimiterPos_3).c_str());
	//	delimiterPos_4 = readLine.find(" ", delimiterPos_3 + 1);
	//	v2 = atoi(readLine.substr(delimiterPos_3, delimiterPos_4).c_str());
	//	this->faces(0, n) = v0;
	//	this->faces(1, n) = v1;
	//	this->faces(2, n) = v2;
	}*/
}

Matrix<double, 3, Dynamic> MyMesh::getVertices()
{
	return this->vertices;
}

Matrix<int, 3, Dynamic> MyMesh::getFaces()
{
	return this->faces;
}

void MyMesh::setVertices(Matrix<double, 3, Dynamic> &v)
{
	this->vertices = v;
	this->numberOfVertices = v.cols();
}

void MyMesh::setFaces(Matrix<int, 3, Dynamic> &f)
{
	this->faces = f;
	this->numberOfFaces = f.cols();
}

int MyMesh::getNumberOfVertices()
{
	return this->numberOfVertices;
}

int MyMesh::getNumberOfFaces()
{
	return this->numberOfFaces;
}

void MyMesh::copyNumberFromMesh( const MyMesh &myMesh)
{
	this->faces = myMesh.faces;
	this->numberOfFaces = myMesh.numberOfFaces;
	this->vertices = myMesh.vertices;
	this->numberOfVertices = myMesh.numberOfVertices;
}

/*
/// @brief      compute area of all faces
/// @details    
/// @param[in]  : 
/// @return     
/// @attention  
*/
VectorXd MyMesh::areaFaces()
{
	VectorXd area = VectorXd::Zero(this->numberOfFaces);
	for (int i = 0; i < this->numberOfFaces; ++i)
	{
		area[i] = 1.0/((this->vertices.col(this->faces(1, i)) - this->vertices.col(this->faces(0, i))).cross(this->vertices.col(this->faces(2, i)) - this->vertices.col(this->faces(0, i)))).norm();
	}
	return area;
}


/*
/// @brief      mean value coordinates in 3D
/// @details
/// @param[in]  p:
/// @return
/// @attention
*/
VectorXd MyMesh::meanValueCoordinates(Vector3d p)
{
	//this->faces = this->fixMesh().faces;   // repair the faces index, so that the all normal vector of faces outward (inward)
	VectorXd mvc_coor = VectorXd::Zero(this->numberOfVertices);

	// get projection vertices on unit sphere
	Matrix<double, 3, Dynamic> pro_vertices;

	pro_vertices.resize(3, this->numberOfVertices);
	for (int i = 0; i < this->numberOfVertices; ++i)
	{
		double r = (this->vertices.col(i) - p).norm();
		pro_vertices(0, i) = (this->vertices(0, i) - p(0)) / r;
		pro_vertices(1, i) = (this->vertices(1, i) - p(1)) / r;
		pro_vertices(2, i) = (this->vertices(2, i) - p(2)) / r;
	}

	Matrix<double, 3, 3> T;
	Vector3d n0, n1, n2;
	double beta0, beta1, beta2;
	for (int i = 0; i < this->numberOfFaces; ++i)
	{
		int i0, i1, i2;		// index of vertices on the current face
		i0 = this->faces(0, i);
		i1 = this->faces(1, i);
		i2 = this->faces(2, i);

		T.col(0) = pro_vertices.col(i0);
		T.col(1) = pro_vertices.col(i1);
		T.col(2) = pro_vertices.col(i2);

		n0 = pro_vertices.col(i0).cross(pro_vertices.col(i1)).normalized();
		n1 = pro_vertices.col(i1).cross(pro_vertices.col(i2)).normalized();
		n2 = pro_vertices.col(i2).cross(pro_vertices.col(i0)).normalized();

		// angle
		beta0 = acos(pro_vertices.col(i0).dot(pro_vertices.col(i1)));
		beta1 = acos(pro_vertices.col(i1).dot(pro_vertices.col(i2)));
		beta2 = acos(pro_vertices.col(i2).dot(pro_vertices.col(i0)));

		mvc_coor[i0] += (beta1 + beta0 * n0.dot(n1) + beta2 * n2.dot(n1)) / (2 * pro_vertices.col(i0).dot(n1));
		mvc_coor[i1] += (beta2 + beta1 * n1.dot(n2) + beta0 * n0.dot(n2)) / (2 * pro_vertices.col(i1).dot(n2));
		mvc_coor[i2] += (beta0 + beta2 * n2.dot(n0) + beta1 * n1.dot(n0)) / (2 * pro_vertices.col(i2).dot(n0));
	}
	for (int i = 0; i < this->numberOfVertices; ++i)
	{
		mvc_coor[i] /= (this->vertices.col(i) - p).norm();
	}
	mvc_coor /= mvc_coor.sum();
	return mvc_coor;
}


/*
/// @brief      mean value coordinates in 3D
/// @details    
/// @param[in]  p: 
/// @return     
/// @attention  
*/ 
VectorXd MyMesh::meanValueCoordinates(Vector3d p, bool &postiveity)
{
	postiveity = true;
	VectorXd mvc_coor = this->meanValueCoordinates(p);
	postiveity = mvc_coor.minCoeff() >= 0;
	return mvc_coor;
}

/*
/// @brief      Write the mesh to file
/// @details    
/// @param[in]  outputFileName: 
/// @return     
/// @attention  
*/
bool MyMesh::writeMeshToFile(string outputFileName)
{
	try
	{
		ofstream outfile(outputFileName);
		outfile << "OFF" << endl;
		outfile << this->numberOfVertices << " " << this->numberOfFaces << " 0" << endl;
		for (int i = 0; i < this->numberOfVertices; ++i)
		{
			outfile << this->vertices(0, i) << " " << this->vertices(1, i) << " " << this->vertices(2, i) << endl;
		}
		for (int i = 0; i < this->numberOfFaces; ++i)
		{
			outfile << "3 " << this->faces(0, i) << " " << this->faces(1, i) << " " << this->faces(2, i) << endl;
		}
		outfile.close();
		return true;
	}
	catch (const std::exception&)
	{
		return false;
	}
}

/*
/// @brief      Repair the model so that the normal direction of all faces are consistent
/// @details    Using a tool: meshfix.exe
/// @param[in]  : 
/// @return     
/// @attention  
*/
bool MyMesh::fixMesh()
{
	try
	{
		this->writeMeshToFile("temp.off");
		system("meshfix temp.off output.off");
		copyNumberFromMesh(MyMesh("output.off"));
		return true;
	}
	catch (const std::exception&)
	{
		return false;
	}
}

/*
/// @brief      to do sqrt(3) subdivision for polyhedron
/// @details    new face vertex is v, the topology rule is same as sqrt(3)-subdivision scheme.
/// @param[in]  v: 
/// @return     
/// @attention  
*/
MyMesh MyMesh::sqrt3subdivision(Matrix<double, 3, Dynamic> &newv)
{
	Matrix<double, 3, Dynamic> v;
	v.resize(3, this->numberOfVertices + this->numberOfFaces);
	v.leftCols(this->numberOfVertices) = this->vertices;
	v.rightCols(this->numberOfFaces) = newv;

	Matrix<int, 3, Dynamic> f = -MatrixXi::Ones(3, 6 * this->numberOfFaces);
	Matrix<int, 3, Dynamic> tf = -MatrixXi::Ones(3, 3 * this->numberOfFaces);
	Matrix<int, 3, Dynamic> flag_face = MatrixXi::Zero(3, this->numberOfFaces);

	MatrixXi face_ring_each = this->computeFaceRing();
	vector<int> index;
	for (int i = 0; i < this->numberOfFaces; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (!flag_face(j, i))  // if j-th edge of current triangle don't conection
			{
				f.col(6 * i + 2 * j) = Vector3i(this->numberOfVertices + i, this->faces((j + 1) % 3, i), this->numberOfVertices + face_ring_each(j, i));
				index.push_back(6 * i + 2 * j);
				f.col(6 * i + 2 * j + 1) = Vector3i(this->numberOfVertices + i, this->numberOfVertices + face_ring_each(j, i), this->faces((j + 2) % 3, i));
				index.push_back(6 * i + 2 * j + 1);
				for (int k = 0; k < 3; ++k)
				{
					if (face_ring_each(k, face_ring_each(j, i)) == i)
					{
						flag_face(k, face_ring_each(j, i)) = 1;
						break;
					}
				}
			}
		}
	}
	int j = 0;
	for (vector<int>::iterator iter = index.begin(); iter != index.end(); ++iter)
	{
		tf.col(j++) = (f.col(*iter));
	}
	return MyMesh(v,tf);
}

/*
/// @brief      compute ring of each face
/// @details    
/// @param[in]  : 
/// @return     
/// @attention  
*/
MatrixXi MyMesh::computeFaceRing()
{
	MatrixXi edge = -MatrixXi::Ones(2, 3 * this->numberOfFaces);
	edge.leftCols(this->numberOfFaces) = this->faces.bottomRows(2);
	edge.rightCols(this->numberOfFaces) = this->faces.topRows(2);
	edge.middleCols(this->numberOfFaces, this->numberOfFaces).row(0) = this->faces.row(0);
	edge.middleCols(this->numberOfFaces, this->numberOfFaces).row(1) = this->faces.row(2);
	int temp;
	for (int i = 0; i < edge.cols(); ++i)
	{
		if (edge(0, i) > edge(1, i))
		{
			temp = edge(0, i);
			edge(0, i) = edge(1, i);
			edge(1, i) = temp;
		}
	}
	MatrixXi face_ring_each = MatrixXi::Zero(3, this->numberOfFaces);
	int x = -1, y = -1, z = -1;
	for (int i = 0; i < this->numberOfFaces; ++i)
	{
		for (int j = 0; j < edge.cols(); ++j)
		{
			if (j % this->numberOfFaces == i) continue;
			if (edge.col(i) == edge.col(j))
			{
				x = j % this->numberOfFaces;
			}
			if (edge.col(i + this->numberOfFaces) == edge.col(j))
			{
				y = j % this->numberOfFaces;
			}
			if (edge.col(i + 2 * this->numberOfFaces) == edge.col(j))
			{
				z = j % this->numberOfFaces;
			}
		}
		face_ring_each.col(i) = Vector3i(x, y, z);
		x = 0; y = 0; z = 0;
	}
	return face_ring_each;
}


/*
/// @brief      find all adjacent faces of the v_i
/// @details    
/// @param[in]  int indexOfVertex: 
/// @return     
/// @attention  
*/
vector<int> MyMesh::findAdjacentFaces(int indexOfVertex)
{
	vector<int> adjacentFaces;
	for (int i = 0; i < this->numberOfFaces; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			if (this->faces(j,i)==indexOfVertex)
			{
				adjacentFaces.push_back(i);
				break;
			}
		}
	}
	return adjacentFaces;
}

/*
/// @brief      distance from p to face_i
/// @details    
/// @param[in]  p: 
/// @return     
/// @attention  
*/
VectorXd MyMesh::distancePointFace(Vector3d p)
{
	VectorXd dist = VectorXd::Zero(this->numberOfFaces);
	VectorXd area = this->areaFaces();
	Matrix3d T = Matrix3d::Zero(3, 3);
	for (int i = 0; i < this->numberOfFaces; ++i)
	{
		T.col(0) = this->vertices.col(this->faces(0, i)) - p;
		T.col(1) = this->vertices.col(this->faces(1, i)) - p;
		T.col(2) = this->vertices.col(this->faces(2, i)) - p;
		dist[i] = T.determinant() / 2 / area[i];
	}
	return dist;
}

MyMesh::~MyMesh()
{
}
