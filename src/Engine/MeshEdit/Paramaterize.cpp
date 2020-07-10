#include <Engine/MeshEdit/Paramaterize.h>
#include <Engine/MeshEdit/CParameterize.h>
#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>
#include <Engine/Scene/CmptGeometry.h>
#include <Engine/Scene/CmptMaterial.h>
#include <Engine/Material/BSDF_Frostbite.h>
#include <Engine/Scene/SObj.h>

#include <chrono>
#include <ctime> 

//#include <Basic/Math.h>
#define PI 3.14159

using namespace Ubpa;
using namespace Eigen;
using namespace std;

Paramaterize::Paramaterize(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_flaten, int mode, int pType)
	:CParameterize(triMeshObj, triMesh)
{
	Init(triMeshObj, triMesh, is_flaten, mode, pType);
}

void Paramaterize::Clear() {
	boundary_V.clear();
	coeff.clear();
	reshaped_boundary_V_pos.clear();
//	//texCoor.clear();
	//if(is_boundary_array != nullptr)
		//delete []is_boundary_array;
}

Paramaterize::~Paramaterize()
{
	Clear();
}

bool Paramaterize::Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex, int mode, int pType)
{
	Clear();
	this->mode = (mode == 0) ? Circle_boundary : Square_boundary;
	this->pType = (pType == 0) ? Naive : Cotan;
	this->is_tex = is_tex;

	//initial matrix
	mat_A = SparseMatrix<double>(nV, nV);
	mat_A.setZero();
	mat_bx = MatrixXd(nV, 3);
	mat_bx.setZero();
	is_boundary_array = new size_t[nV](); // initialize to be all 0 

	return true;
}

bool Paramaterize::Run() {
	//minSurf.Run();
	
	// First, get boundary information
	GetBoundaryInfo();

	// Second, set inner vertix information
	SetInnerVCoefficient();

	// Third, reshape boundary
	ReshapeBoundary(this->mode);

	// Fourth, after boundary is reshaped, set the coefficients of boundary vertixes. 
	SetBoundaryVCoefficient();

	// Then, solve equation   
	Minimize(this->pType);

	//boundaryObj->GetComponent<CmptGeometry>()->primitive = GenMesh(this->adj_V);

	// Finally, half-edge structure -> triangle mesh, end. 
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	//triMeshObj = TriMesh::New(indice, positions);
	//this->triMesh->Init(indice, positions);
	if (this->is_tex)
	{
		this->triMesh->Update(texCoor);
	}
	else 
	{
		this->triMesh->Update(positions);
	}
	return true;
}

/* 
	Used to get infomation of boundary:
	1. get the length of boundary edges to boundary reshape --> store in mat_nE_boundary_len
	2. get the boundary vertixes. --> store in boundary_V
*/
void Paramaterize::GetBoundaryInfo()
{
	auto boundary_HE = heMesh->Boundaries()[0];
	size_t nE_boundary = boundary_HE.size();
	MatrixXd mat_nE_boundary_len(nE_boundary, 1);
	for (size_t i = 0; i < nE_boundary; i++)
	{
		V* v1 = boundary_HE[i]->Origin();
		V* v2 = boundary_HE[i]->Pair()->Origin();
		mat_nE_boundary_len(i, 0) = pointf3::distance(v1->pos.cast_to<pointf3>(), v2->pos.cast_to<pointf3>());
		boundary_V.push_back(v1);  // the angle of v1 with index j will be mat_nE_boundary_len[j-1];
	}
	this->mat_nE_boundary_len = mat_nE_boundary_len;
}

/*
	Used to set the coefficient of inner vertixes.
	The reason why to split the coefficient of inner vertixes and boundary vertixes is:
	The boundary vertixes should not be reshaped before the time when use Cotan way to get coefficient.
*/
void Paramaterize::SetInnerVCoefficient()
{
	size_t nV = heMesh->NumVertices();

	for (size_t i = 0; i < boundary_V.size(); i++)
	{
		size_t idx = heMesh->Index(boundary_V[i]);
		is_boundary_array[idx] = 1; //  1: boundary v, 0: inner v 
	}

	for (size_t i = 0; i < nV; i++)
	{
		// if not boundary vertix 
		if (is_boundary_array[i] != 1) {
			auto v = heMesh->Vertices()[i];
			size_t idx = heMesh->Index(v);
			vector<V*> adjV = v->AdjVertices(); // get adjcent v

			// 1. Naive way: only use degree
			if (Naive == pType)
			{
				for (size_t j = 0; j < adjV.size(); j++)
				{
					size_t adj_idx = heMesh->Index(adjV[j]);
					coeff.push_back(Eigen::Triplet<double>(idx, adj_idx, -1));
				}
				size_t degree = v->Degree();
				coeff.push_back(Eigen::Triplet<double>(idx, idx, degree));
			}

			// 2. Cotan: more sophisticated way to use tan 
			else if (Cotan == pType)
			{
				MatrixXd adj_coeff(nV, 1);
				adj_coeff.setZero();

				// cout << adjV.size() << endl; 
				for (size_t j = 0; j < adjV.size(); j++)
				{
					V* adj_v = adjV[j];
					V* adj_v1 = adjV[(j + 1) % adjV.size()];
					V* adj_v2;
					if (0 == j)
						adj_v2 = adjV[adjV.size() - 1];
					else
						adj_v2 = adjV[(j - 1) % adjV.size()];

					size_t adj_idx = heMesh->Index(adj_v);
					double dist = pointf3::distance(adj_v->pos.cast_to<pointf3>(), v->pos.cast_to<pointf3>()); // distance between adj_v and v
					auto cos1 = vecf3::cos_theta((adj_v->pos - adj_v1->pos), (v->pos - adj_v1->pos));
					auto cos2 = vecf3::cos_theta((adj_v->pos - adj_v2->pos), (v->pos - adj_v2->pos));
					auto sin1 = sqrt(1 - cos1 * cos1);
					auto sin2 = sqrt(1 - cos2 * cos2);
					//auto sin1 = vecf3::sin_theta((adj_v->pos - v->pos), (adj_v1->pos - adj_v->pos));
					//auto sin2 = vecf3::sin_theta((adj_v->pos - v->pos), (adj_v2->pos - adj_v->pos));
					adj_coeff(adj_idx, 0) = (sin1 / (cos1 + 1) + sin2 / (cos2 + 1)) / dist;
				}
				//adj_coeff(adj_idx, 0) = (cos1 / sin1 + cos2 / sin2 );

				adj_coeff = adj_coeff / adj_coeff.sum();
				for (size_t j = 0; j < adjV.size(); j++)
				{
					size_t adj_idx = heMesh->Index(adjV[j]);
					coeff.push_back(Eigen::Triplet<double>(idx, adj_idx, adj_coeff(adj_idx, 0)));
				}
				coeff.push_back(Eigen::Triplet<double>(idx, idx, -1)); // minus sum 
			}
			mat_bx(idx, 0) = mat_bx(idx, 1) = mat_bx(idx, 2) = 0;
		}
	}

}

/*
	Used to set the coefficient of inner vertixes.
*/
void Paramaterize::SetBoundaryVCoefficient()
{
	size_t nV = heMesh->NumVertices();

	// boundary vertix: coefficient be 1 
	for (size_t i = 0; i < reshaped_boundary_V_pos.size(); i++)
	{
		auto v = get<0>(reshaped_boundary_V_pos[i]);
		size_t idx = heMesh->Index(v);
		size_t degree = v->Degree();

		coeff.push_back(Eigen::Triplet<double>(idx, idx, 1));
		is_boundary_array[idx] = 1;

		// Note: correct idx is important to mat_bx and mat_A
		mat_bx(idx, 0) = get<1>(reshaped_boundary_V_pos[i])[0];
		mat_bx(idx, 1) = get<1>(reshaped_boundary_V_pos[i])[1];
		mat_bx(idx, 2) = get<1>(reshaped_boundary_V_pos[i])[2];
	}
}

/*
	Used to reshape the boundary vertixes. 
*/
void Paramaterize::ReshapeBoundary(Boundary_Shape mode = Circle_boundary)
{
	// First, get the angle of boundary edges according to the length of boundary edges.
	mat_nE_boundary_len = (mat_nE_boundary_len / mat_nE_boundary_len.sum()) * 2 * PI;
	for (size_t i = 1; i < mat_nE_boundary_len.size(); i++)
	{
		mat_nE_boundary_len(i,0) += mat_nE_boundary_len(i - 1,0);
	}
	
	// Second, map boundary vertixes to a circle with radius be 0.5 according to their angles
	if (Circle_boundary == mode)
	{
		for (size_t i = 0; i < boundary_V.size(); i++)
		{
			reshaped_boundary_V_pos.push_back(tuple<V*, pointf3>(
				boundary_V[i], 
				pointf3(2 * cos(mat_nE_boundary_len(i, 0)), 
						2 * sin(mat_nE_boundary_len(i, 0))
						,0)
				)
			);
			/*boundary_V[i]->pos[0] = 0.5* cos(mat_nE_boundary_len(i, 0)); 
			boundary_V[i]->pos[1] = 0.5* sin(mat_nE_boundary_len(i, 0));
			boundary_V[i]->pos[2] = 0;*/
		}
	}
	//map boundary vertixes to a square with width be 1 according to their angles
	else if (Square_boundary == mode)
	{
		for (size_t i = 0; i < boundary_V.size(); i++)
		{
			//boundary_V[i]->pos[2] = 0;
			// right
			if (mat_nE_boundary_len(i, 0) < PI / 4 || mat_nE_boundary_len(i, 0) >= 7 * PI / 4)
			{
				/*boundary_V[i]->pos[0] = 0.5;
				boundary_V[i]->pos[1] = 0.5 * tan(mat_nE_boundary_len(i, 0));*/
				reshaped_boundary_V_pos.push_back(tuple<V*, pointf3>(
					boundary_V[i],
					pointf3(1,
						0.5 * tan(mat_nE_boundary_len(i, 0)) + 0.5
						, 0)
					)
				);
			}
			// up
			else if (mat_nE_boundary_len(i, 0) >= PI / 4 && mat_nE_boundary_len(i, 0) < 3 * PI / 4)
			{
				//boundary_V[i]->pos[0] = 0.5 * tan(PI/2 - mat_nE_boundary_len(i, 0));
				//boundary_V[i]->pos[1] = 0.5;
				reshaped_boundary_V_pos.push_back(tuple<V*, pointf3>(
					boundary_V[i],
					pointf3(0.5 * tan(PI / 2 - mat_nE_boundary_len(i, 0)) + 0.5,
						1
						, 0)
					)
				);
			}
			// left
			else if (mat_nE_boundary_len(i, 0) >= 3*PI / 4 && mat_nE_boundary_len(i, 0) < 5 * PI / 4)
			{
				//boundary_V[i]->pos[0] = -0.5;
				//boundary_V[i]->pos[1] = 0.5 * tan(PI - mat_nE_boundary_len(i, 0));
				reshaped_boundary_V_pos.push_back(tuple<V*, pointf3>(
					boundary_V[i],
					pointf3(0,
						0.5 * tan(PI - mat_nE_boundary_len(i, 0)) + 0.5
						, 0)
					)
				);
			}
			// down 
			else if (mat_nE_boundary_len(i, 0) >= 5*PI / 4 && mat_nE_boundary_len(i, 0) < 7 * PI / 4)
			{
				//boundary_V[i]->pos[0] = 0.5 * tan(mat_nE_boundary_len(i, 0) - 3*PI / 2);
				//boundary_V[i]->pos[1] = -0.5;
				reshaped_boundary_V_pos.push_back(tuple<V*, pointf3>(
					boundary_V[i],
					pointf3(0.5 * tan(mat_nE_boundary_len(i, 0) - 3 * PI / 2) + 0.5,
						 0
						, 0)
					)
				);
			}
			//boundary_V[i]->pos[0] += 0.5;
			//boundary_V[i]->pos[1] += 0.5;
		}
	}
}

/*
	Used to solve the equation and get result. 
*/
void Paramaterize::Minimize(Parameterize_Type pType = Naive)
{
	size_t nV = heMesh->NumVertices();
	mat_A.setFromTriplets(coeff.begin(), coeff.end());
	//cout << MatrixXd(mat_A) << endl; // DEBUG

	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	 
	solver.analyzePattern(mat_A); 
	solver.factorize(mat_A);
	solver.compute(mat_A);
	MatrixXd min_surf_V;
	//cout << mat_bx << endl;  // DEBUG 

	min_surf_V = solver.solve(mat_bx);
	//cout << min_surf_V << endl; // DEBUG 
	// update the position of Vertices in heMesh 
	for (size_t i = 0; i < nV; i++)
	{
		if (! is_tex) { // if need to flaten to a plane
			heMesh->Vertices()[i]->pos[0] = (abs(min_surf_V(i, 0)) > 1e-4 ) ? min_surf_V(i, 0): 0;
			heMesh->Vertices()[i]->pos[1] = (abs(min_surf_V(i, 1)) > 1e-4 ) ? min_surf_V(i, 1) : 0;
			heMesh->Vertices()[i]->pos[2] = (abs(min_surf_V(i, 2)) > 1e-4 ) ? min_surf_V(i, 2) : 0;
		}
		texCoor.push_back(pointf2(min_surf_V(i, 0), min_surf_V(i, 1)));
	}
}
