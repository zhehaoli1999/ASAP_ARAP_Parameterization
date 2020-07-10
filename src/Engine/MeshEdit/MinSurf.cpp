#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>
#include <Engine/Scene/CmptGeometry.h>
#include <Engine/Scene/CmptMaterial.h>
#include <Engine/Material/BSDF_Frostbite.h>
#include <Engine/Scene/SObj.h>

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include <chrono>
#include <ctime>    

using namespace Ubpa;
using namespace std;
using namespace Eigen;

MinSurf::MinSurf(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMeshObj, triMesh);
}

void MinSurf::Clear() {
	heMesh->Clear();
	triMeshObj = nullptr;
	boundaryObj = nullptr;
	triMesh = nullptr;
	boundary_V.clear();
}

Ptr<SObj> MinSurf::GetBoundaryObj()
{
	return this->boundaryObj;
}

bool MinSurf::Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh) {
	Clear();

	if (triMeshObj == nullptr)
		return true;

	//auto triMesh = CastTo<TriMesh>(triMeshObj->GetComponent<CmptGeometry>()->primitive);

	if (!triMesh || triMesh->GetType() == TriMesh::INVALID) 
	{
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}

	this->triMesh = triMesh;

	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<pointf3>();
	}

	// add boundary object 
	this->triMeshObj = triMeshObj;
	string name = "$[MinSurf]boundaryObj";
	for (auto child : triMeshObj->GetChildren()) {
		if (child->name == name)
			boundaryObj = child;
	}
	if (!boundaryObj) {
		boundaryObj = SObj::New(triMeshObj, name);
		boundaryObj->AddComponent<CmptGeometry>(nullptr);
		boundaryObj->AddComponent<CmptMaterial>(BSDF_Frostbite::New(rgbf{ 1.f,0.f,0.f }));
	}
	return true;
}

/* This function is used to get the boundary of the 3D mesh and store the boundary in this->boundary_V */
void MinSurf::GetBoundary()
{
	//vector<V*> boundary_V;
	auto boundary_HE = heMesh->Boundaries()[0];
	size_t nE_boundary = boundary_HE.size();
	for (size_t i = 0; i < nE_boundary; i++)
	{
		V* v1 = boundary_HE[i]->Origin();
		this->boundary_V.push_back(v1);
	}
}

bool MinSurf::Run() {
	if (heMesh->IsEmpty() || !triMeshObj) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}

	// first, get the boundary vertice of the 3d object 
	GetBoundary();

	// Second, show the boundary in engine 
	boundaryObj->GetComponent<CmptGeometry>()->primitive = GenMesh(this->boundary_V);
	
	// Third, build the sparse equation and solve it 
	auto start = std::chrono::system_clock::now();
	Minimize();
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "minimize time:" << elapsed_seconds.count() << "s" << endl;

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
	this->triMesh->Init(indice, positions);

	return true;
}

/* This function is used to generate the mesh for boundary to clearly show boundary in Engine. */
Ptr<TriMesh> MinSurf::GenMesh(const std::vector<V*>& boundary) {
	if ( 0 == boundary.size()) 
	{
		return nullptr;
	}
	
	vector<pointf3> positions;
	positions.reserve(2 * boundary.size());
	vector<unsigned> indices;
	indices.reserve(6 * boundary.size());
	vecf3 off0(0, 0.01f, 0);
	vecf3 off1(0, -0.01f, 0);
	
	positions.push_back(boundary[0]->pos + off0);
	positions.push_back(boundary[0]->pos + off1);

	for (size_t i = 1; i < boundary.size(); i++) {
		positions.push_back(boundary[i]->pos + off0);
		positions.push_back(boundary[i]->pos + off1);

		size_t curIndices[6] = {
			2 * i - 2,2 * i - 1,2 * i,
			2 * i + 1,2 * i,2 * i - 1
		};
		for (auto idx : curIndices)
			indices.push_back(static_cast<unsigned>(idx));
	}
	auto boundaryMesh = TriMesh::New(indices, positions);
	cout << "Gen Mesh" << endl;
	return boundaryMesh;

}

/* This function is used to build the sparse linear equation and solve it. */
void MinSurf::Minimize() {
	// first, get the boundary vertice of the 3d object 
	// Second, show the boundary in engine 
	// Third, build the sparse equation  
	vector<pointf3> boundary_Vpos;
	for (size_t i = 0; i < boundary_V.size(); i++)
	{
		boundary_Vpos.push_back(boundary_V[i]->pos.cast_to<pointf3>());  // store the boundary position in boundary_Vpos. 
	}

	size_t nV = heMesh->NumVertices();
	
	size_t* is_boundary_array = new size_t[nV](); // initialize to be all 0 

	SparseMatrix<double> mat_A(nV, nV);
	mat_A.setZero(); // set zeros initially 
	std::vector<Eigen::Triplet<double> > coeff;
	MatrixXd mat_bx(nV, 3);
	mat_bx.setZero();

	// boundary vertix: coefficient be 1 
	for (size_t i = 0; i < boundary_Vpos.size(); i++)
	{
		size_t idx = heMesh->Index(boundary_V[i]);
		coeff.push_back(Eigen::Triplet<double>(idx, idx, 1));
		is_boundary_array[idx] = 1;

		// Note: correct idx is important to mat_bx and mat_A
		mat_bx(idx, 0) = boundary_Vpos[i][0];
		mat_bx(idx, 1) = boundary_Vpos[i][1];
		mat_bx(idx, 2) = boundary_Vpos[i][2];
	}
	for (size_t i = 0; i < nV; i++)
	{
		// if not boundary vertix 
		if (is_boundary_array[i] != 1) { 
			auto v = heMesh->Vertices()[i];
			size_t idx = heMesh->Index(v);
			vector<V*> adjV = v->AdjVertices(); // get adjcent v
			for (size_t j = 0; j < adjV.size(); j++)
			{
				size_t adj_idx = heMesh->Index(adjV[j]);
				coeff.push_back(Eigen::Triplet<double>(idx, adj_idx, -1));
			}
			size_t degree = v->Degree();
			coeff.push_back(Eigen::Triplet<double>(idx, idx, degree));

			mat_bx(idx, 0) = mat_bx(idx, 1) = mat_bx(idx, 2) = 0;
		}
	}
	mat_A.setFromTriplets(coeff.begin(), coeff.end());

	//cout << MatrixXd(mat_A) << endl; // DEBUG

	//Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > LLT;
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	solver.compute(mat_A);
	MatrixXd min_surf_V;
	min_surf_V = solver.solve(mat_bx); 
	
	//cout << mat_bx << endl; // DEBUG
	//cout << min_surf_V << endl; // DEBUG

	// update the position of Vertices in heMesh 
	for (size_t i = 0; i < nV; i++)
	{
		heMesh->Vertices()[i]->pos[0] = min_surf_V(i, 0);
		heMesh->Vertices()[i]->pos[1] = min_surf_V(i, 1);
		heMesh->Vertices()[i]->pos[2] = min_surf_V(i, 2);
	}

	/*cout << "WARNING::MinSurf::Minimize:" << endl
		<< "\t" << "not implemented" << endl;*/
}
