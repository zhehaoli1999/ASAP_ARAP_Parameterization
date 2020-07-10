#include<Engine/MeshEdit/CParameterize.h>
#include <Engine/Primitive/TriMesh.h>
#include <Engine/Scene/CmptGeometry.h>
#include <Engine/Scene/CmptMaterial.h>
#include <Engine/Material/BSDF_Frostbite.h>
#include <Engine/Scene/SObj.h>

using namespace Ubpa;
using namespace std;

CParameterize::CParameterize(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMeshObj, triMesh);
}

bool CParameterize::Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh)
{
	Clear();

	if (triMeshObj == nullptr)
		return true;

	if (!triMesh || triMesh->GetType() == TriMesh::INVALID)
	{
		printf("ERROR::Paramater::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}
	
	this->triMesh = triMesh;

	// init half-edge structure
	this->nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	this->nT = triMesh->GetTriangles().size();
	triangles.reserve(nT);
	for (auto triangle : triMesh->GetTriangles()) {
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	}
	
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::Paramater::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<pointf3>();
	}

	// add object 
	this->triMeshObj = triMeshObj;
	return true;
}

void CParameterize::Clear()
{
	//
	heMesh->Clear();
	triMesh = nullptr;
	triMeshObj = nullptr;
	texCoor.clear();
}

bool CParameterize::Run()
{
	CongruentMapping2D();
	return true;
}

void CParameterize::CongruentMapping2D()
{
	this->points2d.clear();
	this->cot2d.clear();
	for (auto triangle : this->heMesh->Polygons())
	{
		if (triangle != nullptr) {
			map<V*, pointf3> map_v_2d;

			auto v0 = triangle->BoundaryVertice()[0];
			auto v1 = triangle->BoundaryVertice()[1];
			auto v2 = triangle->BoundaryVertice()[2];

			map_v_2d[v0] = pointf3(0, 0, 0); // set new v0 = (0,0,0)

			// get 2d coordinates
			double dist01 = pointf3::distance(v0->pos.cast_to<pointf3>(), v1->pos.cast_to<pointf3>());
			double dist02 = pointf3::distance(v0->pos.cast_to<pointf3>(), v2->pos.cast_to<pointf3>());
			double cos01_02 = vecf3::cos_theta((v1->pos - v0->pos), (v2->pos - v0->pos)); // angle of edge 01 and 02
			double sin01_02 = sqrt(1 - cos01_02 * cos01_02);

			map_v_2d[v1] = pointf3(dist01, 0, 0);   // set new v1 = (dist01,0,0)
			map_v_2d[v2] = pointf3(dist02 * cos01_02, dist02 * sin01_02, 0);  // set new v2 = (d02 * cos, d02 * sin, 0)

			this->points2d.push_back(map_v_2d); // triangle in points2d has the same sequence with heMesh->triangles()

			// get cot
			map<V*, double> cotan_2d;
			for (size_t i = 0; i < 3; i++)
			{
				double cos = vecf3::cos_theta((triangle->BoundaryVertice()[i]->pos - triangle->BoundaryVertice()[(i + 2) % 3]->pos),
					(triangle->BoundaryVertice()[(i + 1) % 3]->pos - triangle->BoundaryVertice()[(i + 2) % 3]->pos));
				double sin = sqrt(1 - cos * cos);
				cotan_2d[triangle->BoundaryVertice()[(i + 2) % 3]] = (cos / sin); // store the angle of edge i, i+1. idx is V_{(i+2) %3}
			}
			assert(cotan_2d.size() == 3);
			this->cot2d.push_back(cotan_2d);
		}
	}
	
	assert(this->points2d.size() == nT);
	assert(this->cot2d.size() == nT);

}

double CParameterize::getCotan(int tri_idx, V* v)
{
	assert(this->cot2d.size() == nT);
	try
	{
		auto cotan = cot2d[tri_idx].at(v);
		return cotan;
	}
	catch (const std::exception & e)
	{
		cout << "[Error] Cannot find V in map. ";
	}
}
double CParameterize::getCotan(V* v0, V* v1, V* v2)
{
	double cos = vecf3::cos_theta((v0->pos - v2->pos), (v1->pos - v2->pos));
	double sin = sqrt(1 - cos * cos);
	return cos / sin;
}

//map<V*,pointf3> CParameterize::getMappedV(size_t idx)
//{
//	//vector<pointf3> v_2d;
//	//v_2d.push_back(pointf3(0, 0, 0)); // set new v0 = (0,0,0)
//
//	//double dist01 = pointf3::distance(v0->pos.cast_to<pointf3>(), v1->pos.cast_to<pointf3>());
//	//double dist02 = pointf3::distance(v0->pos.cast_to<pointf3>(), v2->pos.cast_to<pointf3>());
//	//double cos01_02 = vecf3::cos_theta((v1->pos - v0->pos), (v2->pos - v0->pos)); // angle of edge 01 and 02
//	//double sin01_02 = sqrt(1 - cos01_02 * cos01_02);
//
//	//v_2d.push_back(pointf3(dist01, 0, 0));   // set new v1 = (dist01,0,0)
//	//v_2d.push_back(pointf3(dist02 * cos01_02, dist02 * sin01_02, 0));  // set new v2 = (d02 * cos, d02 * sin, 0)
//	map<V*, pointf3> v_2d = points2d[idx];
//	return v_2d;
//}
