#include <Engine/MeshEdit/CParameterize.h>
#include <Engine/MeshEdit/ASAP.h>
#include <Engine/Primitive/TriMesh.h>
#include <Engine/Scene/CmptGeometry.h>
#include <Engine/Scene/CmptMaterial.h>
#include <Engine/Material/BSDF_Frostbite.h>

using namespace Ubpa;
using namespace Eigen; 

ASAP::ASAP(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex)
	: CParameterize(triMeshObj, triMesh)
{
	Init(triMeshObj, triMesh,is_tex);
}

void ASAP::Clear()
{
	coeff.clear();
}

bool ASAP::Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex)
{
	Clear();
	this->is_tex = is_tex;

	mat_A = SparseMatrix<double>(2*(nV +nT), 2*(nV + nT));
	mat_A.setZero();
	mat_b = MatrixXd(2*(nV+nT), 1);
	mat_b.setZero();
	return true;
}

bool ASAP::Run()
{
	// step 1. mapping from 3D to 2D
	CongruentMapping2D();
	// step 2. Choose two boundary points, to make their distance as far as possible, one idx is [0], the other is HE_boundaries.size() / 2. 
	auto HE_boundaries = heMesh->Boundaries()[0];
	auto he1 = HE_boundaries[0];
	auto he2 = HE_boundaries[HE_boundaries.size() / 2];
	auto v1 = he1->Origin();
	auto v2 = he2->Origin();

	// step 3. set one anchor point (0,0), the other (1,1), and set coefficient of Ux, Uy
	setCoefficient_u(heMesh->Index(v1), pointf3(0,0,0),
			heMesh->Index(v2), pointf3(1, 1, 0));
	// step 4. set the coefficient of Lt
	setCoefficient_Lt();
	// step 5. solve equation.
	solve();

	// Finally, half-edge structure -> triangle mesh, end. 
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nT);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	if (this->is_tex)
	{
		assert(texCoor.size() == nV);
		this->triMesh->Update(texCoor);
	}
	else
	{
		this->triMesh->Update(positions);
	}
	return true;
}

void ASAP::setCoeff_ux(size_t idx, size_t adj_idx, double coeffcient)
{
	coeff.push_back(Eigen::Triplet<double>(idx, adj_idx, coeffcient));
}

void ASAP::setCoeff_uy(size_t idx, size_t adj_idx, double coeffcient)
{
	coeff.push_back(Eigen::Triplet<double>(nV + idx, adj_idx, coeffcient));
}

void ASAP::setCoeff_a(size_t idx, size_t adj_idx, double coeffcient)
{
	coeff.push_back(Eigen::Triplet<double>(2* nV + idx, adj_idx, coeffcient));
}

void ASAP::setCoeff_b(size_t idx, size_t adj_idx, double coeffcient)
{
	coeff.push_back(Eigen::Triplet<double>(2 * nV + nT + idx, adj_idx, coeffcient));
}

void ASAP::setCoefficient_u(/*size_t idx1, pointf3 pos1,*/ size_t idx2, pointf3 pos2, size_t idx3, pointf3 pos3)
{
	// Two Anchor points 
	coeff.push_back(Eigen::Triplet<double>(idx2, idx2, 1));
	coeff.push_back(Eigen::Triplet<double>(nV + idx2, nV + idx2, 1));
	mat_b(idx2) = pos2[0];  // x
	mat_b(nV + idx2) = pos2[1]; // y

	coeff.push_back(Eigen::Triplet<double>(idx3, idx3, 1));
	coeff.push_back(Eigen::Triplet<double>(nV + idx3, nV + idx3, 1));
	mat_b(idx3) = pos3[0];  // x
	mat_b(nV + idx3) = pos3[1]; // y

	// Other non-anchor points
	for (size_t i = 0; i < nV; i++)
	{
		if (/*i != idx1 &&*/ i != idx2  && i !=idx3 ) {
			auto v = heMesh->Vertices()[i];
			double cotan_sum = 0.0;

			for (auto adj_v : v->AdjVertices())
			{
				size_t adj_idx = heMesh->Index(adj_v);

				auto e = v->EdgeWith(adj_v);
				auto he1 = e->HalfEdge();
				auto he2 = e->HalfEdge()->Pair();
				
				double cot1 = 0.0;
				double cot2 = 0.0;

				// to get congruently mapped x
				auto triangle1 = he1->Polygon();
				if (triangle1 != nullptr)
				{
					auto tri_v1 = he1->Next()->End(); // get vertix of adjacent triangle 
					assert(heMesh->Index(tri_v1) != heMesh->Index(v)
						&& heMesh->Index(tri_v1) != heMesh->Index(adj_v));

					size_t tri1_idx = heMesh->Index(triangle1);
					cot1 = getCotan(tri1_idx, tri_v1);

					map<V*, pointf3> mapped_v1 = this->points2d[tri1_idx]; // congruent mapping of triangle1 
					// x
					setCoeff_ux(i, 2 * nV + tri1_idx, -cot1 * (mapped_v1[v][0] - mapped_v1[adj_v][0])); // a_ij
					setCoeff_ux(i, 2 * nV + nT + tri1_idx, -cot1 * (mapped_v1[v][1] - mapped_v1[adj_v][1])); // b_ij
				
					setCoeff_uy(i, 2 * nV + nT + tri1_idx, cot1 * (mapped_v1[v][0] - mapped_v1[adj_v][0])); // b_ji
					setCoeff_uy(i, 2 * nV + tri1_idx, -cot1 * (mapped_v1[v][1] - mapped_v1[adj_v][1])); // a_ji
				}
				auto triangle2 = he2->Polygon();

				if (triangle2 != nullptr)
				{
					auto tri_v2 = he2->Next()->End();
					assert(heMesh->Index(tri_v2) != heMesh->Index(v)
						&& heMesh->Index(tri_v2) != heMesh->Index(adj_v));

					size_t tri2_idx = heMesh->Index(triangle2);
					cot2 = getCotan(tri2_idx, tri_v2);
					map<V*, pointf3> mapped_v2 = this->points2d[tri2_idx]; // congruent mapping of triangle2
					setCoeff_ux(i, 2 * nV + tri2_idx, -cot2 * (mapped_v2[v][0] - mapped_v2[adj_v][0])); // a_ij
					setCoeff_ux(i, 2 * nV + nT + tri2_idx, -cot2 * (mapped_v2[v][1] - mapped_v2[adj_v][1])); // b_ij

					setCoeff_uy(i, 2 * nV + nT + tri2_idx, cot2 * (mapped_v2[v][0] - mapped_v2[adj_v][0])); // b_ji
					setCoeff_uy(i, 2 * nV + tri2_idx, -cot2 * (mapped_v2[v][1] - mapped_v2[adj_v][1])); // a_ji
				}
				cotan_sum += (cot1 + cot2);
				setCoeff_ux(i, adj_idx, -(cot1 + cot2));
 				setCoeff_uy(i, nV + adj_idx, -(cot1 + cot2));
			}
			setCoeff_ux(i, i, cotan_sum);
			setCoeff_uy(i, nV + i, cotan_sum);
		}
	}
}

void ASAP::setCoefficient_Lt()
{
	for (auto triangle : heMesh->Polygons())
	{
		if (triangle != nullptr)
		{
			auto tri_idx = heMesh->Index(triangle);
			auto vec_V = triangle->BoundaryVertice();
			
			size_t v0_idx = heMesh->Index(vec_V[0]);
			size_t v1_idx = heMesh->Index(vec_V[1]);
			size_t v2_idx = heMesh->Index(vec_V[2]);

			assert(vec_V.size() == 3);
			double cot0 = getCotan(tri_idx, vec_V[2]);
			double cot1 = getCotan(tri_idx, vec_V[0]);
			double cot2 = getCotan(tri_idx, vec_V[1]);
			map<V*, pointf3> mapped_v = this->points2d[tri_idx];
			auto map_v0 = mapped_v[vec_V[0]];
			auto map_v1 = mapped_v[vec_V[1]];
			auto map_v2 = mapped_v[vec_V[2]];
			double deltaX0 = map_v0[0] - map_v1[0];
			double deltaX1 = map_v1[0] - map_v2[0];
			double deltaX2 = map_v2[0] - map_v0[0];
			
			double deltaY0 = map_v0[1] - map_v1[1];
			double deltaY1 = map_v1[1] - map_v2[1];
			double deltaY2 = map_v2[1] - map_v0[1];
			
			setCoeff_a(tri_idx, 2*nV+tri_idx, - (cot0 * (deltaX0* deltaX0 + deltaY0 * deltaY0)
												+ cot1 * (deltaX1 * deltaX1 + deltaY1 * deltaY1)
												+ cot2 * (deltaX2 * deltaX2 + deltaY2 * deltaY2)));
			// x
			setCoeff_a(tri_idx, v0_idx, cot0 * deltaX0 - cot2 * deltaX2);
			setCoeff_a(tri_idx, v1_idx, cot1 * deltaX1 - cot0 * deltaX0);
			setCoeff_a(tri_idx, v2_idx, cot2 * deltaX2 - cot1 * deltaX1);
			//y 
			setCoeff_a(tri_idx, nV+v0_idx, cot0 * deltaY0 - cot2 * deltaY2);
			setCoeff_a(tri_idx, nV+v1_idx, cot1 * deltaY1 - cot0 * deltaY0);
			setCoeff_a(tri_idx, nV+v2_idx, cot2 * deltaY2 - cot1 * deltaY1);


			setCoeff_b(tri_idx, 2 * nV + nT + tri_idx, -(cot0 * (deltaX0 * deltaX0 + deltaY0 * deltaY0)
											+ cot1 * (deltaX1 * deltaX1 + deltaY1 * deltaY1)
											+ cot2 * (deltaX2 * deltaX2 + deltaY2 * deltaY2)));
			// x
			setCoeff_b(tri_idx, v0_idx, cot0 * deltaY0 - cot2 * deltaY2);
			setCoeff_b(tri_idx, v1_idx, cot1 * deltaY1 - cot0 * deltaY0);
			setCoeff_b(tri_idx, v2_idx, cot2 * deltaY2 - cot1 * deltaY1);
			//// y
			setCoeff_b(tri_idx, nV + v0_idx, -(cot0 * deltaX0 - cot2 * deltaX2));
			setCoeff_b(tri_idx, nV + v1_idx, -(cot1 * deltaX1 - cot0 * deltaX0));
			setCoeff_b(tri_idx, nV + v2_idx, -(cot2 * deltaX2 - cot1 * deltaX1));
		}
	}
}

void ASAP::solve()
{
	// pre-factorize 
	mat_A.setFromTriplets(coeff.begin(), coeff.end());
	mat_A.makeCompressed();
	solver.compute(mat_A.transpose()*mat_A);
	//cout << MatrixXd(mat_A) << endl; // DEBUG

	solution = solver.solve(mat_A.transpose()*mat_b);
	//cout <<"b:" << endl << mat_b << endl << endl; // DEBUG
	//cout << "solu:" <<endl << solution << endl<<endl; // DEBUG

	for (size_t i = 0; i < nV; i++)
	{
		if (!is_tex) 
		{
			heMesh->Vertices()[i]->pos[0] = solution(i);
			heMesh->Vertices()[i]->pos[1] = solution(nV + i);
			heMesh->Vertices()[i]->pos[2] = 0.0;
		}
		texCoor.push_back(pointf2(solution(i), solution(nV + i)));
	}
}