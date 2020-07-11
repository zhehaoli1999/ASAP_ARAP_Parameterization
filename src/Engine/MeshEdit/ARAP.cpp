#include <Engine/MeshEdit/CParameterize.h>
#include <Engine/MeshEdit/ARAP.h>
#include <Engine/Primitive/TriMesh.h>
#include <Engine/Scene/CmptGeometry.h>
#include <Engine/Scene/CmptMaterial.h>
#include <Engine/Material/BSDF_Frostbite.h>
#include <Eigen/SVD>

using namespace Ubpa;
using namespace Eigen;

ARAP::ARAP(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex)
	: Paramaterize(triMeshObj, triMesh, is_tex, 0, 1)/* ASAP(triMeshObj, triMesh, true)*/ // Circle boundary and Naive mode
{
	this->is_tex = is_tex;
	Init(triMeshObj, triMesh);
}

void ARAP::Clear()
{
	ARAP_coeff.clear();
	Lt_array.clear();
	triangle_points.clear();
}

bool ARAP::Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh)
{
	ARAP_mat_A = SparseMatrix<double>(nV , nV);
	ARAP_mat_A.setZero();
	ARAP_mat_b = MatrixXd(nV, 2);
	ARAP_mat_b.setZero();
	//Initialize 1.map 3D triangles to 2d 
	CongruentMapping2D();

	// Initialize 2. Select anchor points
	auto triangle = heMesh->Polygons().back();
	auto v1 = triangle->BoundaryVertice()[0];
	anchor_v1_idx = heMesh->Index(v1);
	auto v2 = triangle->BoundaryVertice()[1]; 
	anchor_v2_idx = heMesh->Index(v2);
	size_t tri_idx = heMesh->Index(triangle);
	//anchor_pos1 = pointf2(points2d[tri_idx][v1][0], points2d[tri_idx][v1][1]);
	//anchor_pos2 = pointf2(points2d[tri_idx][v2][0], points2d[tri_idx][v2][1]);
	anchor_pos2 = pointf2(1, 1);
	// Initialize 3. set coefficients of A 
	setCoefficientA(anchor_v1_idx, anchor_v2_idx);

	//Initialize 4. get initial parameterization 
	this->Run();

	return true;
}

bool ARAP::RunARAP(int iter_n, double error_threshold = 0.01, int debug=5)
{
	// set debug mode
	if (debug < 50) this->is_debug = false;
	else this->is_debug = true;

	// do local/global iteration
	cout <<"iter:" << iter_n << endl;

	for (int i = 0; i < iter_n; i++)
	{
		localupdate();
		double max_error = globalupdate();
		cout << "iter:" << i << " error: " << max_error << endl;
		if (max_error < error_threshold) break;
	}

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
		this->triMesh->Update(texCoor);
	}
	else
	{
		this->triMesh->Update(positions);
	}
	return true;
}

void ARAP::localupdate()
{
	getTrianglePoints();  // get the newest iteration result from trangle_points
	getLt();
}

double ARAP::globalupdate()
{
	// first, set b
	setb(anchor_v1_idx, anchor_pos1, anchor_v2_idx, anchor_pos2);
	double max_error = -1 ; 
	// Then solve the equation
	ARAP_solution = ARAP_solver.solve(ARAP_mat_A.transpose() * ARAP_mat_b);

	if (is_debug) {
		//cout << MatrixXd(ARAP_mat_A) << endl << endl; // DEBUG
		//cout << "b:" << endl << ARAP_mat_b << endl << endl;
		//cout << "solu:" << endl << ARAP_solution << endl << endl;
	}
	

	// End, update textCoordinate
	this->texCoor.clear();
	for (size_t i = 0; i < nV; i++)
	{
		double dist = pointf3::distance(pointf3(ARAP_solution(i, 0), ARAP_solution(i, 1), 0.0), heMesh->Vertices()[i]->pos.cast_to<pointf3>());
		if (max_error < 0) 
		{
			max_error = dist;
		}
		else if(max_error < dist)
		{
			max_error = dist;
		}
		
		heMesh->Vertices()[i]->pos[0] = ARAP_solution(i, 0);
		heMesh->Vertices()[i]->pos[1] = ARAP_solution(i, 1);
		heMesh->Vertices()[i]->pos[2] = 0.0;
		// update tex coordinates
		texCoor.push_back(pointf2(ARAP_solution(i, 0), ARAP_solution(i, 1)));
	}

	return max_error;
}

/* Get newest u per iteration */
void ARAP::getTrianglePoints()
{
	triangle_points.clear();	// remember to clear
	for (auto triangle : heMesh->Polygons())
	{
		if(triangle != nullptr)
		{
			this->triangle_points.push_back(triangle->BoundaryVertice());
		}
	}
}

void ARAP::getLt()
{
	assert(triangle_points.size() == nT);
	Lt_array.clear();

	for (size_t t = 0; t < nT; t++)
	{
		auto vec_u = triangle_points[t];
		auto mapped_u = points2d[t];
		Matrix2d St; 
		St.setZero();

		// get St
		for (int i = 0; i < 3; i++)
		{
			V* u0 = vec_u[i];
			V* u1 = vec_u[(i+1)%3]; // % ?
			MatrixXd delta_u(2, 1);
			delta_u <<
				u0->pos[0] - u1->pos[0], u0->pos[1] - u1->pos[1];
			
			pointf3 x0 = mapped_u[u0];
			pointf3 x1 = mapped_u[u1];
			MatrixXd delta_x(2, 1);
			delta_x <<
				x0[0] - x1[0], x0[1] - x1[1];

			double cot = getCotan(t, vec_u[(i + 2) % 3]);

			St += cot * delta_u * delta_x.transpose();
		}

		// Do SVD Composition on St
		JacobiSVD<MatrixXd> svd(St, ComputeThinU | ComputeThinV);

		Matrix2d Lt = svd.matrixU() * svd.matrixV().transpose(); // Lt = U * V^T
		/*cout << svd.matrixV().transpose() << endl << endl;*/ // DEBUG
		//if (Lt.determinant() < 0 )
		//{
		//	if (is_debug) {
		//		//cout << "before" << Lt.determinant() << endl;
		//	}
		//	Matrix2d newV;
		//	newV <<
		//		svd.matrixV().transpose()(0, 0), svd.matrixV().transpose()(0, 1),
		//		-(svd.matrixV().transpose()(1, 0)), -(svd.matrixV().transpose()(1, 1));
		//	if (is_debug) {
		//		//cout << newV << endl;
		//	}
		//	Lt = svd.matrixU() * newV;
		//	assert(Lt.determinant() > 0);
		//}

		if (is_debug)
		{
			cout << Lt.determinant() << endl;
			//cout << svd.singularValues() << endl << endl;
		}

		//cout << Lt << endl; // DEBUG
		
		//Lt << 1, 0, 0, 1; // DEBUG
 
		Lt_array.push_back(Lt);
		
	}
	assert(Lt_array.size() == nT);
}

void ARAP::setb(size_t idx1, pointf2 pos1, size_t idx2, pointf2 pos2)
{
	ARAP_mat_b.setZero();

	// two anchor points
	//ARAP_mat_b(idx1, 0) = pos1[0];  // x
	//ARAP_mat_b(idx1, 1) = pos1[1]; // y

	ARAP_mat_b(idx2, 0) = pos2[0];  // x
	ARAP_mat_b(idx2, 1) = pos2[1]; // y

	// Other non-anchor points
	for (size_t i = 0; i < nV; i++)
	{
		if (/*i != idx1 && */i != idx2) {
			auto v = heMesh->Vertices()[i];
			MatrixXd b(2, 1);
			b.setZero();
			//cout << "v idx: " << heMesh->Index(v) << endl<<endl;
			// traverse adjecent vertrices
			for (auto adj_v : v->AdjVertices())
			{
				size_t adj_idx = heMesh->Index(adj_v);  // get index of adj_v
				// get adjacent triangles
				auto e = v->EdgeWith(adj_v);
				auto he1 = e->HalfEdge();				
				auto he2 = e->HalfEdge()->Pair();

				// get Lt and cot
				auto triangle1 = he1->Polygon();
				if (triangle1 != nullptr)
				{
					auto tri_v1 = he1->Next()->End();
					assert(heMesh->Index(tri_v1) != heMesh->Index(v)
						&& heMesh->Index(tri_v1) != heMesh->Index(adj_v));

					size_t tri_idx = heMesh->Index(triangle1);
					double cot1 = getCotan(tri_idx, tri_v1);
					map<V*, pointf3> mapped_v = this->points2d[tri_idx]; // congruent mapping of triangle1 
					MatrixXd Lt = Lt_array[tri_idx];
					MatrixXd delta_x(2, 1);
					delta_x <<
						mapped_v[v][0] - mapped_v[adj_v][0],
						mapped_v[v][1] - mapped_v[adj_v][1];
					/*cout << "cot1:" << cot1 << endl;
					cout << "mapped_v[v]:" << mapped_v[v][0] <<"," << mapped_v[v][1] << endl;
					cout << "mapped_[adj_v]:" << mapped_v[adj_v][0] << "," << mapped_v[adj_v][1] << endl;
					cout << "cot1 * Lt * delta_x:" << cot1 * Lt * delta_x << endl;
					cout << endl;*/
					
					b += cot1 * Lt * delta_x;
				}

				auto triangle2 = he2->Polygon();
				if (triangle2 != nullptr)
				{
					auto tri_v2 = he2->Next()->End();
					assert(heMesh->Index(tri_v2) != heMesh->Index(v)
						&& heMesh->Index(tri_v2) != heMesh->Index(adj_v));

					size_t tri_idx = heMesh->Index(triangle2);
					double cot2 = getCotan(tri_idx, tri_v2);
					map<V*, pointf3> mapped_v = this->points2d[tri_idx]; // congruent mapping of triangle1 
					MatrixXd Lt = Lt_array[tri_idx];
					MatrixXd delta_x(2, 1);
					delta_x <<
						mapped_v[v][0] - mapped_v[adj_v][0],
						mapped_v[v][1] - mapped_v[adj_v][1];

					/*cout << "cot2:" << cot2 << endl;
					cout << "mapped_v[v]:" << mapped_v[v][0] << "," << mapped_v[v][1] << endl;
					cout << "mapped_v[adj_v]:" << mapped_v[adj_v][0] << "," << mapped_v[adj_v][1] << endl;
					cout << "cot2 * Lt * delta_x:" << cot2 * Lt * delta_x << endl;
					cout << endl;*/
					b += cot2 * Lt * delta_x;
				}
			}
			ARAP_mat_b(i, 0) = b(0); // x
			ARAP_mat_b(i, 1) = b(1); // y 
		}
	}
}


void ARAP::setCoefficientA(size_t idx1, size_t idx2)
{
	// Two Anchor points 
	//ARAP_coeff.push_back(Eigen::Triplet<double>(idx1, idx1, 1));
	ARAP_coeff.push_back(Eigen::Triplet<double>(idx2, idx2, 1));

	// Other non-anchor points
	for (size_t i = 0; i < nV; i++)
	{
		if (/*i != idx1 &&*/ i != idx2) {
			auto v = heMesh->Vertices()[i];
			double cotan_sum = 0.0;

			// traverse adjacent vertrices
			for (auto adj_v : v->AdjVertices())
			{
				size_t adj_idx = heMesh->Index(adj_v);

				// To get cotan 
				auto e = v->EdgeWith(adj_v);
				auto he1 = e->HalfEdge();

				double cot1 = 0.0;
				if (he1->Polygon() != nullptr)
				{
					auto tri1_idx = heMesh->Index(he1->Polygon()); // get index of adjacent triangle
					try 
					{
						auto tri_v1 = he1->Next()->End(); // get vertix of adjacent triangle 
						assert(heMesh->Index(tri_v1) != heMesh->Index(v)
							&& heMesh->Index(tri_v1) != heMesh->Index(adj_v));

						cot1 = getCotan(tri1_idx, tri_v1);
					}
					catch (const std::exception& e)
					{
						cout << "cannot find cot 1 in map. " << endl;
					}
				}
				auto he2 = e->HalfEdge()->Pair();

				double cot2 = 0.0;
				if (he2->Polygon() != nullptr)
				{
					auto tri2_idx = heMesh->Index(he2->Polygon());
					try
					{
						auto tri_v2 = he2->Next()->End();
						assert(heMesh->Index(tri_v2) != heMesh->Index(v)
							&& heMesh->Index(tri_v2) != heMesh->Index(adj_v));
						cot2 = getCotan(tri2_idx, tri_v2);
					}
					catch (const std::exception & e)
					{
						cout << "cannot find cot 2 in map. " << endl;
					}
				}

				ARAP_coeff.push_back(Eigen::Triplet<double>(i, adj_idx, -(cot1 + cot2)));

				cotan_sum += (cot1 + cot2);
			}

			ARAP_coeff.push_back(Eigen::Triplet<double>(i, i, cotan_sum));
		}
	}
	ARAP_mat_A.setFromTriplets(ARAP_coeff.begin(), ARAP_coeff.end());
	// pre-computation
	ARAP_mat_A.makeCompressed();
	ARAP_solver.compute(ARAP_mat_A.transpose() * ARAP_mat_A);
}