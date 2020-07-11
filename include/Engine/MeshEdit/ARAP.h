#pragma once
#include <Engine/MeshEdit/CParameterize.h>
#include <Engine/MeshEdit/Paramaterize.h>
#include <Engine/MeshEdit/ASAP.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include<map>

using namespace Eigen;
using namespace std;

namespace Ubpa {
	class ASAP;

	class ARAP : public Paramaterize/* public ASAP*/
	{
	public:
		ARAP(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex);
	public:
		static const Ptr<ARAP> New(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex) {
			return Ubpa::New<ARAP>(triMeshObj, triMesh, is_tex);
		}
	public:
		bool RunARAP(int iter_n, double error, int debug);
		void localupdate();
		double globalupdate();
		void Clear();
		bool Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh);

		/*
			Used to set the coefficient for sparse linear system and prefactorize
			Param: size_t idx1, pointf3 pos1, size_t idx2, pointf3 pos2: Constrain 2 anchor vertixes for solving equation
		*/
		void setCoefficientA(size_t idx1, size_t idx2);

		void setb(size_t idx1, pointf2 pos1, size_t idx2, pointf2 pos2);

		/* Used to solve sparse linear system and set the new coordinates */

		/* Used to get Lt (t = 1, 2,..., nT), and store Lt into Lt_array according to the sequence of heMesh->triangle() */
		void getLt();
		
		/* Used to get newest u per iteration */
		void getTrianglePoints();
	
	public:
		bool is_tex;
		bool is_debug;
	private:
		SparseMatrix<double> ARAP_mat_A;
		std::vector<Eigen::Triplet<double> > ARAP_coeff;
		MatrixXd ARAP_mat_b;
		MatrixXd ARAP_solution;
		SimplicialCholesky<SparseMatrix<double> > ARAP_solver;
		//Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > ARAP_solver;
		std::vector<Matrix2d> Lt_array;
		std::vector<std::vector<V*>> triangle_points;
	
		size_t anchor_v1_idx;
		pointf2 anchor_pos1;
		size_t anchor_v2_idx;
		pointf2 anchor_pos2;
	};


}