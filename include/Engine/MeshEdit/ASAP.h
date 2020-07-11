#ifndef ASAP_H
#define ASAP_H

#include <Engine/MeshEdit/CParameterize.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

using namespace Eigen; 
using namespace std;

namespace Ubpa {

	class ASAP : public CParameterize
	{
	public:
		ASAP(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex);
	public:
		static const Ptr<ASAP> New(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex) {
			return Ubpa::New<ASAP>(triMeshObj, triMesh, is_tex);
		}
	public:
		bool Run();
		void Clear();
		bool Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex);
		
		/* 
			Used to set the coefficient for sparse linear system and prefactorize
			Param: size_t idx2, pointf3 pos2: Constrain 2 anchor vertixes for solving equation
		*/
		void setCoefficient_u(/*size_t idx1, pointf3 pos1,*/ size_t idx2, pointf3 pos2, size_t idx3, pointf3 pos3);
		void setCoefficient_Lt();
		void setCoeff_ux(size_t row, size_t col, double coeffcient);
		void setCoeff_uy(size_t row, size_t col, double coeffcient);
		void setCoeff_a(size_t row, size_t col, double coeffcient);
		void setCoeff_b(size_t row, size_t col, double coeffcient);
	
		/* Used to solve sparse linear system and set the new coordinates  */
		void solve();

	private:
		bool is_tex; 
		SparseMatrix<double> mat_A;
		std::vector<Eigen::Triplet<double> > coeff;
		MatrixXd mat_b;
		MatrixXd solution;
		SimplicialCholesky<SparseMatrix<double> > solver;
		//Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	};
	

}

#endif // ASAP_H