#pragma once

#include <Basic/HeapObj.h>
#include <UGM/UGM>
#include <Engine/MeshEdit/Simulate.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace Ubpa 
{
	class FastSimulate : public Simulate
	{
	public:
		FastSimulate(const std::vector<pointf3>& plist, const std::vector<unsigned>& elist);

	public:
		static const Ptr<FastSimulate> New(const std::vector<pointf3>& plist,
			const std::vector<unsigned>& elist) {
			return Ubpa::New<FastSimulate>(plist, elist);
		}
	public:
		void Clear();

		void Init();

		bool Run();

		void FastSimulateOnce(int iter);
		
		void getD(); 
		void getA();
		
		void InitXf(); // Used to initialize Xf, Xn, Xn_
		
		void InitConstant(); // Used to initialize matrix of gravity, mass and stiffness

		void LocalUpdate(); // fix x, and update D
		
		void GlobalUpdate(); // fix D, and solve x;

	private:
		SparseMatrix<double> M; // diagonal matrix of mass 
		SparseMatrix<double> A; // Matrix of string incidence vector
		MatrixXd J;  
		SparseMatrix<double> L;
		MatrixXd D; // matrix of di, i in [1,..., nE]
		SparseMatrix<double> mat_stiff; // diagonal matrix of stiffness coedfficient 
		 
		MatrixXd y;  
		MatrixXd Gravity; // gravity matrix
		SimplicialCholesky<SparseMatrix<double> > solver;

		MatrixXd newX; 
		MatrixXd Xn;  // Xn
		MatrixXd Xn_; // Xn-1 
		MatrixXd Xf;  // X not fix

		size_t n_notfix; // number of not fix points 
		SparseMatrix<double> K; // Xf = K * X
		MatrixXd b;  //  b =  X - Xf
	};
};