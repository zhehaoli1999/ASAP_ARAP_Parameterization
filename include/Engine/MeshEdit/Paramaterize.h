#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <Engine/Scene/SObj.h>
#include<Engine/MeshEdit/MinSurf.h>
#include<Engine/MeshEdit/CParameterize.h>
#include <UGM/UGM>
#include <limits>

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCholesky>

#include <tuple>
using namespace Eigen;

namespace Ubpa {
	class TriMesh;
	class MinSurf;

	// mesh boundary == 1
	class Paramaterize : public CParameterize {
	public:
		enum Boundary_Shape{
			Circle_boundary = 0,
			Square_boundary = 1
		};

		enum Parameterize_Type {
			Naive = 0,
			Cotan = 1
		};

	public:
		Paramaterize(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex, int mode, int pType);
		~Paramaterize();
	public:
		static const Ptr<Paramaterize> New(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex, int mode, int pType) {
			return Ubpa::New<Paramaterize>(triMeshObj, triMesh, is_tex, mode, pType);
		}
	public:
		void Clear();
		bool Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh, bool is_tex,  int mode, int pType);
		bool Run();
		void ReshapeBoundary(Boundary_Shape mode);
		void Minimize(Parameterize_Type pType);
		void SetInnerVCoefficient();
		void SetBoundaryVCoefficient();
		void GetBoundaryInfo();

	protected:
		std::vector<pointf2> texCoor;

	private:
		SparseMatrix<double> mat_A; 
		std::vector<Eigen::Triplet<double> > coeff;
		MatrixXd mat_bx;

		std::vector<V*> boundary_V;
		std::vector<std::tuple<V*, pointf3>> reshaped_boundary_V_pos;
		MatrixXd mat_nE_boundary_len; 
		size_t* is_boundary_array;
		Boundary_Shape mode; 
		Parameterize_Type pType;
		bool is_tex;
	};

}
