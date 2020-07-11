#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <Engine/Scene/SObj.h>
#include <UGM/UGM>
#include <limits>

namespace Ubpa {
	class TriMesh;
	class Paramaterize;

	class MinSurf : public HeapObj {
	public:
		struct V;
		struct E;
		struct P;
		struct V : public TVertex<V, E, P> {
			pointf3 pos;
			V* pre;
			float weight;
		};
		struct E : public TEdge<V, E, P> { };
		struct P :public TPolygon<V, E, P> { };
	
	public:
		MinSurf(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh);
	
	public:
		static const Ptr<MinSurf> New(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh) {
			return Ubpa::New<MinSurf>(triMeshObj, triMesh);
		}
		static Ptr<TriMesh> GenMesh(const std::vector<V*>& boundary);
	
	public:
		// clear cache data
		void Clear();
		Ptr<SObj> GetBoundaryObj();
		//Ptr<vector<V*>> GetBoundary();
		void GetBoundary();
		// init cache data (eg. half-edge structure) for Run()
		bool Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh);

		// call it after Init()
		bool Run();
		

	private:
		// kernel part of the algorithm
		void Minimize();

	private:
		friend class Paramaterize;

		Ptr<SObj> triMeshObj;
		Ptr<SObj> boundaryObj;

		std::vector<V*> boundary_V;

		Ptr<TriMesh> triMesh;
		const Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh
	};
}
