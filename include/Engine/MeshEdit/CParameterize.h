#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <Engine/Scene/SObj.h>
#include <UGM/UGM>
#include <limits>
#include <map>

namespace Ubpa {
	class TriMesh;

	class CParameterize : public HeapObj {
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
		CParameterize(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh);
	public:
		static const Ptr<CParameterize> New(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh) {
			return Ubpa::New<CParameterize>(triMeshObj, triMesh);
		}
	public:
		void Clear();
		bool Init(Ptr<SObj> triMeshObj, Ptr<TriMesh> triMesh);
		bool Run();
		//void GetBoundaryInfo();

		/* Get cotan of angle of edge v1v2 in triangle v1v2v3 */
		double getCotan(V* v1, V* v2, V* v3);

		/* get the angle of edge i, i+1. idx is V_{(i+2) %3} */
		double getCotan(int tri_idx,  V* v);

		/*
			Congruently Map triangles on 3D mesh to 2D
			return: vector of new vertrixes
		*/
		void CongruentMapping2D();

	protected:
		std::vector<pointf2> texCoor;

		Ptr<TriMesh> triMesh;
		Ptr<SObj> triMeshObj;
		const Ptr<HEMesh<V>> heMesh;

		std::vector<std::map<V*, pointf3>> points2d; // 2d points mapped from 3D
		std::vector<std::map<V*, double>> cot2d; // cotan value of origin 3D mesh

		size_t nV; // num of vertices
		size_t nT; // num of triangles
	};

}
