#pragma once

#include <Basic/HeapObj.h>
//#include <Engine/Primitive/MassSpring.h>
#include <UGM/UGM>

namespace Ubpa {
	class Simulate : public HeapObj {
	public:
		Simulate(const std::vector<pointf3>& plist, const std::vector<unsigned>& elist);
	public:
		static const Ptr<Simulate> New(const std::vector<pointf3>& plist,
			const std::vector<unsigned> &elist) {
			return Ubpa::New<Simulate>(plist, elist);
		}
	public:
		// clear cache data
		void Clear();

		// init cache data (eg. half-edge structure) for Run()
		bool Init();

		// call it after Init()
		bool Run();
		
		const std::vector<pointf3>& GetPositions() const { return positions; };

		const float GetStiff() { return stiff; };

		void SetStiff(float k) { stiff = k; Init();};
		
		const float GetTimeStep() { return h; };
		
		void SetTimeStep(float k) { h = k; Init();};
		
		std::vector<unsigned>& GetFix() { return this->fixed_id; };
		
		void SetFix(const std::vector<unsigned>& f) { this->fixed_id = f; Init();};
		
		const std::vector<vecf3>& GetVelocity() { return velocity; };
		//void SetVelocity(const std::vector<pointf3>& v) { velocity = v; };

		/* This are functions used to fix points */
		void SetLeftFix(); // fix the points on the far left
		void FixTwoPoints(); // fix upper-left & upper-right
		void SetUpFix(); // fix the points on the upper

	private:
		// kernel part of the algorithm
		void SimulateOnce();

	protected:
		float h = 0.03f;  //²½³¤
		float stiff;
		float air_fric_coeff;
		std::vector<unsigned> fixed_id;  //fixed point id

		//mesh data
		std::vector<unsigned> edgelist;

		//simulation data
		std::vector<pointf3> positions;
		std::vector<vecf3> velocity;
		
		std::vector<double> mass;
		std::vector<vecf3> force;
		std::vector<double> stringOriginLen;
		std::vector<double> stringStiff; 
		std::vector<bool> isFixed; 
		vecf3 g;  // gravity  

		size_t nV;
		size_t nE; 

	};
}
