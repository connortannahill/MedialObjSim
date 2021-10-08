#ifndef KD_POINT_CLOUD_3D_H
#define KD_POINT_CLOUD_3D_H

#include "MassSpring3D.h"

#include <stdlib.h>
#include <vector>

// class MassSpring3D;

using namespace std;
using namespace mass_spring;

class KDPointCloud3D {
public:
    KDPointCloud3D();
	void resetCloud();
	void addMSS(mass_spring::MassSpring3D &mss);

    // Must return the number of data points
	inline size_t kdtree_get_point_count() const {
		return points->size();
	};

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
		if (dim == 0) {
			return points->at(idx)->x;
		} else if (dim == 1) {
			return points->at(idx)->y;
		} else {
			return points->at(idx)->z;
		}
	};

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const {
		return false;
	};


    ~KDPointCloud3D();
// protected:
    vector<massPoint3D*> *points;
    
};

#endif