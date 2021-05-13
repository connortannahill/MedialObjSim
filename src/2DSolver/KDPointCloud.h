#ifndef KD_POINT_CLOUD_H
#define KD_POINT_CLOUD_H

#include "./MassSpring2D.h"

#include <stdlib.h>
#include <vector>

using namespace std;
using namespace mass_spring;

class KDPointCloud {
public:
    KDPointCloud();
	void resetCloud();
	void addMSS(MassSpring2D &mss);

    // Must return the number of data points
	inline size_t kdtree_get_point_count() const {
		return points->size();
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	// inline double kdtree_get_pt(const size_t idx, const size_t dim) const;
	inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
		if (dim == 0) {
			return points->at(idx)->x;
		} else {
			return points->at(idx)->y;
		}
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const {
		return false;
	}

    ~KDPointCloud();
    vector<massPoint2D*> *points;
};

#endif