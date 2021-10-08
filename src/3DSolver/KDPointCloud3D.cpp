#include "KDPointCloud3D.h"
#include <stdlib.h>
#include <vector>
#include "MassSpring3D.h"

using namespace std;
using namespace mass_spring;

KDPointCloud3D::KDPointCloud3D() {
    points = new vector<massPoint3D*>();
}

void KDPointCloud3D::resetCloud() {
    points->clear();
}

void KDPointCloud3D::addMSS(MassSpring3D &mss) {
    // Append the current MSS to the nodeList
    for (auto bId = mss.boundaryNodeIdList->begin(); bId != mss.boundaryNodeIdList->end(); ++bId) {
        points->push_back(&mss.pntList->at(*bId));
    }
}

KDPointCloud3D::~KDPointCloud3D() {
    delete points;
}