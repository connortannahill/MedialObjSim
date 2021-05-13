#include "KDPointCloud.h"
#include <stdlib.h>
#include <vector>
#include "./MassSpring2D.h"

using namespace std;
using namespace mass_spring;

KDPointCloud::KDPointCloud() {
    points = new vector<massPoint2D*>();
}

void KDPointCloud::resetCloud() {
    points->clear();
}

void KDPointCloud::addMSS(MassSpring2D &mss) {
    // Append the current MSS to the nodeList
    for (auto bId = mss.boundaryNodeIdList->begin(); bId != mss.boundaryNodeIdList->end(); ++bId) {
        points->push_back(&mss.pntList->at(*bId));
    }
}

KDPointCloud::~KDPointCloud() {
    delete points;
}