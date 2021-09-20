#include "ObjectSeeder.h"
#include "Boundary.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>
#include <iostream>
#include "../Utils/SolidParams.h"
#include <ctime>

using namespace std;

Eigen::Vector2d getRandom() {
    Eigen::Vector2d rand;
    double dummy = std::rand();
    dummy += 1.0;

    rand(0) = ((double) std::rand() / (RAND_MAX));
    rand(1) = ((double) std::rand() / (RAND_MAX));
    return rand;
}

Eigen::Vector2d getCandVec(Eigen::Vector2d &za, Eigen::Vector2d &zb) {
    Eigen::Vector2d rand(getRandom());
    Eigen::Vector2d cand(za + rand.cwiseProduct(zb - za));
    return cand;
}

double distProj(Eigen::Vector2d x1, Eigen::Vector2d x2, Eigen::Vector2d c) {
    double projLen = (c-x1).dot(x2-x1)/(x2 - x1).squaredNorm();
    return (projLen*(x2 - x1) - (c - x1)).norm();
}

bool intersectsBoundary(Eigen::Vector2d &circ, double r, double eps, Boundary &boundary) {
    // Boundary vectors
    Eigen::Vector2d x1;
    Eigen::Vector2d x2;

    // Bot boundary
    x1 << boundary.getXmin(), boundary.getYmin();
    x2 << boundary.getXmax(), boundary.getYmin();

    if (distProj(x1, x2, circ) < r+eps)
        return true;
    
    // Top boundary
    x1 << boundary.getXmin(), boundary.getYmax();
    x2 << boundary.getXmax(), boundary.getYmax();

    if (distProj(x1, x2, circ) < r+eps)
        return true;

    // Left boundary
    x1 << boundary.getXmin(), boundary.getYmin();
    x2 << boundary.getXmin(), boundary.getYmax();

    if (distProj(x1, x2, circ) < r+eps)
        return true;

    // Right boundary
    x1 << boundary.getXmax(), boundary.getYmin();
    x2 << boundary.getXmax(), boundary.getYmax();

    if (distProj(x1, x2, circ) < r+eps)
        return true;
    
    // Indicate no intersection found
    return false;
}

std::vector<SolidObject> ObjectSeeder::randomInitialize(int nObjs,
        SolidObject::ObjectType objType, double r, double eps,
        Boundary &boundary, double (*shapeFun)(double,double,SolidParams&),
        SolidParams &params) {
    
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    
    std::vector<SolidObject> objs;

    // Vector of D+1 dimension vectors, vector of form [X, Y, R]^T
    std::vector<Eigen::Vector2d> circs;
    Eigen::Vector2d circCand;

    Eigen::Vector2d za;
    za << boundary.getXmin()+(r+eps), boundary.getYmin()+(r+eps);

    Eigen::Vector2d zb;
    zb << boundary.getXmax()-(r+eps), boundary.getYmax()-(r+eps);

    bool found;

    for (int i = 0; i < nObjs; i++) {
        found = false;

        if (i == 0) {
            // Random between [0, 1]
            circs.push_back(getCandVec(za, zb));
            continue;
        }

        while (!found) {
            circCand = getCandVec(za, zb);

            found = true;

            // Iterate through existing circles, checking for overlap, and/or intersection with the boundary
            for (auto circ = circs.begin(); circ != circs.end(); ++circ) {
                // Overlap detected
                if ((*circ - circCand).norm() < 2*r + eps) {
                    found = false;
                    break;
                }

                // Boundary intersection detected
                if (intersectsBoundary(circCand, r, eps, boundary)) {
                    found = false;
                    break;
                }
            }
        }

        // If circle location has been found, append
        circs.push_back(circCand);
    }

    // Now, create all of the solid objects
    SolidParams initParams(params);

    for (int i = 0; i < nObjs; i++) {
        initParams.addParam("cx", circs.at(i)(0));
        initParams.addParam("cy", circs.at(i)(1));

        SolidObject obj(0.0, 0.0, objType, shapeFun, initParams);
        objs.push_back(obj);
    }

    std::cout << "About to return objs, circs:" << std::endl;
    for (int i = 0; i < circs.size(); i++) {
        std::cout << circs.at(i).transpose() << std::endl;
    }

    return objs;

}
