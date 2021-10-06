#include "SolidObject3D.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include "../Utils/SolidParams.h"

using namespace std;

/**
 * Define the copy constructor
*/
SolidObject3D::SolidObject3D(const SolidObject3D &obj) {
    this->shapeFun = obj.shapeFun;
    this->params = obj.params;
    this->mass = obj.mass;
    this->density = obj.density;
    this->objType = obj.objType;

    this->E = obj.E;
    this->eta = obj.eta;

    this->u0 = obj.u0;
    this->v0 = obj.v0;
    this->w0 = obj.w0;

}

SolidObject3D::SolidObject3D(double u0, double v0, double w0,
                SolidObject3D::ObjectType objType,
                double (*shapeFun)(double,double,double,SolidParams&),
                SolidParams &params) {
    // Assign the instance variables, move on.
    this->shapeFun = shapeFun;
    this->params = params;

    params.getParam("mass", this->mass);
    params.getParam("density", this->density);

    this->u0 = u0;
    this->v0 = v0;
    this->w0 = w0;

    this->objType = objType;

    if (this->objType == ObjectType::DEFORMABLE) {
        params.getParam("E", this->E);
        params.getParam("eta", this->eta);
    } else {
        this->E = 0.0;
        this->eta = 0.0;
    }
}

double SolidObject3D::shapeEval(double x, double y, double z) {
    return this->shapeFun(x, y, z, this->params);
}


// double SolidObject3D::getCx() {
//     return this->cx;
// }
// double SolidObject3D::getCy() {
//     return this->cy;
// }
// double SolidObject3D::getCz() {
//     return this-cz;
// }

double SolidObject3D::getMass() {
    return this->mass;
}

double SolidObject3D::getDensity() {
    return this->density;
}

double SolidObject3D::getU0() {
    return this->u0;
}

double SolidObject3D::getV0() {
    return this->v0;
}

double SolidObject3D::getW0() {
    return this->w0;
}

double SolidObject3D::getE() {
    return this->E;
}

double SolidObject3D::getEta() {
    return this->eta;
}

SolidObject3D::ObjectType SolidObject3D::getObjType() {
    return this->objType;
}