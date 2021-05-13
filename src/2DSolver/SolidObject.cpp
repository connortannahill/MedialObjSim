#include "SolidObject.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include "../Utils/SolidParams.h"

using namespace std;

/**
 * Define the copy constructor
*/
SolidObject::SolidObject(const SolidObject &obj) {
    this->shapeFun = obj.shapeFun;
    this->params = obj.params;
    this->mass = obj.mass;
    this->density = obj.density;
    this->objType = obj.objType;

    this->E = obj.E;
    this->eta = obj.eta;

    this->u0 = obj.u0;
    this->v0 = obj.v0;
}

SolidObject::SolidObject(double u0, double v0, SolidObject::ObjectType objType,
                            double (*shapeFun)(double,double,SolidParams&),
                            SolidParams &params) {
    // Assign the instance variables, move on.
    this->shapeFun = shapeFun;
    this->params = params;

    params.getParam("mass", this->mass);
    params.getParam("density", this->density);

    this->objType = objType;

    if (this->objType == ObjectType::DEFORMABLE) {
        params.getParam("E", this->E);
        params.getParam("eta", this->eta);
    } else {
        this->E = 0;
        this->eta = 0;
    }

    this->u0 = u0;
    this->v0 = v0;
}

double SolidObject::shapeEval(double x, double y) {
    return this->shapeFun(x, y, this->params);
}

double SolidObject::getMass() {
    return this->mass;
}

double SolidObject::getDensity() {
    return this->density;
}

double SolidObject::getU0() {
    return this->u0;
}

double SolidObject::getV0() {
    return this->v0;
}

double SolidObject::getE() {
    return this->E;
}

double SolidObject::getEta() {
    return this->eta;
}

SolidObject::ObjectType SolidObject::getObjType() {
    return objType;
}