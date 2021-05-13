#ifndef SOLID_OBJECT_3D_H
#define SOLID_OBJECT_3D_H

#include <unordered_map>
#include <string>
#include "../Utils/SolidParams.h"

using namespace std;

/**
 * Class that is used to define an object to be placed in the pool
*/
class SolidObject3D {
public:
    /* Enumeration to handle the different kinds of objects we consider */
    enum ObjectType {
        RIGID,
        DEFORMABLE,
        STATIC
    };

    /* Base class constructor */
    SolidObject3D(double u0, double v0, double w0, ObjectType objType,
                    double (*shapeFun)(double,double,double,SolidParams&),
                    SolidParams &params);
    SolidObject3D(const SolidObject3D &obj);

    /* Wrapped evaluation of the shape function */
    double shapeEval(double x, double y, double z);

    /* Getters & setters */
    double getMass();
    double getDensity();
    double getU0();
    double getV0();
    double getW0();
    double getE();
    double getEta();
    ObjectType getObjType();
private:
    /* Implicit function defining the level set */
    double (*shapeFun)(double,double,double,SolidParams&);

    /* Params for the level set */
    SolidParams params;

    /* Variables relevant to keeping track of information over the course of the fluid sim */
    double mass;
    double density;
    double E;
    double eta;
    double u0;
    double v0;
    double w0;
    ObjectType objType;
};

#endif