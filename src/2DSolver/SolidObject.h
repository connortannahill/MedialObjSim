#ifndef SOLID_OBJECT_H
#define SOLID_OBJECT_H

#include <string>
#include "../Utils/SolidParams.h"

/**
 * Class that is used to define an object to be placed in the pool
*/
class SolidObject {
public:
    /* Enumeration to handle the different kinds of objects we consider */
    enum ObjectType {
        RIGID,
        DEFORMABLE,
        STATIC
    };

    /* Base class constructor */
    SolidObject(double u0, double v0, ObjectType objType,
                double (*shapeFun)(double,double,SolidParams&),
                SolidParams &params);
    SolidObject(const SolidObject &obj);


    /* Wrapped evaluation of the shape function */
    double shapeEval(double x, double y);

    /* Getters & setters */
    double getMass();
    double getDensity();
    double getE();
    double getEta();
    double getU0();
    double getV0();
    ObjectType getObjType();
private:
    /* Implicit function defining the level set */
    double (*shapeFun)(double,double,SolidParams&);

    /* Params for the level set */
    SolidParams params;

    /* Variables relevant to keeping track of information over the course of the fluid sim */
    double mass;
    double density;
    double E;
    double eta;
    double u0;
    double v0;
    ObjectType objType;
};

#endif