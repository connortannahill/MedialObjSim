#ifndef OBJECT_SEEDER_H
#define OBJECT_SEEDER_H
#include <string>
#include <vector>
#include <unordered_map>
#include "SolidObject.h"
#include "Boundary.h"

using namespace std;

class ObjectSeeder {
public:
    std::vector<SolidObject> randomInitialize(int nObjs,
            SolidObject::ObjectType objType, double r, double eps,
            Boundary &boundary, double (*shapeFun)(double,double,SolidParams&),
            SolidParams &params);
};

#endif