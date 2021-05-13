#ifndef DIFF_3D_COMPARATOR_H
#define DIFF_3D_COMPARATOR_H

#include "../Utils/SimUtilities.h"

using namespace simstructs;

class Diff3DComparator {
public:
    int operator() (const struct surfaceDiff& p1, const struct surfaceDiff& p2);
};

#endif