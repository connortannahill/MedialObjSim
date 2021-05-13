#include "Diff3DComparator.h"
#include "../Utils/SimUtilities.h"
#include <math.h>

using namespace simstructs;
using namespace simutils;

int Diff3DComparator::operator() (const struct surfaceDiff& p1, const struct surfaceDiff& p2) {
    return p1.val > p2.val;
}