#include "GridValComparator.h"
#include "GridVal2D.h"
#include "GridVal3D.h"

int GridValComparator::operator() (const GridVal2D& p1, const GridVal2D& p2) {
    return p1.getVal() > p2.getVal();
}

int GridValComparator::operator() (const GridVal3D& p1, const GridVal3D& p2) {
    return p1.getVal() > p2.getVal();
}
