#ifndef GRID_VAL_COMPARATOR_H
#define GRID_VAL_COMPARATOR_H

#include "GridVal2D.h"
#include "GridVal3D.h"

class GridValComparator {
public:
    int operator() (const GridVal2D& p1, const GridVal2D& p2);
    int operator() (const GridVal3D& p1, const GridVal3D& p2);
};

#endif