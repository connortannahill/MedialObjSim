#include "GridVal3D.h"

GridVal3D::GridVal3D(int x, int y, int z, double val) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->val = val;
}

int GridVal3D::getX() const {
    return this->x;
}

int GridVal3D::getY() const {
    return this->y;
}

int GridVal3D::getZ() const {
    return this->z;
}

double GridVal3D::getVal() const {
    return this->val;
}