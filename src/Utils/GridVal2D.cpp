#include "GridVal2D.h"

GridVal2D::GridVal2D(int x, int y, double val) {
    this->x = x;
    this->y = y;
    this->val = val;
}

int GridVal2D::getX() const {
    return this->x;
}

int GridVal2D::getY() const {
    return this->y;
}

double GridVal2D::getVal() const {
    return this->val;
}