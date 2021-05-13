#ifndef GRID_VAL_3D_H
#define GRID_VAL_3D_H

class GridVal3D {
public:
    GridVal3D(int x, int y, int z, double val);
    int getX() const;
    int getY() const;
    int getZ() const;
    double getVal() const;
private:
    int x, y, z;
    double val;
};

#endif