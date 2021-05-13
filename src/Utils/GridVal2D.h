#ifndef GRID_VAL_2D_H
#define GRID_VAL_2D_H

class GridVal2D {
public:
    GridVal2D(int x, int y, double val);
    int getX() const;
    int getY() const;
    double getVal() const;
private:
    int x, y;
    double val;
};

#endif