#ifndef SIM_PARAMS_3D_H
#define SIM_PARAMS_3D_H

class SimParams3D {
protected:
    bool reSet=false;
    bool mssNxSet=false;
    bool mssNySet=false;
    bool mssNzSet=false;
    bool nxSet=false;
    bool nySet=false;
    bool nzSet=false;
    bool muSet=false;
    bool repulseModeSet=false;
    bool repulseDistSet=false;
    bool collisionStiffnessSet=false;
    bool collisionDistSet=false;
public:
    void setRe(double Re);
    void setMssNx(int mssNx);
    void setMssNy(int mssNy);
    void setMssNz(int mssNy);
    void setNx(int nx);
    void setNy(int ny);
    void setNz(int ny);
    void setDtFix(double dtFix);
    void setMu(double mu);
    void setUseEno(bool useEno);
    bool useEno = false;

    double Re;
    bool dtFixSet=false;
    int mssNx;
    int mssNy;
    int mssNz;
    int nx;
    int ny;
    int nz;
    double dtFix;
    int methodOrd = 1;
    int updateMode = 2;
    int elementMode = 0;
    double mu;
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    double u0 = 0.0;
    double v0 = 0.0;
    double w0 = 0.0;

    bool checkParams();

    double admmTol = 1e-6;

    // This variable determines whether repulsive forces are to be used to keep objects apart and which
    // kind should be used
    //  0: no attempt is made to handle collisions
    //  1: use fixed repulsion vector to move entire objects. Collision avoidance, not handling
    //  2: trivial collision handling. Inaccurate, but applies a repulsive forces for points detected in radius search.
    //     requires the use of kd trees.
    int repulseMode=-1;
    double repulseDist=-1;
    double collisionStiffness;
    double collisionDist;
    void setRepulseMode(int mode);
    void setRepulseDist(double dist);
    void setCollisionStiffness(double stiffness);
    void setCollisionDist(double colDist);
    void setAdmmTol(double admmTol);
    void setUpdateMode(int mode);
    void setElementMode(int mode);
    void setGx(double gx);
    void setGy(double gy);
    void setGz(double gz);
    double getGx();
    double getGy();
    double getGz();

    void setU0(double u0);
    void setV0(double v0);
    void setW0(double w0);
    double getU0();
    double getV0();
    double getW0();
};

#endif