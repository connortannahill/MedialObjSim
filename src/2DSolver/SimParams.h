#ifndef SIM_PARAMS_H
#define SIM_PARAMS_H

class SimParams {
protected:
    bool reSet=false;
    bool mssNxSet=false;
    bool mssNySet=false;
    bool nxSet=false;
    bool nySet=false;
    bool muSet=false;
    bool repulseModeSet=false;
    bool repulseDistSet=false;
    bool collisionStiffnessSet=false;
    bool collisionDistSet=false;
public:
    void setRe(double Re);
    void setMssNx(int mssNx);
    void setMssNy(int mssNy);
    void setNx(int nx);
    void setNy(int ny);
    void setDtFix(double dtFix);
    void setMu(double mu);
    void setUseEno(bool use);
    bool useEno = false;

    double Re;
    bool dtFixSet=false;
    int mssNx;
    int mssNy;
    int nx;
    int ny;
    int methodOrd = 1;
    double dtFix;
    double mu;
    int updateMode = 2;
    int elementMode = 0;
    double admmTol = 1e-6;
    double h;
    double gx = 0.0;
    double gy = 0.0;

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
    double getGx();
    double getGy();
    double getH();


    bool checkParams();
};

#endif