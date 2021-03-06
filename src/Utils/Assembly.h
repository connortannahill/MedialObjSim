#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <Eigen/Sparse>
#include <unordered_map>

using namespace std;

class Assembly {
public:
    Assembly();
    // ~Assembly();
    Eigen::SparseMatrix<double> *M;
    Eigen::SparseMatrix<double> *D;
    Eigen::SparseMatrix<double> *W;

    int getNPnts();
    int getD();
    void setNPnts(int nPnts);
    void setD(int d);
    // virtual double U(int nodeId) = 0;
    // virtual double nablaU(int nodeId) = 0;
    virtual void prox(double dt, Eigen::VectorXd &x, Eigen::VectorXd &DXpU, Eigen::VectorXd &z) = 0;
    virtual void updateAfterStep(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x) = 0;
    virtual void copyX(Eigen::VectorXd &tar) = 0;
    virtual void predictX(double dt, Eigen::VectorXd &xPrev, Eigen::VectorXd &x, Eigen::VectorXd &xBar) = 0;
    void printDiff();
    virtual ~Assembly();
protected:
    virtual void buildMassMatrix(Eigen::VectorXd &m);
    virtual void buildDMatrix();
    virtual void buildWMatrix(double w);
    bool nPntsSet;
    bool dSet;
    bool dAlloc;
    bool wAlloc;
    bool mAlloc;
    double w;
    int nPnts;
    int d;
     

   
};

#endif