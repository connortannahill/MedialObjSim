#include "SimParams3D.h"
#include <iostream>

using namespace std;

void SimParams3D::setRe(double Re) {
    this->Re = Re;
    reSet = true;
}

void SimParams3D::setMssNx(int mssNx) {
    this->mssNx = mssNx;
    mssNxSet = true;
}

void SimParams3D::setMssNy(int mssNy) {
    this->mssNy = mssNy;
    mssNySet = true;
}

void SimParams3D::setMssNz(int mssNz) {
    this->mssNz = mssNz;
    mssNySet = true;
}

void SimParams3D::setGx(double gx) {
    this->gx = gx;
}
void SimParams3D::setGy(double gy) {
    this->gy = gy;
}
void SimParams3D::setGz(double gz) {
    this->gz = gz;
}

double SimParams3D::getGx() {
    return this->gx;
}
double SimParams3D::getGy() {
    return this->gy;
}
double SimParams3D::getGz() {
    return this->gz;
}

void SimParams3D::setU0(double u0) {
    this->u0 = u0;
}
void SimParams3D::setV0(double v0) {
    this->v0 = v0;
}
void SimParams3D::setW0(double w0) {
    this->w0 = w0;
}

double SimParams3D::getU0() {
    return this->u0;
}
double SimParams3D::getV0() {
    return this->v0;
}
double SimParams3D::getW0() {
    return this->w0;
}

void SimParams3D::setUseEno(bool useEno) {
    this->useEno = useEno;
}

void SimParams3D::setNx(int nx) {
    this->nx = nx;

    if (!mssNxSet) {
        mssNx = nx;
    }

    nxSet = true;
}

void SimParams3D::setNy(int ny) {
    this->ny = ny;

    if (!mssNySet) {
        mssNy = ny;
    }

    nySet = true;
}

void SimParams3D::setNz(int nz) {
    this->nz = nz;

    if (!mssNzSet) {
        mssNz = nz;
    }

    nzSet = true;
}

void SimParams3D::setDtFix(double dtFix) {
    this->dtFix = dtFix;
    dtFixSet = true;
}

void SimParams3D::setMu(double mu) {
    this->mu = mu;
    muSet = true;
}

bool SimParams3D::checkParams() {
    cout << nxSet << endl;
    cout << nySet << endl;
    cout << nzSet << endl;
    cout << reSet << endl;
    cout << muSet << endl;
    cout << repulseModeSet << endl;
    cout << repulseDistSet << endl;
    cout << collisionStiffnessSet << endl;
    cout << collisionDistSet << endl;
    return nxSet && nySet && nzSet && reSet && muSet && repulseModeSet
            && collisionStiffnessSet && collisionDistSet;
}

void SimParams3D::setRepulseMode(int mode) {
    this->repulseMode = mode;
    this->repulseModeSet = true;
}

void SimParams3D::setRepulseDist(double dist) {
    this->repulseDist = dist;
    this->repulseDistSet = true;
}

void SimParams3D::setCollisionStiffness(double stiffness) {
    this->collisionDist = stiffness;
    this->collisionStiffnessSet = true;
}

void SimParams3D::setCollisionDist(double colDist) {
    this->collisionDist = colDist;
    this->collisionDistSet = true;
}

void SimParams3D::setAdmmTol(double admmTol) {
    this->admmTol = admmTol;
}

void SimParams3D::setUpdateMode(int mode) {
    this->updateMode = mode;
}

void SimParams3D::setElementMode(int mode) {
    this->elementMode = mode;
}
