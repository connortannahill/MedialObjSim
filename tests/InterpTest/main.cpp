
#include <iostream>
#include "./src/SimUtilities.h"

using namespace std;

double fun(double x) {
    // return x*x*x + 3*x*x + 2;
    return 2*x + 1;
}

double DfunDx(double x) {
    return 3*x*x + 6*x;
}

int main() {
    // Testing the interpolation routines

    // Degree of the peicewise polynomial
    int n = 1;

    // Get the Barycentric weights
    double bw[n+1];
    simutils::barycentricInterp(n, bw);

    cout << "The bary weights: " << endl;
    for (int i = 0; i < n+1; i++) {
        cout << bw[i] << ", ";
    }
    cout << endl;

    // Create mesh to evaluate the solution
    // double xa = -1;
    double xa = 0;
    double xb = 1;
    double width = xb - xa;
    int N = 10;
    double x[N+1];
    simutils::linspace(xa, xb, N+1, x);

    // Evaluate the soln
    double y[N+1];
    for (int i = 0; i < N+1; i++) {
        y[i] = fun(x[i]);
    }

    // Evaluate the exact solution at midpoints
    double xmid[N];
    double ymid[N];
    for (int i = 0; i < N; i++) {
        xmid[i] = (x[i+1] + x[i])/2.0;
        ymid[i] = fun(xmid[i]);
    }

    // Evaluate the interpolating polynomial at the midpoint
    double yeval[N];
    double ypnts[2];
    for (int i = 0; i < N; i++) {
        ypnts[0] = y[i];
        ypnts[1] = y[i+1];

        yeval[i] = simutils::barycentricInterv(xmid[i], x[i], x[i+1], width, n, bw, ypnts);
    }

    cout << "h from mesh" << (x[1] - x[0]) << endl;
    cout << "h pred" << 1/(n+1) << endl;

    cout << "Printing out the values:" << endl;
    cout << "x:" << endl;
    for (int i = 0; i < N+1; i++) {
        cout << x[i] << ", ";
    }
    cout << "\n";

    cout << "xmid:" << endl;
    for (int i = 0; i < N; i++) {
        cout << xmid[i] << ", ";
    }
    cout << "\n";


    cout << "yexact:" << endl;

    for (int i = 0; i < N; i++) {
        cout << ymid[i] << ", ";
    }
    cout << "\n";

    cout << "yinterp:" << endl;

    for (int i = 0; i < N; i++) {
        cout << yeval[i] << ", ";
    }
    cout << "\n";
    
    

}