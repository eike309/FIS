#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include <fstream>
#include "myMath.h"
#include "myHelperFunctions.h"
using namespace std;


vector<double> CG(bool sym, vector<double>JM, vector<double>VM, vector<double> x_0, vector<double> b)
{
  ofstream myfile;
  myfile.open ("CG_output.txt");

    vector<double> x = x_0;
    vector<double> r = residualCalc(sym,JM,VM,x,b);
    vector<double> r1;
    vector<double> r_0 =r;
    vector<double> p = r;
    vector<double> p1;
    //vector<double> r1;
    double r_0normSQ = DotProd(r,r); // norm squared
    double r_normSQ = r_0normSQ;
    int m = 0;


    while(VecNorm(r) >= 1.0e-8 * VecNorm(r_0))
    {
        double dotprod_r = DotProd(r,r);
        vector<double> z = MatVecProd(sym,JM,VM,p);

        double alpha = dotprod_r / DotProd(z,p);

        x = VecAdd(x,VecScalarMult(p,alpha));

        r1 = VecSub(r, VecScalarMult(z,-alpha));

        double dotprod_r1 = DotProd(r1,r1);
        double beta = dotprod_r1 / dotprod_r;

        p1 = VecAdd(r,VecScalarMult(p,beta));
        p = p1;

        myfile << m << "  " << VecNorm(r) <<'\n';

        r = r1;
        m++;


    }
    cout << "Converged at iteration m: " << m << endl;
    cout << "With norm of residual = " << VecNorm(r) << endl;
    myfile.close();
    return x;
}





int main()
{

    /*
    // Testcase
    vector<double> JM = {5,  7,  8,  9,  2,  3,  1,  1};
    vector<double> VM = {25, 18, 11, 0, 15, -5, 15, -5};
    vector<double> x;
    vector<double> x_0 = {0,0,0};
    vector<double> b = {35, 33, 6};
    bool sym = 1;
    // Testcase
     */


   auto[JM, VM , sym] = readMatrixMSR("cg_test_msr.txt");
   cout << "sym: " << sym;
    int size = JM[0]-2;
    vector<double> x(size,1); // prescription of solution vector

    int m_guess = 200000;
    vector<double> b = MatVecProd(sym,JM,VM,x);
    vector<double> x_0(size,0);

    vector<double> solution = CG(sym, JM, VM, x_0, b);


}
