#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <tuple>
#include <cmath>
#include "myMath.h"
#include "myHelperFunctions.h"
using namespace std;







// no seperate GetKrylov function, all in the same j loop
// precond: 0-> no preconditioning, 1 -> Jacobi, 2 -> Gauss-Seidel, 3 -> ILU;
tuple< vector<double> , double,double> GMRES(bool sym,vector<double> JM, vector<double> VM, vector<double> x_0, vector<double> b , double m, int precond, double r_0normGlob)
{

    double rho;
    vector<double> x = x_0;
    vector<double> Diag; // will contain diagonal elements (for Jacobi preconditioner)
    vector<vector<double>> lowTri; // will contain lower triangle of A
    vector<vector<double>> Hm  = CreateMatrix(m+1,m);
    vector<vector<double>> Vm1;

    vector<double> r_0 = residualCalc(sym,JM, VM , x_0 ,b);
// ------------------------------------Preconditioning ------------------------------

        if(precond == 1)
        {
            cout << "Jacobi preconditioner applied!" << '\n';
            Diag = VecSlice(VM,0,JM[0]-3);
            r_0 = SolveLinearDiag(Diag,r_0);
        }

        if(precond == 2)
        {
            cout << "Gauss Seidel preconditioner applied!" << '\n';
            lowTri = extractLowTri(JM,VM);
            r_0 = forSub(lowTri,r_0);
        }


// ---------------------------------------------------------------------------------


    double r_0norm = VecNorm(r_0);
    //cout << "r_0 norm: " << r_0norm << '\n';
    //cout << "Tolerance = " << 1.0e-8 * r_0norm << '\n';
    vector<double> g (r_0.size(),0);
    g[0]= r_0norm;

    r_0 = VecNormalize(r_0);
    Vm1.push_back(r_0);



    double alpha;
    vector<double> c;
    vector<double> s;
    double m_til; // j at which solution converges

    for(int j = 0; j<m; j++)
    {
        vector<double> w;
        vector<double> z = MatVecProd(sym,JM,VM,Vm1[j]);
// ------------------------------------Preconditioning ------------------------------

        if(precond == 0)
        {
            w = z;
        }

        else if(precond == 1)
        {
            w = SolveLinearDiag(Diag,z);
        }
        else if(precond == 2)
        {
            w = forSub(lowTri,z);
        }

// -----------------------------------------------------------------------------------

        for(int i = 0; i<=j; i++)
        {
            vector<double> v_i = Vm1[i];

            double Hm_ij = DotProd(v_i,w);
            Hm[i][j] = Hm_ij;
            v_i = VecScalarMult(v_i,Hm_ij);
            w = VecSub(w,v_i);
        }
        double Hm_j1j = VecNorm(w);
        Hm[j+1][j] = Hm_j1j;
        w = VecElemDiv(w,Hm_j1j);
        Vm1.push_back(w);

        for(int k = 1; k<=j ; k++)
        {
            double tmp1 = c[k-1]*Hm[k-1][j] + s[k-1]*Hm[k][j];
            double tmp2 = c[k-1]*Hm[k][j] - s[k-1]*Hm[k-1][j];
            Hm[k-1][j] = tmp1;
            Hm[k][j] =  tmp2;
        }


        //cout << "After iteration j: " << j << " , Hm looks like: " << endl;
        //PrintMatrix(Hm);
        alpha = sqrt( pow(Hm[j][j],2.0) + pow(Hm[j+1][j],2.0));
        c.push_back(Hm[j][j] / alpha);
        s.push_back(Hm[j+1][j] / alpha);
        Hm[j][j] = alpha;
        Hm[j+1][j] = 0;

        g[j+1] = -s[j]*g[j];
        g[j] = c[j]*g[j];


        if (fabs(g[j+1]) <= 1.0e-8 * r_0normGlob)
        {
            cout << "Converged at iteration j: " << j << endl;
            cout << "With residual = " << g[j+1] << endl;
            m_til = j;
            rho = fabs(g[j+1]);

            return{x , rho, m_til};
            //exit(0);
            //break;
        }
        else
        {
            m_til = m;

        }
        rho = fabs(g[j+1]);
        cout << rho << '\n';
    }

    if(m_til == m)
    {
    Hm = MatrixTruncate(Hm, m_til-1); // this is Rm now
    g = VecSlice(g,0,m_til-1);  // "smaller system jxj"
    vector<double> y_star = backSub(Hm,g);

    vector<vector<double>> Vm1_T = MatrixTrans(Vm1); // need to transpose because I was pushing vectors into it as coloumns

    for(int i = 0; i<Vm1_T.size(); i++) // truncate matrix Vm1
    {
        Vm1_T[i] = VecSlice(Vm1_T[i],0,m-1);
    }

    vector<double> y_starVM = MatVecProd_simple(Vm1_T,y_star); // = Vm*y_star
    x = VecAdd(x,y_starVM);
    }

    return{x , rho, m_til};
}


vector<double> restarted_GMRES(bool sym,vector<double> JM, vector<double> VM , vector<double> x_0, vector<double> b, double m)
{
    cout << "Want to apply preconditioning? 0-> no preconditioning, 1 -> Jacobi, 2 -> Gauss-Seidel" << '\n';
    int precond;
    cin >> precond;

    double m_total = 0;
    vector<double> res = residualCalc(sym,JM, VM , x_0 ,b);

    double rho = VecNorm(res);
    double r_0normGlob = rho;

    vector<double> x = x_0;
    double tol = 1.0e-8 * r_0normGlob;
    cout << "Tolerance: " << tol << '\n';

    while(rho>=tol)
    {
        auto [xm,rho1, m_till] = GMRES(sym,JM,VM,x,b,m,precond,r_0normGlob);
        x = xm;
        rho = rho1;
        m_total += m_till;
    }
    cout << "m_total =" << m_total << endl;
    return x;
}





int main()
{

    double m = 600;
    auto[JM,VM,sym] = readMatrixMSR("gmres_test_msr.txt");
    int size = JM[0]-2;
    vector<double> x(size,1); // prescription of solution vector

    vector<double> b = MatVecProd(sym,JM,VM,x);
    vector<double> x_0(size,0);

    vector<double> x_result = restarted_GMRES(sym,JM,VM,x_0,b,m);


    return 0;
}
