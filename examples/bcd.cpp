

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <functional>
using namespace std;
using namespace std::placeholders;

#include "Solvers/BCD.hpp"

///
/// example
/// 
///     min  0.5*||Ax-b||_2^2 + tau*||x||_1
/// 
/// here  f(x) = 0.5*||Ax-b||_2^2
/// and   r(x) = ||x||_1
///
/// 

/// f(x)
// parameter of f(x): A & b
struct FParam
{
    Matrix A;
    Vector b;
};
// evaluation of f(x)
Float fEval(const Vector& x, const FParam& fpar)
{
    Vector r = fpar.A*x - fpar.b;
    Float val = r.norm();
    return 0.5 * val * val;
}
// gradient of f(x)
Vector fGrad(const Vector& x, const FParam& fpar)
{
    Matrix aat = fpar.A.transMulti(fpar.A);
    Vector g = aat*x - fpar.A.transpose() * fpar.b;
    return g;
}

/// tool functions
// generate sparse random x of size sz with nnz nonzeros
void generateSparseRandomX(Vector& x, int sz, int nnz);
// inf norm of a vector
Float vectorInfNorm(const Vector& x);

int main(int argc, char *argv[])
{
    /// generate parameter for f
    // matrix dimention
    int nRow = 15;
    int nCol = 75;
    FParam fpar;
    fpar.A.setRandom(nRow, nCol);
    for(int i = 0; i < fpar.A.cols(); ++i)
        fpar.A.col(i).normalize();
    // x0, the ground truth
    Vector x0;
    Float spdensity = 0.05;
    int nnz = max(static_cast<int>(nCol*spdensity), 5);
    generateSparseRandomX(x0, nCol, nnz);   
    // b, the rhs
    fpar.b = fpar.A * x0;
    // bind parameter to function f
    auto f_eval2 = std::bind(fEval, _1, fpar);
    auto f_grad2 = std::bind(fGrad, _1, fpar);
    
    /// set param for r
    vector<FuncRInfo> ri{FuncRInfo(FuncRInfo::Norm1, nCol)};
    
    /// set tau
    Float tau = 0.8;       // tau should be tuned for specific problems
    
    /// set solver param
    BCDSolverParam sol_par;
    sol_par.setF(f_eval2, f_grad2);
    sol_par.setR(ri);
    sol_par.setTau(tau);
    
    /// set solver option
    BCDSolverOption sol_opt;
    sol_opt.MaxIter = 20000;
    sol_opt.XTol    = 1e-6;
    
    /// set solver and solve
    BCDSolver bcd_sol(sol_par, sol_opt);
    bcd_sol.solve();
    
    /// outpout
    cout.unsetf(ios::fixed);
    cout.unsetf(ios::scientific);
    cout << "x0: " << endl;
    cout << x0 << endl;
    cout << "result: " << endl;
    cout << bcd_sol.result() << endl;
    
}


/**************************************************
 * 
 *  sub routines
 * 
 * 
 */


Float vectorInfNorm(const Vector& x)
{
    auto v = x[0];
    for(const auto& xi:x)
    {
        if(v < xi)  v = xi;
        if(v < -xi) v = -xi;
    }
    return v;
}

void generateSparseRandomX(Vector& x, int sz, int nnz)
{
    x.resize(sz);
    x.setZeros();
    // rand engine
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    // rand index
    vector<int> idx(sz);
    for(int i = 0; i < sz; ++i) idx[i] = i;
    shuffle(idx.begin(), idx.end(), generator);
    // rand value
    uniform_real_distribution<Float> distribution(0.0,1.0);
    auto dice = std::bind(distribution, generator);
    for(int i = 0; i < nnz; ++i)
        x[idx[i]] = dice();
}
