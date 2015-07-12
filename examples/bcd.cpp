

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
struct f_param
{
    Matrix A;
    Vector b;
};
// evaluation of f(x)
Float f_eval(const Vector& x, const f_param& fpar)
{
    Vector r = fpar.A*x - fpar.b;
    Float val = r.norm();
    return 0.5 * val * val;
}
// gradient of f(x)
Vector f_grad(const Vector& x, const f_param& fpar)
{
    Matrix aat = fpar.A.transMulti(fpar.A);
    Vector g = aat*x - fpar.A.transpose() * fpar.b;
    return g;
}

/// tool functions
// generate sparse random x of size sz with nnz nonzeros
void gen_sp_rand_x(Vector& x, int sz, int nnz);
// inf norm of a vector
Float vec_inf_norm(const Vector& x);

int main(int argc, char *argv[])
{
    /// generate parameter for f
    // matrix dimention
    int nrow = 15;
    int ncol = 75;
    f_param fpar;
    fpar.A.setRandom(nrow, ncol);
    for(int i = 0; i < fpar.A.cols(); ++i)
        fpar.A.col(i).normalize();
    // x0, the ground truth
    Vector x0;
    Float spdensity = 0.05;
    int nnz = max(static_cast<int>(ncol*spdensity), 5);
    gen_sp_rand_x(x0, ncol, nnz);   
    // b, the rhs
    fpar.b = fpar.A * x0;
    // bind parameter to function f
    auto f_eval2 = std::bind(f_eval, _1, fpar);
    auto f_grad2 = std::bind(f_grad, _1, fpar);
    
    /// set param for r
    vector<FuncRInfo> ri{FuncRInfo(FuncRInfo::Norm1, ncol)};
    
    /// set tau
    Float tau = 0.8;       // tau should be tuned for specific problems
    
    /// set solver param
    BCDSolverParam sol_par;
    sol_par.set_f(f_eval2, f_grad2);
    sol_par.set_r(ri);
    sol_par.set_tau(tau);
    
    /// set solver option
    BCDSolverOption sol_opt;
    sol_opt.MaxIter = 20000;
    sol_opt.xTol    = 1e-6;
    
    /// set solver and solve
    BCDSolver bcdsol(sol_par, sol_opt);
    bcdsol.solve();
    
    /// outpout
    cout.unsetf(ios::fixed);
    cout.unsetf(ios::scientific);
    cout << "x0: " << endl;
    cout << x0 << endl;
    cout << "result: " << endl;
    cout << bcdsol.result() << endl;
    
}


/**************************************************
 * 
 *  sub routines
 * 
 * 
 */


Float vec_inf_norm(const Vector& x)
{
    auto v = x[0];
    for(const auto& xi:x)
    {
        if(v < xi)  v = xi;
        if(v < -xi) v = -xi;
    }
    return v;
}

void gen_sp_rand_x(Vector& x, int sz, int nnz)
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
