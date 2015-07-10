

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>
using namespace std;

#include "Solvers/BCD.hpp"


class MyFuncF : public FuncF
{
public:
    Float eval(const Vector& x)
    {
        Vector r = A*x - b;
        Float val = r.norm();
        return 0.5 * val * val;
    }

    Vector grad(const Vector& x)
    {
        Matrix aat = A.transMulti(A);
        Vector g = aat*x - A.transpose() * b;
        return g;
    }
    
    Matrix A;
    Vector b;
};

// generate sparse random x of size sz with nnz nonzeros
void gen_sp_rand_x(Vector& x, int sz, int nnz);
// inf norm of a vector
Float vec_inf_norm(const Vector& x);

int main(int argc, char *argv[])
{
    
    /// dim
    int nrow = 15;
    int ncol = 75;
    
    
    /// set f, the differentiable part
    MyFuncF f;
    
    f.A.setRandom(nrow, ncol);
    for(int i = 0; i < f.A.cols(); ++i)
        f.A.col(i).normalize();
    
    // x0, the ground truth
    Vector x0;
    Float spdensity = 0.05;
    int nnz = max(static_cast<int>(ncol*spdensity), 5);
    gen_sp_rand_x(x0, ncol, nnz);
    
    // b, the rhs
    f.b = f.A * x0;
    
    
    /// set r_i, the non-smooth part
    int dim = f.A.cols();
    vector<FuncR> ri{FuncR(FuncR::Norm1, dim)};
    
    /// set par
    BCDSolverParam par;
    // set tau
    par.tau = 0.8;      // tau should be tuned for specific problems
    par.maxIter = 20000;
    par.objTol = 1e-6;
    par.xTol = 1e-6;
    
    //
    BCDSolver sol(par, f, ri);
    Vector x_init(dim);
    x_init.setZeros();
    sol.init_x(x_init);
    sol.solve();
    cout << "-----------------------------------" << endl;
    cout << "x0 is " << x0 << endl;
    
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
    // rand id
    vector<int> idx(sz);
    for(int i = 0; i < sz; ++i) idx[i] = i;
    shuffle(idx.begin(), idx.end(), generator);
    // rand val
    uniform_real_distribution<Float> distribution(0.0,1.0);
    auto dice = std::bind(distribution, generator);
    for(int i = 0; i < nnz; ++i)
        x[idx[i]] = dice();
}
