

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include "Solvers/BCD.hpp"


//class MyFuncF : public FuncF
//{
//public:
//    Float eval(const Vector& x)
//    {
//        Float val = 0;
//        int dim = x.size();
//        for(int i = 0; i < dim; ++i)
//        {
//            val += (x[i]-5.0) * (x[i]-5.0);
//        }
//        return val;
//    }

//    Vector grad(const Vector& x)
//    {
//        int dim = x.size();
//        Vector g(dim);
//        for(int i = 0; i < dim; ++i)
//        {
//            g[i] = 2.0 * (x[i]-5.0);
//        }
//        return g;
//    }
//};

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
    
    
    void loadA(const string& fn)
    {
        Float val;
        ifstream ifs(fn.c_str());
        int m, n;
        ifs >> m >> n;
        cout << "m : " << m << ", n : " << n << endl;
        A.resize(m,n);
        for(int i = 0; i < m; ++i)
        {
            for(int j = 0; j < n; ++j)
            {
                ifs >> val;
                A(i,j) = val;
            }
        }
        ifs.close();
//        cout << "load A ok" << endl;
//        cout << "A(m,n)" << A(m-1,n-1) << endl;
        
//        // colume normalize
//        auto ncol = A.cols();
//        auto nrow = A.rows();
//        for(int ic = 0; ic < ncol; ++ ic)
//        {
//            A.col(ic).normalize();
//            Float cc = 0;
//            for(int ir = 0; ir < nrow; ++ir)
//            {
//                auto v = A(ir, ic);
//                cc += v * v;
//            }
//            cc = sqrt(cc);
//            Float cc2 = 1.0 / cc;
//            for(int ir = 0; ir < nrow; ++ir)
//                A(ir,ic) *= cc2;
//        }
//        cout << "Normalized" << endl;
        
    }
    
//    Float get_lambdaMax()
//    {
//        Vector v = A.transpose() * b;
//        Float vm = fabs(v(0));
//        for(int i = 1; i < v.size(); ++i)
//        {
//            Float vi = fabs(v(i));
//            if(vm < vi)
//                vm = vi;
//        }
//        return vm;
//    }

    void loadb(const string& fn)
    {
        Float val;
        ifstream ifs(fn.c_str());
        int n;
        ifs >> n;
        b.resize(n);
        for(int i = 0; i < n; ++i)
        {
            ifs >> val;
            b(i) = val;
        }
        ifs.close();
//        cout << "load b ok" << endl;
//        cout << "b(n)" << b(n-1) << endl;
    }
    
    const Matrix& getA()const{return A;}
        
    Matrix A;
    Vector b;
};

void put_vec(const Vector& x);
void load_vec(const string &fn, Vector& x);
Float load_tau();
void set_vec(Vector& x, Float val);


int main(int argc, char *argv[])
{
    // set f, the differential part
    MyFuncF f;
    f.loadA("A.txt");
    f.loadb("b.txt");
    int dim = f.getA().cols();
    
    
    // set r_i, the non-smooth part
    vector<FuncR> ri{FuncR(FuncR::Norm1, dim)};
    
    
    // set par    
    BCDSolverParam par;
    par.tau = load_tau();
    par.maxIter = 80000;
    par.objTol = 1e-6;
    par.xTol = 1e-6;
    
    //
    BCDSolver sol(par, f, ri);
    Vector x0(dim);
    set_vec(x0,0);
    sol.init_x(x0);
    sol.solve();
    auto x_re = sol.result();
    put_vec(x_re);
        
}

void set_vec(Vector& x, Float val)
{
    for(auto& t:x) t = val;
}

void put_vec(const Vector& x)
{
    ofstream ofs("x_re.txt");
    ofs << x.size() << endl;
    for(int i = 0; i < x.size(); ++i)
        ofs << x(i) << endl;
    ofs.close();
}

void load_vec(const string& fn, Vector& x)
{
    Float val;
    ifstream ifs(fn.c_str());
    int n;
    ifs >> n;
    x.resize(n);
    for(int i = 0; i < n; ++i)
    {
        ifs >> val;
        x(i) = val;
    }
    ifs.close();
}

Float load_tau()
{
    Float tau;
    ifstream ifs("tau.txt");
    ifs >> tau;
    ifs.close();
    return tau;
}
